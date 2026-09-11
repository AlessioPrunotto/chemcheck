"""Train/test split and leakage checks.

Covers 2D-identity leakage, analog leakage (Tc>=threshold), scaffold overlap,
near-dup across splits, and suspiciously-easy-split heuristics.
"""
from __future__ import annotations

from collections import defaultdict

from ..models import Finding, Severity

MAX_ROWS = 500


def _split_groups(records):
    tr = [r for r in records if r.valid and (r.split or "").lower() == "train"]
    te = [r for r in records if r.valid and (r.split or "").lower() in ("test", "valid", "validation")]
    return tr, te


def _no_split_info():
    return Finding(check_id="split_info", severity=Severity.INFO,
                   title="no split information", count=0, rate=0.0,
                   affected_rows=[], total_affected=0, examples=[],
                   recommendation="Pass --split-col or two files (train.csv test.csv) to enable leakage checks.",
                   details="Split/leakage checks skipped.")


def check_scaffold_overlap(records, ctx):
    tr, te = _split_groups(records)
    if not tr or not te:
        return None
    tr_scaf = {r.scaffold for r in tr if r.scaffold}
    overlap_rows = [r.row_id for r in te if r.scaffold in tr_scaf]
    if not overlap_rows:
        return None
    rate = len(overlap_rows) / max(1, len(te))
    ex = [{"test_row": r.row_id, "scaffold": (r.scaffold or "")[:80]} for r in te
          if r.scaffold in tr_scaf][:ctx["max_examples"]]
    return Finding(check_id="scaffold_overlap", severity=Severity.WARNING,
                   title=f"scaffold split overlap ({rate * 100:.1f}% of test scaffolds seen in train)",
                   count=len(overlap_rows), rate=len(overlap_rows) / max(1, len(records)),
                   affected_rows=overlap_rows[:MAX_ROWS], total_affected=len(overlap_rows),
                   examples=ex,
                   recommendation="Use a true scaffold split (or DataSAIL/UMAP clustering) so test scaffolds "
                                  "are unseen; overlapping scaffolds reward memorization.",
                   details=f"{len(overlap_rows)}/{len(te)} test molecules share a Bemis–Murcko scaffold with train.")


def check_identity_leakage(records, ctx):
    tr, te = _split_groups(records)
    if not tr or not te:
        return None
    # stereo-agnostic identity: connectivity (parent) match across splits
    tr_keys = defaultdict(list)
    for r in tr:
        if r.connectivity_smi:
            tr_keys[r.connectivity_smi].append(r.row_id)
    pairs = []
    leaked_test = set()
    for r in te:
        if r.connectivity_smi in tr_keys:
            leaked_test.add(r.row_id)
            pairs.append({"train_rows": tr_keys[r.connectivity_smi][:4],
                          "test_row": r.row_id, "connectivity": (r.connectivity_smi or "")[:100]})
    if not pairs:
        return None
    return Finding(check_id="identity_leakage", severity=Severity.ERROR,
                   title=f"2D-identity leakage ({len(leaked_test)} test molecules seen in train)",
                   count=len(leaked_test), rate=len(leaked_test) / max(1, len(records)),
                   affected_rows=sorted(leaked_test)[:MAX_ROWS], total_affected=len(leaked_test),
                   examples=pairs[:ctx["max_examples"]],
                   recommendation="Remove leaked test molecules or move them to train; identical structures "
                                  "must never span splits.",
                   details="Stereo-agnostic match (connectivity SMILES) across train/test.")


def _cross_similarities(tr, te, thresh, cap_pairs=3000, max_test=3000):
    from rdkit import DataStructs
    te_use = te if len(te) <= max_test else te[:: max(1, len(te) // max_test)][:max_test]
    tr_fps = [r.fp for r in tr]
    pairs = []
    per_test_max = []
    for r in te_use:
        if r.fp is None:
            per_test_max.append((r, 0.0, None))
            continue
        sims = DataStructs.BulkTanimotoSimilarity(r.fp, tr_fps)
        best_i = max(range(len(sims)), key=lambda i: sims[i]) if sims else None
        best = float(sims[best_i]) if best_i is not None else 0.0
        per_test_max.append((r, best, tr[best_i] if best_i is not None else None))
    for r, best, mate in per_test_max:
        if best >= thresh and mate is not None:
            pairs.append((mate, r, best))
    pairs.sort(key=lambda t: -t[2])
    return pairs[:cap_pairs], per_test_max, len(te_use), len(te)


def check_analog_leakage(records, ctx):
    tr, te = _split_groups(records)
    if not tr or not te or not any(r.fp is not None for r in tr):
        return None
    thresh = float(ctx.get("analog_thresh", 0.6))
    pairs, per_test_max, n_used, n_te = _cross_similarities(
        [r for r in tr if r.fp is not None], [r for r in te if r.fp is not None], thresh)
    if not pairs:
        return None
    leaked = sorted({b.row_id for _, b, _ in pairs})
    import statistics
    med = statistics.median([s for _, s, _ in per_test_max]) if per_test_max else 0.0
    examples = []
    for mate, r, s in pairs[:ctx["max_examples"]]:
        examples.append({
            "train_row": mate.row_id, "test_row": r.row_id,
            "tanimoto": round(s, 3),
            "train_smiles": (mate.canon_smi or "")[:100],
            "test_smiles": (r.canon_smi or "")[:100],
            "shared_scaffold": (mate.scaffold if mate.scaffold == r.scaffold else None),
            "recommendation": "move test molecule to train OR use scaffold/cluster split",
        })
    note = "" if n_used == n_te else f" (estimated on {n_used}/{n_te} sampled test molecules)"
    return Finding(check_id="analog_leakage", severity=Severity.WARNING,
                   title=f"analog leakage: {len(leaked)} test molecules with Tc≥{thresh} to train{note}",
                   count=len(leaked), rate=len(leaked) / max(1, len(records)),
                   affected_rows=leaked[:MAX_ROWS], total_affected=len(leaked),
                   examples=examples,
                   recommendation="Inspect the worst pairs below; regroup analog series into the same split. "
                                  "Median max test→train similarity is "
                                  f"{med:.2f} — above ~0.6 the split is easier than it looks.",
                   details=f"ECFP4/Morgan Tc≥{thresh} cross-split pairs: {len(pairs)} enumerated (capped).")


def check_cross_near_duplicates(records, ctx):
    tr, te = _split_groups(records)
    if not tr or not te:
        return None
    thresh = float(ctx.get("near_dup_thresh", 0.95))
    pairs, _, n_used, n_te = _cross_similarities(
        [r for r in tr if r.fp is not None], [r for r in te if r.fp is not None], thresh)
    if not pairs:
        return None
    leaked = sorted({b.row_id for _, b, _ in pairs} | {a.row_id for a, _, _ in pairs})
    examples = [{"train_row": a.row_id, "test_row": b.row_id, "tanimoto": round(s, 3)}
                for a, b, s in pairs[:ctx["max_examples"]]]
    note = "" if n_used == n_te else f" (estimated on {n_used}/{n_te} sampled test molecules)"
    return Finding(check_id="cross_split_near_duplicates", severity=Severity.ERROR,
                   title=f"near-duplicates across splits (Tc≥{thresh}){note}",
                   count=len(leaked), rate=len(leaked) / max(1, len(records)),
                   affected_rows=leaked[:MAX_ROWS], total_affected=len(leaked),
                   examples=examples,
                   recommendation="Near-duplicates must live in one split only — move the test copy to train "
                                  "or drop it; otherwise metrics measure memorization.",
                   details=f"{len(pairs)} cross-split pairs above threshold.")


def check_easy_split(records, ctx):
    tr, te = _split_groups(records)
    if not tr or not te:
        return None
    _, per_test_max, _, _ = _cross_similarities(
        [r for r in tr if r.fp is not None], [r for r in te if r.fp is not None], 2.0)
    import statistics
    sims = [s for _, s, _ in per_test_max]
    if not sims:
        return None
    med = statistics.median(sims)
    mean = statistics.mean(sims)
    if med < 0.55:
        return None
    worst = sorted(per_test_max, key=lambda t: -t[1])[:ctx["max_examples"]]
    ex = [{"test_row": r.row_id, "max_train_sim": round(s, 3)} for r, s, _ in worst]
    return Finding(check_id="suspiciously_easy_split", severity=Severity.WARNING,
                   title=f"suspiciously easy split (median max test→train Tc={med:.2f})",
                   count=len([s for s in sims if s >= 0.6]),
                   rate=len([s for s in sims if s >= 0.6]) / max(1, len(sims)),
                   affected_rows=[r.row_id for r, _, _ in worst][:MAX_ROWS],
                   total_affected=len([s for s in sims if s >= 0.6]),
                   examples=ex,
                   recommendation="Prefer scaffold, Butina-cluster, UMAP-cluster, or time-based splits for a "
                                  "realistic generalization estimate; report random-split numbers as optimistic.",
                   details=f"mean={mean:.2f} median={med:.2f} over {len(sims)} test molecules.")


def _has_any_split(records) -> bool:
    return any((r.split or "") != "" and r.split is not None for r in records)


CHECKS = [check_scaffold_overlap, check_identity_leakage, check_analog_leakage,
          check_cross_near_duplicates, check_easy_split]
