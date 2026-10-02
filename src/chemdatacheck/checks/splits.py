"""Train/test split and leakage checks.

Covers 2D-identity leakage, analog leakage (Tc>=threshold), scaffold overlap,
near-dup across splits, and suspiciously-easy-split heuristics.
"""
from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Any

from ..models import Finding, Severity

if TYPE_CHECKING:
    from ..molecules import MoleculeRecord

MAX_ROWS = 500
DEFAULT_TRAIN_VALUES = {"train", "training", "fit"}
DEFAULT_TEST_VALUES = {"test", "testing", "valid", "validation", "val", "eval", "evaluation", "holdout"}


def _clean_split_value(value: Any) -> str:
    """Normalize a split value for case-insensitive comparison."""
    return str(value).strip().casefold()


def _split_value_sets(
    records: list[MoleculeRecord], ctx: dict[str, Any]
) -> tuple[set[str], set[str], set[str]]:
    """Resolve train/test values and return any values left unassigned.

    When only one side is configured, all other observed split values are
    assigned to the opposite side. This supports fold columns naturally:
    ``--test-value 0`` treats every other fold as training data.
    """
    observed = {_clean_split_value(r.split) for r in records if r.split is not None
                and _clean_split_value(r.split)}
    requested_train = ctx.get("train_values")
    requested_test = ctx.get("test_values")
    train = {_clean_split_value(v) for v in requested_train or []}
    test = {_clean_split_value(v) for v in requested_test or []}
    if requested_train is None and requested_test is None:
        train = observed & DEFAULT_TRAIN_VALUES
        test = observed & DEFAULT_TEST_VALUES
    elif requested_train is None:
        train = observed - test
    elif requested_test is None:
        test = observed - train
    return train, test, observed - train - test


def _split_groups(
    records: list[MoleculeRecord], ctx: dict[str, Any]
) -> tuple[list[MoleculeRecord], list[MoleculeRecord]]:
    """Split valid records into train and test/validation groups.

    Args:
        records: Per-molecule records.

    Returns:
        Tuple of (train records, test/validation records).
    """
    train_values, test_values, _ = _split_value_sets(records, ctx)
    tr = [r for r in records if r.valid and _clean_split_value(r.split) in train_values]
    te = [r for r in records if r.valid and _clean_split_value(r.split) in test_values]
    return tr, te


def _no_split_info() -> Finding:
    """Build the informational finding used when split checks are skipped.

    Returns:
        Finding explaining how to enable leakage checks.
    """
    return Finding(check_id="split_info", severity=Severity.INFO,
                   title="no split information", count=0, rate=0.0,
                   affected_rows=[], total_affected=0, examples=[],
                   recommendation="Pass --split-col or two files (train.csv test.csv) to enable leakage checks.",
                   details="Split/leakage checks skipped.")


def check_split_configuration(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Report missing, incomplete, or partially ignored split configuration."""
    observed = sorted({_clean_split_value(r.split) for r in records if r.split is not None
                       and _clean_split_value(r.split)})
    if not observed:
        return _no_split_info()
    train_values, test_values, ignored = _split_value_sets(records, ctx)
    tr, te = _split_groups(records, ctx)
    if tr and te and not ignored:
        return None
    examples = [{"observed_values": observed,
                 "train_values": sorted(train_values),
                 "test_values": sorted(test_values),
                 "ignored_values": sorted(ignored)}]
    if not train_values or not test_values or not tr or not te:
        title = "split values could not be resolved into non-empty train and test groups"
        recommendation = ("Pass --train-value/--test-value explicitly. For k-fold data, "
                          "pass the held-out fold with --test-value; all other folds become train.")
    else:
        title = f"{len(ignored)} split value(s) ignored by leakage checks"
        recommendation = ("Map every intended split with --train-value or --test-value, or remove rows "
                          "that should not participate in leakage checks.")
    ignored_rows = [r.row_id for r in records if _clean_split_value(r.split) in ignored]
    return Finding(check_id="split_configuration", severity=Severity.WARNING,
                   title=title, count=len(ignored_rows),
                   rate=len(ignored_rows) / max(1, len(records)),
                   affected_rows=ignored_rows[:MAX_ROWS], total_affected=len(ignored_rows), examples=examples,
                   recommendation=recommendation,
                   details="Leakage results are incomplete until split values are resolved.")


def check_scaffold_overlap(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag test molecules whose scaffold was seen in train.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for scaffold overlap, or None if no split or no overlap.
    """
    tr, te = _split_groups(records, ctx)
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


def check_identity_leakage(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag test molecules identical (connectivity) to train molecules.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Error finding for 2D-identity leakage, or None if none found.
    """
    tr, te = _split_groups(records, ctx)
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


def _cross_similarity_profile(
    records: list[MoleculeRecord], ctx: dict[str, Any], max_test: int = 3000,
) -> tuple[
    list[tuple[MoleculeRecord, float, MoleculeRecord | None]],
    int,
    int,
]:
    """Compute and cache each test molecule's closest training neighbor."""
    if "_cross_similarity_profile" in ctx:
        return ctx["_cross_similarity_profile"]
    from rdkit import DataStructs
    tr, te = _split_groups(records, ctx)
    tr = [r for r in tr if r.fp is not None]
    te = [r for r in te if r.fp is not None]
    te_use = te if len(te) <= max_test else te[:: max(1, len(te) // max_test)][:max_test]
    tr_fps = [r.fp for r in tr]
    per_test_max = []
    for r in te_use:
        sims = DataStructs.BulkTanimotoSimilarity(r.fp, tr_fps)
        best_i = max(range(len(sims)), key=lambda i: sims[i]) if sims else None
        best = float(sims[best_i]) if best_i is not None else 0.0
        per_test_max.append((r, best, tr[best_i] if best_i is not None else None))
    result = (per_test_max, len(te_use), len(te))
    ctx["_cross_similarity_profile"] = result
    if len(te_use) < len(te):
        ctx.setdefault("approximations", []).append({
            "check_ids": ["analog_leakage", "cross_split_near_duplicates", "suspiciously_easy_split"],
            "method": "test_sampling", "test_limit": max_test, "test_records": len(te)})
    return result


def _pairs_above(
    per_test_max: list[tuple[MoleculeRecord, float, MoleculeRecord | None]],
    threshold: float,
    cap: int = 3000,
) -> list[tuple[MoleculeRecord, MoleculeRecord, float]]:
    """Filter a cross-split similarity profile at a threshold."""
    pairs = [(mate, test, similarity) for test, similarity, mate in per_test_max
             if similarity >= threshold and mate is not None]
    pairs.sort(key=lambda t: -t[2])
    return pairs[:cap]


def check_analog_leakage(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag test molecules with a close analog in train.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples` and optional
            `analog_thresh`.

    Returns:
        Finding for analog leakage, or None if no split or no pairs.
    """
    thresh = float(ctx.get("analog_thresh", 0.6))
    per_test_max, n_used, n_te = _cross_similarity_profile(records, ctx)
    pairs = _pairs_above(per_test_max, thresh)
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
                   details=f"ECFP4/Morgan Tc≥{thresh} cross-split pairs: {len(pairs)} enumerated (capped).",
                   metadata={"approximate": n_used < n_te, "method": "test_sampling" if n_used < n_te else "all_test",
                             "test_records_scored": n_used, "test_records_total": n_te})


def check_cross_near_duplicates(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag near-duplicate pairs spanning train and test.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples` and optional
            `near_dup_thresh`.

    Returns:
        Error finding for cross-split near-duplicates, or None if none found.
    """
    thresh = float(ctx.get("near_dup_thresh", 0.95))
    per_test_max, n_used, n_te = _cross_similarity_profile(records, ctx)
    pairs = _pairs_above(per_test_max, thresh)
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
                   details=f"{len(pairs)} cross-split pairs above threshold.",
                   metadata={"approximate": n_used < n_te, "method": "test_sampling" if n_used < n_te else "all_test",
                             "test_records_scored": n_used, "test_records_total": n_te})


def check_easy_split(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag splits where test molecules are unusually close to train.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for suspiciously easy splits, or None if the split looks hard.
    """
    per_test_max, n_used, n_te = _cross_similarity_profile(records, ctx)
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
                   details=f"mean={mean:.2f} median={med:.2f} over {len(sims)} test molecules.",
                   metadata={"approximate": n_used < n_te, "method": "test_sampling" if n_used < n_te else "all_test",
                             "test_records_scored": n_used, "test_records_total": n_te})
CHECKS = [check_split_configuration, check_scaffold_overlap, check_identity_leakage, check_analog_leakage,
          check_cross_near_duplicates, check_easy_split]
