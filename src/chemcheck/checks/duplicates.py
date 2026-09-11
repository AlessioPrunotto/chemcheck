"""Dataset-level duplicate checks: exact → canonical → stereo → salt → tautomer → near."""
from __future__ import annotations

from collections import defaultdict

from ..models import Finding, Severity

MAX_ROWS = 500


def _group_finding(check_id, severity, title, groups: dict, n_total, max_examples, rec, details=""):
    # groups: key -> list[row_id]; keep only len>1
    dup_groups = {k: v for k, v in groups.items() if len(v) > 1}
    if not dup_groups:
        return None
    affected = sorted({rid for v in dup_groups.values() for rid in v})
    examples = [{"key": str(k)[:100], "rows": v[:8], "n": len(v)}
                for k, v in sorted(dup_groups.items(), key=lambda kv: -len(kv[1]))[:max_examples]]
    return Finding(check_id=check_id, severity=severity, title=title,
                   count=len(affected), rate=len(affected) / max(1, n_total),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=examples, recommendation=rec,
                   details=details or f"{len(dup_groups)} duplicate groups.")


def check_exact_duplicates(records, ctx):
    g = defaultdict(list)
    for r in records:
        g[r.raw_smiles.strip()].append(r.row_id)
    return _group_finding("exact_duplicates", Severity.WARNING, "exact duplicate rows",
                          g, len(records), ctx["max_examples"],
                          "Deduplicate raw SMILES first — exact dups inflate N and leak across random splits.")


def check_canonical_duplicates(records, ctx):
    g = defaultdict(list)
    for r in records:
        if r.valid and r.canon_smi:
            g[r.canon_smi].append(r.row_id)
    return _group_finding("canonical_duplicates", Severity.WARNING,
                          "duplicate structures (canonical SMILES)", g, len(records),
                          ctx["max_examples"],
                          "Deduplicate on canonical isomeric SMILES; keep one row per structure (aggregate labels).")


def check_stereo_collisions(records, ctx):
    g = defaultdict(set)   # connectivity -> set of distinct isomeric smi
    members = defaultdict(list)
    for r in records:
        if r.valid and r.connectivity_smi and r.canon_smi:
            g[r.connectivity_smi].add(r.canon_smi)
            members[r.connectivity_smi].append(r.row_id)
    colliding = {k: v for k, v in members.items() if len(g[k]) > 1}
    if not colliding:
        return None
    affected = sorted({rid for v in colliding.values() for rid in v})
    examples = [{"connectivity": str(k)[:100], "rows": v[:8],
                 "distinct_stereoforms": len(g[k])} for k, v in
                sorted(colliding.items(), key=lambda kv: -len(kv[1]))[:ctx["max_examples"]]]
    return Finding(check_id="stereochemical_collisions", severity=Severity.WARNING,
                   title="stereochemical collisions", count=len(affected),
                   rate=len(affected) / max(1, len(records)),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=examples,
                   recommendation="Same connectivity with different stereo = different data points only if "
                                  "labels are stereo-dependent. Otherwise collapse to one form and note ambiguity.",
                   details=f"{len(colliding)} connectivity groups with >1 stereoisomer.")


def check_salt_duplicates(records, ctx):
    g = defaultdict(list)
    for r in records:
        if r.valid and r.parent_smi:
            g[r.parent_smi].append(r.row_id)
    # only interesting where parents collide but full structures differ
    full = defaultdict(set)
    for r in records:
        if r.valid and r.parent_smi and r.canon_smi:
            full[r.parent_smi].add(r.canon_smi)
    salt_groups = {k: v for k, v in g.items() if len(v) > 1 and len(full[k]) > 1}
    if not salt_groups:
        return None
    affected = sorted({rid for v in salt_groups.values() for rid in v})
    examples = [{"parent": str(k)[:100], "rows": v[:8]} for k, v in
                sorted(salt_groups.items(), key=lambda kv: -len(kv[1]))[:ctx["max_examples"]]]
    return Finding(check_id="salt_duplicates", severity=Severity.WARNING,
                   title="salts inconsistently represented", count=len(affected),
                   rate=len(affected) / max(1, len(records)),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=examples,
                   recommendation="Apply one salt policy (e.g. strip to parent) and re-run; mixed salt forms "
                                  "fragment SAR and leak across splits.",
                   details=f"{len(salt_groups)} parent structures with multiple salt/solvate forms.")


def check_tautomer_duplicates(records, ctx):
    g = defaultdict(list)
    for r in records:
        if r.valid and r.tautomer_smi:
            g[r.tautomer_smi].append(r.row_id)
    canon = defaultdict(set)
    for r in records:
        if r.valid and r.tautomer_smi and r.canon_smi:
            canon[r.tautomer_smi].add(r.canon_smi)
    tg = {k: v for k, v in g.items() if len(v) > 1 and len(canon[k]) > 1}
    if not tg:
        return None
    affected = sorted({rid for v in tg.values() for rid in v})
    examples = [{"canonical_tautomer": str(k)[:100], "rows": v[:8]} for k, v in
                sorted(tg.items(), key=lambda kv: -len(kv[1]))[:ctx["max_examples"]]]
    return Finding(check_id="tautomer_duplicates", severity=Severity.INFO,
                   title="tautomer duplicates", count=len(affected),
                   rate=len(affected) / max(1, len(records)),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=examples,
                   recommendation="Canonicalize tautomers before dedup/splitting "
                                  "so keto/enol pairs are not treated as independent.",
                   details=f"{len(tg)} tautomer groups with multiple input forms.")


def _valid_with_fp(records):
    return [r for r in records if r.valid and r.fp is not None]


def check_near_duplicates(records, ctx):
    from rdkit import DataStructs
    thresh = float(ctx.get("near_dup_thresh", 0.95))
    valid = _valid_with_fp(records)
    n = len(valid)
    if n < 2:
        return None
    pairs: list[tuple] = []
    # Full comparison for small sets; popcount-windowed for large sets (thin but scalable).
    if n <= 4000:
        fps = [r.fp for r in valid]
        for i in range(n):
            sims = DataStructs.BulkTanimotoSimilarity(valid[i].fp, fps[i + 1:])
            for j, s in enumerate(sims):
                if s >= thresh:
                    pairs.append((valid[i].row_id, valid[i + j + 1].row_id, round(float(s), 3)))
                    if len(pairs) >= 5000:
                        break
            if len(pairs) >= 5000:
                break
    else:
        order = sorted(valid, key=lambda r: r.fp.GetNumOnBits())
        window = 400
        for i in range(n):
            cand = order[i + 1:i + 1 + window]
            if not cand:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(order[i].fp, [c.fp for c in cand])
            for c, s in zip(cand, sims):
                if s >= thresh:
                    pairs.append((order[i].row_id, c.row_id, round(float(s), 3)))
                    if len(pairs) >= 5000:
                        break
            if len(pairs) >= 5000:
                break
    if not pairs:
        return None
    affected = sorted({p[0] for p in pairs} | {p[1] for p in pairs})
    examples = [{"row_a": a, "row_b": b, "tanimoto": s} for a, b, s in pairs[:ctx["max_examples"]]]
    f = Finding(check_id="near_duplicates", severity=Severity.WARNING,
                title=f"near-duplicates (Morgan Tanimoto ≥ {thresh})", count=len(affected),
                rate=len(affected) / max(1, len(records)),
                affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                examples=examples,
                recommendation="Cluster near-duplicates and keep one representative per cluster (or put the "
                               "whole cluster in one split) — otherwise random splits look deceptively easy.",
                details=f"{len(pairs)}+ pairs above threshold (capped at 5000 enumerated).")
    return f


CHECKS = [check_exact_duplicates, check_canonical_duplicates, check_stereo_collisions,
          check_salt_duplicates, check_tautomer_duplicates, check_near_duplicates]
