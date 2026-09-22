"""Chemical-space checks: elements, functional groups, rings, bias, domain gaps."""
from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING, Any

from ..models import Finding, Severity

if TYPE_CHECKING:
    from ..molecules import MoleculeRecord

MAX_ROWS = 500
COMMON_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}

FG_PANEL = [
    ("nitro", "[N+](=O)[O-]"),
    ("peroxide", "OO"),
    ("azide", "[N-]=[N+]=N"),
    ("epoxide", "C1OC1"),
    ("aldehyde", "[CH]=O"),
    ("thiourea", "NC(=S)N"),
    ("phosphonate", "P(=O)(O)O"),
    ("boronic_acid", "B(O)O"),
    ("hydrazine", "NN"),
    ("isocyanate", "N=C=O"),
    ("acyl_halide", "C(=O)[F,Cl,Br,I]"),
    ("metal", "[Na,K,Ca,Mg,Fe,Zn,Cu,Pd,Pt,Li,Al,Sn]"),
]

PAINS_LIKE = [
    ("catechol", "c1ccc(O)c(O)c1"),
    ("quinone", "C1=CC(=O)C=CC1=O"),
    ("rhodanine_like", "S=C1NCC(=O)S1"),
    ("azo", "cN=Nc"),
]


def check_rare_elements(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag molecules containing elements outside the common organic set.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Informational finding for rare elements, or None if none found.
    """
    bad = [r for r in records if r.valid and (r.elements - COMMON_ELEMENTS)]
    if not bad:
        return None
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": (r.canon_smi or "")[:80],
           "rare": sorted(r.elements - COMMON_ELEMENTS)} for r in bad][:ctx["max_examples"]]
    return Finding(check_id="rare_elements", severity=Severity.INFO,
                   title=f"rare elements in {len(bad)} molecules",
                   count=len(bad), rate=len(bad) / max(1, len(records)),
                   affected_rows=rows[:MAX_ROWS], total_affected=len(bad),
                   examples=ex,
                   recommendation="Confirm rare elements are intentional (organometallics, counterions); most "
                                  "descriptors/embeddings handle them poorly — consider a separate subset.",
                   details="")


def check_functional_groups(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag rare (<0.5%) or PAINS-like functional groups.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Informational finding for functional groups, or None if none flagged.
    """
    try:
        from rdkit import Chem
    except Exception:
        return None
    valid = [r for r in records if r.valid and r.mol is not None]
    n = max(1, len(valid))
    patt_hits = []
    for name, smarts in FG_PANEL + PAINS_LIKE:
        try:
            patt = Chem.MolFromSmarts(smarts)
            hits = [r.row_id for r in valid if r.mol.HasSubstructMatch(patt)]
        except Exception:
            hits = []
        patt_hits.append((name, hits))
    rare = [(name, hits) for name, hits in patt_hits if 0 < len(hits) / n < 0.005]
    pains = [(name, hits) for name, hits in patt_hits if name in dict(PAINS_LIKE) and hits]
    examples = ([{"group": name, "count": len(h), "rows": h[:5]} for name, h in rare]
                + [{"group": name + " (PAINS-like)", "count": len(h), "rows": h[:5]} for name, h in pains])
    flagged_rows = sorted({rid for _, h in (rare + pains) for rid in h})
    if not flagged_rows:
        return None
    prev = {name: round(len(h) / n, 4) for name, h in patt_hits if h}
    return Finding(check_id="rare_functional_groups", severity=Severity.INFO,
                   title="rare or promiscuous functional groups",
                   count=len(flagged_rows), rate=len(flagged_rows) / max(1, len(records)),
                   affected_rows=flagged_rows[:MAX_ROWS], total_affected=len(flagged_rows),
                   examples=(examples[:ctx["max_examples"]]),
                   recommendation="Rare groups are outside the applicability domain for most models; PAINS-like "
                                  "hits deserve assay-interference review, not silent training.",
                   details=f"Prevalence among valid: {prev}.")


def check_ring_systems(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag macrocycles and singleton scaffolds.

    Args:
        records: Per-molecule records.
        ctx: Check context (unused beyond record access).

    Returns:
        Finding for unusual ring systems, or None if none found.
    """
    valid = [r for r in records if r.valid]
    macro = [r.row_id for r in valid if r.max_ring_size >= 8]
    scaf_counts = Counter(r.scaffold for r in valid if r.scaffold)
    singletons = sum(1 for c in scaf_counts.values() if c == 1)
    rows = sorted(set(macro))
    examples = []
    if macro:
        examples.append({"kind": "macrocycle (ring>=8)", "count": len(macro),
                         "rows": macro[:5]})
    if len(valid) >= 10:
        examples.append({"kind": "singleton scaffolds", "count": singletons,
                         "scaffolds": [s[:60] for s, c in scaf_counts.items() if c == 1][:5]})
    if not macro and (singletons == 0 or len(valid) < 10):
        return None
    sev = Severity.INFO if len(macro) < 0.05 * max(1, len(valid)) else Severity.WARNING
    n_single = singletons if len(valid) >= 10 else 0
    return Finding(check_id="unusual_ring_systems", severity=sev,
                   title="unusual ring systems",
                   count=len(rows) + n_single, rate=(len(rows) + n_single) / max(1, len(records)),
                   affected_rows=rows[:MAX_ROWS], total_affected=len(rows),
                   examples=examples,
                   recommendation="Macrocycles and singleton scaffolds are generalization tests, not training "
                                  "bulk — keep them, but evaluate them as a separate slice.",
                   details=f"{singletons} singleton scaffolds; {len(macro)} macrocycles.")


def check_representation_bias(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag datasets dominated by a few scaffolds.

    Args:
        records: Per-molecule records.
        ctx: Check context (unused beyond record access).

    Returns:
        Finding when top-5 scaffolds cover >=30% of data, else None.
    """
    valid = [r for r in records if r.valid and r.scaffold]
    if len(valid) < 30:
        return None
    scaf_counts = Counter(r.scaffold for r in valid)
    top5 = scaf_counts.most_common(5)
    top5_cover = sum(c for _, c in top5) / max(1, len(valid))
    if top5_cover < 0.3:
        return None
    ex = [{"scaffold": s[:80], "count": c} for s, c in top5]
    return Finding(check_id="representation_bias", severity=Severity.WARNING,
                   title=f"representation bias: top-5 scaffolds cover {top5_cover * 100:.1f}% of data",
                   count=sum(c for _, c in top5), rate=top5_cover,
                   affected_rows=[r.row_id for r in valid if r.scaffold in dict(top5)][:MAX_ROWS],
                   total_affected=sum(c for _, c in top5),
                   examples=ex,
                   recommendation="Rebalance or report per-scaffold performance; a model can score well by "
                                  "memorizing 5 chemotypes while failing everywhere else.",
                   details="")


def check_applicability_gaps(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag molecules isolated in chemical space (nearest-neighbor Tc < 0.3).

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Informational finding for isolated molecules, or None if none found.
    """
    from rdkit import DataStructs
    valid = [r for r in records if r.valid and r.fp is not None]
    n = len(valid)
    if n < 10:
        return None
    # For large sets, sample only to find a superset of possible gaps, then
    # verify those candidates against every molecule. This retains scalability
    # without turning reference-sampling misses into false alerts.
    ref = valid if n <= 3000 else valid[:: max(1, n // 3000)][:3000]
    ref_fps = [r.fp for r in ref]
    candidates = []
    for r in valid:
        sims = DataStructs.BulkTanimotoSimilarity(r.fp, ref_fps)
        # exclude self-match (=1.0) by taking second max
        top2 = sorted(sims, reverse=True)[:2]
        nn = top2[1] if len(top2) > 1 and top2[0] >= 0.999 else (top2[0] if top2 else 0.0)
        if nn < 0.3:
            candidates.append((r, float(nn)))
    if n > 3000:
        all_fps = [r.fp for r in valid]
        positions = {id(r): index for index, r in enumerate(valid)}
        sparse = []
        for r, _ in candidates:
            sims = DataStructs.BulkTanimotoSimilarity(r.fp, all_fps)
            sims[positions[id(r)]] = -1.0
            nn = max(sims, default=0.0)
            if nn < 0.3:
                sparse.append((r.row_id, round(float(nn), 3)))
    else:
        sparse = [(r.row_id, round(nn, 3)) for r, nn in candidates]
    if not sparse:
        return None
    sparse.sort(key=lambda t: t[1])
    return Finding(check_id="applicability_gaps", severity=Severity.INFO,
                   title=f"{len(sparse)} molecules isolated in chemical space (nearest-neighbor Tc<0.3)",
                   count=len(sparse), rate=len(sparse) / max(1, len(records)),
                   affected_rows=[rid for rid, _ in sparse][:MAX_ROWS], total_affected=len(sparse),
                   examples=[{"row": rid, "max_sim_to_rest": s} for rid, s in sparse[:ctx["max_examples"]]],
                   recommendation="Treat isolated molecules as domain-gap probes: expect poor predictions there; "
                                  "consider acquiring analogs or flagging them at inference time.",
                   details="", metadata={"approximate": False,
                                         "method": ("reference_prefilter_exact_refinement"
                                                    if n > 3000 else "all_pairs"),
                                         "records": n, "reference_records": len(ref)})


CHECKS = [check_rare_elements, check_functional_groups, check_ring_systems,
          check_representation_bias, check_applicability_gaps]
