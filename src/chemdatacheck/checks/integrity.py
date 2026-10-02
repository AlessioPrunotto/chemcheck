"""Chemical-integrity checks (per-molecule chemistry validity)."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..models import Finding, Severity

if TYPE_CHECKING:
    from ..molecules import MoleculeRecord

MAX_ROWS = 200


def _mk(
    check_id: str,
    severity: Severity,
    title: str,
    rows: list[str],
    examples: list[dict[str, Any]],
    rec: str,
    n_total: int,
    details: str = "",
) -> Finding | None:
    """Build a finding, or None when no rows are affected.

    Args:
        check_id: Stable machine-readable check identifier.
        severity: Severity level for the finding.
        title: Short human-readable title.
        rows: Affected row ids.
        examples: Example payloads illustrating the issue.
        rec: Remediation recommendation.
        n_total: Total number of records (for rate computation).
        details: Extra details about the finding.

    Returns:
        A finding, or None if `rows` is empty.
    """
    if not rows:
        return None
    return Finding(check_id=check_id, severity=severity, title=title,
                   count=len(rows), rate=len(rows) / max(1, n_total),
                   affected_rows=rows[:MAX_ROWS], total_affected=len(rows),
                   examples=examples, recommendation=rec, details=details)


def check_invalid_smiles(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag rows that failed RDKit parsing/sanitization.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for invalid SMILES, or None if all rows are valid.
    """
    rows = [r.row_id for r in records if not r.valid]
    ex = [{"row": r.row_id, "smiles": r.raw_smiles[:80], "error": r.sanitize_error[:160]}
          for r in records if not r.valid][:ctx["max_examples"]]
    return _mk("invalid_smiles", Severity.ERROR, "invalid SMILES", rows, ex,
               "Drop or fix these rows; they cannot be parsed and will silently shrink training data.",
               len(records))


def check_valence(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag invalid rows whose error mentions valence.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for impossible valence states, or None if none found.
    """
    rows = [r.row_id for r in records if not r.valid and "valence" in (r.sanitize_error or "").lower()]
    ex = [{"row": r.row_id, "smiles": r.raw_smiles[:80], "error": r.sanitize_error[:160]}
          for r in records if r.row_id in set(rows)][:ctx["max_examples"]]
    return _mk("valence_error", Severity.ERROR, "impossible valence states", rows, ex,
               "Inspect hypervalent atoms (e.g. pentavalent carbon); usually a drawing/SMILES error.",
               len(records))


def check_aromaticity(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag invalid rows with kekulization/aromaticity errors.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for aromaticity inconsistencies, or None if none found.
    """
    rows = [r.row_id for r in records
            if not r.valid and ("kekul" in (r.sanitize_error or "").lower()
                                or "aromat" in (r.sanitize_error or "").lower())]
    ex = [{"row": r.row_id, "smiles": r.raw_smiles[:80], "error": r.sanitize_error[:160]}
          for r in records if r.row_id in set(rows)][:ctx["max_examples"]]
    return _mk("aromaticity_suspect", Severity.WARNING, "aromaticity/kekulization inconsistencies",
               rows, ex, "Re-kekulize with RDKit or redraw aromatic systems; check lowercase aromatic symbols.",
               len(records))


def check_impossible_charge(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules with suspicious formal charges.

    Flags molecules with |formal charge| > 4 or carbon atoms with
    |charge| > 1.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for suspicious charges, or None if none found.
    """
    bad = []
    for r in records:
        if not r.valid:
            continue
        if abs(r.formal_charge) > 4:
            bad.append(r)
            continue
        # carbon with |charge|>1 is almost always wrong
        try:
            for a in r.mol.GetAtoms():
                if a.GetSymbol() == "C" and abs(a.GetFormalCharge()) > 1:
                    bad.append(r)
                    break
        except Exception:
            pass
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi, "formal_charge": r.formal_charge} for r in bad][:ctx["max_examples"]]
    return _mk("impossible_charge", Severity.WARNING, "impossible or suspicious charges", rows, ex,
               "Neutralize salts/solvents or correct formal charges; verify nitro/quaternary N patterns.",
               len(records))


def check_disconnected(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules with multiple disconnected fragments.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for salts/mixtures, or None if none found.
    """
    bad = [r for r in records if r.valid and r.n_fragments > 1]
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi, "n_fragments": r.n_fragments} for r in bad][:ctx["max_examples"]]
    return _mk("disconnected_components", Severity.WARNING, "disconnected components (salts/mixtures)",
               rows, ex, "Decide a salt policy: keep parent only (desalt) or encode counterions "
                         "explicitly — be consistent.",
               len(records))


def check_isotopes(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules containing isotopically labeled atoms.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Informational finding for isotopes, or None if none found.
    """
    bad = [r for r in records if r.valid and r.has_isotope]
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi} for r in bad][:ctx["max_examples"]]
    return _mk("isotope_flag", Severity.INFO, "isotopically labeled atoms", rows, ex,
               "Usually harmless, but confirm labels are intentional — fingerprints may treat them as distinct.",
               len(records))


def check_radicals(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules with radical electrons.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for radicals, or None if none found.
    """
    bad = [r for r in records if r.valid and r.has_radical]
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi} for r in bad][:ctx["max_examples"]]
    return _mk("radical_flag", Severity.WARNING, "radical electrons present", rows, ex,
               "Radicals are rarely valid training data for property models; verify or remove.",
               len(records))


def check_unspecified_stereo(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules with unassigned tetrahedral centers.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Finding for unspecified stereocenters, or None if none found.
    """
    bad = [r for r in records if r.valid and r.n_stereo_unassigned > 0]
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi,
           "unassigned_centers": r.n_stereo_unassigned} for r in bad][:ctx["max_examples"]]
    f = _mk("unspecified_stereo", Severity.WARNING, "molecules with unspecified stereocenters",
            rows, ex, "Assign stereochemistry where known; otherwise record that stereo is unknown "
                      "and use stereo-agnostic evaluation to avoid inflated scores.",
            len(records))
    if f is not None:
        f.details = f"{len(rows)} molecules have ≥1 unassigned tetrahedral center."
    return f


def check_tautomer_ambiguity(records: list[MoleculeRecord], ctx: dict[str, Any]) -> Finding | None:
    """Flag valid molecules with multiple enumerable tautomers.

    Args:
        records: Per-molecule records.
        ctx: Check context with `max_examples`.

    Returns:
        Informational finding for tautomer ambiguity, or None if none found.
    """
    bad = [r for r in records if r.valid and (r.n_tautomers or 1) > 1]
    rows = [r.row_id for r in bad]
    ex = [{"row": r.row_id, "smiles": r.canon_smi, "n_tautomers": r.n_tautomers,
           "canonical_tautomer": r.tautomer_smi} for r in bad][:ctx["max_examples"]]
    return _mk("tautomer_ambiguity", Severity.INFO, "tautomer ambiguity", rows, ex,
               "Canonicalize tautomers (RDKit TautomerEnumerator) before dedup/split so keto/enol forms match.",
               len(records))


CHECKS = [check_invalid_smiles, check_valence, check_aromaticity,
          check_impossible_charge, check_disconnected, check_isotopes,
          check_radicals, check_unspecified_stereo, check_tautomer_ambiguity]
