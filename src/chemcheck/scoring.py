"""Dataset quality scoring: transparent weighted deductions → 0-100."""
from __future__ import annotations

import math
from typing import Any

from .models import Finding

WEIGHTS = {
    # chemical integrity
    "invalid_smiles": 8.0,
    "valence_error": 7.0,
    "aromaticity_suspect": 3.0,
    "impossible_charge": 3.0,
    "disconnected_components": 2.0,
    "isotope_flag": 0.5,
    "radical_flag": 4.0,
    "unspecified_stereo": 2.5,
    "tautomer_ambiguity": 1.0,
    # duplicates
    "exact_duplicates": 3.0,
    "canonical_duplicates": 4.0,
    "stereochemical_collisions": 4.0,
    "salt_duplicates": 3.0,
    "tautomer_duplicates": 1.5,
    "near_duplicates": 4.0,
    # splits / leakage
    "scaffold_overlap": 5.0,
    "identity_leakage": 10.0,
    "analog_leakage": 7.0,
    "cross_split_near_duplicates": 9.0,
    "suspiciously_easy_split": 4.0,
    # ML
    "duplicated_measurements": 3.0,
    "conflicting_measurements": 8.0,
    "target_distribution_shift": 4.0,
    "label_outliers": 1.0,
    "target_leakage_split_predicts_label": 6.0,
    # chemical space
    "rare_elements": 1.0,
    "rare_functional_groups": 1.0,
    "unusual_ring_systems": 1.5,
    "representation_bias": 3.0,
    "applicability_gaps": 1.5,
}


def _factor(rate: float, count: int) -> float:
    """Map prevalence rate to a [0, 1] deduction factor.

    Saturating curve: even small counts deduct something, large rates cap out.
    factor = 1 - exp(-(rate*12 + (0.15 if count>0 else 0)))

    Args:
        rate: Fraction of rows affected.
        count: Number of affected rows.

    Returns:
        Deduction factor between 0.0 and 1.0.
    """
    if count <= 0:
        return 0.0
    return 1.0 - math.exp(-(rate * 12.0 + 0.15))


def score_findings(findings: list[Finding], n_total: int) -> tuple[int, list[dict[str, Any]]]:
    """Score findings with transparent weighted deductions.

    Args:
        findings: Findings to score (`split_info` is skipped).
        n_total: Total number of records (unused beyond signature symmetry).

    Returns:
        Tuple of (score from 0 to 100, per-check deduction breakdown
        sorted by deduction descending).
    """
    breakdown = []
    total_deduction = 0.0
    for f in findings:
        if f.check_id == "split_info":
            continue
        w = WEIGHTS.get(f.check_id, 2.0)
        # target shift / leakage flags are dataset-level: use fixed partial factor
        if f.check_id in ("target_distribution_shift", "target_leakage_split_predicts_label",
                          "suspiciously_easy_split") and f.rate == 0.0 and f.count > 0:
            fac = 0.6
        else:
            fac = _factor(f.rate, f.count)
        ded = round(w * fac, 2)
        total_deduction += ded
        breakdown.append({"check_id": f.check_id, "title": f.title,
                          "severity": f.severity.value, "count": f.count,
                          "weight": w, "deduction": ded})
    breakdown.sort(key=lambda d: -d["deduction"])
    score = max(0, min(100, round(100 - total_deduction)))
    return score, breakdown
