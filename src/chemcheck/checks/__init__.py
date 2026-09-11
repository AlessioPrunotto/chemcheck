"""Check registry: each check returns a Finding or None."""
from __future__ import annotations

from typing import TYPE_CHECKING

from . import duplicates, integrity, ml, space, splits

if TYPE_CHECKING:
    import pandas as pd

    from ..models import Finding
    from ..molecules import MoleculeRecord

ALL_MODULES = [integrity, duplicates, splits, ml, space]


def run_all(records: list[MoleculeRecord], df: pd.DataFrame, max_examples: int = 20,
            near_dup_thresh: float = 0.95, analog_thresh: float = 0.6) -> list[Finding]:
    """Run every registered check over records.

    Checks never crash an audit: exceptions are caught and reported
    as informational findings.

    Args:
        records: Per-molecule records built by `build_records`.
        df: Normalized input table (passed through in check context).
        max_examples: Maximum examples stored per finding.
        near_dup_thresh: Tanimoto threshold for near-duplicates.
        analog_thresh: Tanimoto threshold for analog leakage.

    Returns:
        List of findings (checks returning None are skipped).
    """
    findings = []
    ctx = {"df": df, "max_examples": max_examples,
           "near_dup_thresh": near_dup_thresh, "analog_thresh": analog_thresh}
    for mod in ALL_MODULES:
        for fn in getattr(mod, "CHECKS", []):
            try:
                f = fn(records, ctx)
            except Exception as e:  # checks must never crash an audit
                from ..models import Finding, Severity
                f = Finding(check_id=getattr(fn, "__name__", "check"),
                            severity=Severity.INFO,
                            title=f"check crashed: {e}", count=0,
                            details="Internal check error; please report.")
            if f is not None:
                findings.append(f)
    return findings
