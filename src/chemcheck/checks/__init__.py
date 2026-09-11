"""Check registry: each check returns a Finding or None."""
from __future__ import annotations

from . import duplicates, integrity, ml, space, splits

ALL_MODULES = [integrity, duplicates, splits, ml, space]


def run_all(records, df, max_examples: int = 20, near_dup_thresh: float = 0.95,
            analog_thresh: float = 0.6):
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
