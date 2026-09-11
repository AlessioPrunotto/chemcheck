"""Audit orchestration + reporters (terminal / JSON / HTML / JUnit)."""
from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from typing import Any

from .checks import run_all
from .io import load_table
from .models import AuditReport, Severity
from .molecules import build_records
from .scoring import score_findings


def audit(paths: list[str], smiles_col: str | None = None, id_col: str | None = None,
          label_col: str | None = None, split_col: str | None = None,
          max_examples: int = 20, near_dup_thresh: float = 0.95,
          analog_thresh: float = 0.6) -> AuditReport:
    """Audit a molecular dataset and return a quality report.

    Loads input file(s), builds molecule records, runs all checks,
    scores findings, and assembles report metadata.

    Args:
        paths: One dataset path, or two paths treated as train/test.
        smiles_col: SMILES column name, or None to auto-detect.
        id_col: ID column name, or None to auto-detect.
        label_col: Label column name, or None for unlabeled data.
        split_col: Split column name, or None to auto-detect.
        max_examples: Maximum examples stored per finding.
        near_dup_thresh: Tanimoto threshold for near-duplicates.
        analog_thresh: Tanimoto threshold for analog leakage.

    Returns:
        Populated audit report with findings ordered by severity.
    """
    df = load_table(paths, smiles_col, id_col, label_col, split_col)
    records = build_records(df)
    n_total = len(records)
    n_valid = sum(1 for r in records if r.valid)
    findings = run_all(records, df, max_examples, near_dup_thresh, analog_thresh)
    # order: errors, warnings, infos; within group by count desc
    order = {Severity.ERROR: 0, Severity.WARNING: 1, Severity.INFO: 2}
    findings.sort(key=lambda f: (order[f.severity], -f.count))
    score, breakdown = score_findings(findings, n_total)
    n_unassigned = sum(r.n_stereo_unassigned for r in records if r.valid)
    n_centers = sum(r.n_stereo_defined + r.n_stereo_unassigned for r in records if r.valid)
    meta = {
        "source": df.attrs.get("source"),
        "smiles_col": df.attrs.get("smiles_col"),
        "split_col": df.attrs.get("split_col"),
        "label_col": df.attrs.get("label_col") or label_col,
        "n_total": n_total, "n_valid": n_valid, "n_invalid": n_total - n_valid,
        "frac_unspecified_stereo": (n_unassigned / n_centers) if n_centers else 0.0,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    return AuditReport(n_total=n_total, n_valid=n_valid, n_invalid=n_total - n_valid,
                       findings=findings, score=score, score_breakdown=breakdown, meta=meta)


# ---------- reporters ----------

def _icon(sev: Severity) -> str:
    """Return a glyph for a severity level.

    Args:
        sev: Severity to render.

    Returns:
        Unicode icon string.
    """
    return {"error": "✕", "warning": "⚠", "info": "ℹ"}[sev.value]


def to_dict(report: AuditReport) -> dict[str, Any]:
    """Convert a report to a JSON-serializable dictionary.

    Args:
        report: Audit report to serialize.

    Returns:
        Dictionary with meta, summary, score breakdown, and findings.
    """
    return {
        "meta": report.meta,
        "summary": {"n_total": report.n_total, "n_valid": report.n_valid,
                    "n_invalid": report.n_invalid, "score": report.score},
        "score_breakdown": report.score_breakdown,
        "findings": [
            {"check_id": f.check_id, "severity": f.severity.value, "title": f.title,
             "count": f.count, "rate": round(f.rate, 4), "total_affected": f.total_affected,
             "affected_rows": f.affected_rows[:200], "examples": f.examples[:20],
             "recommendation": f.recommendation, "details": f.details}
            for f in report.findings
        ],
    }


def render_terminal(report: AuditReport) -> str:
    """Render a report as plain-text terminal output.

    Args:
        report: Audit report to render.

    Returns:
        Multi-line human-readable report string.
    """
    L = []
    L.append("")
    L.append(f"✓ {report.n_valid:,} valid molecules / ✕ {report.n_invalid:,} invalid "
             f"(n={report.n_total:,})")
    if not report.findings:
        L.append("✓ no issues found")
    for f in report.findings:
        if f.check_id == "split_info" and f.count == 0 and not f.examples:
            L.append(f"ℹ {f.title} — {f.recommendation}")
            continue
        pct = f" ({f.rate * 100:.1f}%)" if f.rate else ""
        L.append(f"{_icon(f.severity)} [{f.severity.value}] {f.title}: "
                 f"{f.count:,}{pct}  [{f.check_id}]")
        for ex in f.examples[:5]:
            L.append(f"    e.g. {json.dumps(ex)[:220]}")
        if f.recommendation:
            L.append(f"    → {f.recommendation[:240]}")
        if f.details:
            L.append(f"    ({f.details[:200]})")
    L.append("")
    L.append(f"Dataset quality score: {report.score}/100")
    if report.score_breakdown[:5]:
        top = ", ".join(f"{d['check_id']} -{d['deduction']}" for d in report.score_breakdown[:5])
        L.append(f"  top deductions: {top}")
    L.append("")
    return "\n".join(L)


def render_json(report: AuditReport) -> str:
    """Render a report as indented JSON.

    Args:
        report: Audit report to render.

    Returns:
        JSON string of `to_dict(report)`.
    """
    return json.dumps(to_dict(report), indent=2)


def render_html(report: AuditReport) -> str:
    """Render a report as a standalone HTML page.

    Args:
        report: Audit report to render.

    Returns:
        HTML document string.
    """
    rows = "\n".join(
        f"<tr><td>{_icon(f.severity)} {f.severity.value}</td><td>{html.escape(f.title)}</td>"
        f"<td>{f.count}</td><td>{html.escape(f.recommendation[:300])}</td>"
        f"<td><pre>{html.escape(json.dumps(f.examples[:3], indent=1)[:1500])}</pre></td></tr>"
        for f in report.findings)
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>chemcheck report — score {report.score}/100</title>
<style>body{{font-family:system-ui,sans-serif;margin:2em}}table{{border-collapse:collapse;width:100%}}
td,th{{border:1px solid #ccc;padding:6px;vertical-align:top}}pre{{white-space:pre-wrap}}</style>
</head><body><h1>chemcheck report — score {report.score}/100</h1>
<p>{report.n_valid}/{report.n_total} valid molecules. Source: {html.escape(str(report.meta.get('source')))}</p>
<table><tr><th>severity</th><th>finding</th><th>n</th><th>recommendation</th><th>examples</th></tr>
{rows}</table></body></html>"""


def render_junit(report: AuditReport) -> str:
    """Render a report as JUnit XML for CI integration.

    Args:
        report: Audit report to render.

    Returns:
        JUnit XML string with one testcase per finding.
    """
    cases = []
    for f in report.findings:
        status = "failed" if f.severity == Severity.ERROR else "warning"
        body = html.escape(f"{f.title} (n={f.count}). {f.recommendation} {f.details}"[:2000])
        if f.severity == Severity.ERROR:
            cases.append(f'  <testcase classname="chemcheck" name="{f.check_id}">'
                         f"<failure message=\"{html.escape(f.title)}\">{body}</failure></testcase>")
        else:
            cases.append(f'  <testcase classname="chemcheck" name="{f.check_id}">'
                         f"<system-out>{status}: {body}</system-out></testcase>")
    if not cases:
        cases.append('  <testcase classname="chemcheck" name="all_checks_passed"/>')
    n_fail = sum(1 for f in report.findings if f.severity == Severity.ERROR)
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n<testsuite name="chemcheck" '
            f'tests="{len(cases)}" failures="{n_fail}">\n' + "\n".join(cases) + "\n</testsuite>\n")
