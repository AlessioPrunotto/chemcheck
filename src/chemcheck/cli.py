"""CLI: chemcheck dataset.csv [test.csv] [options]. Exit codes pytest-like."""
from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from typing import TYPE_CHECKING

from . import __version__
from .report import audit, render_html, render_json, render_junit, render_terminal

if TYPE_CHECKING:
    from .models import AuditReport


def _unit_interval(value: str) -> float:
    """Parse a floating-point CLI value constrained to [0, 1]."""
    parsed = float(value)
    if not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return parsed


def _positive_int(value: str) -> int:
    """Parse a strictly positive integer CLI value."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _non_negative_float(value: str) -> float:
    """Parse a non-negative floating-point CLI value."""
    parsed = float(value)
    if parsed < 0.0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    """Build the chemcheck argument parser.

    Returns:
        Configured argument parser for the CLI.
    """
    p = argparse.ArgumentParser(
        prog="chemcheck",
        description="pytest for molecular datasets — sanity checks for chemistry ML datasets.")
    p.add_argument("inputs", nargs="+",
                   help="one dataset file (CSV/TSV/SDF/Parquet/Excel/JSONL) or two files as train/test")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--smiles-col", default=None)
    p.add_argument("--id-col", default=None)
    p.add_argument("--label-col", default=None)
    p.add_argument("--split-col", default=None)
    p.add_argument("--train-value", action="append", dest="train_values", default=None,
                   help="split value treated as train; repeat for multiple values")
    p.add_argument("--test-value", action="append", dest="test_values", default=None,
                   help="split value treated as test; repeat for multiple values (one fold is sufficient)")
    p.add_argument("--format", dest="fmt", default="terminal",
                   choices=["terminal", "json", "html", "junit"])
    p.add_argument("--output", "-o", default=None, help="write report to file instead of stdout")
    p.add_argument("--fail-on", default="warning",
                   choices=["error", "warning", "info", "never"],
                   help="minimum severity that triggers a nonzero exit code")
    p.add_argument("--max-examples", type=_positive_int, default=20)
    p.add_argument("--near-dup-thresh", type=_unit_interval, default=0.95)
    p.add_argument("--analog-thresh", type=_unit_interval, default=0.6)
    p.add_argument("--conflict-thresh", type=_non_negative_float, default=1.0,
                   help="numeric label spread considered conflicting (default: 1.0)")
    p.add_argument("--quiet", "-q", action="store_true")
    return p


def exit_code_for(report: AuditReport, fail_on: str) -> int:
    """Map report findings to a pytest-like exit code.

    Args:
        report: Audit report to evaluate.
        fail_on: Minimum severity triggering nonzero exit
            ("error", "warning", "info", or "never").

    Returns:
        0 for pass, 1 for warnings, 2 for errors.
    """
    levels = {"never": 99, "info": 0, "warning": 1, "error": 2}
    threshold = levels[fail_on]
    if threshold == 99:
        return 0
    worst = -1
    for f in report.findings:
        lv = {"info": 0, "warning": 1, "error": 2}[f.severity.value]
        worst = max(worst, lv)
    if worst < threshold:
        return 0
    return 2 if worst == 2 else 1


def main(argv: Sequence[str] | None = None) -> int:
    """Run the chemcheck CLI.

    Args:
        argv: Argument list excluding the program name, or None to use
            `sys.argv`.

    Returns:
        Process exit code (0 pass, 1 warnings, 2 errors).
    """
    args = build_parser().parse_args(argv)
    try:
        report = audit(args.inputs, smiles_col=args.smiles_col, id_col=args.id_col,
                       label_col=args.label_col, split_col=args.split_col,
                       max_examples=args.max_examples,
                       near_dup_thresh=args.near_dup_thresh,
                       analog_thresh=args.analog_thresh,
                       train_values=args.train_values,
                       test_values=args.test_values,
                       conflict_thresh=args.conflict_thresh)
    except ImportError as e:
        print(f"chemcheck error: {e}", file=sys.stderr)
        return 2
    except (FileNotFoundError, ValueError) as e:
        print(f"chemcheck error: {e}", file=sys.stderr)
        return 2
    if args.fmt == "json":
        out = render_json(report)
    elif args.fmt == "html":
        out = render_html(report)
    elif args.fmt == "junit":
        out = render_junit(report)
    else:
        # terminal: try rich coloring if available, else plain
        out = render_terminal(report)
        try:
            from rich.console import Console  # optional
            Console().print(out)
            out = None
        except Exception:
            pass
    if out is not None:
        if args.output:
            with open(args.output, "w") as fh:
                fh.write(out)
        else:
            print(out)
    elif args.output:
        # rich path + output file: write plain too
        with open(args.output, "w") as fh:
            fh.write(render_terminal(report))
    code = exit_code_for(report, args.fail_on)
    if not args.quiet and args.fmt == "terminal":
        sev = "clean" if code == 0 else ("warnings" if code == 1 else "errors")
        print(f"(exit {code}: {sev}, --fail-on {args.fail_on})")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
