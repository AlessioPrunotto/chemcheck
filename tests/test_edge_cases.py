"""Boundary behavior that should remain stable across releases."""
from __future__ import annotations

import argparse

import pandas as pd
import pytest

from chemdatacheck.checks import integrity, run_all
from chemdatacheck.cli import _non_negative_float, _positive_int, _unit_interval, main
from chemdatacheck.io import _guess_smiles_column, load_table
from chemdatacheck.models import AuditReport, Finding, Severity
from chemdatacheck.report import audit


@pytest.mark.parametrize("dtype", ["object", "string"])
def test_smiles_heuristic_accepts_string_dtypes(dtype):
    frame = pd.DataFrame({
        "measurement": [1.0, 2.0],
        "fold": pd.Series(["train", "test"], dtype=dtype),
        "compound": pd.Series([" CCO ", "CCN"], dtype=dtype),
    })
    assert _guess_smiles_column(frame) == "compound"


def test_smiles_heuristic_and_generated_ids(tmp_path):
    path = tmp_path / "unusual.csv"
    path.write_text("compound,fold\n CCO ,train\nCCN,test\n")
    frame = load_table([str(path)])
    assert frame.attrs["smiles_col"] == "compound"
    assert frame.attrs["split_col"] == "fold"
    assert frame["_smiles"].tolist() == ["CCO", "CCN"]
    assert frame["_id"].tolist() == ["row_0", "row_1"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"smiles_col": "missing"}, "Could not detect SMILES"),
        ({"id_col": "missing"}, "--id-col"),
        ({"label_col": "missing"}, "--label-col"),
        ({"split_col": "missing"}, "--split-col"),
    ],
)
def test_explicit_missing_columns_are_rejected(tmp_path, kwargs, message):
    path = tmp_path / "molecules.csv"
    path.write_text("smiles,id\nCCO,a\n")
    with pytest.raises(ValueError, match=message):
        load_table([str(path)], **kwargs)


def test_unsupported_format_and_input_count_are_rejected(tmp_path):
    path = tmp_path / "molecules.xyz"
    path.write_text("smiles\nCCO\n")
    with pytest.raises(ValueError, match="Unsupported input format"):
        load_table([str(path)])
    with pytest.raises(ValueError, match="Pass 1 dataset file"):
        load_table([])
    with pytest.raises(ValueError, match="Pass 1 dataset file"):
        load_table(["a.csv", "b.csv", "c.csv"])


def test_header_only_dataset_is_a_valid_empty_audit(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("id,smiles\n")
    report = audit([str(path)])
    assert (report.n_total, report.n_valid, report.n_invalid) == (0, 0, 0)
    assert report.score == 100
    assert {finding.check_id for finding in report.findings} == {"split_info"}


@pytest.mark.parametrize(
    ("parser", "value"),
    [(_unit_interval, "nan"), (_unit_interval, "1.01"),
     (_positive_int, "0"), (_non_negative_float, "-0.1")],
)
def test_cli_numeric_parsers_reject_invalid_bounds(parser, value):
    with pytest.raises(argparse.ArgumentTypeError):
        parser(value)


def test_cli_error_and_file_report_formats(tmp_path, capsys):
    missing = tmp_path / "missing.csv"
    assert main([str(missing), "-q"]) == 2
    assert "chemdatacheck error" in capsys.readouterr().err

    data = tmp_path / "molecules.csv"
    data.write_text("id,smiles\na,CCO\nb,CCN\n")
    for fmt, marker in (("html", "<!doctype html>"), ("junit", "<testsuite")):
        output = tmp_path / f"report.{fmt}"
        assert main([str(data), "--format", fmt, "--output", str(output),
                     "--fail-on", "never", "-q"]) == 0
        assert marker in output.read_text()


def test_categorical_and_numeric_target_shift_paths(tmp_path):
    categorical = tmp_path / "categorical.csv"
    categorical_rows = [f"t{i},C{'C' * (i + 1)},train,A" for i in range(10)]
    categorical_rows += [f"v{i},N{'C' * (i + 1)},test,B" for i in range(5)]
    categorical.write_text("id,smiles,split,label\n" + "\n".join(categorical_rows) + "\n")
    cat_report = audit([str(categorical)], label_col="label")
    cat_shift = next(f for f in cat_report.findings if f.check_id == "target_distribution_shift")
    assert "prevalence shift" in cat_shift.title

    numeric = tmp_path / "numeric.csv"
    numeric_rows = [f"t{i},C{'C' * (i + 1)},train,{i / 100}" for i in range(10)]
    numeric_rows += [f"v{i},N{'C' * (i + 1)},test,{10 + i / 100}" for i in range(5)]
    numeric.write_text("id,smiles,split,label\n" + "\n".join(numeric_rows) + "\n")
    num_report = audit([str(numeric)], label_col="label")
    ids = {finding.check_id for finding in num_report.findings}
    assert {"target_distribution_shift", "target_leakage_split_predicts_label"} <= ids


def test_report_filters_and_exit_never():
    warning = Finding("warning", Severity.WARNING, "warning", 1)
    error = Finding("error", Severity.ERROR, "error", 1)
    report = AuditReport(1, 1, 0, [warning, error], 50, [])
    assert report.errors() == [error]
    assert report.warnings() == [warning]


def test_check_failure_is_isolated(monkeypatch):
    def explode(records, ctx):
        raise RuntimeError("deliberate failure")

    monkeypatch.setattr(integrity, "CHECKS", [explode])
    findings = run_all([], pd.DataFrame())
    assert findings[0].check_id == "explode"
    assert "check crashed" in findings[0].title
