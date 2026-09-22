"""Split/leakage, ML, space, scoring, CLI."""
import json

import pytest

from chemcheck.checks import splits
from chemcheck.cli import main as cli_main
from chemcheck.io import load_table
from chemcheck.models import Finding, Severity
from chemcheck.molecules import build_records
from chemcheck.report import audit
from chemcheck.scoring import score_findings


def _audit(tmp_path, text, **kw):
    p = tmp_path / "d.csv"
    p.write_text(text)
    return audit([str(p)], **kw)


def test_identity_and_scaffold_leakage(tmp_path):
    text = ("id,smiles,split\n"
            "a1,c1ccccc1,train\na2,CCO,train\n"
            "b1,c1ccccc1,test\nb2,CCOCC,test\n")
    rep = _audit(tmp_path, text)
    ids = {f.check_id for f in rep.findings}
    assert "identity_leakage" in ids  # benzene spans splits


def test_analog_leakage_reports_pair(tmp_path):
    text = ("id,smiles,split\n"
            "a1,CCCCCCCCO,train\n"
            "b1,CCCCCCCCCO,test\n")
    rep = _audit(tmp_path, text, analog_thresh=0.5)
    leak = next((f for f in rep.findings if f.check_id == "analog_leakage"), None)
    assert leak is not None and leak.examples
    assert {"train_row", "test_row", "tanimoto"} <= set(leak.examples[0])


def test_cross_split_similarity_profile_is_cached(tmp_path, monkeypatch):
    path = tmp_path / "split.csv"
    path.write_text("id,smiles,split\na,CCO,train\nb,CCN,train\nc,CCC,test\nd,CCCl,test\n")
    records = build_records(load_table([str(path)]))
    ctx = {"train_values": None, "test_values": None, "approximations": []}
    from rdkit import DataStructs
    original = DataStructs.BulkTanimotoSimilarity
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(DataStructs, "BulkTanimotoSimilarity", counted)
    first = splits._cross_similarity_profile(records, ctx, max_test=1)
    second = splits._cross_similarity_profile(records, ctx, max_test=1)
    assert first is second
    assert calls == 1
    assert ctx["approximations"][0]["method"] == "test_sampling"


def test_custom_fold_mapping_enables_leakage_checks(tmp_path):
    text = ("id,smiles,fold\n"
            "a,CCO,0\n"
            "b,CCO,1\n")
    rep = _audit(tmp_path, text, test_values=["1"])
    ids = {f.check_id for f in rep.findings}
    assert "identity_leakage" in ids
    assert "split_configuration" not in ids
    assert rep.meta["settings"]["resolved_train_values"] == ["0"]
    assert rep.meta["settings"]["resolved_test_values"] == ["1"]


def test_unrecognized_split_values_are_not_silent(tmp_path):
    text = ("id,smiles,split\n"
            "a,CCO,calibration\n"
            "b,CCO,production\n")
    rep = _audit(tmp_path, text)
    finding = next(f for f in rep.findings if f.check_id == "split_configuration")
    assert finding.severity == Severity.WARNING
    assert finding.examples[0]["observed_values"] == ["calibration", "production"]


def test_no_split_is_reported_as_info(tmp_path):
    rep = _audit(tmp_path, "id,smiles\na,CCO\n")
    finding = next(f for f in rep.findings if f.check_id == "split_info")
    assert finding.severity == Severity.INFO


def test_conflicting_and_duplicated_measurements(tmp_path):
    text = ("id,smiles,activity\n"
            "a,CCO,1.0\nb,CCO,1.0\nc,c1ccccc1O,2.0\nd,c1ccccc1O,8.0\n")
    rep = _audit(tmp_path, text, label_col="activity")
    ids = {f.check_id for f in rep.findings}
    assert "duplicated_measurements" in ids
    assert "conflicting_measurements" in ids


def test_numeric_conflict_threshold_is_configurable(tmp_path):
    text = "id,smiles,activity\na,CCO,1.0\nb,CCO,1.5\n"
    default = _audit(tmp_path, text, label_col="activity")
    strict = _audit(tmp_path, text, label_col="activity", conflict_thresh=0.25)
    assert "conflicting_measurements" not in {f.check_id for f in default.findings}
    assert "conflicting_measurements" in {f.check_id for f in strict.findings}
    assert strict.meta["settings"]["conflict_thresh"] == 0.25


def test_binary_labels_are_categorical_not_numeric_outliers(tmp_path):
    rows = [f"r{i},CC{'C' * i},{1 if i == 20 else 0}" for i in range(1, 21)]
    rep = _audit(tmp_path, "id,smiles,activity\n" + "\n".join(rows) + "\n",
                 label_col="activity")
    assert "label_outliers" not in {f.check_id for f in rep.findings}


def test_conflicting_binary_labels_are_detected(tmp_path):
    rep = _audit(tmp_path, "id,smiles,activity\na,CCO,0\nb,CCO,1\n",
                 label_col="activity")
    assert "conflicting_measurements" in {f.check_id for f in rep.findings}


def test_space_and_score(tmp_path):
    p = tmp_path / "x.csv"
    import shutil
    shutil.copy("tests/fixtures/demo.csv", p)
    rep = audit([str(p)], label_col="activity", split_col="split")
    assert 0 <= rep.score <= 100
    assert rep.n_total == 20
    assert any(f.severity.value == "error" for f in rep.findings)
    # breakdown sums roughly to 100 - score
    ded = sum(d["deduction"] for d in rep.score_breakdown)
    assert abs((100 - rep.score) - ded) < 5.0


def test_score_caps_overlapping_duplicate_deductions():
    findings = [
        Finding(check_id=check_id, severity=Severity.WARNING, title=check_id,
                count=100, rate=1.0)
        for check_id in ("exact_duplicates", "canonical_duplicates", "near_duplicates")
    ]
    _, breakdown = score_findings(findings, 100)
    assert sum(item["deduction"] for item in breakdown) <= 6.01
    assert all(item["overlap_adjusted"] for item in breakdown)


def test_cli_exit_codes_and_json(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("id,smiles\na,CCO\nb,c1ccccc1\n")
    assert cli_main([str(p), "--fail-on", "error", "-q"]) == 0
    bad = tmp_path / "bad.csv"
    bad.write_text("id,smiles\na,not_a_smiles!!\n")
    assert cli_main([str(bad), "--fail-on", "error", "-q"]) == 2
    out = tmp_path / "rep.json"
    assert cli_main([str(p), "--format", "json", "--output", str(out), "-q"]) == 0
    data = json.loads(out.read_text())
    assert data["summary"]["n_total"] == 2 and "score" in data["summary"]
    assert data["meta"]["chemcheck_version"] == "0.1.0"
    assert data["meta"]["rdkit_version"]
    assert len(next(iter(data["meta"]["input_sha256"].values()))) == 64
    assert data["meta"]["settings"]["near_dup_thresh"] == 0.95


def test_audit_rejects_invalid_settings(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("id,smiles\na,CCO\n")
    with pytest.raises(ValueError, match="between 0 and 1"):
        audit([str(p)], analog_thresh=1.1)
    with pytest.raises(ValueError, match="both train and test"):
        audit([str(p)], train_values=["same"], test_values=["SAME"])
    with pytest.raises(ValueError, match="non-negative"):
        audit([str(p)], conflict_thresh=-0.1)
