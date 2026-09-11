"""Split/leakage, ML, space, scoring, CLI."""
import json

from chemcheck.io import load_table
from chemcheck.molecules import build_records
from chemcheck.checks import splits, ml, space
from chemcheck.report import audit
from chemcheck.cli import main as cli_main


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


def test_conflicting_and_duplicated_measurements(tmp_path):
    text = ("id,smiles,activity\n"
            "a,CCO,1.0\nb,CCO,1.0\nc,c1ccccc1O,2.0\nd,c1ccccc1O,8.0\n")
    rep = _audit(tmp_path, text, label_col="activity")
    ids = {f.check_id for f in rep.findings}
    assert "duplicated_measurements" in ids
    assert "conflicting_measurements" in ids


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
