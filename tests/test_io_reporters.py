"""Input-format and report-renderer coverage."""
import json
from xml.etree import ElementTree

import pandas as pd
import pytest

from chemcheck.io import load_table
from chemcheck.models import AuditReport, Finding, Severity
from chemcheck.report import render_html, render_json, render_junit


@pytest.mark.parametrize("extension", ["tsv", "txt", "jsonl", "json", "parquet", "xlsx"])
def test_tabular_input_formats(tmp_path, extension):
    raw = pd.DataFrame({"id": ["a", "b"], "smiles": ["CCO", "CCN"]})
    path = tmp_path / f"molecules.{extension}"
    if extension in ("tsv", "txt"):
        raw.to_csv(path, sep="\t", index=False)
    elif extension == "jsonl":
        raw.to_json(path, orient="records", lines=True)
    elif extension == "json":
        raw.to_json(path, orient="records")
    elif extension == "parquet":
        raw.to_parquet(path, index=False)
    else:
        raw.to_excel(path, index=False)
    loaded = load_table([str(path)])
    assert loaded["_id"].tolist() == ["a", "b"]
    assert loaded["_smiles"].tolist() == ["CCO", "CCN"]


def test_sdf_input(tmp_path):
    from rdkit import Chem
    path = tmp_path / "molecules.sdf"
    writer = Chem.SDWriter(str(path))
    for name, smiles in (("a", "CCO"), ("b", "CCN")):
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", name)
        writer.write(mol)
    writer.close()
    loaded = load_table([str(path)])
    assert loaded["_id"].tolist() == ["a", "b"]
    assert len(loaded) == 2


def test_reporters_escape_content_and_serialize_metadata():
    finding = Finding(check_id="unsafe", severity=Severity.ERROR,
                      title='unsafe <tag> & "quote"', count=1,
                      metadata={"approximate": True, "method": "sample"})
    report = AuditReport(n_total=1, n_valid=1, n_invalid=0, findings=[finding],
                         score=90, score_breakdown=[], meta={"source": "x<&"})
    ElementTree.fromstring(render_junit(report))
    assert "&lt;tag&gt;" in render_html(report)
    data = json.loads(render_json(report))
    assert data["findings"][0]["metadata"]["approximate"] is True
