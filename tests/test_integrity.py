"""Integrity checks."""
from chemcheck.checks import integrity
from chemcheck.io import load_table
from chemcheck.molecules import build_records


def _recs(tmp_path, rows="id,smiles\n"):
    p = tmp_path / "t.csv"
    p.write_text(rows)
    df = load_table([str(p)])
    return build_records(df), {"max_examples": 5}


def test_invalid_smiles(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CCO\nb,not_a_smiles!!\n")
    f = integrity.check_invalid_smiles(recs, ctx)
    assert f is not None and f.count == 1 and "b" in f.affected_rows


def test_valence_error(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CCO\nb,C(C)(C)(C)(C)C\n")
    f = integrity.check_valence(recs, ctx)
    assert f is not None and f.count >= 1


def test_impossible_charge_and_disconnected(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,[C-2]\nb,CCO.[Na+]\n")
    fc = integrity.check_impossible_charge(recs, ctx)
    dc = integrity.check_disconnected(recs, ctx)
    assert fc is not None and "a" in fc.affected_rows
    assert dc is not None and "b" in dc.affected_rows


def test_unspecified_stereo_and_isotope(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CC(Cl)Br\nb,[13C]C\n")
    fs = integrity.check_unspecified_stereo(recs, ctx)
    fi = integrity.check_isotopes(recs, ctx)
    assert fs is not None and "a" in fs.affected_rows
    assert fi is not None and "b" in fi.affected_rows
