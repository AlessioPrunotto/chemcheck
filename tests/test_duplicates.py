"""Duplicate checks."""
from chemdatacheck.checks import duplicates
from chemdatacheck.io import load_table
from chemdatacheck.molecules import build_records


def _recs(tmp_path, rows):
    p = tmp_path / "t.csv"
    p.write_text(rows)
    df = load_table([str(p)])
    return build_records(df), {"max_examples": 5, "near_dup_thresh": 0.5}


def test_exact_and_canonical(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CCO\nb,CCO\nc,not_a_smiles!!\n")
    # a/b share raw + canonical; c is invalid and excluded from canonical dups
    fe = duplicates.check_exact_duplicates(recs, ctx)
    fc = duplicates.check_canonical_duplicates(recs, ctx)
    assert fe is not None and fe.count == 2
    assert fc is not None and fc.count == 2


def test_stereo_collision(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,C[C@H](Cl)Br\nb,C[C@@H](Cl)Br\n")
    f = duplicates.check_stereo_collisions(recs, ctx)
    assert f is not None and f.count == 2


def test_salt_and_tautomer(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CCO\nb,CCO.[Na+]\nc,CC(=O)C\nd,CC(O)=C\n")
    fs = duplicates.check_salt_duplicates(recs, ctx)
    ft = duplicates.check_tautomer_duplicates(recs, ctx)
    assert fs is not None and set(fs.affected_rows) >= {"a", "b"}
    assert ft is not None  # keto/enol collapse


def test_near_duplicates(tmp_path):
    recs, ctx = _recs(tmp_path, "id,smiles\na,CCCCCCCCO\nb,CCCCCCCCCO\n")
    f = duplicates.check_near_duplicates(recs, ctx)
    assert f is not None and f.count == 2
    assert f.metadata["approximate"] is False
