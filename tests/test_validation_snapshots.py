"""Integrity checks for committed external-validation sample definitions."""
from pathlib import Path


def test_database_sample_identifier_snapshots_are_fixed_and_unique():
    sample_dir = Path(__file__).parents[1] / "validation" / "samples"
    for filename in ("chembl37_ids.txt", "pubchem_cids.txt"):
        identifiers = (sample_dir / filename).read_text().splitlines()
        assert len(identifiers) == 5_000
        assert len(set(identifiers)) == 5_000
