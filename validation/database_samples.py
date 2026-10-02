"""Acquire fixed, structure-only samples from the official ChEMBL and PubChem APIs."""
from __future__ import annotations

import csv
import hashlib
import json
import random
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

SAMPLE_SIZE = 5_000
CHEMBL_VERSION = "ChEMBL_37"
CHEMBL_RELEASE_DATE = "2026-05-01"
CHEMBL_STRUCTURE_COUNT = 2_897_819
CHEMBL_PAGE_SIZE = 500
CHEMBL_OFFSETS = tuple(
    round(index * (CHEMBL_STRUCTURE_COUNT - CHEMBL_PAGE_SIZE) / 9)
    for index in range(10)
)
PUBCHEM_SEED = 20_260_924
PUBCHEM_LIVE_COMPOUND_COUNT = 124_663_543
PUBCHEM_CID_RANGE_MAX = PUBCHEM_LIVE_COMPOUND_COUNT
PUBCHEM_POPULATION_DATE = "2026-09-06"
USER_AGENT = "chemdatacheck-scientific-validation/0.1 (+https://github.com/AlessioPrunotto/chemdatacheck)"


@dataclass(frozen=True)
class DatabaseSample:
    name: str
    filename: str
    ids_filename: str
    source_url: str
    expected_sha256: str | None = None


DATABASE_SAMPLES = {
    "chembl37_sample": DatabaseSample(
        "chembl37_sample",
        "chembl37_sample.csv",
        "chembl37_ids.txt",
        "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
        "50790ba94e76376528efb48d60944679b202e5e693ce7ca66c5e50926513f081",
    ),
    "pubchem_sample": DatabaseSample(
        "pubchem_sample",
        "pubchem_sample.csv",
        "pubchem_cids.txt",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/SMILES,InChIKey/JSON",
        "7c3090c52692447825ce449a7f554a6337d3f990945b381b496f1d167efcf718",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_request(url: str, data: bytes | None = None, attempts: int = 6) -> dict[str, Any]:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    for attempt in range(attempts):
        try:
            with urlopen(Request(url, data=data, headers=headers), timeout=120) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError):
            if attempt == attempts - 1:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _write_ids(path: Path, identifiers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(identifiers) + "\n")


def _read_ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def initialize_identifier_lists(sample_dir: Path) -> None:
    """Create stable identifier lists; intended only when defining a new sample."""
    chembl_ids: list[str] = []
    for offset in CHEMBL_OFFSETS:
        query = urlencode({
            "limit": CHEMBL_PAGE_SIZE,
            "offset": offset,
            "molecule_structures__isnull": "false",
            "order_by": "molecule_chembl_id",
            "only": "molecule_chembl_id",
        })
        payload = _json_request(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?{query}")
        chembl_ids.extend(item["molecule_chembl_id"] for item in payload["molecules"])
    if len(chembl_ids) != SAMPLE_SIZE or len(set(chembl_ids)) != SAMPLE_SIZE:
        raise RuntimeError(f"expected {SAMPLE_SIZE} distinct ChEMBL IDs, got {len(set(chembl_ids))}")
    _write_ids(sample_dir / DATABASE_SAMPLES["chembl37_sample"].ids_filename, chembl_ids)

    candidates = random.Random(PUBCHEM_SEED).sample(
        range(1, PUBCHEM_CID_RANGE_MAX + 1), SAMPLE_SIZE * 2)
    selected: list[str] = []
    endpoint = DATABASE_SAMPLES["pubchem_sample"].source_url
    for start in range(0, len(candidates), 100):
        batch = candidates[start:start + 100]
        payload = _json_request(endpoint, urlencode({"cid": ",".join(map(str, batch))}).encode())
        selected.extend(str(item["CID"]) for item in payload.get("PropertyTable", {}).get("Properties", []))
        if len(selected) >= SAMPLE_SIZE:
            break
        time.sleep(0.2)
    selected = selected[:SAMPLE_SIZE]
    if len(selected) != SAMPLE_SIZE or len(set(selected)) != SAMPLE_SIZE:
        raise RuntimeError(f"expected {SAMPLE_SIZE} distinct PubChem CIDs, got {len(set(selected))}")
    _write_ids(sample_dir / DATABASE_SAMPLES["pubchem_sample"].ids_filename, selected)


def _chembl_rows(ids: list[str]) -> list[tuple[str, str]]:
    found: dict[str, str] = {}
    for start in range(0, len(ids), 50):
        batch = ids[start:start + 50]
        query = urlencode({
            "limit": len(batch),
            "molecule_chembl_id__in": ",".join(batch),
            "only": "molecule_chembl_id,molecule_structures",
        })
        payload = _json_request(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?{query}")
        for item in payload["molecules"]:
            structures = item.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            if smiles:
                found[item["molecule_chembl_id"]] = smiles
        time.sleep(0.1)
    missing = [identifier for identifier in ids if identifier not in found]
    if missing:
        raise RuntimeError(f"ChEMBL returned no structure for {len(missing)} fixed IDs: {missing[:5]}")
    return [(identifier, found[identifier]) for identifier in ids]


def _pubchem_rows(ids: list[str]) -> list[tuple[str, str]]:
    found: dict[str, str] = {}
    endpoint = DATABASE_SAMPLES["pubchem_sample"].source_url
    for start in range(0, len(ids), 100):
        batch = ids[start:start + 100]
        payload = _json_request(endpoint, urlencode({"cid": ",".join(batch)}).encode())
        for item in payload.get("PropertyTable", {}).get("Properties", []):
            smiles = item.get("SMILES")
            if smiles:
                found[str(item["CID"])] = smiles
        time.sleep(0.2)
    missing = [identifier for identifier in ids if identifier not in found]
    if missing:
        raise RuntimeError(f"PubChem returned no structure for {len(missing)} fixed CIDs: {missing[:5]}")
    return [(identifier, found[identifier]) for identifier in ids]


def _write_sample(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["id", "smiles"])
        writer.writerows(rows)


def fetch_database_samples(data_dir: Path, sample_dir: Path) -> dict[str, dict[str, Any]]:
    metadata: dict[str, dict[str, Any]] = {}
    for name, spec in DATABASE_SAMPLES.items():
        ids = _read_ids(sample_dir / spec.ids_filename)
        if len(ids) != SAMPLE_SIZE or len(set(ids)) != SAMPLE_SIZE:
            raise RuntimeError(f"{spec.ids_filename} must contain {SAMPLE_SIZE} distinct identifiers")
        path = data_dir / spec.filename
        cached = path.exists() and spec.expected_sha256 is not None and sha256(path) == spec.expected_sha256
        if not cached:
            rows = _chembl_rows(ids) if name == "chembl37_sample" else _pubchem_rows(ids)
            _write_sample(path, rows)
        digest = sha256(path)
        if spec.expected_sha256 is not None and digest != spec.expected_sha256:
            raise RuntimeError(f"checksum mismatch for {name}: {digest}")
        metadata[name] = {
            "url": spec.source_url,
            "sha256": digest,
            "rows": SAMPLE_SIZE,
            "identifier_file": f"validation/samples/{spec.ids_filename}",
            "selection": (
                f"10 evenly spaced pages of {CHEMBL_PAGE_SIZE} structure-bearing records, "
                f"ordered by molecule_chembl_id in {CHEMBL_VERSION}"
                if name == "chembl37_sample"
                else f"seeded live-CID sample (seed {PUBCHEM_SEED}) from 1..{PUBCHEM_CID_RANGE_MAX}"
            ),
        }
    metadata["chembl37_sample"].update({
        "database_version": CHEMBL_VERSION,
        "release_date": CHEMBL_RELEASE_DATE,
        "structure_population": CHEMBL_STRUCTURE_COUNT,
        "status_url": "https://www.ebi.ac.uk/chembl/api/data/status.json",
    })
    metadata["pubchem_sample"].update({
        "population_date": PUBCHEM_POPULATION_DATE,
        "live_compound_count": PUBCHEM_LIVE_COMPOUND_COUNT,
        "cid_range_max": PUBCHEM_CID_RANGE_MAX,
        "population_source": "https://pubchem.ncbi.nlm.nih.gov/docs/statistics",
    })
    (data_dir / "database_samples_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def load_database_metadata(data_dir: Path) -> dict[str, dict[str, Any]]:
    return json.loads((data_dir / "database_samples_metadata.json").read_text())
