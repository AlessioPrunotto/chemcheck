#!/usr/bin/env python3
"""Reproducible scientific validation and scaling study for chemdatacheck.

Downloads public MoleculeNet data and fixed ChEMBL/PubChem samples, runs
natural-data audits, measures the sensitivity of the similarity thresholds,
and injects known defects into real ESOL records to test detection. Raw
third-party datasets stay untracked.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from database_samples import (
    DATABASE_SAMPLES,
    fetch_database_samples,
    initialize_identifier_lists,
    load_database_metadata,
)

from chemdatacheck.io import load_table
from chemdatacheck.molecules import build_records
from chemdatacheck.report import audit


@dataclass(frozen=True)
class Dataset:
    url: str
    filename: str
    sha256: str
    smiles: str
    label: str | None
    identifier: str | None = None
    split: str | None = None
    cohort: str = "ml_benchmark"


DATASETS = {
    "freesolv": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv",
        "freesolv.csv", "ab5895d914ee87cb563bd7b9611e869527bba45bec6b014d34dc495a0f9dcb72",
        "smiles", "expt", "iupac"),
    "delaney": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv",
        "delaney.csv", "8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189",
        "smiles", "measured log solubility in mols per litre", "Compound ID"),
    "bace": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/bace.csv",
        "bace.csv", "f3fb9ce90bada3e2bd6148b0df13f8f8145a357bf87df0dd5b391ede974fc737",
        "mol", "Class", "CID", "Model"),
    "lipophilicity": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv",
        "lipophilicity.csv", "aed41590cb30609d51d8e08ad3ff06495a76e80e211358801f596b10da69bacd",
        "smiles", "exp", "CMPD_CHEMBLID"),
    "tox21": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz",
        "tox21.csv.gz", "45d09792492ce049039dd24aa27b07fc79ce20c573187d4d90bcd178c0c0d360",
        "smiles", "NR-AR", "mol_id"),
    "hiv": Dataset(
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv",
        "hiv.csv", "9ffa7fe57dc86c342627ee1d5255e937e2ab812393c73c4d16c697022f6e1d22",
        "smiles", "HIV_active"),
}

for _name, _sample in DATABASE_SAMPLES.items():
    DATASETS[_name] = Dataset(
        _sample.source_url, _sample.filename, _sample.expected_sha256 or "",
        "smiles", None, "id", cohort="database_sample")

THRESHOLDS = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(data_dir: Path, sample_dir: Path) -> dict[str, dict[str, Any]]:
    data_dir.mkdir(parents=True, exist_ok=True)
    for name, spec in DATASETS.items():
        if spec.cohort == "database_sample":
            continue
        path = data_dir / spec.filename
        if not path.exists() or sha256(path) != spec.sha256:
            print(f"downloading {name}: {spec.url}", file=sys.stderr)
            urllib.request.urlretrieve(spec.url, path)
        actual = sha256(path)
        if actual != spec.sha256:
            raise RuntimeError(f"checksum mismatch for {name}: {actual}")
    return fetch_database_samples(data_dir, sample_dir)


def normalized_frame(name: str, data_dir: Path) -> pd.DataFrame:
    spec = DATASETS[name]
    raw = pd.read_csv(data_dir / spec.filename)
    frame = pd.DataFrame({
        "id": raw[spec.identifier].astype(str) if spec.identifier else raw.index.astype(str),
        "smiles": raw[spec.smiles],
    })
    if spec.label:
        frame["label"] = raw[spec.label]
    if spec.split:
        frame["split"] = raw[spec.split].astype(str).str.casefold()
    else:
        # Deterministic 80/20 row split emulates a common random holdout.
        frame["split"] = np.where(raw.index.to_numpy() % 5 == 0, "test", "train")
    return frame


def peak_rss_mb() -> float:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    return value / (1024 * 1024) if sys.platform == "darwin" else value / 1024


def natural_audits(data_dir: Path, work_dir: Path) -> tuple[list[dict], list[dict]]:
    summaries: list[dict] = []
    findings: list[dict] = []
    for name in DATASETS:
        frame = normalized_frame(name, data_dir)
        path = work_dir / f"{name}.csv"
        frame.to_csv(path, index=False)
        started = time.perf_counter()
        audit_kwargs = {"smiles_col": "smiles", "id_col": "id", "split_col": "split"}
        if "label" in frame:
            audit_kwargs["label_col"] = "label"
        report = audit([str(path)], **audit_kwargs)
        elapsed = time.perf_counter() - started
        summaries.append({
            "dataset": name, "cohort": DATASETS[name].cohort,
            "rows": len(frame), "valid": report.n_valid,
            "invalid": report.n_invalid, "score": report.score,
            "wall_seconds": round(elapsed, 3), "process_peak_rss_mb": round(peak_rss_mb(), 1),
            "finding_count": len(report.findings),
            "approximations": json.dumps(report.meta.get("approximations", []), sort_keys=True),
        })
        for finding in report.findings:
            findings.append({
                "dataset": name, "cohort": DATASETS[name].cohort,
                "check_id": finding.check_id,
                "severity": finding.severity.value, "count": finding.count,
                "rate": finding.rate, "title": finding.title,
                "approximate": bool(finding.metadata.get("approximate", False)),
            })
        print(f"audited {name}: {len(frame):,} rows in {elapsed:.2f}s", file=sys.stderr)
    return summaries, findings


def threshold_sweep(data_dir: Path, work_dir: Path) -> list[dict]:
    from rdkit import DataStructs

    rows: list[dict] = []
    # HIV is covered by the end-to-end benchmark; its full 80/20 NN sweep adds
    # little threshold information while dominating study time.
    for name in ("freesolv", "delaney", "bace", "lipophilicity", "tox21"):
        frame = normalized_frame(name, data_dir)
        path = work_dir / f"threshold-{name}.csv"
        frame.to_csv(path, index=False)
        normalized = load_table([str(path)], smiles_col="smiles", id_col="id",
                                label_col="label", split_col="split")
        records = build_records(normalized)
        train = [r for r in records if r.valid and r.fp is not None and r.split == "train"]
        test = [r for r in records if r.valid and r.fp is not None and r.split in {"test", "valid"}]
        train_fps = [r.fp for r in train]
        maxima = []
        for record in test:
            similarities = DataStructs.BulkTanimotoSimilarity(record.fp, train_fps)
            maxima.append(max(similarities, default=0.0))
        values = np.asarray(maxima, dtype=float)
        base = {
            "dataset": name, "train": len(train), "test": len(test),
            "median_max_tc": round(float(np.median(values)), 4),
            "p90_max_tc": round(float(np.quantile(values, 0.9)), 4),
            "p95_max_tc": round(float(np.quantile(values, 0.95)), 4),
        }
        for threshold in THRESHOLDS:
            rows.append({**base, "threshold": threshold,
                         "test_fraction_at_or_above": round(float(np.mean(values >= threshold)), 6)})
        print(f"swept thresholds for {name}", file=sys.stderr)
    return rows


def approximation_fidelity(data_dir: Path, work_dir: Path) -> list[dict]:
    """Compare the >4k-record approximation paths with exact computation."""
    from rdkit import DataStructs

    frame = normalized_frame("lipophilicity", data_dir)
    path = work_dir / "fidelity-lipophilicity.csv"
    frame.to_csv(path, index=False)
    records = build_records(load_table([str(path)], smiles_col="smiles", id_col="id",
                                       label_col="label", split_col="split"))
    valid = [record for record in records if record.valid and record.fp is not None]
    fps = [record.fp for record in valid]

    started = time.perf_counter()
    exact_pairs: set[tuple[str, str]] = set()
    for i, record in enumerate(valid):
        similarities = DataStructs.BulkTanimotoSimilarity(record.fp, fps[i + 1:])
        exact_pairs.update((record.row_id, valid[i + j + 1].row_id)
                           for j, similarity in enumerate(similarities) if similarity >= 0.95)
    near_exact_seconds = time.perf_counter() - started
    ordered = sorted(valid, key=lambda record: record.fp.GetNumOnBits())
    approximate_pairs: set[tuple[str, str]] = set()
    for i, record in enumerate(ordered):
        candidates = ordered[i + 1:i + 401]
        similarities = DataStructs.BulkTanimotoSimilarity(record.fp,
                                                           [item.fp for item in candidates])
        for candidate, similarity in zip(candidates, similarities):
            if similarity >= 0.95:
                approximate_pairs.add(tuple(sorted((record.row_id, candidate.row_id))))
    exact_pairs = {tuple(sorted(pair)) for pair in exact_pairs}
    shared_pairs = exact_pairs & approximate_pairs

    started = time.perf_counter()
    exact_isolated: set[str] = set()
    for i, record in enumerate(valid):
        similarities = DataStructs.BulkTanimotoSimilarity(record.fp, fps)
        similarities[i] = -1.0
        if max(similarities, default=0.0) < 0.3:
            exact_isolated.add(record.row_id)
    gap_exact_seconds = time.perf_counter() - started
    reference = valid[:: max(1, len(valid) // 3000)][:3000]
    reference_fps = [record.fp for record in reference]
    approximate_isolated: set[str] = set()
    for record in valid:
        similarities = DataStructs.BulkTanimotoSimilarity(record.fp, reference_fps)
        top_two = sorted(similarities, reverse=True)[:2]
        nearest = (top_two[1] if len(top_two) > 1 and top_two[0] >= 0.999
                   else (top_two[0] if top_two else 0.0))
        if nearest < 0.3:
            approximate_isolated.add(record.row_id)
    shared_isolated = exact_isolated & approximate_isolated
    refined_isolated = approximate_isolated & exact_isolated

    def ratio(numerator: int, denominator: int) -> float:
        return round(numerator / denominator, 6) if denominator else 1.0

    return [
        {"dataset": "lipophilicity", "check_id": "near_duplicates",
         "records": len(valid), "exact_positive": len(exact_pairs),
         "approximate_positive": len(approximate_pairs),
         "precision": ratio(len(shared_pairs), len(approximate_pairs)),
         "recall": ratio(len(shared_pairs), len(exact_pairs)),
         "exact_seconds": round(near_exact_seconds, 3)},
        {"dataset": "lipophilicity", "check_id": "applicability_gaps_prefilter",
         "records": len(valid), "exact_positive": len(exact_isolated),
         "approximate_positive": len(approximate_isolated),
         "precision": ratio(len(shared_isolated), len(approximate_isolated)),
         "recall": ratio(len(shared_isolated), len(exact_isolated)),
         "exact_seconds": round(gap_exact_seconds, 3)},
        {"dataset": "lipophilicity", "check_id": "applicability_gaps",
         "records": len(valid), "exact_positive": len(exact_isolated),
         "approximate_positive": len(refined_isolated),
         "precision": ratio(len(refined_isolated), len(refined_isolated)),
         "recall": ratio(len(refined_isolated), len(exact_isolated)),
         "exact_seconds": round(gap_exact_seconds, 3)},
    ]


def _alternative_smiles(smiles: str, mode: str) -> str | None:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if mode == "canonical":
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        for _ in range(20):
            candidate = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False, doRandom=True)
            if candidate != smiles and candidate != canonical:
                return candidate
    elif mode == "stereo":
        for atom in mol.GetAtoms():
            tag = atom.GetChiralTag()
            if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                return Chem.MolToSmiles(mol, isomericSmiles=True)
            if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                return Chem.MolToSmiles(mol, isomericSmiles=True)
    elif mode == "tautomer":
        enum = rdMolStandardize.TautomerEnumerator()
        original = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        for tautomer in enum.Enumerate(mol):
            candidate = Chem.MolToSmiles(tautomer, isomericSmiles=True, canonical=True)
            if candidate != original:
                return candidate
    return None


def _candidate(frame: pd.DataFrame, mode: str) -> tuple[pd.Series, str]:
    for _, row in frame.iterrows():
        alternative = _alternative_smiles(str(row["smiles"]), mode)
        if alternative:
            return row, alternative
    raise RuntimeError(f"no ESOL candidate found for seeded {mode} test")


def seeded_detection(data_dir: Path, work_dir: Path) -> list[dict]:
    base = normalized_frame("delaney", data_dir).iloc[:300].copy()
    base["id"] = [f"base-{i}" for i in range(len(base))]
    base["split"] = np.where(np.arange(len(base)) % 5 == 0, "test", "train")
    scenarios: list[tuple[str, str, pd.DataFrame, str]] = []

    source = base.iloc[1].copy()
    source["split"] = "train"
    base.iloc[1, base.columns.get_loc("split")] = "train"

    def add(name: str, expected: str, original: pd.Series, **changes: Any) -> None:
        injected = original.copy()
        injected["id"] = f"seed-{name}"
        for key, value in changes.items():
            injected[key] = value
        scenarios.append((name, expected, pd.concat([base, injected.to_frame().T], ignore_index=True),
                          str(injected["id"])))

    add("invalid", "invalid_smiles", source, smiles="C1(CC")
    add("exact_duplicate", "exact_duplicates", source)
    canonical_row, canonical_alt = _candidate(base, "canonical")
    add("canonical_duplicate", "canonical_duplicates", canonical_row, smiles=canonical_alt)
    add("salt_duplicate", "salt_duplicates", source, smiles=f'{source["smiles"]}.[Na+]')
    # ESOL's first 300 published SMILES contain no assigned stereocentre, so
    # use an explicitly stereochemical BACE record and include both forms.
    stereo_row, stereo_alt = _candidate(normalized_frame("bace", data_dir), "stereo")
    stereo_original = stereo_row.copy()
    stereo_original["id"] = "seed-stereo-original"
    stereo_original["split"] = "train"
    stereo_mutant = stereo_original.copy()
    stereo_mutant["id"] = "seed-stereo_collision"
    stereo_mutant["smiles"] = stereo_alt
    stereo_frame = pd.concat(
        [base, stereo_original.to_frame().T, stereo_mutant.to_frame().T], ignore_index=True)
    scenarios.append(("stereo_collision", "stereochemical_collisions", stereo_frame,
                      str(stereo_mutant["id"])))
    tautomer_row, tautomer_alt = _candidate(base, "tautomer")
    add("tautomer_duplicate", "tautomer_duplicates", tautomer_row, smiles=tautomer_alt)
    add("repeated_measurement", "duplicated_measurements", source)
    add("conflicting_measurement", "conflicting_measurements", source,
        label=float(source["label"]) + 2.0)
    add("identity_leakage", "identity_leakage", source, split="test")

    results = []
    for name, expected, frame, injected_id in scenarios:
        path = work_dir / f"seed-{name}.csv"
        frame.to_csv(path, index=False)
        report = audit([str(path)], smiles_col="smiles", id_col="id",
                       label_col="label", split_col="split")
        finding = next((item for item in report.findings if item.check_id == expected), None)
        detected = finding is not None and injected_id in finding.affected_rows
        results.append({"scenario": name, "expected_check": expected,
                        "injected_row": injected_id, "detected": detected})
        print(f"seeded {name}: {'detected' if detected else 'MISSED'}", file=sys.stderr)
    return results


def write_results(output_dir: Path, summaries: list[dict], findings: list[dict],
                  thresholds: list[dict], seeded: list[dict], fidelity: list[dict],
                  database_metadata: dict[str, dict[str, Any]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(output_dir / "benchmark.csv", index=False)
    pd.DataFrame(findings).to_csv(output_dir / "natural_findings.csv", index=False)
    pd.DataFrame(thresholds).to_csv(output_dir / "threshold_sweep.csv", index=False)
    pd.DataFrame(seeded).to_csv(output_dir / "seeded_detection.csv", index=False)
    pd.DataFrame(fidelity).to_csv(output_dir / "approximation_fidelity.csv", index=False)
    import rdkit
    metadata = {
        "python": platform.python_version(), "platform": platform.platform(),
        "rdkit": rdkit.__version__, "datasets": {
            name: (database_metadata[name] if spec.cohort == "database_sample"
                   else {"url": spec.url, "sha256": spec.sha256})
            for name, spec in DATASETS.items()
        },
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("validation/data"))
    parser.add_argument("--output-dir", type=Path, default=Path("validation/results"))
    parser.add_argument("--sample-dir", type=Path, default=Path("validation/samples"))
    parser.add_argument("--download", action="store_true")
    parser.add_argument(
        "--initialize-database-samples",
        action="store_true",
        help="replace the committed ChEMBL/PubChem ID snapshots (requires --download)",
    )
    args = parser.parse_args()
    if args.initialize_database_samples and not args.download:
        parser.error("--initialize-database-samples requires --download")
    if args.initialize_database_samples:
        initialize_identifier_lists(args.sample_dir)
    if args.download:
        database_metadata = fetch(args.data_dir, args.sample_dir)
    else:
        database_metadata = load_database_metadata(args.data_dir)
    missing = [name for name, spec in DATASETS.items() if not (args.data_dir / spec.filename).exists()]
    if missing:
        parser.error(f"missing datasets {missing}; run with --download")
    import tempfile
    with tempfile.TemporaryDirectory(prefix="chemdatacheck-validation-") as temporary:
        work_dir = Path(temporary)
        summaries, findings = natural_audits(args.data_dir, work_dir)
        thresholds = threshold_sweep(args.data_dir, work_dir)
        seeded = seeded_detection(args.data_dir, work_dir)
        fidelity = approximation_fidelity(args.data_dir, work_dir)
    write_results(args.output_dir, summaries, findings, thresholds, seeded, fidelity,
                  database_metadata)
    return 0 if all(row["detected"] for row in seeded) else 1


if __name__ == "__main__":
    raise SystemExit(main())
