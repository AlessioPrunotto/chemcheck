# chemcheck — pytest for molecular datasets

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/AlessioPrunotto/chemcheck/actions/workflows/ci.yml/badge.svg)](https://github.com/AlessioPrunotto/chemcheck/actions/workflows/ci.yml)
<!-- Uncomment after the first GitHub release + Zenodo hookup:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

One command to sanity-check a chemistry dataset before you train, publish, or trust it:

```bash
chemcheck dataset.csv
chemcheck train.csv test.csv --label-col activity --format json
chemcheck dataset.sdf --split-col split --fail-on error
```

It reports valid molecules, duplicates, stereochemical collisions, salt/tautomer
ambiguity, split leakage, target shift, and chemical-space bias.
Every warning comes with row IDs, evidence, and an actionable recommendation, plus
a 0–100 dataset quality score.

## Install

```bash
conda install -c conda-forge rdkit pandas scikit-learn scipy pyarrow openpyxl
pip install -e .          # plus `pip install rich` for colored output, `pytest` for dev
```

RDKit is a hard dependency and must come from conda (pip wheels also work where
available). chemcheck fails fast with an install hint if RDKit is missing.

## Inputs

CSV/TSV, SDF/SD, Parquet, Excel, JSON-lines. The SMILES column is autodetected
(`smiles`, `canonical_smiles`, …) or pass `--smiles-col`. Splits via `--split-col`
or by passing two files (`train.csv test.csv`).

## Exit codes (pytest-like)

- `0` clean · `1` warnings · `2` errors. Tune with `--fail-on error|warning|info|never`
  for CI gates. `--format junit` emits JUnit XML for CI dashboards.

## What it checks

- **Chemical integrity:** invalid SMILES, valence errors, aromaticity/kekulization,
  impossible charges, disconnected components (salts/mixtures), isotopes,
  radicals, unspecified stereocenters, tautomer ambiguity.
- **Dataset duplicates:** exact, canonical, stereochemical collisions, salt/solvate
  duplicates, tautomer duplicates, near-duplicates (Morgan Tanimoto).
- **Split leakage:** 2D-identity leakage, analog leakage (Tc ≥ 0.6), near-duplicates
  across splits, scaffold overlap, suspiciously-easy-split heuristic.
- **ML readiness:** duplicated/conflicting measurements, target distribution shift,
  split-predicts-label leakage suspects, label outliers.
- **Chemical space:** rare elements, rare/promiscuous functional groups, unusual
  ring systems, representation bias, applicability-domain gaps.

Every finding carries row IDs, evidence examples, and a recommendation, e.g.
`train_row` ↔ `test_row` pairs with Tanimoto and shared scaffold for leakage.

## Learn

- **Tutorial:** [`docs/tutorial.md`](docs/tutorial.md) — from install to CI-gated
  audits in ~20 minutes.
- **Examples:** [`examples/`](examples/) — three executed notebooks, all data
  generated inline (no downloads):
  - `01_quickstart.ipynb` — CLI + Python API on a deliberately dirty dataset.
  - `02_leakage_splits.ipynb` — identity vs analog vs scaffold leakage, and fixes.
  - `03_curation_case_study.ipynb` — full curation loop (22 → 80) with HTML report.
