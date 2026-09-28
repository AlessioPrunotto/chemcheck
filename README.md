# chemcheck — pytest for molecular datasets

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/AlessioPrunotto/chemcheck/actions/workflows/ci.yml/badge.svg)](https://github.com/AlessioPrunotto/chemcheck/actions/workflows/ci.yml)
<!-- Uncomment after the first GitHub release + Zenodo hookup:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

One command to sanity-check a chemistry dataset before you train, publish, or trust it:

```bash
chemcheck dataset.csv
```

It checks for:
 - valid molecules
 - duplicates
 - stereochemical collisions
 - salt/tautomer
 - ambiguity
 - split leakage
 - target shift
 - chemical-space bias
 - etc.

Every warning comes with row IDs, evidence, and an actionable recommendation, plus a 0–100 dataset quality score.

## Install

```bash
conda install -c conda-forge rdkit pandas scipy pyarrow openpyxl
pip install -e .              # core package
pip install -e ".[pretty]"   # optional colored terminal output
pip install -e ".[dev]"      # contributor tools and coverage reporting
```

RDKit is a hard dependency and must come from conda (pip wheels also work where
available). chemcheck fails fast with an install hint if RDKit is missing.

## Inputs

CSV/TSV, SDF/SD, Parquet, Excel, JSON-lines. The SMILES column is autodetected
(`smiles`, `canonical_smiles`, …) or pass `--smiles-col`. Splits via `--split-col`
or by passing two files (`train.csv test.csv`).

## Exit codes (pytest-like)

- `0` clean
- `1` warnings
- `2` errors

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

Split values named `train`/`test`/`valid` are recognized automatically. For
custom values, repeat `--train-value` and `--test-value` as needed. With k-fold
data, passing only `--test-value FOLD` treats every other fold as training data.
Unresolved or ignored split values produce a warning instead of silently
skipping leakage checks.

The quality score is a transparent prioritization heuristic, not a validated
scientific metric. Related findings are overlap-capped so the same underlying
invalid structure, duplicate, or leakage issue is not fully deducted multiple
times. JSON reports record tool/RDKit versions, settings, resolved columns, and
SHA-256 input hashes for reproducibility.
When a large-dataset check uses sampling or a bounded candidate search, JSON
reports expose it in `meta.approximations` and in the finding's `metadata`.

## Learn

- **Tutorial:** [`docs/tutorial.md`](docs/tutorial.md) — from install to CI-gated
  audits in ~20 minutes.
- **Scientific validation:** [`docs/validation.md`](docs/validation.md) — public
  ML datasets plus fixed ChEMBL/PubChem samples, seeded defects, threshold
  sensitivity, false-positive guidance, and runtime through 41,127 molecules.
- **Release process:** [`docs/releasing.md`](docs/releasing.md) — TestPyPI
  verification and tokenless trusted publishing.
- **Examples:** [`examples/`](examples/) — three executed notebooks, all data
  generated inline (no downloads):
  - `01_quickstart.ipynb` — CLI + Python API on a deliberately dirty dataset.
  - `02_leakage_splits.ipynb` — identity vs analog vs scaffold leakage, and fixes.
  - `03_curation_case_study.ipynb` — full curation loop (39 → 80) with HTML report.
