# chemcheck: a pytest for molecular datasets

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

From a cloned copy of this repository:

```bash
python -m pip install .             # core package
python -m pip install ".[pretty]"  # core package plus colored terminal output
```

Choose one of these commands; they are alternatives, not consecutive steps.
Required dependencies (including RDKit, pandas, SciPy, PyArrow, and openpyxl) are
declared in `pyproject.toml` and installed automatically by pip.

If a compatible RDKit wheel is not available for your Python version or
platform, install RDKit from conda-forge first and then install chemcheck:

```bash
conda install -c conda-forge rdkit
python -m pip install .
```

Contributors who need an editable installation, tests, coverage reporting, or
the notebooks should follow [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Inputs

CSV/TSV, SDF/SD, Parquet, Excel, JSON-lines. The SMILES column is autodetected
(`smiles`, `canonical_smiles`, ...); alternatively, you can pass `--smiles-col`. Splits (training set / test set)
can be passed via a column (`--split-col`) or by passing two files (`train.csv test.csv`).

## Exit codes (pytest-like)

- `0` clean
- `1` warnings
- `2` errors

## What it checks

- **Chemical integrity:** invalid SMILES, valence errors, aromaticity/kekulization,
  impossible charges, disconnected components (salts/mixtures), isotopes,
  radicals, unspecified stereocenters, tautomer ambiguity.
- **Dataset duplicates:** exact, canonical, stereochemical collisions, salt/solvate
  duplicates, tautomer duplicates, near-duplicates (Morgan fingerprints, similarity
  measured with Tanimoto distance).
- **Split leakage:** 2D-identity leakage, analog leakage (Tc ≥ 0.6), near-duplicates
  across splits, scaffold overlap, suspiciously-easy-split heuristic.
- **ML readiness:** duplicated/conflicting measurements, target distribution shift,
  split-predicts-label leakage suspects, label outliers.
- **Chemical space:** rare elements, rare/promiscuous functional groups, unusual
  ring systems, representation bias, applicability-domain gaps.

Every finding carries row IDs, evidence examples, and a recommendation, e.g.
`train_row` ↔ `test_row` pairs with Tanimoto and shared scaffold for leakage.

To run leakage checks, Chemcheck must know which rows are training data and
which are held out for testing or validation. Common split values such as
`train`, `test`, and `valid` are recognized automatically. If your dataset uses
different names, map them explicitly; for example:

```bash
chemcheck dataset.csv --split-col partition \
  --train-value development --test-value external
```

For numbered cross-validation folds, you only need to identify the held-out
fold. This treats fold `0` as test data and every other observed fold as
training data:

```bash
chemcheck dataset.csv --split-col fold --test-value 0
```

If Chemcheck cannot form both a non-empty training group and a non-empty test
group, it reports a warning and does not run the leakage checks. If only some
split values are mapped, it warns that the remaining rows were excluded from
those checks.

The quality score is a prioritization heuristic, not a validated scientific metric.
The quality score starts from 100, and points are subtracted for each bad quality finding.
Related findings are overlap-capped so the same underlying
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
  sensitivity, false-positive guidance, and runtime through 40k+ molecules.
- **Examples:** [`examples/`](examples/) — three executed notebooks, all data
  generated inline (no downloads):
  - `01_quickstart.ipynb` — CLI + Python API on a deliberately dirty dataset.
  - `02_leakage_splits.ipynb` — identity vs analog vs scaffold leakage, and fixes.
  - `03_curation_case_study.ipynb` — full curation loop (39 → 80) with HTML report.
