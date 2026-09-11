# Changelog

All notable changes to chemcheck are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.0] - 2026-09-11

First public release.

### Added
- `chemcheck` CLI: audit CSV/TSV, SDF/SD, Parquet, Excel, JSON-lines datasets
  with SMILES-column autodetection; single-file (`--split-col`) and two-file
  train/test modes; terminal, JSON, HTML, and JUnit reports; pytest-like
  exit codes with `--fail-on`.
- 26 checks across 5 families: chemical integrity (invalid SMILES, valence,
  aromaticity, charges, disconnected components, isotopes, radicals,
  unspecified stereo, tautomers), dataset duplicates (exact, canonical,
  stereochemical, salt, tautomer, near-duplicates), split leakage
  (2D-identity, analog, cross-split near-duplicates, scaffold overlap,
  easy-split heuristic), ML readiness (duplicated/conflicting measurements,
  target shift, split-predicts-label, outliers), chemical space (rare
  elements/groups, rings, representation bias, applicability gaps).
- Transparent 0–100 dataset quality score with per-check deduction breakdown.
- Python API (`chemcheck.report.audit`) plus tutorial (`docs/tutorial.md`)
  and three executed example notebooks (`examples/`).
