# Changelog

All notable changes to ChemDataCheck are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.0] - 2026-10-02

First public release as ChemDataCheck.

### Added
- Explicit `--train-value` / `--test-value` split mapping, including one-option
  holdout-fold mode, plus warnings when split values cannot be interpreted.
- Reproducibility metadata in reports: chemdatacheck/RDKit versions, audit settings,
  resolved columns, and SHA-256 hashes of inputs.
- CLI `--version` and validation for thresholds and example limits.
- Machine-readable approximation metadata at report and finding level.
- Configurable numeric measurement conflicts via `--conflict-thresh` /
  `conflict_thresh`.
- Package authorship, classifier, keyword, and project URL metadata.
- Reproducible scientific-validation harness and committed results covering six
  public datasets (642–41,127 molecules), seeded defects, threshold sensitivity,
  approximation fidelity, expected false positives, and runtime scaling.
- Fixed 5,000-structure ChEMBL 37 and PubChem validation samples, with committed
  identifier lists, official-API acquisition, checksums, and structure-only
  audit results kept distinct from ML-labelled benchmarks.
- Tested minimum versions for runtime and optional dependencies, branch coverage
  reporting, distribution build/install checks, and a Trusted Publishing
  workflow with TestPyPI artifact verification before production release.
- Edge-case coverage for input inference and validation, empty datasets, CLI
  failures and report formats, target-shift modes, report filters, and check
  crash isolation.

### Changed
- Separated user-facing validation guidance from the study's maintainer and
  reviewer notes.
- Related score deductions are capped by overlap family to avoid repeatedly
  penalizing the same invalid-structure, duplicate, split-similarity, or target-shift issue.
- Cohen's d now uses the conventional pooled sample standard deviation.
- Two-file mode preserves resolved-column metadata and permits labels on only one input.
- Removed the unused scikit-learn runtime dependency.
- Cross-split similarity is computed once per audit and reused by analog,
  near-duplicate, and easy-split checks.
- Approximation metadata is now recorded by the check execution path itself,
  eliminating duplicate prediction logic in report assembly.

### Fixed
- Standard JSON arrays are parsed as ordinary JSON rather than being silently
  mis-shaped by the JSON-lines reader.
- Binary 0/1 targets are treated as categorical labels: minority-class rows are
  no longer reported as numeric outliers, and contradictory binary labels for
  one structure are reported as measurement conflicts.
- Large-dataset applicability-gap candidates are now exactly verified against
  the full dataset, eliminating false alerts introduced by reference sampling.

### Removed
- Unused internal RDKit-availability and split-detection helpers.

### Core capabilities
- `chemdatacheck` CLI: audit CSV/TSV, SDF/SD, Parquet, Excel, JSON-lines datasets
  with SMILES-column autodetection; single-file (`--split-col`) and two-file
  train/test modes; terminal, JSON, HTML, and JUnit reports; pytest-like
  exit codes with `--fail-on`.
- 30 dataset checks across 5 families: chemical integrity (invalid SMILES, valence,
  aromaticity, charges, disconnected components, isotopes, radicals,
  unspecified stereo, tautomers), dataset duplicates (exact, canonical,
  stereochemical, salt, tautomer, near-duplicates), split leakage
  (2D-identity, analog, cross-split near-duplicates, scaffold overlap,
  easy-split heuristic), ML readiness (duplicated/conflicting measurements,
  target shift, split-predicts-label, outliers), chemical space (rare
  elements/groups, rings, representation bias, applicability gaps).
- Transparent 0–100 dataset quality score with per-check deduction breakdown.
- Python API (`chemdatacheck.report.audit`) plus tutorial (`docs/tutorial.md`)
  and three executed example notebooks (`examples/`).
