# Validation maintainer notes

This file is for the maintainers/reviewers of chemcheck's scientific
validation. User-facing interpretation belongs in
[`docs/validation.md`](../docs/validation.md).

## Reproduce the published validation

From the repository root, run:

```bash
python validation/run_validation.py --download
```

The run writes machine-readable results under `validation/results/`. Record the
software versions, hardware, date, and whether cached downloads were used when
publishing new benchmark numbers.

The study has four components:

1. seeded-defect sensitivity;
2. characterization of findings in unmodified public data;
3. comparison of scalable similarity paths with exact calculations; and
4. runtime measurement on collections of realistic size.

## Public database snapshots

The ChEMBL and PubChem samples are snapshots, not live queries during ordinary
validation runs.

- Compound identifiers are stored in `validation/samples/`.
- Retrieval URLs, timestamps, HTTP validators, and SHA-256 hashes are stored in
  `validation/data/database_samples_metadata.json`; run-level environment and
  source metadata are stored in `validation/results/metadata.json`.
- Do not pass `--initialize-database-samples` during routine reproduction. It is
  an explicit sample-initialization operation that replaces the committed
  identifier lists.
- Do not silently accept an upstream response whose hash has changed. Inspect
  the difference, decide whether the snapshot should be updated, and document
  the reason.

Keeping the identifiers and hashes versioned makes the study reproducible while
acknowledging that public APIs and database records can change.

## Default-threshold decisions

The current defaults are:

- analog warning: Tanimoto similarity `>= 0.60`;
- near-duplicate review: Tanimoto similarity `>= 0.95`.

The public-dataset sweep supports retaining `0.60` as an early-warning default
and `0.95` as a high-similarity review threshold. Neither value is an identity
test or a universal definition of leakage.

Change a default only after repeating the validation and documenting:

- the fingerprint and molecular representation used;
- the effect on seeded-defect detection;
- the effect on unmodified public datasets and likely false positives;
- approximation fidelity on the large-data path; and
- the migration impact on existing users and CI configurations.

The overall score weights and severity mapping are also policy choices. Validate
them against real curation decisions before presenting the score as a calibrated
measure of dataset quality.

## Findings that changed the implementation

The validation exercise exposed two issues now covered by regression tests:

1. Numeric `0/1` endpoints could be interpreted as regression because pandas
   loaded them as integers. Low-cardinality numeric targets are now treated as
   categorical for target diagnostics, and contradictory binary labels for one
   structure are recognized as conflicts.
2. Reference-sampled applicability gaps could include molecules whose true
   nearest neighbor was omitted from the sample. Candidates are now verified
   against the complete dataset before being reported.

Do not remove or weaken these regression cases when refactoring the checks.

## Approximation checks

Large-data paths must be compared with exact calculations on controlled
fixtures. At minimum, verify:

- recovery of injected close cross-split pairs;
- exact confirmation of candidates returned by a bounded neighbor search;
- agreement on applicability gaps when the reference set is sampled; and
- recording of the selected mode in report metadata.

The current study establishes perfect recovery only for the tested fixtures. It
does not prove perfect recall for every dataset or fingerprint.

## When to rerun the study

Rerun the relevant validation sections when changing:

- check semantics or severities;
- default thresholds or fingerprints;
- standardization or scaffold logic;
- sampling, neighbor-search, or applicability-domain code;
- supported RDKit or dependency versions; or
- report fields used to interpret validation results.

Update the user-facing page and committed result files together when conclusions
or benchmark numbers change.

## Future validation work

Useful extensions include:

- blinded expert review of mixed clean and perturbed records, enabling estimates
  of precision, recall, and inter-rater agreement;
- evaluation with additional fingerprints and RDKit versions;
- larger approximation-fidelity studies across several random seeds; and
- repeated runtime measurements on isolated hardware with peak-memory tracking.

These extensions would strengthen external validity, but they should not be
presented as completed evidence until the corresponding results are available.
