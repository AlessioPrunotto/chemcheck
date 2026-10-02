# ChemDataCheck tutorial

With this tutorial, you will learn how to install chemdatacheck, inspect a dataset,
and configure an automated check that runs whenever the project changes.
This tutorial focuses on short command-line examples.

If you prefer interactive explanations, plots,
and experiments, see the following notebooks in [`examples/`](../examples/).

| notebook | what you learn | time |
|---|---|---|
| `01_quickstart.ipynb` | CLI + Python API on a dirty dataset | 10 min |
| `02_leakage_splits.ipynb` | identity vs analog vs scaffold leakage, and how to fix each | 15 min |
| `03_curation_case_study.ipynb` | full curation loop: dirty file with quality score going from 39 to 80 after cleaning. HTML report included | 20 min |

In the notebooks, all data is generated inline with RDKit, therefore all notebooks can be run offline.

## 1. Install

ChemDataCheck requires Python 3.10 or newer and is tested with Python 3.10–3.12.
Create an isolated environment and install the package from PyPI. On macOS or
Linux:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install "chemdatacheck[pretty]"
```

The `pretty` option adds colored terminal output. You can replace the final
command with `python -m pip install chemdatacheck` if you do not want it. Windows activation
instructions and the Conda fallback for systems without a compatible RDKit
wheel are documented in the [README](../README.md#install).

## 2. Your first audit (60 seconds)

Download the deliberately problematic demo dataset. It contains examples of
many ChemDataCheck findings, including two invalid structures and several
dataset-level problems involving otherwise valid molecules.

```bash
curl -LO https://raw.githubusercontent.com/AlessioPrunotto/chemdatacheck/main/tests/fixtures/demo.csv
chemdatacheck demo.csv --label-col activity --split-col split
```
`label-col`: tells chemdatacheck which column contains the target value (typically activity, pIC50, etc.). Given these values, chemdatacheck can monitor conflicting measurements, unusual values, etc. <br>
`split-col`: the column which splits your data into training set and test set

Your output will look like:

```
✓ 18 valid molecules / ✕ 2 invalid (n=20)
✕ [error] 2D-identity leakage (2 test molecules seen in train): 2  [identity_leakage]
    e.g. {"train_rows": ["m4"], "test_row": "m5", "connectivity": "CC(Cl)Br"}
    → Remove leaked test molecules or move them to train; ...
✕ [error] conflicting measurements (same structure, different labels): 2  [conflicting_measurements]
    e.g. {"structure": "Oc1ccccc1", "rows": ["m14", "m15"], "labels": ["2.5", "7.9"], ...}
    → Curate conflicts: check assay conditions, average with uncertainty, ...
...
Dataset quality score: 39/100
  top deductions: cross_split_near_duplicates -8.61, identity_leakage -7.41, ...
```

If you input two files, these will be automatically read as training set and test set (no `--split-col` needed):

```bash
chemdatacheck train.csv test.csv --label-col activity
```

If your split column does not have explicit values "training", "test", but rather custom values, you can map
them explicitly. For example, you can select which value of `fold` is
associated to the test set, and chemdatacheck will automatically treats every other value as training set:

```bash
chemdatacheck dataset.csv --split-col split --test-value 0
```

Repeat `--train-value` or `--test-value` when multiple values belong to one
side. Unrecognized or ignored split values produce a warning rather than
silently disabling leakage checks.

## 3. Understanding the output

Every finding has the same anatomy:

- **severity:** `error` (must fix: invalid chemistry, leakage, contradictions),
  `warning` (probably fix: duplicates, unspecified stereo, easy splits),
  `info` (know about it: isotopes, tautomers, domain gaps).
- **check_id:** stable machine-readable name (`identity_leakage`,
  `conflicting_measurements`, …). When processing JSON reports automatically,
  use check_id rather than matching the human-readable title.
- **count / rate + affected_rows:** how many rows, and *which* rows (your IDs).
- **examples:** the evidence: molecule pairs, Tanimoto scores, shared scaffolds.
- **recommendation:** what to do about it.
- **score:** 0–100 with a printed deduction breakdown, so any score is auditable.
  Weights live in `src/chemdatacheck/scoring.py` and every deduction scales with
  prevalence: a handful of bad rows costs little, a systemic problem costs a lot.

A warning should not be just a count for the user: it's meant to deliver relevant information.

## 4. Splits and leakage (the important chapter)

Based on the finding, you can have 3 different verdicts. Notebook 02 (`examples/02_leakage_splits.ipynb`) describes an example of one for each:

1. **Identity leakage** (`error`): same connectivity in train and test.
   This is never acceptable, as metrics become memorization scores. Fix: drop or move.
2. **Analog leakage** (`warning`, Morgan Tc ≥ 0.6, tunable via `--analog-thresh`):
   the test molecule is a close cousin of a training molecule. The model is
   interpolating, not generalizing. Fix: keep chemical series in one split.
3. **Scaffold overlap** (`warning`): shared Bemis–Murcko scaffold, dissimilar
   molecules. Context-dependent: fine for lead-optimization claims, fatal for
   "works on new chemotypes" claims. Match the split to the claim.

Rule of thumb: re-audit after every fix. The score should improve every time a fix is performed (notebook 02 goes 71 → 82 → 93).

## 5. Python API cookbook

You can also run chemdatacheck from Python instead of the command line. Pass the
input files as a list; one file is audited as a single dataset, while two files
are treated as training and test data.

```python
from chemdatacheck.report import audit, render_html, render_json

report = audit(["dataset.csv"], label_col="activity", split_col="split")
print(report.score, f"{report.n_valid}/{report.n_total} valid")
```

`audit()` returns an `AuditReport`. Its `findings` list contains the problems
and review signals detected by the individual checks.

### Recipe 1: inspect one type of finding.
If you are particularly interested in one specific type of finding, you can inspect that one alone. For example,
if you are particularly interested in leakage of analogs from the training set to the test set:

```python
analog_leakage = next(
    (finding for finding in report.findings if finding.check_id == "analog_leakage"),
    None,
)

if analog_leakage is None:
    print("No close train/test analogues were found.")
else:
    for example in analog_leakage.examples:
        print(
            "test row", example["test_row"],
            "has similarity", example["tanimoto"],
            "to train row", example["train_row"],
        )
```

### Recipe 2: make similarity checks more sensitive.
You can also decide to increase the sensitivity to similarity checks.

```python
more_sensitive = audit(
    ["train.csv", "test.csv"],
    analog_thresh=0.5,
    near_dup_thresh=0.9,
)
```

The defaults are `0.60` for analogues and `0.95` for near duplicates. Lowering
these thresholds reports more pairs, so this configuration is **more
sensitive**. Raising the same values will report fewer similar pairs.

For repeated numeric measurements, `conflict_thresh` specifies how far apart
the largest and smallest values for the same structure may be before chemdatacheck
reports a conflict. The value uses the same units as the label. Here, repeated
pIC50 measurements are reported when their spread is greater than `0.5`:

```python
assay = audit(["assay.csv"], label_col="pIC50", conflict_thresh=0.5)
```

On large datasets, it would be impractical to compare every molecule with all other molecules. In this case,
similarity checks use deterministic sampling. `report.meta["approximations"]`
records whether this happened and which checks were affected. The JSON report
contains the same information.

### Recipe 3: extract counts from the JSON report.
You can extract other relevant information, such as the number of corrupted SMILES, from the JSON report:

```python
import json

info = json.loads(render_json(report))
counts = {entry["check_id"]: entry["count"] for entry in info["findings"]}
print(counts)
```

This produces a dictionary such as `{"invalid_smiles": 2,
"analog_leakage": 14}`, which can be stored or compared between dataset
versions.

### Recipe 4: fail a pipeline on errors only.
With `--fail-on error`, warnings are allowed: the exit code is `0` when there
are no error-level findings and `2` when at least one error is found. Use
`--fail-on warning` if warnings should also fail the pipeline; in that mode a
warning produces exit code `1` and an error produces exit code `2`.

```python
from chemdatacheck.cli import main as chemdatacheck_main

exit_code = chemdatacheck_main(["dataset.csv", "--fail-on", "error", "--quiet"])
raise SystemExit(exit_code)
```

### Recipe 5: save an HTML quality-control report.

```python
with open("dataset_qc.html", "w", encoding="utf-8") as fh:
    fh.write(render_html(report))
```

The resulting file is a standalone report that can be opened in a browser or
shared with the dataset.

## 6. CI integration

Run ChemDataCheck automatically whenever the repository changes. The following GitHub Actions step fails if dataset.sdf contains any error-level finding and writes a JUnit report that CI tools can display.

```yaml
# .github/workflows/chemdatacheck.yml
- name: Audit dataset
  run: chemdatacheck data/dataset.sdf --split-col split --fail-on error --format junit -o chemdatacheck.xml
```

## 7. FAQ

**Is this a replacement for RDKit / MolVS / the ChEMBL curation pipeline?**
No. chemdatacheck is built *on top* of RDKit. Standardizers fix molecules one at a time;
chemdatacheck reasons across rows and splits (duplicates, leakage, shift).
Use both: standardize first, then audit.

**How big a dataset can it handle?**
~100k molecules on a laptop. Exact checks are linear; pairwise similarity uses
blocking plus deterministic sampling above a few thousand rows, and reports
counts as lower bounds. Thresholds are tunable.

**My score is low. Where do I start?**
Errors first, in this order: `invalid_smiles` → `identity_leakage` /
`cross_split_near_duplicates` → `conflicting_measurements`. Then warnings by
deduction size (`score_breakdown` tells you). Notebook 03 walks a 39 → 80
curation end to end.

**A standardizer corrupted my molecule: now what?**
That happened to us too: RDKit's default salt stripper deletes neutral acetic
acid (carboxylic acids match its salt list). Defensive rule from notebook 03:
check every transformation's output, keep the input on failure, and let the
re-audit confirm the fix.

**How do I cite / reference the report?**
Export `--format json` (exact counts, row IDs, thresholds in one file) or the
HTML report alongside the dataset. Both record the chemdatacheck version
implicitly via the finding set. Pin the version in your environment.
