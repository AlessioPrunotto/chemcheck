# chemcheck tutorial: "a pytest for molecular datasets"

In about 20 minutes, you will learn how to install chemcheck, inspect a dataset,
and configure an automated check that runs whenever the project changes.
The tutorial focuses on short command-line examples. For interactive explanations, plots,
and experiments, see the notebooks in [`examples/`](../examples/).

| notebook | what you learn | time |
|---|---|---|
| `01_quickstart.ipynb` | CLI + Python API on a dirty dataset | 10 min |
| `02_leakage_splits.ipynb` | identity vs analog vs scaffold leakage, and how to fix each | 15 min |
| `03_curation_case_study.ipynb` | full curation loop: dirty file with quality score going from 39 to 80 after cleaning. HTML report included | 20 min |

In the notebooks, all data is generated inline with RDKit, therefore all notebooks can be run offline.

## 1. Install

```bash
conda create -n chemcheck -c conda-forge python=3.11 rdkit
conda activate chemcheck
pip install -e ".[pretty]"  # [pretty] is optional: it adds a colored terminal output
```

## 2. Your first audit (60 seconds)

The repository includes a deliberately problematic demo dataset. It contains examples of many chemcheck findings, including two invalid structures and several dataset-level problems involving otherwise valid molecules.


```bash
chemcheck tests/fixtures/demo.csv --label-col activity --split-col split
```

You get something like:

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

Two files mean train/test mode (no `--split-col` needed):

```bash
chemcheck train.csv test.csv --label-col activity
```

Custom split names can be mapped explicitly. For a fold column, selecting the
held-out fold automatically treats every other observed fold as train:

```bash
chemcheck dataset.csv --split-col fold --test-value 0
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
  Weights live in `src/chemcheck/scoring.py` and every deduction scales with
  prevalence — a handful of bad rows costs little, a systemic problem costs a lot.

A warning should not be just a count for the user: it's meant to deliver relevant information.

## 4. Splits and leakage (the important chapter)

Based on the finding, you can have 3 different verdicts. Notebook 02 describes an example of one for
each:

1. **Identity leakage** (`error`): same connectivity in train and test.
   This is never acceptable, as metrics become memorization scores. Fix: drop or move.
2. **Analog leakage** (`warning`, Morgan Tc ≥ 0.6, tunable via `--analog-thresh`):
   the test molecule is a close cousin of a training molecule. The model is
   interpolating, not generalizing. Fix: keep chemical series in one split.
3. **Scaffold overlap** (`warning`): shared Bemis–Murcko scaffold, dissimilar
   molecules. Context-dependent: fine for lead-optimization claims, fatal for
   "works on new chemotypes" claims. Match the split to the claim.

Rule of thumb: re-audit after every fix. The score should only go in one
direction (notebook 02 goes 71 → 82 → 93).

## 5. Python API cookbook

```python
from chemcheck.report import audit, render_json, render_html, render_terminal

report = audit(["dataset.csv"], label_col="activity", split_col="split")
print(report.score, f"{report.n_valid}/{report.n_total} valid")
```

**Recipe 1 — work with one finding at a time:**

```python
leaks = [f for f in report.findings if f.check_id == "analog_leakage"]
for ex in leaks[0].examples:          # train_row, test_row, tanimoto, ...
    print(ex["test_row"], "is a", ex["tanimoto"], "neighbor of", ex["train_row"])
```

**Recipe 2 — tune sensitivity:**

```python
strict = audit(["train.csv", "test.csv"], analog_thresh=0.5, near_dup_thresh=0.9)
```

For repeated numeric measurements, configure the label spread that counts as
contradictory in your target's units:

```python
assay = audit(["assay.csv"], label_col="pIC50", conflict_thresh=0.5)
```

For large datasets, inspect `report.meta["approximations"]`. It records when a
check used deterministic sampling or a bounded candidate window; the same
information is included per finding in JSON output.

**Recipe 3 — mine the JSON (diff dataset versions in a PR):**

```python
import json
info = json.loads(render_json(report))
{entry["check_id"]: entry["count"] for entry in info["findings"]}
```

**Recipe 4 — fail a pipeline on errors only:**

```python
from chemcheck.cli import main as chemcheck_main
raise SystemExit(chemcheck_main(["dataset.csv", "--fail-on", "error", "-q"]))
# exit codes: 0 clean · 1 warnings · 2 errors
```

**Recipe 5 — export the QC record for a submission:**

```python
with open("dataset_qc.html", "w") as fh:
    fh.write(render_html(report))
```

## 6. CI integration

Run Chemcheck automatically whenever the repository changes. The following GitHub Actions step fails if dataset.sdf contains any error-level finding and writes a JUnit report that CI tools can display.

```yaml
# .github/workflows/chemcheck.yml
- name: Audit dataset
  run: chemcheck data/dataset.sdf --split-col split --fail-on error --format junit -o chemcheck.xml
```

## 7. FAQ

**Is this a replacement for RDKit / MolVS / the ChEMBL curation pipeline?**
No. chemcheck is built *on top* of RDKit. Standardizers fix molecules one at a time;
chemcheck reasons across rows and splits (duplicates, leakage, shift).
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
HTML report alongside the dataset. Both record the chemcheck version
implicitly via the finding set. Pin the version in your environment.
