# chemcheck tutorial — pytest for molecular datasets

This tutorial takes you from zero to a CI-gated dataset audit in ~20 minutes.
It is CLI-first and skimmable; for the interactive versions with plots and
experiments, see the notebooks in [`examples/`](../examples/):

| notebook | what you learn | time |
|---|---|---|
| `01_quickstart.ipynb` | CLI + Python API on a dirty dataset | 10 min |
| `02_leakage_splits.ipynb` | identity vs analog vs scaffold leakage, and how to fix each | 15 min |
| `03_curation_case_study.ipynb` | full curation loop: dirty file → 22 → 80, HTML report included | 20 min |

All notebook data is generated inline with RDKit — everything here runs offline.

## 1. Install

```bash
conda install -c conda-forge rdkit pandas scikit-learn scipy pyarrow openpyxl
pip install -e .          # from the repo root
pip install rich          # optional: colored terminal output
```

RDKit is a hard dependency. If it is missing, chemcheck fails fast with an
install hint instead of a cryptic traceback.

## 2. Your first audit (60 seconds)

The repo ships a deliberately dirty demo dataset — every row is a lesson:

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
Dataset quality score: 22/100
  top deductions: cross_split_near_duplicates -8.61, identity_leakage -7.41, ...
```

Two files mean train/test mode (no `--split-col` needed):

```bash
chemcheck train.csv test.csv --label-col activity
```

## 3. Understanding the output

Every finding has the same anatomy:

- **severity** — `error` (must fix: invalid chemistry, leakage, contradictions),
  `warning` (probably fix: duplicates, unspecified stereo, easy splits),
  `info` (know about it: isotopes, tautomers, domain gaps).
- **check_id** — stable machine-readable name (`identity_leakage`,
  `conflicting_measurements`, …). Gate on these in CI, not on prose.
- **count / rate + affected_rows** — how many rows, and *which* rows (your IDs).
- **examples** — the evidence: molecule pairs, Tanimoto scores, shared scaffolds.
- **recommendation** — what to do about it.
- **score** — 0–100 with a printed deduction breakdown, so any score is auditable.
  Weights live in `src/chemcheck/scoring.py` and every deduction scales with
  prevalence — a handful of bad rows costs little, a systemic problem costs a lot.

A warning is never just a count. That is the whole point of the tool.

## 4. Splits and leakage (the important chapter)

Three different findings, three different verdicts — notebook 02 plants one of
each so you can feel the difference:

1. **Identity leakage** (`error`) — same connectivity in train and test.
   Never acceptable; metrics become memorization scores. Fix: drop or move.
2. **Analog leakage** (`warning`, Morgan Tc ≥ 0.6, tunable via `--analog-thresh`) —
   the test molecule is a close cousin of a training molecule. Your model is
   interpolating, not generalizing. Fix: keep chemical series in one split.
3. **Scaffold overlap** (`warning`) — shared Bemis–Murcko scaffold, dissimilar
   molecules. Context-dependent: fine for lead-optimization claims, fatal for
   "works on new chemotypes" claims. Match the split to the claim.

Rule of thumb: re-audit after every fix. The score should only go in one
direction (notebook 02 goes 55 → 82 → 93).

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

Gate datasets the way pytest gates code. `dataset.sdf` failing on any error,
with JUnit output for the build dashboard:

```yaml
# .github/workflows/chemcheck.yml
- name: Audit dataset
  run: chemcheck data/dataset.sdf --split-col split --fail-on error --format junit -o chemcheck.xml
```

## 7. FAQ

**Is this a replacement for RDKit / MolVS / the ChEMBL curation pipeline?**
No — it is built *on top* of RDKit. Standardizers fix molecules one at a time;
chemcheck reasons across rows and splits (duplicates, leakage, shift).
Use both: standardize first, then audit.

**How big a dataset can it handle?**
~100k molecules on a laptop. Exact checks are linear; pairwise similarity uses
blocking plus deterministic sampling above a few thousand rows, and reports
counts as lower bounds. Thresholds are tunable.

**My score is low. Where do I start?**
Errors first, in this order: `invalid_smiles` → `identity_leakage` /
`cross_split_near_duplicates` → `conflicting_measurements`. Then warnings by
deduction size (`score_breakdown` tells you). Notebook 03 walks a 22 → 80
curation end to end.

**A standardizer corrupted my molecule — now what?**
That happened to us too: RDKit's default salt stripper deletes neutral acetic
acid (carboxylic acids match its salt list). Defensive rule from notebook 03:
check every transformation's output, keep the input on failure, and let the
re-audit confirm the fix.

**How do I cite / reference the report?**
Export `--format json` (exact counts, row IDs, thresholds in one file) or the
HTML report alongside the dataset. Both record the chemcheck version
implicitly via the finding set — pin the version in your environment.
