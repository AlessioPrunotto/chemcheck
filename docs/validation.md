# Scientific validation

This page is for chemcheck users who want to understand what the checks have
been tested on, how to interpret their results, and where the evidence stops.
Instructions and decision records for reproducing or updating the study are in
the [maintainer notes](../validation/MAINTAINERS.md).

The study evaluates chemcheck as an **audit and triage tool**. It does not claim
that every finding is a chemical error, or that the quality score is a validated
scientific endpoint. The aim is to show that representative known defects can
be detected, observe the checks on real public data, measure threshold
sensitivity, and characterize runtime and large-dataset behavior.

The committed results were generated on 2026-09-24 with Python 3.10.19, RDKit
2025.09.4, and macOS 15.7.4 on arm64. Dataset URLs, SHA-256 hashes, and complete
machine-readable results are in [`../validation/results/`](../validation/results/).

## What was tested

Six public datasets distributed through the DeepChem/MoleculeNet object store
and two fixed chemical-database samples were used. MoleculeNet describes the
provenance and intended tasks of its benchmarks
([Wu et al., 2018](https://doi.org/10.1039/C7SC02664A)). ESOL is the
aqueous-solubility dataset introduced by
[Delaney (2004)](https://pubmed.ncbi.nlm.nih.gov/15154768/), and FreeSolv is a
curated hydration-free-energy collection
([Mobley and Guthrie, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4113415/)).

The database cohorts contain 5,000 structures each. They are fixed,
checksum-pinned samples from the official
[ChEMBL data service](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services)
and [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), so the
published results can be reproduced even if the live databases later change.
The exact identifiers are in
[`../validation/samples/`](../validation/samples/). Sampling and refresh rules
are documented separately for maintainers.

| Dataset | Rows | Task used | Split used |
|---|---:|---|---|
| FreeSolv | 642 | hydration free energy | deterministic 80/20 row holdout |
| ESOL/Delaney | 1,128 | aqueous solubility | deterministic 80/20 row holdout |
| BACE | 1,513 | binary inhibition class | supplied `Model` train/test/valid values |
| Lipophilicity | 4,200 | experimental logD | deterministic 80/20 row holdout |
| Tox21 | 7,831 | NR-AR binary endpoint | deterministic 80/20 row holdout |
| HIV | 41,127 | binary activity | deterministic 80/20 row holdout |
| ChEMBL 37 sample | 5,000 | structure-only audit | deterministic 80/20 row holdout |
| PubChem sample | 5,000 | structure-only audit | deterministic 80/20 row holdout |

The deterministic row holdouts deliberately emulate a common random-split
workflow. They are useful for diagnosing overlap, but they do not constitute
proper benchmark splits. BACE is also unusual in this snapshot: 203 rows are
marked train and 1,310 test/valid, so its rates should not be compared directly
with the 80/20 datasets.

The study asks four practical questions:

1. What does the complete default audit report on unmodified public data?
2. Does the intended check find each of nine deliberately inserted defects?
3. How sensitive are train/test similarity results to thresholds between 0.40
   and 0.99?
4. Do the faster large-dataset paths agree with exact computation?

## Results

### Known-defect sensitivity

All 9/9 inserted cases were detected: invalid SMILES, exact and canonical
duplicates, inconsistent salt forms, stereochemical and tautomer collisions,
repeated and conflicting measurements, and identity leakage across a split.
This is a focused sensitivity test, hence not a claim of 100% sensitivity across all
30 checks or all chemical representations.

### Findings on unmodified public data

| Dataset | Invalid | Identity leakage | Cross-split Tc ≥ 0.95 | Test Tc ≥ 0.60 | Scaffold overlap |
|---|---:|---:|---:|---:|---:|
| FreeSolv | 0 | 2 | 14 | 51 | 127 |
| ESOL | 0 | 6 | 26 | 87 | 176 |
| BACE | 0 | 9 | 24 | 819 | 498 |
| Lipophilicity | 0 | 29 | 75 | 526 | 447 |
| Tox21 | 8 | 36 | 205 | 726 | 1,197 |
| HIV | 7 | 0 | 155 sampled | 1,606 sampled | 5,214 |
| ChEMBL 37 sample | 0 | 5 | 49 | 327 | 316 |
| PubChem sample | 0 | 0 | 0 | 10 | 250 |

The HIV similarity counts are based on a deterministic sample of 3,000 of
8,225 test rows and are labeled approximate in the report. Identity and
scaffold checks still use all eligible records.

For ChEMBL and PubChem, the train/test column was created only to ask what a
naive row holdout of each sample would look like. The reported overlap is not
"leakage in ChEMBL" or "leakage in PubChem." Likewise, zero invalid structures
is expected from standardized Molecule/Compound endpoints and says nothing
about rejected or unstandardized depositor records. PubChem Compound represents
unique standardized structures, so its lack of duplicate findings is also not
evidence that PubChem Substance records are duplicate-free.

The structure-only cohorts also exposed patterns that are uncommon or invisible
in the small ML benchmarks:

| Dataset | Unspecified stereo | Components | Radicals | Rare elements | Applicability gaps |
|---|---:|---:|---:|---:|---:|
| ChEMBL 37 sample | 1,096 | 196 | 2 | 62 | 790 |
| PubChem sample | 1,649 | 220 | 40 | 133 | 1,391 |

These results demonstrate why the output must be interpreted as evidence, not
as an automatic reject list. Even established benchmark datasets contain
invalid structures, repeated chemistry, salts, and structurally easy random
splits. Conversely, a high count can describe the intended domain rather than
bad curation.

The six ML-benchmark quality scores ranged from 63 to 74; the structure-only
ChEMBL and PubChem samples scored 77 and 87. Those higher values are not evidence
that the databases are intrinsically cleaner: label-dependent checks do not
apply to them, and both APIs expose standardized structures. The score is a
prioritization heuristic, not a scientific ranking or publication gate.

### How to interpret the default similarity thresholds

Chemcheck already applies these defaults; users do not need to set them for a
first audit:

- `0.60` is the default warning threshold for a close train/test analogue.
- `0.95` is the default review threshold for a possible near duplicate.

These are review thresholds, not universal chemical laws. A flagged pair is a
reason to inspect the compounds and the purpose of the split, not automatic
proof that the dataset is wrong.

| Dataset | Median maximum test→train Tc | Fraction at Tc ≥ 0.60 | Fraction at Tc ≥ 0.95 |
|---|---:|---:|---:|
| FreeSolv | 0.556 | 39.5% | 5.4% |
| ESOL | 0.556 | 38.5% | 6.2% |
| BACE | 0.659 | 62.5% | 0.9% |
| Lipophilicity | 0.694 | 62.6% | 4.5% |
| Tox21 | 0.577 | 46.4% | 7.2% |

The 0.60 analog threshold is therefore a deliberately sensitive split-difficulty
alert: it flags roughly 38–63% of held-out molecules in these workflows. It
must not be read as a universal boundary between valid and leaked chemistry.
The 0.95 near-duplicate threshold is much more selective, flagging roughly
1–7% here, but it is still fingerprint- and dataset-dependent. Published work
likewise finds that no generally applicable similarity threshold reliably
implies a biological relationship
([Vogt and Bajorath, 2017](https://doi.org/10.1002/minf.201600131)) and that the
significance of a Tanimoto score depends on the representation, query, and
database size
([Baldi and Nasr, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2914517/)).

Practical guidance:

- Start with the defaults and inspect representative flagged pairs.
- Read `0.60` as an early warning that the test set may be unusually similar to
  the training set.
- Read `0.95` as a strong prompt for review, not as an identity test.
- For consequential studies, report the full nearest-neighbor distribution and
  check whether the conclusion changes at nearby thresholds.
- Override the defaults with `--analog-thresh` or `--near-dup-thresh` only when
  the project has a documented reason to use a different fingerprint,
  molecular representation, or chemical-series policy. Do not tune a threshold
  merely to make a finished split look clean.

Scaffold separation is not sufficient by itself. Recent large evaluations show
that scaffold splits can still leave structurally similar molecules across the
boundary and can misrepresent out-of-distribution performance
([Guo et al., 2024](https://arxiv.org/abs/2406.00873);
[Kretschmer et al., 2025](https://doi.org/10.1038/s41467-024-55462-w)).

### Reliability on large datasets

Some expensive checks use bounded searches or sampling to remain practical on
large collections. On all 4,200 Lipophilicity records, the scalable
near-duplicate path found the same 151 Tc ≥ 0.95 pairs as exact all-pairs
computation: precision and recall were both 1.00 in this test. This single
dataset does not prove that the shortcut can never miss a pair, so reports
still identify when this path is used.

The applicability-domain optimization also matched exact computation after
candidate verification: all 200 exact gaps were returned, with no extra
findings. The implementation history and associated regression requirements
are recorded in the maintainer notes.

### Runtime

| Dataset | Rows | Complete audit | Process high-water RSS by completion |
|---|---:|---:|---:|
| FreeSolv | 642 | 3.9 s | 187 MB |
| ESOL | 1,128 | 3.0 s | 196 MB |
| BACE | 1,513 | 10.8 s | 273 MB |
| Lipophilicity | 4,200 | 53.9 s | 273 MB |
| Tox21 | 7,831 | 42.3 s | 273 MB |
| HIV | 41,127 | 608.2 s | 950 MB |
| ChEMBL 37 sample | 5,000 | 70.7 s | 950 MB |
| PubChem sample | 5,000 | 44.5 s | 950 MB |

Runtime is not monotonic in row count because molecular complexity, tautomer
enumeration, and triggered checks differ. Memory values are process-wide
high-water marks from a sequential run, not isolated per-dataset peaks.

The practical conclusion is that full audits through roughly 8,000 molecules
are interactive on this machine. A 41,000-molecule audit is a batch job (about
10.1 minutes here). For CI, use a representative
subset on every change and schedule the full collection audit separately.

## Expected false positives and context-dependent alerts

“False positive” here means “the pattern is real, but interpreting it as a
dataset defect would be wrong.” Chemcheck intentionally reports several such
review signals.

| Finding | Legitimate explanation | Appropriate follow-up |
|---|---|---|
| Exact/canonical duplicate | independent experimental replicates or different assay conditions | retain provenance; aggregate only when conditions are comparable |
| Conflicting measurement | real assay variability, different protocol, pH, temperature, endpoint, or units | reconcile metadata before dropping either value |
| Stereochemical collision | enantiomers/diastereomers are intentionally distinct | confirm that labels and identifiers preserve stereo |
| Salt/disconnected components | tested formulation, ion pair, mixture, or organometallic is intentional | apply a parent policy only if it matches the scientific question |
| Tautomer ambiguity | multiple theoretical tautomers are normal and do not prove inconsistent input | review only when representation affects deduplication or modeling |
| Unspecified stereo | racemate, unknown stereochemistry, or deliberately stereo-agnostic assay | annotate rather than invent a configuration |
| Analog/scaffold overlap | random split or interpolation performance is the intended experiment | describe the split honestly; do not call the dataset corrupt |
| Target shift | deliberate temporal, prospective, or extrapolation test set | keep it, but report it as distribution shift |
| Label outlier | a genuine extreme or assay floor/ceiling | verify units and provenance; avoid automatic deletion |
| Rare/PAINS-like group | valid active chemistry or intentional library coverage | use orthogonal assay evidence; a substructure alert is not proof |
| Applicability gap | intentionally difficult probe of extrapolation | retain as a separate evaluation slice |
| Representation bias | focused chemical series or an acyclic-heavy task is intentional | report domain limits and per-series performance |

The PAINS-like patterns in chemcheck are especially narrow heuristic alerts.
The literature explicitly warns that filters can label genuine compounds as
artifacts and should prompt experimental follow-up rather than automatic
exclusion
([Aldrich et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5364449/);
[Capuzzi et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411023/)).

## What this validation does not prove

- ChEMBL Molecule and PubChem Compound expose standardized database structures,
  not raw depositor records. These samples broaden structural diversity but do
  not validate assay reconciliation, depositor-level standardization conflicts,
  or the full messiness of ChEMBL activities and PubChem Substances.
- The 5,000-record database cohorts are fixed coverage samples, not probability
  samples. ChEMBL uses evenly spaced positions in one ordering; PubChem samples
  a bounded CID range and omits unavailable CIDs. Their finding rates must not
  be extrapolated to either complete database.
- The database cohorts have no coherent endpoint label. They validate
  structure-dependent checks and runtime only; target, measurement-conflict,
  and label-shift checks remain evaluated on the ML datasets.
- The seeded study covers nine high-priority checks, not every check or every
  possible molecular notation.
- There is no independent expert adjudication set labeling each natural finding
  as correct, benign, or erroneous. Natural finding counts measure prevalence,
  not precision.
- Only one fingerprint definition (2048-bit radius-2 Morgan) and one RDKit
  version were studied.
- Approximate near-duplicate fidelity was tested on one 4,200-record dataset.
  Larger and chemically different collections may expose misses.
- Deterministic row splits are not repeated random seeds. Threshold fractions
  should be accompanied by split uncertainty in a model-comparison study.
- Runtime is one machine/run and is a reference measurement, not a service-level
  guarantee.
- The score weights and severity levels remain expert heuristics rather than a
  calibrated measure of dataset quality.

Chemcheck should therefore be described as an evidence-producing audit tool,
not an automated curator. A domain expert must still decide whether a finding
is an error, an intentional feature of the dataset, or uncertain.
