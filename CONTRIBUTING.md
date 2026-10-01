# Contributing to chemcheck

Thanks for stopping by — issues and pull requests are welcome.

## Setup

Chemcheck requires Python 3.10 or newer and is tested with Python 3.10–3.12.
From a cloned copy of the repository, create an isolated environment and
install chemcheck in editable mode with the testing and notebook dependencies:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,notebooks]"
```

On Windows, use the environment creation and activation commands in the
[README](README.md#install), then run the same editable pip installation.

This command also installs the core dependencies declared in `pyproject.toml`,
including RDKit, pandas, SciPy, PyArrow, and openpyxl. The editable installation
means that changes to the source code take effect without reinstalling the
package.

If a compatible RDKit wheel is not available for your Python version or
platform, use a Conda environment instead:

```bash
conda create -n chemcheck-dev -c conda-forge python=3.11 rdkit pip
conda activate chemcheck-dev
python -m pip install -e ".[dev,notebooks]"
```

## Workflow

1. Open an issue first for anything non-trivial, so design can be discussed.
2. Keep changes focused; one PR per concern.
3. Must-haves before a PR is merged:
   - Tests and branch coverage green
     (`python -m pytest -q --cov=chemcheck --cov-report=term-missing`).
   - `ruff check` clean (config in `pyproject.toml`).
   - If you touched `src/`, add or update a test.
   - If you changed findings, scoring, or CLI output, re-execute the affected
     notebooks: `jupyter nbconvert --to notebook --execute --inplace examples/*.ipynb`
     and commit the refreshed outputs.
   - If you changed check semantics, default thresholds, or a large-dataset
     execution path, rerun `python validation/run_validation.py --download`
     and commit the refreshed result tables and interpretation. Follow the
     process and snapshot rules in `validation/MAINTAINERS.md`.
4. Follow the existing finding contract: every new check returns row IDs,
   evidence examples, and an actionable recommendation — never a bare count.

## Reporting bugs

Please include: `chemcheck --version` (or commit hash), a minimal CSV/SDF
reproducing the issue, the full report (`--format json`), and what you
expected instead.

## Code of Conduct

Be kind and professional. See `CODE_OF_CONDUCT.md`.
