# Contributing to chemcheck

Thanks for stopping by — issues and pull requests are welcome.

## Setup

```bash
conda install -c conda-forge rdkit pandas scikit-learn scipy pyarrow openpyxl
pip install -e ".[dev,notebooks]"
```

## Workflow

1. Open an issue first for anything non-trivial, so design can be discussed.
2. Keep changes focused; one PR per concern.
3. Must-haves before a PR is merged:
   - `pytest` green (`python -m pytest tests -q`).
   - `ruff check` clean (config in `pyproject.toml`).
   - If you touched `src/`, add or update a test.
   - If you changed findings, scoring, or CLI output, re-execute the affected
     notebooks: `jupyter nbconvert --to notebook --execute --inplace examples/*.ipynb`
     and commit the refreshed outputs.
4. Follow the existing finding contract: every new check returns row IDs,
   evidence examples, and an actionable recommendation — never a bare count.

## Reporting bugs

Please include: `chemcheck --version` (or commit hash), a minimal CSV/SDF
reproducing the issue, the full report (`--format json`), and what you
expected instead.

## Code of Conduct

Be kind and professional. See `CODE_OF_CONDUCT.md`.
