# Release process

ChemDataCheck publishes without long-lived API tokens. GitHub Actions obtains a
short-lived OpenID Connect credential for each upload and the package index
accepts it only for the configured repository, workflow, and environment.
This follows the
[PyPA trusted-publishing guide](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
and [PyPI's trusted-publisher security model](https://docs.pypi.org/trusted-publishers/security-model/).

## One-time trusted-publisher setup

Create GitHub environments named `testpypi` and `pypi`. Require manual approval
for the `pypi` environment.

Register a trusted publisher in both package indexes with these exact values.
If the project does not exist on an index yet, use that index's pending-publisher
form instead:

| Setting | Value |
|---|---|
| PyPI project | `chemdatacheck` |
| GitHub owner | `AlessioPrunotto` |
| GitHub repository | `chemdatacheck` |
| Workflow | `publish.yml` |
| TestPyPI environment | `testpypi` |
| PyPI environment | `pypi` |

TestPyPI uses a separate account and publisher configuration from PyPI. No
`PYPI_API_TOKEN` or `TEST_PYPI_API_TOKEN` repository secret is required.

## TestPyPI verification

1. Set a version that has not already been uploaded with different content.
2. Run **Publish Python distribution** manually from GitHub Actions.
3. The workflow builds one wheel and source archive in an unprivileged job,
   publishes them to TestPyPI, downloads both files back, verifies that they are
   byte-for-byte identical to the locally built artifacts, installs the wheel
   with normal PyPI dependencies, and runs version smoke tests.

TestPyPI does not permit replacing a file for an existing version. The workflow
allows an existing TestPyPI upload only when all downloaded artifacts are
identical; otherwise the byte comparison fails.

## Production release

1. Update `pyproject.toml`, `src/chemdatacheck/__init__.py`, `CITATION.cff`, and the
   changelog to the same version.
2. Run the full local checks and confirm the main CI workflow is green.
3. Tag the release as `vX.Y.Z` and create a GitHub release from that tag. The
   workflow rejects a tag that does not match `project.version`.
4. The exact artifacts are published and verified on TestPyPI first.
5. After approval in the protected `pypi` environment, those verified artifacts
   are uploaded to PyPI with a signed publication attestation.

The production job runs only for a published GitHub release; a manual workflow
run can publish to TestPyPI but never to PyPI.
