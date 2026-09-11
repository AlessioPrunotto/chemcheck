"""Shared dataclasses for findings and reports."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class Severity(str, Enum):
    """Severity level for a finding.

    Attributes:
        ERROR: Blocking issue that should fail an audit.
        WARNING: Suspicious issue worth reviewing.
        INFO: Informational note, not a failure on its own.
    """

    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass
class Finding:
    """A single audit finding produced by a check.

    Attributes:
        check_id: Stable machine-readable check identifier.
        severity: Severity level of the finding.
        title: Short human-readable title.
        count: Number of affected rows (after truncation logic).
        rate: Fraction of total rows affected (`count / n_total`).
        affected_rows: Truncated row ids affected by this finding.
        total_affected: Full count before truncation.
        examples: Example payloads illustrating the issue.
        recommendation: Suggested remediation for the user.
        details: Extra details or statistics about the finding.
    """

    check_id: str
    severity: Severity
    title: str
    count: int
    rate: float = 0.0  # count / n_total
    affected_rows: list[str] = field(default_factory=list)  # row ids (truncated)
    total_affected: int = 0  # full count before truncation
    examples: list[dict[str, Any]] = field(default_factory=list)
    recommendation: str = ""
    details: str = ""


@dataclass
class CheckResult:
    """Result wrapper for a single check.

    Attributes:
        check_id: Stable machine-readable check identifier.
        finding: The finding, or None if the check passed.
    """

    check_id: str
    finding: Finding | None = None  # None => passed, no finding


@dataclass
class AuditReport:
    """Aggregated result of auditing a molecular dataset.

    Attributes:
        n_total: Total number of input rows.
        n_valid: Number of rows with valid chemistry.
        n_invalid: Number of rows with invalid chemistry.
        findings: All findings, ordered errors/warnings/infos.
        score: Dataset quality score from 0 to 100.
        score_breakdown: Per-check score deductions.
        meta: Free-form metadata (source, columns, timestamp).
    """

    n_total: int
    n_valid: int
    n_invalid: int
    findings: list[Finding]
    score: int
    score_breakdown: list[dict[str, Any]]
    meta: dict[str, Any] = field(default_factory=dict)

    def errors(self) -> list[Finding]:
        """Return findings with error severity.

        Returns:
            List of error findings.
        """
        return [f for f in self.findings if f.severity == Severity.ERROR]

    def warnings(self) -> list[Finding]:
        """Return findings with warning severity.

        Returns:
            List of warning findings.
        """
        return [f for f in self.findings if f.severity == Severity.WARNING]
