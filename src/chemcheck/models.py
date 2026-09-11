"""Shared dataclasses for findings and reports."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass
class Finding:
    check_id: str
    severity: Severity
    title: str
    count: int
    rate: float = 0.0  # count / n_total
    affected_rows: list = field(default_factory=list)  # row ids (truncated)
    total_affected: int = 0  # full count before truncation
    examples: list[dict] = field(default_factory=list)
    recommendation: str = ""
    details: str = ""


@dataclass
class CheckResult:
    check_id: str
    finding: Finding | None = None  # None => passed, no finding


@dataclass
class AuditReport:
    n_total: int
    n_valid: int
    n_invalid: int
    findings: list[Finding]
    score: int
    score_breakdown: list[dict]
    meta: dict = field(default_factory=dict)

    def errors(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == Severity.ERROR]

    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == Severity.WARNING]
