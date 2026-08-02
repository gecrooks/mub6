"""Structured, fail-closed outcomes for numerical certificates.

Boolean success is not enough for this project: an experimental shortcut
must not silently flow into theorem-grade reporting.  These small value
objects make the evidence grade and its dependencies part of the result.
"""

from dataclasses import dataclass, field
from enum import IntEnum
from types import MappingProxyType


class CertificateGrade(IntEnum):
    """Strength of the weakest ingredient supporting a result."""

    EXPERIMENTAL = 0
    SAMPLED_BOUND = 1
    RIGOROUS = 2


@dataclass(frozen=True)
class Evidence:
    name: str
    grade: CertificateGrade
    detail: str = ""

    def __post_init__(self):
        if not self.name.strip():
            raise ValueError("evidence name must not be empty")
        object.__setattr__(self, "grade", CertificateGrade(self.grade))


@dataclass(frozen=True)
class CertificateResult:
    ok: bool
    grade: CertificateGrade
    dependencies: tuple[Evidence, ...]
    reason: str = ""
    metadata: object = field(default_factory=dict)

    def __post_init__(self):
        grade = CertificateGrade(self.grade)
        dependencies = tuple(self.dependencies)
        if not dependencies:
            raise ValueError("a certificate result must declare its evidence")
        weakest = min(item.grade for item in dependencies)
        if grade > weakest:
            raise ValueError(
                f"result grade {grade.name} exceeds weakest dependency "
                f"{weakest.name}"
            )
        object.__setattr__(self, "grade", grade)
        object.__setattr__(self, "dependencies", dependencies)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @classmethod
    def from_evidence(
        cls,
        ok,
        dependencies,
        *,
        reason="",
        ceiling=CertificateGrade.RIGOROUS,
        metadata=None,
    ):
        """Construct a result whose grade is derived, never asserted."""
        dependencies = tuple(dependencies)
        if not dependencies:
            raise ValueError("a certificate result must declare its evidence")
        grade = min(
            CertificateGrade(ceiling),
            *(item.grade for item in dependencies),
        )
        return cls(
            bool(ok),
            grade,
            dependencies,
            reason,
            {} if metadata is None else metadata,
        )

    def accepted_at(self, required=CertificateGrade.RIGOROUS):
        """Whether this successful result meets a caller's evidence policy."""
        return self.ok and self.grade >= CertificateGrade(required)

    def __bool__(self):
        """Preserve fail-closed compatibility with legacy Boolean callers."""
        return self.ok

    def as_dict(self):
        return {
            "ok": self.ok,
            "grade": self.grade.name,
            "reason": self.reason,
            "dependencies": [
                {
                    "name": item.name,
                    "grade": item.grade.name,
                    "detail": item.detail,
                }
                for item in self.dependencies
            ],
            "metadata": dict(self.metadata),
        }


def combine_results(name, results, *, ceiling=CertificateGrade.RIGOROUS):
    """Combine stages, inheriting failure and the weakest stage grade."""
    results = tuple(results)
    if not results:
        raise ValueError("cannot combine an empty result sequence")
    evidence = tuple(
        Evidence(name=f"{name}:{i}", grade=result.grade, detail=result.reason)
        for i, result in enumerate(results)
    )
    failed = [result.reason for result in results if not result.ok and result.reason]
    return CertificateResult.from_evidence(
        all(result.ok for result in results),
        evidence,
        reason="; ".join(failed),
        ceiling=ceiling,
        metadata={"stage_count": len(results)},
    )
