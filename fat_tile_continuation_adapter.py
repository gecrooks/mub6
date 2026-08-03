"""Fail-closed adapter boundary from fat-tile runs to continuation v2."""

from dataclasses import dataclass

from certificate_result import CertificateGrade, CertificateResult
from continuation_artifact import ContinuationCoverageArtifact
from coverage_contract import ParentCoverageArtifact


MISSING_PARENT = (
    "complete ParentCoverageArtifact (global sweep, parameter box, zones)",
    "stable root_index and tube/fold kind for every coverage zone",
)
MISSING_CONTINUATION = (
    "regular roots: certified S[5,3] and symmetric Q[5,3,3]",
    "regular roots: residual, Jacobian-remainder, and tube-radius bounds",
    "fold roots: certified oracle frame and all split octants when used",
)


class ContinuationAdapterError(RuntimeError):
    def __init__(self, report):
        self.report = report
        super().__init__("; ".join(report.missing) or report.reason)


@dataclass(frozen=True)
class ContinuationAdapterReport:
    ready: bool
    source_ok: bool
    source_grade: CertificateGrade
    required_grade: CertificateGrade
    root_count: int | None
    missing: tuple[str, ...]
    reason: str
    artifact_id: str | None = None

    def as_dict(self):
        return {
            "ready": self.ready,
            "source_ok": self.source_ok,
            "source_grade": self.source_grade.name,
            "required_grade": self.required_grade.name,
            "root_count": self.root_count,
            "missing": list(self.missing),
            "reason": self.reason,
            "artifact_id": self.artifact_id,
        }


def _metadata(result):
    if not isinstance(result, CertificateResult):
        raise TypeError("fat-tile adapter requires CertificateResult")
    return dict(result.metadata)


def continuation_readiness(result, *, parent=None, continuations=None,
                           required_grade=CertificateGrade.RIGOROUS):
    """Report whether supplied fat-tile components can form v2.

    The source result's overall grade is reported but is not substituted for
    the artifact's own coverage and continuation grades.  A failed pair leg
    does not invalidate independently complete coverage components.
    """
    metadata = _metadata(result)
    required = CertificateGrade(required_grade)
    missing = []
    if parent is None:
        missing.extend(MISSING_PARENT)
    elif not isinstance(parent, ParentCoverageArtifact) or not parent.complete:
        missing.append("complete ParentCoverageArtifact")
    if continuations is None:
        missing.extend(MISSING_CONTINUATION)
    root_count = metadata.get("n_roots")
    if missing:
        return ContinuationAdapterReport(
            False, result.ok, result.grade, required, root_count,
            tuple(missing), "fat-tile output lacks v2 proof components",
        )
    try:
        artifact = ContinuationCoverageArtifact(parent, tuple(continuations))
    except (TypeError, ValueError, KeyError) as error:
        return ContinuationAdapterReport(
            False, result.ok, result.grade, required, root_count,
            (str(error),), "invalid v2 proof components",
        )
    if root_count is not None and int(root_count) != len(
            artifact.parent.witness.zones):
        return ContinuationAdapterReport(
            False, result.ok, result.grade, required, root_count,
            ("fat-tile root count does not match coverage zones",),
            "source/artifact root mismatch",
        )
    if artifact.grade < required:
        return ContinuationAdapterReport(
            False, result.ok, result.grade, required, root_count,
            (f"artifact grade {artifact.grade.name} is below "
             f"{required.name}",), "artifact grade rejected",
            artifact.artifact_id,
        )
    return ContinuationAdapterReport(
        True, result.ok, result.grade, required, root_count, (),
        "v2 components validated; child Krawczyk verification still required",
        artifact.artifact_id,
    )


def adapt_fat_tile_result(result, *, parent=None, continuations=None,
                          required_grade=CertificateGrade.RIGOROUS):
    """Build v2 or raise with a complete machine-readable readiness report."""
    report = continuation_readiness(
        result, parent=parent, continuations=continuations,
        required_grade=required_grade,
    )
    if not report.ready:
        raise ContinuationAdapterError(report)
    return ContinuationCoverageArtifact(parent, tuple(continuations))
