"""Content-addressed tile artifacts and rigorous ledger binding checks."""

from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile

from ball_coverage_artifact import BallCoverageArtifact, SCHEMA as BALL_SCHEMA
from certificate_result import CertificateGrade
from ledger_bits import decode_box_record, float_bits


@dataclass(frozen=True)
class ArtifactBindingReport:
    accepted: bool
    reason: str
    artifact_schema: str | None
    artifact_id: str | None


class TileArtifactStore:
    """Immutable JSON objects keyed by their verified SHA-256 artifact ID."""

    def __init__(self, directory):
        self.directory = Path(directory)

    def path_for(self, artifact_id):
        if not isinstance(artifact_id, str) or len(artifact_id) != 64 \
                or any(char not in "0123456789abcdef"
                       for char in artifact_id):
            raise ValueError("artifact id must be 64 lowercase hex digits")
        return self.directory / f"{artifact_id}.json"

    def save(self, artifact):
        value = artifact.as_dict()
        artifact_id = artifact.artifact_id
        value = {**value, "artifact_id": artifact_id}
        self.directory.mkdir(parents=True, exist_ok=True)
        target = self.path_for(artifact_id)
        if target.exists():
            loaded = self.load_raw(artifact_id)
            if loaded != value:
                raise ValueError("artifact id collision")
            return target
        handle, temporary = tempfile.mkstemp(
            prefix=f".{artifact_id}.", suffix=".tmp", dir=self.directory
        )
        try:
            with os.fdopen(handle, "w") as stream:
                json.dump(value, stream, sort_keys=True,
                          separators=(",", ":"), allow_nan=False)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return target

    def load_raw(self, artifact_id):
        with self.path_for(artifact_id).open() as stream:
            return json.load(stream)


def bind_tile_artifact(record, artifact):
    """Bind a ledger record to a full tile artifact after exact box check."""
    beta, widths, _interval, _token = decode_box_record(record)
    if not isinstance(artifact, BallCoverageArtifact):
        raise TypeError("only full ball-coverage artifacts certify tiles")
    if tuple(float_bits(x) for x in beta) != tuple(
            float_bits(x) for x in artifact.parameter_center):
        raise ValueError("ledger beta does not match artifact parameter box")
    if tuple(float_bits(x) for x in widths) != tuple(
            float_bits(x) for x in artifact.parameter_half_widths):
        raise ValueError("ledger hv does not match artifact parameter box")
    return {**record, "artifact_schema": BALL_SCHEMA,
            "artifact_id": artifact.artifact_id}


def verify_tile_binding(record, store,
                        required_grade=CertificateGrade.RIGOROUS):
    schema = record.get("artifact_schema") if isinstance(record, dict) else None
    artifact_id = record.get("artifact_id") if isinstance(record, dict) else None
    if schema is None or artifact_id is None:
        return ArtifactBindingReport(False, "missing tile artifact reference",
                                     schema, artifact_id)
    if schema != BALL_SCHEMA:
        return ArtifactBindingReport(False, "schema is not a full tile artifact",
                                     schema, artifact_id)
    if store is None:
        return ArtifactBindingReport(False, "artifact store unavailable",
                                     schema, artifact_id)
    try:
        value = store.load_raw(artifact_id)
        artifact = BallCoverageArtifact.from_dict(value)
        if artifact.artifact_id != artifact_id:
            raise ValueError("artifact digest mismatch")
        if not artifact.complete:
            raise ValueError("artifact is incomplete")
        if artifact.grade < CertificateGrade(required_grade):
            raise ValueError(f"artifact grade {artifact.grade.name} is below "
                             f"{CertificateGrade(required_grade).name}")
        failures = artifact.check()
        if failures:
            raise ValueError("artifact kernel check failed: " + "; ".join(failures))
        bind_tile_artifact(record, artifact)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        return ArtifactBindingReport(False, str(error), schema, artifact_id)
    return ArtifactBindingReport(True, "artifact verified", schema, artifact_id)


def verified_resume_records(records, store,
                            required_grade=CertificateGrade.RIGOROUS):
    """Return verified records plus one binding report per input record."""
    verified = []
    reports = []
    for record in records:
        report = verify_tile_binding(record, store, required_grade)
        reports.append(report)
        if report.accepted:
            verified.append(record)
    return verified, tuple(reports)
