"""Persist one coverage parent and run rigorously restricted child tiles."""

import argparse
import itertools
import json
import os
from pathlib import Path
import tempfile
import time

import numpy as np

from campaign_coverage import grade_accepted
from certificate_result import CertificateGrade
from coverage_contract import ParentCoverageArtifact


class CoverageArtifactStore:
    """Content-addressed JSON store with digest verification on every load."""

    def __init__(self, directory):
        self.directory = Path(directory)

    def path_for(self, artifact_id):
        if len(artifact_id) != 64 or any(
                char not in "0123456789abcdef" for char in artifact_id):
            raise ValueError("invalid coverage artifact id")
        return self.directory / f"{artifact_id}.json"

    def save(self, artifact):
        if not isinstance(artifact, ParentCoverageArtifact):
            raise TypeError("expected ParentCoverageArtifact")
        if not artifact.complete:
            raise ValueError("cannot persist incomplete coverage artifact")
        self.directory.mkdir(parents=True, exist_ok=True)
        target = self.path_for(artifact.artifact_id)
        if target.exists():
            existing = self.load(artifact.artifact_id)
            if existing != artifact:
                raise ValueError("artifact id collision")
            return target
        handle, temporary = tempfile.mkstemp(
            prefix=f".{artifact.artifact_id}.", suffix=".tmp",
            dir=self.directory,
        )
        try:
            with os.fdopen(handle, "w") as stream:
                json.dump(artifact.as_dict(), stream, sort_keys=True)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return target

    def load(self, artifact_id):
        with self.path_for(artifact_id).open() as stream:
            value = json.load(stream)
        artifact = ParentCoverageArtifact.from_dict(value)
        if artifact.artifact_id != artifact_id:
            raise ValueError("artifact filename does not match its digest")
        return artifact


def certify_and_store_parent(beta, half_widths, store, *, certifier=None,
                             required_grade=CertificateGrade.RIGOROUS,
                             verbose=True):
    """Run one global sweep and persist its independently graded artifact."""
    if certifier is None:
        from parametric import certify_tile
        certifier = certify_tile
    run = certifier(beta, half_widths, verbose=verbose, use_certified=True)
    artifact = run.get("coverage_artifact") if isinstance(run, dict) else None
    if artifact is None:
        raise RuntimeError("parent certification produced no coverage artifact")
    grade = artifact.witness.evidence().grade
    if not artifact.complete:
        raise RuntimeError("parent coverage artifact is incomplete")
    if grade < CertificateGrade(required_grade):
        raise RuntimeError(
            f"parent coverage grade {grade.name} is below "
            f"{CertificateGrade(required_grade).name}"
        )
    path = store.save(artifact)
    return {"artifact": artifact, "path": path, "coverage_grade": grade,
            "parent_result_ok": bool(run.get("ok", False))}


def _outward_interval(center, half_width):
    return (float(np.nextafter(center - half_width, -np.inf)),
            float(np.nextafter(center + half_width, np.inf)))


def _axis_children(center, parent_half_width, child_half_width):
    values = (center, parent_half_width, child_half_width)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("box dimensions must be finite")
    if parent_half_width < 0 or child_half_width <= 0:
        raise ValueError("parent width must be nonnegative; child positive")
    if child_half_width > parent_half_width:
        raise ValueError("child half-width exceeds parent half-width")
    if child_half_width == parent_half_width:
        return (float(center),)
    count = int(np.ceil(parent_half_width / child_half_width))
    first = center - parent_half_width + child_half_width
    last = center + parent_half_width - child_half_width
    candidates = [float(value) for value in np.linspace(first, last, count)]
    parent_lo, parent_hi = _outward_interval(center, parent_half_width)
    adjusted = []
    for candidate in candidates:
        # Algebraically boundary-aligned centers can land one ulp outside
        # after serialization/arithmetic. Move inward, never outward.
        for _ in range(8):
            lo, hi = _outward_interval(candidate, child_half_width)
            if lo >= parent_lo and hi <= parent_hi:
                break
            candidate = float(np.nextafter(candidate, center))
        else:
            raise RuntimeError("cannot represent contained child interval")
        adjusted.append(candidate)
    return tuple(adjusted)


def partition_children(artifact, child_half_widths):
    """Return a validated Cartesian closed cover inside ``artifact``.

    The children may overlap.  Validation fails unless each represented
    floating-point child is contained and every parent axis is covered
    continuously from boundary to boundary.
    """
    parent = artifact.witness
    widths = tuple(float(value) for value in child_half_widths)
    if len(widths) != len(parent.parameter_center):
        raise ValueError("child dimension does not match parent")
    axes = tuple(
        _axis_children(center, parent_width, child_width)
        for center, parent_width, child_width in zip(
            parent.parameter_center, parent.parameter_half_widths, widths
        )
    )
    for axis, (centers, pc, ph, ch) in enumerate(zip(
            axes, parent.parameter_center, parent.parameter_half_widths,
            widths)):
        intervals = sorted(_outward_interval(center, ch)
                           for center in centers)
        parent_lo, parent_hi = _outward_interval(pc, ph)
        if intervals[0][0] > parent_lo or intervals[-1][1] < parent_hi:
            raise RuntimeError(f"child partition misses axis {axis} boundary")
        if any(right[0] > left[1]
               for left, right in zip(intervals, intervals[1:])):
            raise RuntimeError(f"child partition has axis {axis} gap")
    children = []
    for center in itertools.product(*axes):
        restricted = artifact.restrict(center, widths)
        if not restricted.complete:
            raise RuntimeError("generated child protrudes beyond parent")
        children.append((tuple(center), widths))
    return tuple(children)


def child_key(center, half_widths):
    return (tuple(float(value) for value in center),
            tuple(float(value) for value in half_widths))


def completed_child_keys(records, artifact_id,
                         required_grade=CertificateGrade.RIGOROUS):
    """Accepted resume keys bound to exactly one parent artifact."""
    completed = set()
    for record in records:
        if not isinstance(record, dict) \
                or record.get("mode") != "parent-coverage-child" \
                or record.get("coverage_artifact_id") != artifact_id \
                or not grade_accepted(record, required_grade):
            continue
        try:
            completed.add(child_key(record["beta"], record["hv"]))
        except (KeyError, TypeError, ValueError):
            continue
    return completed


def run_children(artifact, child_half_width, ledger, *, prior_records=(),
                 required_grade=CertificateGrade.RIGOROUS,
                 certify_child=None, verbose=True):
    """Certify isotropic children, appending artifact-bound JSONL records."""
    if certify_child is None:
        from rigor_tile import fully_rigorous_signed_tile
        certify_child = fully_rigorous_signed_tile
    h = float(child_half_width)
    children = partition_children(artifact, (h, h, h))
    done = completed_child_keys(
        prior_records, artifact.artifact_id, required_grade
    )
    written = skipped = accepted = 0
    for center, widths in children:
        if child_key(center, widths) in done:
            skipped += 1
            continue
        started = time.time()
        result = certify_child(
            center[0] - h, center[0] + h, center[1], center[2],
            hf=h, hf3=h, coverage_artifact=artifact, verbose=False,
        )
        record = {
            "mode": "parent-coverage-child",
            "coverage_artifact_id": artifact.artifact_id,
            "beta": list(center),
            "hv": list(widths),
            "ok": bool(result.ok),
            "grade": result.grade.name,
            "seconds": time.time() - started,
            "reason": result.reason,
            "result": result.as_dict(),
        }
        ledger.write(json.dumps(record, sort_keys=True) + "\n")
        ledger.flush()
        written += 1
        accepted += grade_accepted(record, required_grade)
        if verbose:
            print(f"child {center}: {'OK' if record['ok'] else 'FAIL'}"
                  f"[{record['grade']}]", flush=True)
    return {"children": len(children), "written": written,
            "skipped": skipped, "accepted": accepted,
            "artifact_id": artifact.artifact_id}


def _read_records(path):
    if not Path(path).exists():
        return []
    records = []
    with open(path) as stream:
        for line in stream:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_store")
    parser.add_argument("artifact_id")
    parser.add_argument("--child-h", type=float, required=True)
    parser.add_argument("--ledger", default="parent_children.jsonl")
    parser.add_argument("--required-grade",
                        choices=[grade.name for grade in CertificateGrade],
                        default=CertificateGrade.RIGOROUS.name)
    args = parser.parse_args()
    artifact = CoverageArtifactStore(args.artifact_store).load(
        args.artifact_id
    )
    prior = _read_records(args.ledger)
    with open(args.ledger, "a") as ledger:
        summary = run_children(
            artifact, args.child_h, ledger, prior_records=prior,
            required_grade=CertificateGrade[args.required_grade],
        )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
