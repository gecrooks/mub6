"""Fail-closed closure contract for survivor boxes.

This module deliberately contains no root finder or clustering logic.  Those
operations may suggest which certified collection zone to try, but they do
not prove coverage.  A survivor box is closed only when the *whole* box is
componentwise inside one certified root zone, or when an independent sweep
has certified that box excluded.

All phase differences are computed on the torus.  ``zone_radii`` are proof
inputs supplied by the certificate that constructed the root structures;
this module never invents a universal localization radius.
"""

from dataclasses import dataclass
import hashlib
import json
from enum import Enum

import numpy as np

from certificate_result import CertificateGrade, Evidence


class BoxDisposition(Enum):
    """How a survivor box was closed."""

    COLLECTED = "collected"
    EXCLUDED = "excluded"
    UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class CoverageReport:
    dispositions: tuple[BoxDisposition, ...]
    root_indices: tuple[int | None, ...]

    @property
    def complete(self):
        """True only when every input box has an explicit proof path."""
        return BoxDisposition.UNRESOLVED not in self.dispositions

    @property
    def unresolved_indices(self):
        return tuple(
            i for i, status in enumerate(self.dispositions)
            if status is BoxDisposition.UNRESOLVED
        )


@dataclass(frozen=True)
class RootCoverageZone:
    """One certified destination of the global coverage sweep.

    ``guard_radii`` describe the region collected by the global sweep;
    ``collected_reach`` records what it actually collected.  Regular roots
    additionally expose the final Krawczyk ``enclosure_radii``.  Fold
    structures have no rectangular final tube and therefore use ``None``.
    ``handoff_complete`` asserts that the guard-to-structure shell was
    independently closed, rather than assuming the whole guard is a tube.
    """

    kind: str
    center: tuple[float, ...]
    guard_radii: tuple[float, ...]
    collected_reach: tuple[float, ...]
    enclosure_radii: tuple[float, ...] | None
    handoff_complete: bool
    grade: CertificateGrade

    def __post_init__(self):
        if self.kind not in {"tube", "fold", "fold-split"}:
            raise ValueError(f"unknown coverage-zone kind: {self.kind}")
        dimension = len(self.center)
        if dimension == 0:
            raise ValueError("coverage-zone dimension must be positive")
        for name in ("guard_radii", "collected_reach"):
            values = np.asarray(getattr(self, name), dtype=float)
            if values.shape != (dimension,) or not np.isfinite(values).all():
                raise ValueError(f"{name} must be finite and match the center")
            if (values < 0).any():
                raise ValueError(f"{name} must be nonnegative")
        if self.enclosure_radii is not None:
            values = np.asarray(self.enclosure_radii, dtype=float)
            if values.shape != (dimension,) or not np.isfinite(values).all():
                raise ValueError(
                    "enclosure_radii must be finite and match the center"
                )
            if (values < 0).any():
                raise ValueError("enclosure_radii must be nonnegative")
        if self.kind == "tube" and self.enclosure_radii is None:
            raise ValueError("tube zones require final enclosure radii")
        if self.kind != "tube" and self.enclosure_radii is not None:
            raise ValueError("fold zones use their certified oracle, not a tube")
        object.__setattr__(self, "grade", CertificateGrade(self.grade))

    @property
    def complete(self):
        return self.handoff_complete

    def as_dict(self):
        return {
            "kind": self.kind,
            "center": list(self.center),
            "guard_radii": list(self.guard_radii),
            "collected_reach": list(self.collected_reach),
            "enclosure_radii": (
                None if self.enclosure_radii is None
                else list(self.enclosure_radii)
            ),
            "handoff_complete": self.handoff_complete,
            "grade": self.grade.name,
        }

    @classmethod
    def from_dict(cls, value):
        return cls(
            kind=value["kind"],
            center=tuple(value["center"]),
            guard_radii=tuple(value["guard_radii"]),
            collected_reach=tuple(value["collected_reach"]),
            enclosure_radii=(
                None if value["enclosure_radii"] is None
                else tuple(value["enclosure_radii"])
            ),
            handoff_complete=bool(value["handoff_complete"]),
            grade=CertificateGrade[value["grade"]],
        )


@dataclass(frozen=True)
class SweepCoverageWitness:
    """Auditable global exclusion-to-structure coverage witness."""

    zones: tuple[RootCoverageZone, ...]
    global_sweep_complete: bool
    arithmetic_grade: CertificateGrade
    boxes_processed: int
    phantom_count: int = 0
    parameter_center: tuple[float, ...] = ()
    parameter_half_widths: tuple[float, ...] = ()

    def __post_init__(self):
        zones = tuple(self.zones)
        if not zones:
            raise ValueError("coverage witness requires at least one root zone")
        if self.boxes_processed < 0 or self.phantom_count < 0:
            raise ValueError("coverage witness counts must be nonnegative")
        pc = np.asarray(self.parameter_center, dtype=float)
        ph = np.asarray(self.parameter_half_widths, dtype=float)
        if pc.ndim != 1 or len(pc) == 0 or ph.shape != pc.shape:
            raise ValueError("coverage witness requires a parameter box")
        if not np.isfinite(pc).all() or not np.isfinite(ph).all() \
                or (ph < 0).any():
            raise ValueError("coverage witness parameter box is invalid")
        object.__setattr__(self, "zones", zones)
        object.__setattr__(
            self, "arithmetic_grade", CertificateGrade(self.arithmetic_grade)
        )

    @property
    def complete(self):
        return self.global_sweep_complete and all(zone.complete for zone in self.zones)

    def evidence(self):
        grade = min(self.arithmetic_grade, *(zone.grade for zone in self.zones))
        detail = (
            f"global sweep {'closed' if self.global_sweep_complete else 'open'}; "
            f"{len(self.zones)} root structures; "
            f"{self.boxes_processed} boxes; {self.phantom_count} phantoms"
        )
        return Evidence("enumeration-coverage", grade, detail)

    def matches(self, parameter_center, parameter_half_widths, roots,
                *, atol=1e-12):
        """Check that this witness covers the requested box and root pool."""
        center = np.asarray(parameter_center, dtype=float)
        half_widths = np.asarray(parameter_half_widths, dtype=float)
        witness_center = np.asarray(self.parameter_center, dtype=float)
        witness_widths = np.asarray(self.parameter_half_widths, dtype=float)
        if center.shape != witness_center.shape \
                or half_widths.shape != witness_widths.shape:
            return False
        if not np.allclose(center, witness_center, rtol=0.0, atol=atol):
            return False
        if not np.all(witness_widths + atol >= half_widths):
            return False
        roots = np.asarray(roots, dtype=float)
        witness_roots = np.asarray([zone.center for zone in self.zones],
                                   dtype=float)
        if roots.shape != witness_roots.shape:
            return False
        if roots.size == 0:
            return True
        delta = np.abs((witness_roots - roots + np.pi)
                       % (2.0 * np.pi) - np.pi)
        return bool(np.max(delta) <= atol)

    def as_dict(self):
        return {
            "complete": self.complete,
            "global_sweep_complete": self.global_sweep_complete,
            "arithmetic_grade": self.arithmetic_grade.name,
            "grade": self.evidence().grade.name,
            "boxes_processed": self.boxes_processed,
            "phantom_count": self.phantom_count,
            "parameter_center": list(self.parameter_center),
            "parameter_half_widths": list(self.parameter_half_widths),
            "zones": [zone.as_dict() for zone in self.zones],
        }

    def artifact(self):
        """Freeze this sweep witness as a reusable parent artifact."""
        return ParentCoverageArtifact(self)

    @classmethod
    def from_dict(cls, value):
        return cls(
            zones=tuple(RootCoverageZone.from_dict(zone)
                        for zone in value["zones"]),
            global_sweep_complete=bool(value["global_sweep_complete"]),
            arithmetic_grade=CertificateGrade[value["arithmetic_grade"]],
            boxes_processed=int(value["boxes_processed"]),
            phantom_count=int(value.get("phantom_count", 0)),
            parameter_center=tuple(value["parameter_center"]),
            parameter_half_widths=tuple(value["parameter_half_widths"]),
        )


def _closed_box_inside(child_center, child_half_widths,
                       parent_center, parent_half_widths):
    """Fail-closed containment of one ordinary (unwrapped) closed box."""
    cc = np.asarray(child_center, dtype=float)
    ch = np.asarray(child_half_widths, dtype=float)
    pc = np.asarray(parent_center, dtype=float)
    ph = np.asarray(parent_half_widths, dtype=float)
    if cc.ndim != 1 or cc.shape != ch.shape or cc.shape != pc.shape \
            or cc.shape != ph.shape or cc.size == 0:
        return False
    if not all(np.isfinite(x).all() for x in (cc, ch, pc, ph)):
        return False
    if (ch < 0).any() or (ph < 0).any():
        return False
    # Compare outward-rounded endpoints.  This admits no numerical tolerance:
    # both represented boxes are enlarged by one floating-point step before
    # containment is tested.
    child_lo = np.nextafter(cc - ch, -np.inf)
    child_hi = np.nextafter(cc + ch, np.inf)
    parent_lo = np.nextafter(pc - ph, -np.inf)
    parent_hi = np.nextafter(pc + ph, np.inf)
    return bool(np.all(child_lo >= parent_lo) and
                np.all(child_hi <= parent_hi))


@dataclass(frozen=True)
class RestrictedCoverageWitness:
    """A parent sweep witness restricted to one proved-contained child."""

    parent: SweepCoverageWitness
    child_center: tuple[float, ...]
    child_half_widths: tuple[float, ...]
    artifact_id: str

    @property
    def zones(self):
        return self.parent.zones

    @property
    def boxes_processed(self):
        return self.parent.boxes_processed

    @property
    def phantom_count(self):
        return self.parent.phantom_count

    @property
    def complete(self):
        return (self.parent.complete and _closed_box_inside(
            self.child_center, self.child_half_widths,
            self.parent.parameter_center, self.parent.parameter_half_widths
        ))

    def evidence(self):
        parent_ev = self.parent.evidence()
        detail = (f"parent artifact {self.artifact_id[:12]}; closed child "
                  f"box {'contained' if self.complete else 'not contained'}")
        return Evidence("enumeration-coverage:parent-restriction",
                        parent_ev.grade, detail)

    def matches(self, parameter_center, parameter_half_widths, roots):
        """Bind use to this exact child box and the parent's structures."""
        center = np.asarray(parameter_center, dtype=float)
        widths = np.asarray(parameter_half_widths, dtype=float)
        if not np.array_equal(center, np.asarray(self.child_center)) \
                or not np.array_equal(widths,
                                      np.asarray(self.child_half_widths)):
            return False
        roots = np.asarray(roots, dtype=float)
        parent_roots = np.asarray([zone.center for zone in self.zones],
                                  dtype=float)
        return (self.complete and roots.shape == parent_roots.shape
                and np.array_equal(roots, parent_roots))

    def as_dict(self):
        return {
            "artifact_id": self.artifact_id,
            "child_center": list(self.child_center),
            "child_half_widths": list(self.child_half_widths),
            "complete": self.complete,
        }


@dataclass(frozen=True)
class ParentCoverageArtifact:
    """Content-addressed global coverage proof reusable by child boxes."""

    witness: SweepCoverageWitness

    @property
    def artifact_id(self):
        payload = json.dumps(self.witness.as_dict(), sort_keys=True,
                             separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @property
    def complete(self):
        return self.witness.complete

    def restrict(self, child_center, child_half_widths):
        center = tuple(float(x) for x in child_center)
        widths = tuple(float(x) for x in child_half_widths)
        return RestrictedCoverageWitness(
            self.witness, center, widths, self.artifact_id
        )

    def as_dict(self):
        return {
            "schema": "mub6-parent-coverage-v1",
            "artifact_id": self.artifact_id,
            "witness": self.witness.as_dict(),
        }

    @classmethod
    def from_dict(cls, value):
        if value.get("schema") != "mub6-parent-coverage-v1":
            raise ValueError("unknown parent coverage artifact schema")
        artifact = cls(SweepCoverageWitness.from_dict(value["witness"]))
        if value.get("artifact_id") != artifact.artifact_id:
            raise ValueError("parent coverage artifact digest mismatch")
        return artifact


def _box_array(value, name, dimension=None):
    arr = np.asarray(value, dtype=float)
    if arr.size == 0:
        if arr.ndim == 2:
            if dimension is not None and arr.shape[1] != dimension:
                raise ValueError(
                    f"{name} must have dimension {dimension}"
                )
            return arr.reshape(0, arr.shape[1])
        if dimension is None:
            raise ValueError(f"cannot infer {name} dimension from empty input")
        return np.empty((0, dimension), dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.ndim != 2 or (dimension is not None and arr.shape[1] != dimension):
        suffix = "" if dimension is None else f" with dimension {dimension}"
        raise ValueError(f"{name} must be a two-dimensional box array{suffix}")
    if not np.isfinite(arr).all():
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _torus_distance(left, right):
    return np.abs((left - right + np.pi) % (2.0 * np.pi) - np.pi)


def close_survivor_boxes(
    centers,
    half_widths,
    root_centers,
    zone_radii,
    *,
    candidate_roots=None,
    excluded=None,
):
    """Classify boxes using certified inclusion or exclusion evidence.

    ``candidate_roots`` is optional advisory output from enumeration or
    polishing.  When supplied, box ``i`` may be collected only into exactly
    that root's zone; proximity to a different root cannot accidentally
    complete the proof.  Use ``-1`` for a box with no candidate.

    ``excluded`` is a Boolean mask produced by an independent certified
    exclusion sweep.  It is consulted only when collection-zone inclusion
    fails.  With neither proof, the disposition is ``UNRESOLVED``.
    """
    centers = _box_array(centers, "centers")
    dimension = centers.shape[1]
    half_widths = _box_array(half_widths, "half_widths", dimension)
    roots = _box_array(root_centers, "root_centers", dimension)
    radii = _box_array(zone_radii, "zone_radii", dimension)

    n_boxes = len(centers)
    n_roots = len(roots)
    if len(half_widths) != n_boxes:
        raise ValueError("centers and half_widths must have equal length")
    if len(radii) != n_roots:
        raise ValueError("root_centers and zone_radii must have equal length")
    if (half_widths < 0).any() or (radii < 0).any():
        raise ValueError("half-widths and zone radii must be nonnegative")

    if candidate_roots is None:
        candidates = None
    else:
        candidates = np.asarray(candidate_roots, dtype=int).reshape(-1)
        if len(candidates) != n_boxes:
            raise ValueError("candidate_roots must have one entry per box")
        if ((candidates < -1) | (candidates >= n_roots)).any():
            raise ValueError("candidate root index is out of range")

    if excluded is None:
        excluded_mask = np.zeros(n_boxes, dtype=bool)
    else:
        excluded_mask = np.asarray(excluded, dtype=bool).reshape(-1)
        if len(excluded_mask) != n_boxes:
            raise ValueError("excluded must have one entry per box")

    dispositions = []
    associations = []
    for i, (center, width) in enumerate(zip(centers, half_widths)):
        if candidates is None:
            allowed = range(n_roots)
        elif candidates[i] >= 0:
            allowed = (int(candidates[i]),)
        else:
            allowed = ()

        containing = [
            root_i for root_i in allowed
            if (_torus_distance(center, roots[root_i]) + width
                <= radii[root_i]).all()
        ]
        if containing:
            # Prefer the zone with the largest componentwise clearance.
            root_i = max(
                containing,
                key=lambda j: float(np.min(
                    radii[j] - _torus_distance(center, roots[j]) - width
                )),
            )
            dispositions.append(BoxDisposition.COLLECTED)
            associations.append(root_i)
        elif excluded_mask[i]:
            dispositions.append(BoxDisposition.EXCLUDED)
            associations.append(None)
        else:
            dispositions.append(BoxDisposition.UNRESOLVED)
            associations.append(None)

    return CoverageReport(tuple(dispositions), tuple(associations))
