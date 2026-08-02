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
