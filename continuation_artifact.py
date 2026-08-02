"""Versioned, fail-closed schema for parent-to-child root continuation.

This module validates and transports proof payloads; it does not manufacture
their bounds.  In particular, predicted child seeds are only inputs to a
child Krawczyk verification and are never themselves called certified roots.
"""

from dataclasses import dataclass
import hashlib
import json

import numpy as np

from certificate_result import CertificateGrade
from coverage_contract import ParentCoverageArtifact


def _finite_tuple(value, shape, name):
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def _nonnegative_tuple(value, shape, name):
    array = _finite_tuple(value, shape, name)
    if (array < 0).any():
        raise ValueError(f"{name} must be nonnegative")
    return array


@dataclass(frozen=True)
class RegularContinuationPayload:
    root_index: int
    seed: tuple[float, ...]
    linear: tuple[tuple[float, ...], ...]
    quadratic: tuple[tuple[tuple[float, ...], ...], ...]
    residual_bound: float
    jacobian_remainder: tuple[float, ...]
    tube_radii: tuple[float, ...]
    grade: CertificateGrade

    kind = "regular"

    def __post_init__(self):
        if self.root_index < 0:
            raise ValueError("root_index must be nonnegative")
        _finite_tuple(self.seed, (5,), "seed")
        _finite_tuple(self.linear, (5, 3), "linear")
        quadratic = _finite_tuple(self.quadratic, (5, 3, 3), "quadratic")
        if not np.array_equal(quadratic, quadratic.transpose(0, 2, 1)):
            raise ValueError("quadratic matrices must be symmetric")
        if not np.isfinite(self.residual_bound) or self.residual_bound < 0:
            raise ValueError("residual_bound must be finite and nonnegative")
        _nonnegative_tuple(self.jacobian_remainder, (5,),
                           "jacobian_remainder")
        _nonnegative_tuple(self.tube_radii, (5,), "tube_radii")
        object.__setattr__(self, "grade", CertificateGrade(self.grade))

    def predicted_seed(self, child_center, parent_center):
        delta = np.asarray(child_center, dtype=float) \
            - np.asarray(parent_center, dtype=float)
        _finite_tuple(delta, (3,), "child displacement")
        seed = np.asarray(self.seed, dtype=float)
        linear = np.asarray(self.linear, dtype=float)
        quadratic = np.asarray(self.quadratic, dtype=float)
        predicted = seed + linear @ delta \
            + 0.5 * np.einsum("ikl,k,l->i", quadratic, delta, delta)
        return tuple(float(value) for value in predicted)

    def as_dict(self):
        return {
            "kind": self.kind,
            "root_index": self.root_index,
            "seed": list(self.seed),
            "linear": [list(row) for row in self.linear],
            "quadratic": [[list(row) for row in plane]
                          for plane in self.quadratic],
            "residual_bound": self.residual_bound,
            "jacobian_remainder": list(self.jacobian_remainder),
            "tube_radii": list(self.tube_radii),
            "grade": self.grade.name,
        }

    @classmethod
    def from_dict(cls, value):
        return cls(
            root_index=int(value["root_index"]),
            seed=tuple(value["seed"]),
            linear=tuple(tuple(row) for row in value["linear"]),
            quadratic=tuple(tuple(tuple(row) for row in plane)
                            for plane in value["quadratic"]),
            residual_bound=float(value["residual_bound"]),
            jacobian_remainder=tuple(value["jacobian_remainder"]),
            tube_radii=tuple(value["tube_radii"]),
            grade=CertificateGrade[value["grade"]],
        )


@dataclass(frozen=True)
class FoldOraclePayload:
    octant: tuple[int, int, int] | None
    seed: tuple[float, ...]
    tangent: tuple[float, ...]
    transverse_basis: tuple[tuple[float, ...], ...]
    t_grid: tuple[float, ...]
    center_curve: tuple[tuple[float, ...], ...]
    transverse_radii: tuple[float, ...]
    half_length: float
    grade: CertificateGrade

    def __post_init__(self):
        if self.octant is not None and (len(self.octant) != 3 or
                                       any(x not in {-1, 1}
                                           for x in self.octant)):
            raise ValueError("octant must contain three signs")
        _finite_tuple(self.seed, (5,), "fold seed")
        _finite_tuple(self.tangent, (5,), "fold tangent")
        _finite_tuple(self.transverse_basis, (5, 4), "transverse_basis")
        grid = np.asarray(self.t_grid, dtype=float)
        curve = np.asarray(self.center_curve, dtype=float)
        if grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all():
            raise ValueError("t_grid must contain at least two finite points")
        if np.any(np.diff(grid) <= 0):
            raise ValueError("t_grid must be strictly increasing")
        if curve.shape != (len(grid), 4) or not np.isfinite(curve).all():
            raise ValueError("center_curve must have shape (len(t_grid), 4)")
        _nonnegative_tuple(self.transverse_radii, (4,),
                           "transverse_radii")
        if not np.isfinite(self.half_length) or self.half_length <= 0:
            raise ValueError("half_length must be finite and positive")
        if grid[0] > -self.half_length or grid[-1] < self.half_length:
            raise ValueError("t_grid does not span the fold half-length")
        object.__setattr__(self, "grade", CertificateGrade(self.grade))

    def as_dict(self):
        return {
            "octant": None if self.octant is None else list(self.octant),
            "seed": list(self.seed),
            "tangent": list(self.tangent),
            "transverse_basis": [list(row)
                                 for row in self.transverse_basis],
            "t_grid": list(self.t_grid),
            "center_curve": [list(row) for row in self.center_curve],
            "transverse_radii": list(self.transverse_radii),
            "half_length": self.half_length,
            "grade": self.grade.name,
        }

    @classmethod
    def from_dict(cls, value):
        return cls(
            octant=(None if value["octant"] is None else
                    tuple(int(x) for x in value["octant"])),
            seed=tuple(value["seed"]),
            tangent=tuple(value["tangent"]),
            transverse_basis=tuple(tuple(row)
                                   for row in value["transverse_basis"]),
            t_grid=tuple(value["t_grid"]),
            center_curve=tuple(tuple(row) for row in value["center_curve"]),
            transverse_radii=tuple(value["transverse_radii"]),
            half_length=float(value["half_length"]),
            grade=CertificateGrade[value["grade"]],
        )


@dataclass(frozen=True)
class FoldContinuationPayload:
    root_index: int
    zone_kind: str
    oracles: tuple[FoldOraclePayload, ...]

    kind = "fold"

    def __post_init__(self):
        if self.root_index < 0:
            raise ValueError("root_index must be nonnegative")
        if self.zone_kind not in {"fold", "fold-split"}:
            raise ValueError("invalid fold zone kind")
        oracles = tuple(self.oracles)
        if self.zone_kind == "fold":
            if len(oracles) != 1 or oracles[0].octant is not None:
                raise ValueError("fold payload requires one unsplit oracle")
        else:
            expected = {(a, b, c) for a in (-1, 1) for b in (-1, 1)
                        for c in (-1, 1)}
            actual = {oracle.octant for oracle in oracles}
            if len(oracles) != 8 or actual != expected:
                raise ValueError("fold-split payload requires all 8 octants")
        object.__setattr__(self, "oracles", oracles)

    @property
    def grade(self):
        return min(oracle.grade for oracle in self.oracles)

    def as_dict(self):
        return {"kind": self.kind, "root_index": self.root_index,
                "zone_kind": self.zone_kind,
                "oracles": [oracle.as_dict() for oracle in self.oracles]}

    @classmethod
    def from_dict(cls, value):
        return cls(root_index=int(value["root_index"]),
                   zone_kind=value["zone_kind"],
                   oracles=tuple(FoldOraclePayload.from_dict(oracle)
                                 for oracle in value["oracles"]))


def _payload_from_dict(value):
    if value.get("kind") == "regular":
        return RegularContinuationPayload.from_dict(value)
    if value.get("kind") == "fold":
        return FoldContinuationPayload.from_dict(value)
    raise ValueError("unknown continuation payload kind")


@dataclass(frozen=True)
class ChildRootRestriction:
    root_index: int
    zone_kind: str
    proposed_seed: tuple[float, ...]
    grade: CertificateGrade
    requires_krawczyk: bool = True


@dataclass(frozen=True)
class RestrictedContinuationArtifact:
    artifact: "ContinuationCoverageArtifact"
    coverage: object
    roots: tuple[ChildRootRestriction, ...]

    @property
    def ready_for_child_verification(self):
        return self.coverage.complete and len(self.roots) == len(
            self.artifact.parent.witness.zones
        )


@dataclass(frozen=True)
class ContinuationCoverageArtifact:
    """V2 artifact: V1 coverage plus one continuation payload per zone."""

    parent: ParentCoverageArtifact
    continuations: tuple[object, ...]

    def __post_init__(self):
        if not self.parent.complete:
            raise ValueError("v2 artifact requires complete parent coverage")
        payloads = tuple(self.continuations)
        zones = self.parent.witness.zones
        if len(payloads) != len(zones):
            raise ValueError("one continuation payload is required per zone")
        indices = [payload.root_index for payload in payloads]
        if sorted(indices) != list(range(len(zones))):
            raise ValueError("continuation root indices must be a permutation")
        ordered = tuple(sorted(payloads, key=lambda payload: payload.root_index))
        for zone, payload in zip(zones, ordered):
            if zone.kind == "tube":
                if not isinstance(payload, RegularContinuationPayload):
                    raise ValueError("tube zone requires regular continuation")
                if tuple(payload.seed) != tuple(zone.center):
                    raise ValueError("regular seed must equal its zone center")
            elif not isinstance(payload, FoldContinuationPayload) \
                    or payload.zone_kind != zone.kind:
                raise ValueError("fold zone requires matching fold payload")
        object.__setattr__(self, "continuations", ordered)

    @property
    def grade(self):
        return min(self.parent.witness.evidence().grade,
                   *(payload.grade for payload in self.continuations))

    @property
    def artifact_id(self):
        value = self._unsigned_dict()
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                             allow_nan=False).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def _unsigned_dict(self):
        return {
            "schema": "mub6-parent-continuation-v2",
            "parent": self.parent.as_dict(),
            "continuations": [payload.as_dict()
                              for payload in self.continuations],
        }

    def as_dict(self):
        return {**self._unsigned_dict(), "artifact_id": self.artifact_id,
                "grade": self.grade.name}

    @classmethod
    def from_dict(cls, value):
        if value.get("schema") != "mub6-parent-continuation-v2":
            raise ValueError("unknown continuation artifact schema")
        artifact = cls(
            ParentCoverageArtifact.from_dict(value["parent"]),
            tuple(_payload_from_dict(item)
                  for item in value["continuations"]),
        )
        if value.get("artifact_id") != artifact.artifact_id:
            raise ValueError("continuation artifact digest mismatch")
        if value.get("grade") != artifact.grade.name:
            raise ValueError("continuation artifact grade mismatch")
        return artifact

    def restrict(self, child_center, child_half_widths):
        coverage = self.parent.restrict(child_center, child_half_widths)
        roots = []
        parent_center = self.parent.witness.parameter_center
        for zone, payload in zip(self.parent.witness.zones,
                                 self.continuations):
            proposed = (payload.predicted_seed(child_center, parent_center)
                        if isinstance(payload, RegularContinuationPayload)
                        else tuple(zone.center))
            roots.append(ChildRootRestriction(
                payload.root_index, zone.kind, proposed, payload.grade, True
            ))
        return RestrictedContinuationArtifact(self, coverage, tuple(roots))
