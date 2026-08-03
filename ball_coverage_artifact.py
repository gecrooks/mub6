"""Ball-coverage artifact: schema ``mub6-ball-coverage-v1``.

The ball-vertex route (PROOF_SKELETON section 8) proves a different
statement than continuation v2: no per-root S/Q continuation, no
tube radii — the proof payload is (a) the box-wise association of
every surviving sweep hull cell to a ball, (b) the certified
blanket drift bound feeding the pair budget, (c) the interval
phasor pair layer, and (d) an exhibited coloring.  Per
FAT_TILE_V2_GAP.md this receives its own versioned artifact rather
than dummy-filled continuation-v2 fields.

Serialization is Lean-re-checker compatible (PROOF_SKELETON
section 9): every float in the canonical dict is an IEEE-754
binary64 bit pattern (16 hex digits, big-endian), never a decimal
string; the artifact digest is taken over that canonical form, so
digest and semantics can never disagree by a parse ulp.  ``pretty``
renders decimals for humans only.

The stored drift value is a *claim*: ``rederive=True`` obliges the
re-checker to recompute BU from the named derivation and confirm
``bu_max`` is an upper bound.  The sweep is likewise a *replay
spec* plus an explicit exception list — surviving-cell geometry is
reproduced by deterministic replay, not stored (uncovered cells,
which break completeness, are always stored explicitly).

``check()`` is the Python model of the Lean kernel: it re-verifies
everything derivable from the payload alone (pair layer replay,
coloring properness, budget arithmetic) and fails closed.
"""

import hashlib
import json
import math
import struct
from dataclasses import dataclass

from certificate_result import CertificateGrade, Evidence

SCHEMA = "mub6-ball-coverage-v1"


def float_hex(x):
    """Canonical IEEE-754 binary64 bit pattern (16 hex, big-endian)."""
    x = float(x)
    if not math.isfinite(x):
        raise ValueError("artifact floats must be finite")
    return struct.pack(">d", x).hex()


def hex_float(s):
    if not isinstance(s, str) or len(s) != 16:
        raise ValueError(f"not a binary64 bit pattern: {s!r}")
    x = struct.unpack(">d", bytes.fromhex(s))[0]
    if not math.isfinite(x):
        raise ValueError("artifact floats must be finite")
    return x


def _hex_vec(v):
    return [float_hex(x) for x in v]


def _vec_hex(v):
    return tuple(hex_float(x) for x in v)


@dataclass(frozen=True)
class DriftBoundClaim:
    """Certified blanket phase drift over the parameter box.

    ``bu_max`` is the claimed sup of the certified per-unit rate
    over the box; ``drift = bu_max * h_max * 3`` (l1 over the
    3-box).  ``rederive`` marks the value as a claim the re-checker
    must reproduce from ``derivation`` — it is never trusted.
    """

    derivation: str
    beta: tuple
    half_widths: tuple
    bu_max: float
    drift: float
    rederive: bool = True

    def __post_init__(self):
        if not self.derivation.strip():
            raise ValueError("drift claim requires a derivation reference")
        if len(self.beta) != 3 or len(self.half_widths) != 3:
            raise ValueError("drift claim requires a 3-dimensional box")
        h_max = max(float(h) for h in self.half_widths)
        if any(float(h) < 0 for h in self.half_widths):
            raise ValueError("half-widths must be nonnegative")
        if self.bu_max < 0 or self.drift < 0:
            raise ValueError("drift bounds must be nonnegative")
        # The stored drift may exceed the product (outward rounding)
        # but must never undercut it.
        if self.drift < self.bu_max * h_max * 3.0 * (1.0 - 1e-12):
            raise ValueError("stored drift undercuts bu_max * h * 3")

    def as_dict(self):
        return {
            "derivation": self.derivation,
            "beta": _hex_vec(self.beta),
            "half_widths": _hex_vec(self.half_widths),
            "bu_max": float_hex(self.bu_max),
            "drift": float_hex(self.drift),
            "rederive": self.rederive,
        }

    @classmethod
    def from_dict(cls, value):
        return cls(
            derivation=value["derivation"],
            beta=_vec_hex(value["beta"]),
            half_widths=_vec_hex(value["half_widths"]),
            bu_max=hex_float(value["bu_max"]),
            drift=hex_float(value["drift"]),
            rederive=bool(value["rederive"]),
        )


@dataclass(frozen=True)
class SweepReplaySpec:
    """Deterministic replay spec + outcome of the coverage sweep.

    Surviving-cell geometry is not stored: the re-checker replays
    ``sweep`` with the stored knobs (the split rule is deterministic)
    and re-runs the box-wise ball-inclusion test.  Exceptions —
    cells that failed inclusion — are stored explicitly as
    ``uncovered`` (center, radii) pairs; completeness requires the
    list be empty AND the resume-frontier accounting be closed
    (``frontier_complete``, the coverage-verifier's ledger).
    """

    sweep: str
    wmin: float
    cell: float
    tax_derivation: str
    boxes_swept: int
    hull_cells: int
    uncovered: tuple = ()
    frontier_complete: bool = False
    arithmetic_grade: CertificateGrade = CertificateGrade.SAMPLED_BOUND

    def __post_init__(self):
        if not self.sweep.strip() or not self.tax_derivation.strip():
            raise ValueError("sweep spec requires sweep and tax references")
        if self.wmin <= 0 or self.cell <= 0:
            raise ValueError("sweep knobs must be positive")
        if self.boxes_swept < 0 or self.hull_cells < 0:
            raise ValueError("sweep counts must be nonnegative")
        for cen, rad in self.uncovered:
            if len(cen) != 5 or len(rad) != 5:
                raise ValueError("uncovered cells live in the 5-torus")
            if any(float(r) < 0 for r in rad):
                raise ValueError("uncovered radii must be nonnegative")
        object.__setattr__(
            self, "arithmetic_grade", CertificateGrade(self.arithmetic_grade)
        )

    @property
    def complete(self):
        return not self.uncovered and self.frontier_complete

    @property
    def grade(self):
        if not self.complete:
            return min(self.arithmetic_grade, CertificateGrade.SAMPLED_BOUND)
        return self.arithmetic_grade

    def as_dict(self):
        return {
            "sweep": self.sweep,
            "wmin": float_hex(self.wmin),
            "cell": float_hex(self.cell),
            "tax_derivation": self.tax_derivation,
            "boxes_swept": self.boxes_swept,
            "hull_cells": self.hull_cells,
            "uncovered": [[_hex_vec(c), _hex_vec(r)]
                          for c, r in self.uncovered],
            "frontier_complete": self.frontier_complete,
            "arithmetic_grade": self.arithmetic_grade.name,
        }

    @classmethod
    def from_dict(cls, value):
        return cls(
            sweep=value["sweep"],
            wmin=hex_float(value["wmin"]),
            cell=hex_float(value["cell"]),
            tax_derivation=value["tax_derivation"],
            boxes_swept=int(value["boxes_swept"]),
            hull_cells=int(value["hull_cells"]),
            uncovered=tuple((_vec_hex(c), _vec_hex(r))
                            for c, r in value["uncovered"]),
            frontier_complete=bool(value["frontier_complete"]),
            arithmetic_grade=CertificateGrade[value["arithmetic_grade"]],
        )


@dataclass(frozen=True)
class BallCoverageArtifact:
    """Complete proof payload of one coverage-only fat tile."""

    parameter_center: tuple
    parameter_half_widths: tuple
    r_loc: float
    balls: tuple                 # phase 5-vectors (ball centers)
    drift: DriftBoundClaim
    sweep: SweepReplaySpec
    budget_w: float              # per-coordinate pair budget
    nonorth_pairs: tuple         # (i, j, certified |<u_i,u_j>| lower bound)
    coloring: tuple              # one color per ball
    chi_bound: int

    def __post_init__(self):
        if len(self.parameter_center) != 3 \
                or len(self.parameter_half_widths) != 3:
            raise ValueError("artifact requires a 3-dimensional tile box")
        if self.r_loc <= 0:
            raise ValueError("r_loc must be positive")
        if self.r_loc >= 0.886:
            raise ValueError("r_loc exceeds the ball-vertex lemma radius")
        n = len(self.balls)
        if n == 0:
            raise ValueError("artifact requires at least one ball")
        for b in self.balls:
            if len(b) != 5:
                raise ValueError("balls live in the 5-torus")
        if len(self.coloring) != n:
            raise ValueError("coloring must assign a color to every ball")
        if any(int(c) < 0 for c in self.coloring):
            raise ValueError("colors must be nonnegative")
        if self.chi_bound < 1:
            raise ValueError("chi bound must be positive")
        seen = set()
        for i, j, lo in self.nonorth_pairs:
            if not (0 <= i < j < n):
                raise ValueError("pair indices must satisfy 0 <= i < j < n")
            if (i, j) in seen:
                raise ValueError("duplicate pair entry")
            seen.add((i, j))
            if float(lo) <= 0:
                raise ValueError("non-orthogonality bound must be positive")
        if self.budget_w < self.r_loc + self.drift.drift - 1e-15:
            raise ValueError("pair budget undercuts r_loc + drift")

    # ---- grading ---------------------------------------------------

    def evidence(self):
        return (
            Evidence("ball-vertex-lemma", CertificateGrade.RIGOROUS,
                     f"SKELETON s8: r_loc={self.r_loc:g} < 0.886"),
            Evidence("blanket-drift",
                     CertificateGrade.RIGOROUS if self.drift.rederive
                     else CertificateGrade.SAMPLED_BOUND,
                     f"{self.drift.derivation}; rederive obliged"),
            self_sweep_evidence(self.sweep),
            Evidence("pairs", CertificateGrade.RIGOROUS,
                     f"{len(self.nonorth_pairs)} pairs non-orthogonal "
                     f"by interval phasor replay, w={self.budget_w:g}"),
            Evidence("coloring", CertificateGrade.RIGOROUS,
                     f"exhibited, chi <= {self.chi_bound}"),
        )

    @property
    def grade(self):
        return min(ev.grade for ev in self.evidence())

    @property
    def complete(self):
        return self.sweep.complete

    @property
    def ok(self):
        return self.complete and self.chi_bound <= 5

    # ---- kernel model ----------------------------------------------

    def check(self):
        """Re-verify everything derivable from the payload alone.

        Python model of the Lean kernel: pair-layer replay from the
        stored balls and budget, coloring properness on the
        complement graph, budget arithmetic.  Returns a list of
        failure strings; empty means the payload is self-consistent
        (the drift claim and sweep replay still need their own
        re-derivation — this checks the closed part).
        """
        from fat_tile import _pair_lo
        failures = []
        if self.budget_w < self.r_loc + self.drift.drift - 1e-15:
            failures.append("budget_w < r_loc + drift")
        n = len(self.balls)
        claimed = {(i, j) for i, j, _lo in self.nonorth_pairs}
        for i in range(n):
            for j in range(i + 1, n):
                lo = _pair_lo(self.balls[i], self.balls[j], self.budget_w)
                if (i, j) in claimed:
                    if lo <= 0.0:
                        failures.append(f"pair ({i},{j}) replay refutes "
                                        f"stored bound")
                elif lo > 0.0:
                    # adjacency in the stored graph is a superset of
                    # the true possibly-orthogonal graph: sound for
                    # the coloring, so only report, never fail
                    pass
        for i, j, lo in self.nonorth_pairs:
            replay = _pair_lo(self.balls[i], self.balls[j], self.budget_w)
            if replay < float(lo) * (1.0 - 1e-12):
                failures.append(f"pair ({i},{j}) stored bound {lo:g} "
                                f"exceeds replay {replay:g}")
        colors = [int(c) for c in self.coloring]
        adjacent = {(i, j) for i in range(n) for j in range(i + 1, n)
                    if (i, j) not in claimed}
        for i, j in adjacent:
            if colors[i] == colors[j]:
                failures.append(f"coloring not proper on edge ({i},{j})")
        if max(colors) + 1 > self.chi_bound:
            failures.append("coloring uses more colors than chi_bound")
        return failures

    # ---- serialization ---------------------------------------------

    def as_dict(self):
        return {
            "schema": SCHEMA,
            "parameter_center": _hex_vec(self.parameter_center),
            "parameter_half_widths": _hex_vec(self.parameter_half_widths),
            "r_loc": float_hex(self.r_loc),
            "balls": [_hex_vec(b) for b in self.balls],
            "drift": self.drift.as_dict(),
            "sweep": self.sweep.as_dict(),
            "budget_w": float_hex(self.budget_w),
            "nonorth_pairs": [[i, j, float_hex(lo)]
                              for i, j, lo in self.nonorth_pairs],
            "coloring": [int(c) for c in self.coloring],
            "chi_bound": int(self.chi_bound),
            "grade": self.grade.name,
            "complete": self.complete,
        }

    @property
    def artifact_id(self):
        payload = json.dumps(self.as_dict(), sort_keys=True,
                             separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, value):
        if value.get("schema") != SCHEMA:
            raise ValueError("unknown ball-coverage artifact schema")
        artifact = cls(
            parameter_center=_vec_hex(value["parameter_center"]),
            parameter_half_widths=_vec_hex(value["parameter_half_widths"]),
            r_loc=hex_float(value["r_loc"]),
            balls=tuple(_vec_hex(b) for b in value["balls"]),
            drift=DriftBoundClaim.from_dict(value["drift"]),
            sweep=SweepReplaySpec.from_dict(value["sweep"]),
            budget_w=hex_float(value["budget_w"]),
            nonorth_pairs=tuple((int(i), int(j), hex_float(lo))
                                for i, j, lo in value["nonorth_pairs"]),
            coloring=tuple(int(c) for c in value["coloring"]),
            chi_bound=int(value["chi_bound"]),
        )
        stored = value.get("artifact_id")
        if stored is not None and stored != artifact.artifact_id:
            raise ValueError("ball-coverage artifact digest mismatch")
        return artifact

    def pretty(self):
        """Human-readable summary (decimals; NOT the canonical form)."""
        return (f"BallCoverageArtifact[{self.artifact_id[:12]}] "
                f"beta={tuple(round(float(b), 6) for b in self.parameter_center)} "
                f"h={float(max(self.parameter_half_widths)):g} "
                f"balls={len(self.balls)} r_loc={self.r_loc:g} "
                f"w={self.budget_w:.4f} chi<={self.chi_bound} "
                f"grade={self.grade.name} "
                f"{'COMPLETE' if self.complete else 'OPEN'}")


def self_sweep_evidence(sweep):
    detail = (f"{sweep.sweep} wmin={sweep.wmin:g} cell={sweep.cell:g}; "
              f"{sweep.boxes_swept} boxes, {sweep.hull_cells} hull cells, "
              f"{len(sweep.uncovered)} uncovered; frontier "
              f"{'closed' if sweep.frontier_complete else 'OPEN'}")
    return Evidence("coverage-replay", sweep.grade, detail)


@dataclass(frozen=True)
class BallReadinessReport:
    """Honest machine-readable readiness report (fail-closed).

    Mirrors ContinuationAdapterReport but targets ball-coverage-v1.
    Never claims continuation-v2 readiness; the fat tile's own
    CertificateResult grade is reported but never substituted for
    the artifact's independently derived grade (FAT_TILE_V2_GAP.md).
    """

    ready: bool
    source_ok: bool
    source_grade: CertificateGrade
    artifact_grade: CertificateGrade
    required_grade: CertificateGrade
    ball_count: int
    missing: tuple = ()
    reason: str = ""
    artifact_id: str | None = None

    def as_dict(self):
        return {
            "ready": self.ready,
            "source_ok": self.source_ok,
            "source_grade": self.source_grade.name,
            "artifact_grade": self.artifact_grade.name,
            "required_grade": self.required_grade.name,
            "ball_count": self.ball_count,
            "missing": list(self.missing),
            "reason": self.reason,
            "artifact_id": self.artifact_id,
        }


def ball_readiness(result, artifact,
                   required_grade=CertificateGrade.RIGOROUS):
    """Report whether a fat-tile run yields a usable v1 artifact."""
    required = CertificateGrade(required_grade)
    missing = []
    if artifact.sweep.uncovered:
        missing.append(f"{len(artifact.sweep.uncovered)} uncovered "
                       f"hull cells (explicit exception list)")
    if not artifact.sweep.frontier_complete:
        missing.append("resume-frontier accounting (coverage-verifier)")
    if artifact.sweep.arithmetic_grade < required:
        missing.append(f"sweep arithmetic at "
                       f"{artifact.sweep.arithmetic_grade.name}")
    kernel = artifact.check()
    missing.extend(f"kernel: {f}" for f in kernel)
    if artifact.chi_bound > 5:
        missing.append(f"chi bound {artifact.chi_bound} > 5")
    ready = not missing and artifact.grade >= required
    if not ready and not missing:
        missing.append(f"artifact grade {artifact.grade.name} below "
                       f"{required.name}")
    return BallReadinessReport(
        ready, result.ok, result.grade, artifact.grade, required,
        len(artifact.balls), tuple(missing),
        "ball-coverage-v1 validated" if ready
        else "ball-coverage payload incomplete",
        artifact.artifact_id,
    )
