# Ball-coverage artifact v1 — schema proposal

*Proposed 2026-08-02 in response to FAT_TILE_V2_GAP.md. For joint
review before any freeze; nothing here is load-bearing yet. The
Python reference implementation is `ball_coverage_artifact.py`
(tests: `test_ball_artifact.py`); `fat_tile.py` emits it via
`fat_tile(..., artifact_out=dict)`.*

## Position

FAT_TILE_V2_GAP.md is right: the ball route proves a different
statement than continuation v2 and gets its **own versioned
artifact** (`mub6-ball-coverage-v1`), never dummy-filled v2 fields.
The gap doc's grade rule is adopted verbatim: the fat tile's
CertificateResult grade is *reported* in the readiness path but
never substituted for the artifact's independently derived grade.

## Proof payload

A tile box `(beta, h)` is closed by:

1. **ball-vertex lemma** (SKELETON §8, global, `r_loc < 0.886`):
   each ball contributes ≤ 1 clique vertex;
2. **coverage**: every surviving hull cell of the deterministic
   sweep fits box-wise inside some `ball(root, r_loc)`;
3. **blanket drift**: certified `bu_max = sup |beta-unit rate|`
   over the box gives per-coordinate budget
   `w = r_loc + bu_max·h·3`;
4. **pair layer**: interval phasor sums at budget `w` certify
   listed pairs non-orthogonal;
5. **coloring**: exhibited proper coloring of the complement
   ("possibly-orthogonal") graph with `chi ≤ 5`.

## Design rules (the two that matter)

**Claims are re-derived, never trusted.** `drift.bu_max` is stored
with `rederive: true` — the re-checker recomputes it from the named
derivation and checks the stored value is an upper bound. Same for
the sweep tax constants (`tax_derivation`).

**Replay spec + explicit exceptions, not stored geometry.** The
sweep's surviving-cell set is reproduced by deterministic replay of
`(sweep, wmin, cell)`; we do not store 10⁴–10⁵ hull cells per
tile. What *is* stored explicitly is every exception: `uncovered`
cells (each breaks completeness) — so an artifact is checkable as
complete without replay, but rigorous only after replay.
*Open: the guarded-sweep measurement now in flight decides whether
survivor counts at fat h are small enough that storing the full
hull-to-ball association is cheaper than replay. The schema admits
either; v1 as implemented stores exceptions only.*

## Fields (canonical dict)

| field | type | notes |
|---|---|---|
| `schema` | str | `"mub6-ball-coverage-v1"` |
| `parameter_center`, `parameter_half_widths` | 3× f64 | tile box |
| `r_loc` | f64 | must be < 0.886 (validated) |
| `balls` | n × 5× f64 | phase 5-vectors |
| `drift` | claim | `derivation`, `beta`, `half_widths`, `bu_max`, `drift`, `rederive` |
| `sweep` | replay spec | `sweep`, `wmin`, `cell`, `tax_derivation`, `boxes_swept`, `hull_cells`, `uncovered[]`, `frontier_complete`, `arithmetic_grade` |
| `budget_w` | f64 | ≥ `r_loc + drift` (validated) |
| `nonorth_pairs` | (i, j, lo)[] | certified lower bounds, i < j |
| `coloring` | n × int | exhibited witness |
| `chi_bound` | int | claim; kernel re-checks |
| `grade`, `complete` | derived | weakest dependency; sweep closure |

`artifact_id` = SHA-256 over the sorted canonical JSON (same
convention as `ParentCoverageArtifact`).

## Serialization: floats are bit patterns

Per the standing Lean-re-checker constraint (SKELETON §9): every
float in the canonical dict is the IEEE-754 binary64 **bit pattern**
(16 lowercase hex digits, big-endian), never a decimal string.
Digest and semantics therefore cannot disagree by a parse ulp, and
the Lean kernel parses floats with zero transcendental or
decimal-rounding trust. `pretty()` renders decimals for humans;
decimals never enter the digest.

**Ask:** adopt the same convention in `coverage_contract` /
`continuation_artifact` serialization before any schema freezes —
`as_dict` currently emits Python floats (decimal JSON), which the
Lean kernel would have to re-parse with rounding trust.
`ball_coverage_artifact.float_hex`/`hex_float` are importable as-is.

## Self-check (Lean kernel model)

`BallCoverageArtifact.check()` re-verifies everything derivable
from the payload alone, fail-closed: pair-layer replay from stored
balls at `budget_w` (refutes inflated claims), coloring properness
on the complement graph, `chi_bound` consistency, budget
arithmetic. Left as external obligations: drift re-derivation,
sweep replay, and **frontier accounting** — the emitter hard-codes
`frontier_complete: false` until it can consume your resume-frontier
ledger, so every artifact honestly reports OPEN today.

## Readiness path

`ball_readiness(result, artifact)` mirrors
`ContinuationAdapterReport` (fail-closed, machine-readable
`missing[]`), targeting v1 — it never claims continuation-v2
readiness. Current honest output on a real run: `ready: false`,
missing = frontier accounting + sweep arithmetic grade.

## Open questions for the freeze

1. `frontier_complete`: should this embed (a reference to) your
   frontier ledger artifact id rather than a bare bool, so the
   claim is content-addressed like everything else? (My preference:
   yes — `frontier_ledger: <artifact_id> | null`, v1.1.)
2. Association encoding: exceptions-only vs. full hull-to-ball
   list — decide on the guarded-probe numbers.
3. Should `SweepReplaySpec.arithmetic_grade` be promoted RIGOROUS
   by the replay itself (kernel-side) rather than asserted by the
   emitter? (My preference: yes; the emitter's value is then only
   advisory.)
