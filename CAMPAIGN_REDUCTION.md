# Reducing the rigorous campaign cost

Current external estimate (2026-08-02): **8,000–15,000 A100-hours**.
Treat this as a provisional upper estimate, not a settled campaign price.
Recent composed certificates already pass at `h=2.5e-4` at both the generic
reference and the former valley failure, while older pricing notes use
smaller pair-layer ceilings in some of the same territory.

The governing rule is that an optimization may reduce theorem cost only if
its saved work remains covered by explicit `RIGOROUS` evidence. Search,
scheduling, and prioritization may be heuristic; ledger acceptance may not.

## Ranked reduction program

### 1. Re-measure the rigorous half-width distribution

Run a stratified dyadic census over generic bulk, valleys, ordinary walls,
wall intersections, and near-corner shells. Record the largest rigorous
half-width, coverage boxes, root-structure mix, chromatic bound, and time.
Integrate the resulting local density instead of pricing the full domain at
the worst observed depth.

Tools: `rigorous_census.py`, `rigorous_census_report.py`.

### 2. Reuse a coarse coverage witness across pair subtiles

The 3–5 million-box torus sweep is the demonstrated cost center. A parent
coverage witness at roughly `h=3e-4` proves that every root over the parent
box belongs to one of its certified tubes/fold structures. Fine signed-pair
children should inherit those indexed structures and recompute only their
pair bounds and coloring.

At a generic pair width of `2.5e-5`, one parent contains about `12^3 = 1728`
isotropic children. Repeating global enumeration coverage in each child
throws away that factor. Required implementation: export each regular
zone's certified Q-curve and enclosure, and export a child-restriction map
for fold oracles. Child association must be proved from the parent
structure, not inferred by nearest-root polishing.

`campaign_geometry.parent_reuse_comparison` separates the coverage and pair
legs. The saving approaches the number of children per parent only when
coverage dominates; pair work remains and must be chained/batched separately.

The first reuse contract is implemented by
`coverage_contract.ParentCoverageArtifact`. `certify_tile` exports a
content-addressed artifact, and `fully_rigorous_signed_tile` accepts it only
when the child's complete closed parameter box is componentwise contained in
the parent box. The restricted witness is bound to the parent's root
structures. Child tube containment, pair bounds, and coloring are recomputed;
none of those claims are inherited from the parent. Serialized artifacts are
digest-checked on load, and an incomplete parent or protruding child fails
closed.

`parent_coverage_campaign.py` supplies the orchestration layer: an atomic
content-addressed artifact store, a validated Cartesian closed-box partition,
artifact-bound JSONL child records, and grade-aware resume. A record from a
different parent digest cannot suppress work, even when its child coordinates
are identical. Parent creation grades and stores the coverage leg separately,
so a valid global sweep remains reusable even if that parent's pair coloring
failed.

The first real paired measurement is recorded in
`PARENT_REUSE_BENCHMARK.md`. It found that box restriction is sound but parent
root coordinates alone do not initialize offset child tubes: three of three
reused children failed all 48 tubes, while their independent counterparts
certified. Reuse economics are therefore blocked on exporting certified root
continuation data; no speedup projection from the failed path is accepted.

### 3. Use anisotropic boxes aligned with special strata

The fine scale near an ordinary wall is normal to a codimension-one
surface. Do not charge that scale in both tangential directions. A normal
width of `1e-6` and tangential widths of `2.5e-5` already reduce wall count
by about `625x` relative to `1e-6` cubes. Use two fine directions only near
codimension-two intersections. Treat exact surfaces with lower-dimensional
or per-sheet certificates.

`campaign_geometry.py` makes this arithmetic explicit. For a representative
`pi/2 x pi x 2.5e-5` wall band, symmetry factor 32, and 6 seconds/tile:

- `(1e-6,1e-6,1e-6)` cubes: about 501 billion tiles;
- tangential `2.5e-5`: about 802 million tiles (`625x` reduction, still far
  too large);
- tangential `3e-4`: about 5.57 million tiles and 9,281 A100-hours.

Thus anisotropy at the old pair width is insufficient. The tangential axes
must remain near the coarse coverage width, and the normal band should use a
dyadic onion whose width scales with distance from the wall. This calculation
also explains the present 8,000–15,000-hour estimate and identifies batching
the per-tile pair work as the next multiplier after geometry.

Using the measured wall-band count above, changing only the per-child pair
time from 6 seconds to the demonstrated warm-chain 1.59 seconds gives about
2,460 A100-hours. A batched 0.5-second pair leg gives about 773 hours. These
are targets to measure, not certificate claims.

### 4. Subdivide only graph-critical children and directions

Attempt the largest pair box first. If coloring fails, identify the roots
and pair rows responsible for the sixth color. Subdivide only directions
that materially change those rows, while inheriting already-certified pair
bounds. Most measured graphs have chromatic bound 2, so uniform refinement
discards large margins.

Tooling: `coloring_witness.py` provides deterministic DSATUR witnesses,
fail-closed search budgets, critical vertices, and conflict-edge directional
uncertainty scores. It is intentionally not wired into `certify_tile` while
the mean-value-rate validation is in progress.

### 5. Stop contracting roots that are irrelevant to a 5-color witness

The proof needs an exhibited proper coloring, not the tightest possible
tube for every root. Tube-less roots can remain conservative full-edge
vertices. Prioritize structures and pair rows on the coloring obstruction
and stop once a rigorous 5-coloring is available.

### 6. Batch and warm-track non-sweep work

Reuse parent enumerations, warm-track root centers, and batch Newton solves,
Jacobians, Taylor models, and pair bounds across children. Submit many
children per GPU launch. Measure GPU utilization: enumeration and Python
launch overhead should not consume billed A100 time.

### 7. Time-box a structural replacement for the pair campaign

Repeated chromatic bound 2 suggests a persistent bipartition or a small
family of analytic overlap inequalities. Spend one or two focused days
searching for a stratum-wise root labeling and symbolic pair lower bounds.
The upside is elimination of most pair tiles; the numerical campaign stays
as the fallback.

## Immediate measurement plan

1. Run 20–50 stratified points on a repeated dyadic `--h` ladder.
2. Compute the local tile-density integral, not a worst-case flat count.
3. Measure how many pair children fit in each successful coverage parent.
4. Survey anisotropic normal/tangential ladders on wall and intersection
   points.
5. Report A100-hours separately for global coverage, root continuation,
   pair bounds, and orchestration overhead.

## Acceptance and assembly safeguards

- Campaign resume uses only the grade-filtered interval component connected
  to the domain start (`campaign_coverage.py`).
- Disconnected successful islands cannot advance the frontier.
- Parent-witness reuse must bind the parameter box and root structures;
  polishing may suggest associations but may not certify them.
- Every cost survey must retain failed and downgraded records.
