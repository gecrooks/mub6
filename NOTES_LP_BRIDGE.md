# The 1/6-quantization and the Matolcsi–Weiner LP bound: hunting a bridge

Research note, first pass (2026-07-29). Status: conjecture + attack
routes; nothing here is proved beyond the certified computations cited.

## 1. The observed 1/6-quantization (what the machine keeps finding)

Three independent obstructions, all landing EXACTLY on 1/d = 1/6:

O1 (Grassl/F6 wall). At the Fourier family F(1/6, 0)-adjacent
    quadruple configurations, the certified minimum of the worst
    unbiasedness defect over the constellation is 1/6 - o(1) with the
    limit value 1/6 exact (Result 8 era; the "quadruple wall").
O2 (X-triple extension gap). For every certified triple {I, H, K}
    across ALL of X6^(2) (B-arc, Dita arc, Szollosi region): the
    best-possible fourth-basis vector misses MU-ness by a defect whose
    observed infimum over the family is again 1/6-quantized — the gap
    never closes and never dips below the 1/6 scale (Results 19-24).
O3 (base-pair deviations). The per-pair Gram deviations |<a|b>|^2-1/6
    at optimal 4-basis constellations concentrate at exact multiples
    of 1/6 (early Results 2-3 numerics; the {5^3, 1} constellation).

The suspicious constant: 1/6 = 1/d is the SAME number that appears in
the Matolcsi–Weiner LP relaxation as the unbiasedness value itself —
suggesting the failures are not soft-analysis accidents but the
integer-programming gap of a hidden LP.

## 2. The Matolcsi–Weiner LP structure (what is known)

MW (and Jaming et al.) attack MUB existence via harmonic analysis /
linear programming: a complete set of d+1 MUBs in C^d yields, through
the Fourier-side positivity constraints on the associated (d^2, d)-
frame difference sets, a feasibility LP. For d = 6 the LP itself is
FEASIBLE (hence no obstruction), which is why the approach stalls:
the LP relaxes the true problem — its feasible point is not known to
come from an actual MUB set. The known integrality-style gap: the LP
value sits at the rational point where the true combinatorial object
would need integer multiplicities, and 1/d is the quantum of the
relaxation's mass assignments.

## 3. The conjectured bridge

**Conjecture (bridge).** The certified 1/6-walls O1-O3 are the
strong-duality face of the MW LP: for any near-quadruple
configuration in the Karlsson-reachable stratum, there is a feasible
MW dual certificate whose objective equals (worst defect) - 1/6, so
the defect cannot drop below 1/6 anywhere the dual stays feasible.
Equivalently: the variational landscape's global gap over K6^(3) is
an LP-duality gap, and the tile campaign is computing, pointwise,
what one dual certificate would give uniformly.

If true, the payoff is decisive: ONE dual function (a trigonometric
polynomial with positivity constraints — a finite object, checkable
by the same interval machinery we already have) replaces ~1e10 tiles.
The campaign becomes the fallback, not the plan.

## 4. Attack routes

R-A (dual reconstruction from data). We hold certified defect
    landscapes at many beta. Fit a candidate MW dual (few dozen
    Fourier coefficients on Z_6 x Z_6-type support) to reproduce the
    1/6 floor across the certified points, THEN try to verify its
    positivity constraints exactly (interval arithmetic — the
    machinery exists). Falsifiable fast: if no low-degree dual fits
    even numerically, the bridge in this form is dead.
R-B (weighted-defect LP per triple). For a FIXED certified triple
    {I, B, K}: the fourth-vector defect min is itself an eigenvalue-
    like quantity; write its LP relaxation explicitly (variables:
    |<u, h_k>|^2 masses; constraints: frame identities; objective:
    max deviation). Solve at a few certified triples and compare the
    LP optimum to the certified 1/6 gap. Equality at even one triple
    = the bridge exists locally; the question becomes uniformity.
R-C (literature triangulation). The needed dual may be implicit in
    Jaming–Matolcsi–Mora–Ruzsa–Weiner's d=6 computations (they
    verified LP feasibility for COMPLETE sets; the QUADRUPLE
    sub-problem's LP may be infeasible — nobody seems to have run
    the LP for 4 bases restricted to an H2-reducible member, which
    is exactly our theorem's hypothesis class).

R-B is the cheapest decisive test — BUT (first-pass finding, this
session): the naive PSD relaxation is TRIVIALLY feasible: P = I/6
satisfies diag(P) = 1/6, b^H P b = 1/6, k^H P k = 1/6 exactly, for
ANY triple. So t*_SDP = 0 always, and the bridge cannot be a generic
rank-relaxation — the maximally-mixed cheat must be excluded by
STRUCTURE, exactly as Delsarte/MW positivity does for group
Hadamards.

R-C' (the sharpened direction). The literature's Delsarte LP needs a
group structure Karlsson matrices lack — but H2-reducibility IS a
structure: an exact Z_2-grading (2x2 block form, the Moebius/cocycle
data of Karlsson's parametrization). Conjecture refinement: the
correct LP lives on the Z_2-graded autocorrelation — positive-
definite functions on the graded object whose even part sees the F2
blocks and whose odd part carries the (theta, phi, lam) phases as
continuous characters. The 1/d quantum then appears as the Delsarte
dual's mass unit on the 6-point orbit, matching O1-O3. Concrete
first computation: write the graded autocorrelation constraints for
a vector unbiased to {I, H(beta)} with H2-reducible H, and check
whether the QUADRUPLE sub-problem's graded LP is INFEASIBLE at one
certified beta — infeasibility of the 4-basis graded LP at a single
point, verified in interval arithmetic, would be the bridgehead: a
new-mathematics kill of the campaign's cost problem.

## 4.5 First experimental cut (2026-07-29): the 1/6 is a CLASS GAP

Probe at the group point F6: the 48 MU vectors of {I, F6} have
pairwise |<v,w>|^2 organized into sharp spectral classes
(counts/1128): 0 x300, 1/36 x18, 2/36 x72, 3/36 x72, 5/36 x72,
**1/6 x216**, 10/36 x288, 1/3 x90 — a few classes carry small
Bjorck-phase irrational offsets, but the 0 and 1/6 classes are
exact. Reading: a fourth basis needs a 6-clique in the ORTHOGONAL
class; the certified {5^3,1} wall-at-1/6 is the best sixth vector
landing in the 1/6 (mutually-unbiased) class instead — **the wall is
the spectral gap between the 0-class and the 1/6-class of a
near-rigid few-class overlap system**. This is association-scheme
shaped: exactly what Delsarte LP formalizes. The deformation
program: track the class structure as beta leaves F6 into K6^(3)
(classes deform continuously; the certified walls along the family
are conjecturally this same gap), and write the graded LP whose dual
face pins the gap at 1/d. Next experiments: (a) clique analysis of
the 0-class graph at F6 (does it top out at 5 + why); (b) the same
spectrum at a certified Karlsson beta — do the classes survive with
deformed values and does the 0-to-MU gap stay >= 1/6 - o(1)?

## 4.6 Second cut: clique structure and the two-regime picture

(a) F6 0-class max clique = 6 — as it must be: the 6-cliques ARE the
partner bases of the known triples. The wall therefore lives one
level up: a quadruple needs TWO disjoint 6-cliques with all 36
cross-pairs in the 1/6 class; the class spectrum forces an expelled
cross-pair to 0 or 1/3 — both at distance EXACTLY 1/6 from the MU
class. The certified defect-1/6 is the class-neighbor distance in a
rigid overlap scheme, at the two-basis level.

(b) Generic Karlsson point (certified reference beta): the
structure collapses — 24 orthogonal pairs (vs 300), max clique 2,
only 98 pairs near 1/6. Matches the tile machine's ubiquitous
"clique <= 2" colorings.

**Two-regime consequence (potentially bigger than the bridge):**
precision is only needed near the group/F-strata, where
near-quadruples almost assemble and the 1/6 scheme-gap operates —
the natural domain of exact Delsarte/character arguments plus a
certified perturbative collar. The generic bulk needs only "max
orthogonal set <= 5" — clique STARVATION, provable with far coarser
and cheaper certificates than 1/6-precision tiles. The campaign
cost model assumed uniform fine precision; a two-regime proof
(scheme rigidity near strata, coarse starvation certificates
elsewhere) could collapse it by orders of magnitude. Next: map how
the 0-class population decays moving from F6 into the family
(where does clique 6 -> 2 happen), and price a coarse
starvation-certificate tile.

## 4.7 The collar is thin (decay map, one ray)

The Fourier locus sits ON the theta=0 branch face of the Karlsson
chart at beta_F ~ (0, 2.0412, pi) — the chart's degeneracy, the
branch-margin hot zone, and the rigidity zone are THE SAME PLACE.
Walking off the face along theta at (., 2.0412, pi):

  theta:   0.001  0.003  0.01  0.03  0.1  0.3  0.8  1.5
  clique:      6      4     4     2    2    2    2    2

The 6-clique (candidate-basis) structure dies within theta ~
0.001-0.003; fully generic by 0.03. If this holds across loci and
directions (CAVEAT: one ray, one locus, multistart enumeration —
needs the same census at other F/Bjorck/S6 strata and transverse
directions), the two-regime split is: a collar of a few percent of
the domain needing fine machinery (or the scheme argument), and a
>= 97% bulk needing only clique-starvation certificates, whose
open design question is root-set control over FAT tiles (h ~
0.01-0.05) — the overlap bounds themselves survive such widths
easily (drift O(1) vs margins 0.16). Cost impact if it works:
bulk tile count ~1e5-1e6 instead of 1e10; the campaign becomes a
$1e3-1e4 object plus a collar theorem.

## 5. What would falsify the bridge

- R-B LPs with optimum strictly below the certified gap (relaxation
  too weak; the walls are genuinely nonlinear).
- A certified point with defect wall NOT at a 1/6 multiple (none
  found in ~30 results, but the survey's h<3e-4 territory has not
  been wall-measured).

## 4.8 Census verdict (60 points, census_2026-07-29.log)

The rigidity zone is the ENTIRE theta=0 branch face (clique-6
regions extend across large 2-D swaths of (phi, lam) at theta=0.002
— the face is the (near-)group-Hadamard locus, plausibly the
Fourier family as the chart boundary), but it is THIN transversally:
every ray tested collapses to clique 2 by theta ~ 0.017. Dangerous
volume: face x [0, ~0.02] < 1% of the domain. The bulk never
exceeds clique 4 (only 6 kills). Even on the face, clique 6 alone
is benign (the known triples); quadruple danger needs two cross-MU
6-cliques simultaneously. TWO-REGIME PLAN CONFIRMED at enumeration
grade. Next: fat-tile starvation certificate design (bulk), and the
collar handled by anisotropic thin-in-theta tiles + the s^2 branch
chart (both already built) or the scheme argument.

## 4.9 Fat-sweep probe status (inconclusive — OOM)

First feasibility probes of the bulk starvation sweep (starve.py, no
collection machinery) died by OOM at h in {0.003, 0.01} — the LIFO
stack balloons where the fine-tile sweeps held 450 MB, i.e. the
fat-tax exclusion profile differs qualitatively from the analytical
model (exclusion is NOT biting at the expected widths somewhere).
Hypotheses: missing coarse root-collection (every root/trench
neighborhood grinds to the floor with no oracle to swallow it), or
a genuine stack-growth pathology of the plain sweep at fat taxes.
NEXT: instrument fat_sweep (per-level counts, exclusion fraction,
RSS), find the explosion level, then either add coarse blob-oracles
(collect first, certify blob-pairs after) or stream survivors to
disk. The two-regime economics stand or fall on this number.

## 4.10 Bulk feasibility MEASURED (streaming wmin scan, h = 3e-3)

Survivor volume fraction at the bulk test point: 1.3e-2 (wmin 0.1)
-> 5.8e-4 (0.05) -> 4.1e-5 (0.025) — ~x20 per halving, i.e. ~w^4:
thin shells around 1-DIM TRENCH CURVES (codim 4 in the 5-torus).
Sweep cost trivial: 23M boxes / 16 s single-core CPU at the finest
rung (the earlier OOM was probe survivor-hoarding, not math).
Economics: ~7e7 bulk tiles at h=3e-3, sub-second GPU sweeps ->
~$1-3k for the whole bulk; h=0.01 if tolerable -> ~$100.
REMAINING DESIGN: segment the surviving filaments (~0.2 rad pieces,
overlap drift ~0.1 vs margins 0.4), certify segment-graph clique
<= 5 via interval overlap bounds (census: actual bulk cliques <= 4).
Then: bulk = starvation tiles, collar = existing fine machinery on
<1% volume, and the full-domain theorem prices at ~$10^3.

## 4.11 Union starvation certificate: NEGATIVE result

The cheap bulk certificate (union-over-beta survivor segments +
pairwise Lipschitz overlap bounds + coloring) FAILS decisively:
at the first bulk point, 82k segments, 1.9e8 possible-orthogonality
edges, coloring bound 14 >> 5 (min self-coherence fine, 0.877).
Cause: the union relaxation discards the fixed-beta correlation the
census exploits — filaments sweep continuous vector families whose
mutual overlaps traverse near-zero broadly, even though at every
FIXED beta the discrete root set has clique <= 4. Lesson recorded:
bulk starvation needs beta-CORRELATED structure (per-branch root
counting + same-beta cross-branch bounds), which is fine-tile-like
machinery — OR the bulk simply runs the existing tile certificates
at their post-Result-34 ceiling (7e-4..1e-3), a 13-37x saving over
the 3e-4 flat plan rather than 1000x. Current honest cost picture:
bulk at ~1e-3 + thin collar fine => ~$30-80k full domain; the
$1e3 figure required the union certificate that just died. The
beta-correlated blob design remains open (genuinely promising, but
research, not engineering).

## 4.12 Beta-correlated blob certificate: DESIGN

The union failure (4.11) dictates the architecture: pair bounds must
compare roots AT THE SAME beta. Design ("fat-tube tile"):

1. COARSE sweep (16 s, wmin 0.025) -> survivor filaments; connected
   components = branches S_1..S_m (~50-100).
2. Branch classes, each with a beta-localization statement:
   R (regular, sigma_min >~ 0.1): parametric Krawczyk with fat box
     (rho ~ 0.05-0.1) — CONTRACTS at h = 3e-3 (q ~ ||J^-1||(RJ(rho)
     + L_HJ rate h) ~ 0.2): theta_i(beta) unique in the box for all
     beta. Localization ~ rho.
   D (deep but continuable, sigma_min ~ 0.01-0.1): 1-param beta-
     continuation (bnrates frame-split machinery, 3 params): motion
     ~ rate*h ~ 0.015 — localization excellent IF sigma_5 clears the
     J-drift over the motion (the known frontier).
   W (wild: near-fold in beta over the box): dip-localization along
     the trench (valley machinery) or per-branch union fallback.
3. Same-beta pair matrix: f_ij(beta) >= f_ij(beta0) - |grad f|_cert
   * h*sqrt3 - (loc_i + loc_j)/6 - SLOP, where f_ij(beta0) is the
   center-root overlap (the tile machine's O0) and |grad f| from
   certified root-motion rates (du/dbeta = -J^-1 dg/dbeta, already
   built). Overlap tolerance is LOOSE (~0.1-0.2 theta), which is
   why fat localization suffices.
4. Coloring of {f_lower <= 0} + <=1-vertex-per-branch-at-fixed-beta
   (from uniqueness in R/D; W branches contribute their dip count).
   Need <= 5 colors; census says actual cliques <= 4 in the bulk.

Cost model: per tile = coarse sweep (seconds) + m Krawczyk/
continuation certificates (ms each) + m^2/2 pair bounds (vectorized)
— SECONDS per tile instead of 30-100 s, and at h = 3e-3 instead of
3e-4: the fine sweep (the actual cost center) disappears from the
bulk entirely. If W-branch handling holds, bulk ~ 7e7 tiles x
seconds-CPU ~ $1-3k on CPU fleets alone; even h = 1e-3 with fat
tubes beats the $30-80k fallback by ~10x. The h-ceiling becomes
W-branch economics, not sweep economics.

RISKS: W-branch fraction in the bulk (census suggests deep roots
are common but folds-in-beta over 3e-3 boxes may be rare); branch
connectivity at wmin 0.025 (components may merge through junctions,
weakening per-branch uniqueness — split at junctions by sigma
profile). Next build: branch decomposition + class R fat-Krawczyk
end-to-end at one bulk point.

## 4.13 Graded LP: the naive grading fails at exactly 1/sqrt2

Projection test at the certified reference beta: the unbiasedness
kernels h_k h_k^H have off-graded Frobenius fraction 0.7071 under
the plain (a - a') grading — exactly half the mass. Signature: the
H2-blocks are GAUGE-TWISTED F2's (internal diagonal phases from
Karlsson's parametrization), so the correct object is a gauged /
cocycle grading. NEXT COMPUTATION (decisive, afternoon-scale):
optimize diagonal gauges d in U(1)^6 (and/or a ZZ_2 cocycle) to
minimize the joint off-graded mass of all six kernels; if a gauge
zeroes it, assemble the LP in that gauge and test quadruple
infeasibility at one certified beta; if no gauge exists, this
bridge route is falsified (the two-regime machine path stands
regardless).

## 4.14 Graded LP: FALSIFIED in all simple forms

Gauge probes at three beta (bulk + near-face): phase gauge U(1)^5
reduces off-graded mass 0.50 -> 0.327; full U(2)^3 block gauge adds
NOTHING (identical plateau); the plateau is family-wide and near the
face equals 1/3 to 7e-6. Z3-circulant gradings on both tensor
factorizations also fail (0.40, 0.50). Conclusion: Karlsson kernels
admit no plain group-difference grading in any gauge — an
independent structural confirmation of the literature's "LP cannot
reach K6(3)". The exact gauge-invariant 1/3 is recorded as an open
curiosity (graded sector is exactly half the Hermitian space; naive
off-mass 1/2; the gauge recovers exactly 1/6 more — why 1/3 is
irreducible deserves a short algebra session someday).
SURVIVING LP SCOPE: the collar — ON the branch face the Hadamards
are group-type and the classical Delsarte/scheme machinery applies;
the LP is a collar tool, the machine covers the bulk. Coherent with
the two-regime split; the bulk's fate rests on the beta-correlated
blob certificate (4.12).

## 4.15 Stage-1 branch build: failure diagnosis and the corrected design

Static parametric Krawczyk FAILS for fat Karlsson tiles — root cause:
4.12's contraction estimate used the B-N arc rate (~0.6); Karlsson
entry rates are ~10x hotter, so root MOTION over the beta box
(S h ~ 0.03-0.3) rivals any contractible box and the H-ball
perturbation (L_HG * entry-drift ~ 0.06-0.2) swamps sigma*rho.
Also cluster_suspects(link) fragments the filaments (451k clusters
from 1.4M boxes) — branch decomposition needs true union-find over
touching survivor boxes.

CORRECTED DESIGN — FAT Q-TUBES (all machinery exists):
per branch, localize the root as theta(beta) in theta0 + S dbeta +-
(TM curve residual + pad): residual scales h^3 (measured 2.2e-8 at
3e-4 => ~2e-5 at 3e-3 — negligible); pair bounds between branches
i,j pay only ||(S_i - S_j) row||_1 * h / 6 (~0.001-0.03 vs 0.4
margins). Deep roots (S undefined/huge) stay WILD pending stage 2.
Build: union-find branches -> root_data2 per branch ->
certified_curve_residual at hv=(3e-3)^3 -> parallelepiped pair
bounds -> coloring + wild count. The open empirical question is the
wild fraction per bulk tile.

## 4.16 Roots-first fat tubes: correlated bounds WORK; architecture
##      converges back to certify_tile-in-coarse-mode

At bulk point 1 (9 s total): 54 roots (4k survivor polishes), 37
S-sloped tubes, pair-bound COLORING = 2 — matching the census's
clique <= 2: the beta-correlated machinery works and is nearly free.
Wild (deep-sigma) roots: 17/54 = 31% — the open front; the crude
bound colors + wild = 19 is useless until wilds get their stage-2
treatment (they do not actually clique). Survivor "coverage" by
tubes is the wrong criterion (survivors are un-excluded territory,
not candidate roots); the right completeness semantics is the fine
machine's: sweep-with-collection. CONCLUSION: the fat bulk
certificate = certify_tile in COARSE MODE (first-order taxes +
coarse wmin + S-fat tube collection + valley windows for wilds) —
retrofit, not a parallel codebase. Next session: coarse-mode
certify_tile at h = 1e-3..3e-3 on the survey points; the pre-
Result-34 1e-3 wall deserves a rematch with the new taxes.

## 4.17 Coarse-mode certify_tile at h = 1e-3: the architecture stands

With stage-B wmin lifted to 5e-3 (Result-34 taxes active): the tile
certifies except 45,572 stuck boxes = vol frac 4e-11 (!) clustering
into 48 WILD BLOBS whose polished roots all have sigma in
0.03-0.10 — precisely the R7 valley-window class. No budget
blowout (the pre-R34 "1e-3 impossible" is overturned in coarse
mode). Coarse-mode bulk tile = coarse sweep + S-fat tubes
(coloring 2) + ~48 valley windows. NEXT (first experiment of next
session): valley windows at fat h (iso 1e-3 -> 3e-3) — how many of
the 48 certify with existing R7 floors, and what does the per-tile
cost settle at; then the h=3e-3 version (7e7 tiles) prices the
bulk for real.

## 4.18 The coarse-mode tile closes end-to-end at h = 1e-3

The 48 wild blobs -> 48 distinct valley roots -> R7 valley windows
at hv = (1e-3)^3: 47 fully consistent certified floors, 1 sampled-
mode fallback, 0 hard failures, 69 s total (~1.4 s/window). The
bulk tile architecture is thus COMPLETE at 1e-3: coarse sweep +
S-fat correlated tubes (colors 2) + valley windows for the deep
class. Remaining for the bulk pricing: the same assembly at
h = 3e-3 (7e7 tiles for the domain) — valley count and window
viability at 3e-3 set the per-tile cost; then assemble the single
"coarse_certify_tile" entry point and re-survey. The road from the
survey's $1M through Result 34 and the two-regime program now ends
at a measured, working per-tile recipe.

## 4.19 The valley-h wall is real beta-motion; recursive splitting
##      is the architecture; fast-mover fraction is the price

First-order valley taxes (retrofitted into floors + T-gate) and
floors-first acceptance (the env-vs-far-tax gate is vestigial in
coarse mode — R33's dissolution one level up; implemented) do NOT
move the 3e-3 wall: probe of the binding valley (sv =
[.32,.32,.23,.071,.054], floors -0.06 at the ends) shows its trench
genuinely translates ~0.05 rad across the beta box — more than its
wall height. Not tax conservatism: real structure motion. Treatment:
RECURSIVE per-valley beta-splitting to each valley's working-h
(octant split generalized; fine work concentrates on fast movers
only). Bulk pricing now reduces to ONE measurable: the fast-mover
fraction per fat tile (valley 11 was the first failure of ~50; if
~5-10/50 need 1e-3-depth recursion, coarse tiles at 3e-3 cost
~3-5x less valley work than all-1e-3, times chain amortization).
NEXT SESSION: (1) count fast movers per tile at 3e-3 across the
survey points; (2) implement recursive valley split; (3) valley
chain-amortization; then re-price.

## 4.20 Fast-mover census (3 points, h = 3e-3): recursion alone
##      is not enough; amortization is load-bearing

Fast fractions 42/48, 22/56, 29/56 (88/39/52%) — heterogeneous;
tube-class (sigma > 0.12) EMPTY at all three points: bulk roots are
soft nearly everywhere. Pricing: recursive splitting buys only
~1-2.5x over all-1e-3 valley work. THE remaining levers, in order:
(1) CHAIN AMORTIZATION of valley windows along beta (windows vary
smoothly; profile once, re-verify cheap — the cache.py pattern,
expected 5-20x, now load-bearing); (2) window batching (30k tiny
solves/tile in python loops, ~10x constant factor); (3) census cost
itself says the per-window 1.4 s must drop for any plan. Bulk cost
today: ~1e-3-equivalent valley work ~ $30k-scale before levers;
with (1)+(2) plausibly $2-6k. Next session: implement valley window
chain-reverification and measure its amortization factor — that
number decides the final price.

## 4.21 Night of engineering: 60x windows; recursion depth = 1 level

Batched valley machinery (foldbatch.py): floors 40x (251->6 ms,
sound regardless of Newton convergence — Fres is measured), block-
warm batched profile, lite corners: 110 ms/window vs 6.5-12.6 s
(~60-110x). Parity 11/13 on the check set (2 marginal windows
become recursion candidates — cost, not soundness). Wide census
(6 fresh points): 89% fast at 3e-3 with heavy variance — BUT
recursion-depth measurement at the worst point (0/52 fat): 11/12
fast movers certify at ONE octant level (8 subs at 1.5e-3); one
root resists to 3.75e-4 (wild-file). True fast multiplier = 8.
PRICING NOW: ~60 s CPU per un-amortized 3e-3 tile x 7e7 tiles ~
$35-58k; the chain lever (windows: skip T-escalation; enumeration:
warm-tracked roots; rates: incremental) at its expected 10-15x
lands $3-6k. Chain mode is the ONE unmeasured lever left.

## 4.22 coarse_certify_tile v1: COMPLETE fat tiles certify at 3e-3

The assembled entry point (coarse_tile.py: enumerate -> batched
windows with one-level recursion -> pair coloring WITH wild vertices
included at blanket 0.75 theta-l1 localization -> coarse sweep with
blob-coverage-vs-root-list semantics) CERTIFIES complete tiles at
h = 3e-3, including the census's worst point (0/52 fat valleys;
colors 3) and a mixed point (colors 4). Coverage exact (0 uncovered
blobs). 47-99 s/tile un-amortized CPU (sweep 30-70 s -> ~2 s on
GPU; windows 13-25 s; enum 4 s). The key conceptual close: wild
roots need no special certificate — they are colored VERTICES with
fat localization, and the overlap-class margins absorb them.
Remaining to the campaign: robustness across more points (running),
chain mode (enum warm-tracking + window T-caching, the 10-15x),
collar tiles, driver integration, and the coarse-mode lemma prose
in PROOF_SKELETON. Pricing at v1 measured numbers: ~$30-60k raw,
$3-6k with chains — the original tweet's number, reached honestly.

## 4.23 End of night: coarse tile 5/6 census points, measured wilds

Per-wild MEASURED localization (sub-beta root spread from the
recursion polishes, x2 + pad) replaced the 0.75 blanket: the
colors-6 failure dropped to colors 2 — wild roots barely MOVE
(their windows fail on wall height, not motion); the blanket
overcharged ~5x. Census scoreboard for coarse_certify_tile v1 at
h = 3e-3: **5 of 6 points CERTIFIED end-to-end** (colors 2-5,
coverage exact everywhere). The sixth has ~10x fatter filaments and
OOMs the in-RAM survivor collection at 8e8 boxes — needs the
streaming-collection engineering (chain-mode list), not new math;
also flag it for the tile-cost distribution tail. Timings/tile:
46-208 s un-amortized CPU (sweep 30-164 s dominates; GPU 41x on it;
windows 13-40 s; enum 4 s). Chain mode remains the one unmeasured
lever between the measured $30-60k raw and the $3-6k target.

## 4.24 CENSUS CLOSED 6/6 — the coarse tile is real

The fat-filament holdout CERTIFIES at the adaptive rung (h = 1.5e-3,
colors 2, coverage exact) — exactly the campaign driver's designed
response; no new machinery needed. FINAL SCOREBOARD, coarse_certify_
tile v1 across all six census points: 6/6 CERTIFIED (five at 3e-3,
one at 1.5e-3), colors 2-5, coverage exact everywhere, 46-208 s
un-amortized CPU per tile. The night's engineering: 60-110x valley
windows (batched floors/profile/lite corners), one-level recursion
(x8, measured), measured per-wild localization (colors 6->2), wilds
as colored vertices, streaming hull sweep (written; wiring proven at
1.5e-3), adaptive h as the fat-filament response. Remaining to the
campaign: chain mode (the 10-15x, sole unmeasured lever), collar
tiles, driver wiring, coarse-mode lemma prose, and the certified-
grade pass over tonight's sampled-grade pieces (wild locs, batched
anchors are definitional — mostly bookkeeping). Pricing: $30-60k
raw, $3-6k with chains, at 3e-3/adaptive.

## 4.25 First certified coarse chain (anchor + 6 steps, all OK)

coarse_chain.py: anchor (160 s) then six 1.6h steps at 30-67 s, ALL
CERTIFIED (colors 2-3, coverage 0 — the w^4 whisker-tail blobs now
get a refinement fallback: mini-sweep at a finer floor; exclusion is
exclusion). The engineering ladder inside: FIRST-ORDER CACHE (signed
margin gradient D1 per cached box; E + D1.db - curv dist^2) cut
re-verify failures 17x; PERSISTENT BLOBS (clustered once, steps
re-associate reps + refinement fallback) removed the per-step stuck
re-sweep. Cache degradation measured: failures 378k -> 3.0M across
6 steps => re-anchor cadence ~5-6. Net amortization ~3x (window
stage ~22-40 s/step is now binding; sub-window caching is the next
~2x). Honest campaign math at today's state: ~50 s/tile chained ->
$25-40k CPU; the remaining gap to $3-6k = window sub-caching + GPU
sweeps + enumeration warm-tracking, all mechanical. Next: collar
validation.

## 4.26 The collar is THICKENED LAYER-3, not coarse tiles

Collar probes near the face FAIL at colors 7 with mass wilds — and
correctly so: the face carries genuine TRIPLES whose partner bases
are 6-cliques, so the Layer-1 tile claim ("no MU triple through
{I, H(beta)}" — deliberately stronger than the theorem) is FALSE
there and no tile-side tuning can fix it. The collar statement must
be the QUADRUPLE level: no fourth basis extending the strata's
triples — exactly what the Layer-3 family certificates prove along
the arcs (B-arc, Dita, Szollosi: done, gap-free). THE COLLAR =
LAYER-3 THICKENED TRANSVERSALLY: extend each triple-carrying
stratum's arc certificate by the transverse beta-directions (the
s^2 branch chart + the bnrates certified K-drift substrate are its
ingredients, both already built/measured). Anatomy of the full
theorem, now explicit:
  bulk (>= 99%):  coarse no-triple tiles (certified chain, working);
  collar (< 1%):  thickened Layer-3 no-fourth-basis certificates
                  around the triple strata (new build, ingredients
                  ready);
  seam:           the tile claim's validity boundary = where triples
                  appear = detectable by the tiles themselves
                  (they fail loudly exactly there).
This also explains the sixth census point's appetite and closes the
conceptual design of the two-regime proof.

## 4.27 Shift close: chain validated 3 points / 18 steps, pricing
##      consolidated

Multi-point chains: 3 anchors + 14 steps, ALL CERTIFIED (colors
2-3, coverage 0 everywhere); class-cached windows lifted the
hardest point from 3.2x to 4.4x; amortization 4.0-5.0x across
territories. MEASURED per-tile at 3e-3, re-anchor cadence ~5:
~45-55 s CPU today -> 7e7 bulk tiles ~ $30-50k spot CPU. Named
mechanical multipliers not yet applied: GPU anchor sweeps (the
100-130 s sweep inside each 134-160 s anchor -> ~4 s at the
measured 41x), transverse-line enum reuse, recursion sub-window
caching — realistic post-stack bulk: **$5-15k**. (The earlier
$3-6k assumed steeper window amortization than measured; numbers
move as measurements land — that is the point of measuring.)
Collar (thickened Layer-3): the remaining new build, unpriced
until designed. The bulk machine, as of this shift: architecture
closed, certified end-to-end, chained, and priced from
measurements at every stage.

## 4.28 FINAL ARCHITECTURE: triples live ONLY on the face

Existence ladder at (theta, 2.0412, pi): 48 MU vectors but ZERO
orthonormal bases among them at theta = 0.001, 0.005, 0.02, 0.05,
0.1 — triples exist only AT theta = 0 (the group-Hadamard face).
Therefore the no-triple claim is TRUE throughout the open collar;
the colors-7 probe failures were certifiability (fat drift bounds
vs near-clique margins), not truth. Thickened Layer-3 is NOT needed
for the collar. The theorem's final anatomy, all on existing
tooling:
  1. BULK: coarse no-triple tiles (chained, certified, priced).
  2. COLLAR (0 < theta <~ 0.02): dyadic theta-scaled anisotropic
     no-triple tiles (precision ~ margins, which shrink toward the
     face; log-many dyadic levels x in-face tiling).
  3. FACE (theta = 0): a 2-parameter Hadamard family = a Layer-3
     no-fourth-basis walk, exactly the (complete) Szollosi-region
     pattern; also the locus where classical LP/Delsarte applies.
Remaining measurement: margin-vs-theta scaling in the collar (sets
the dyadic tile count); remaining builds: collar tiler + face walk
— both compositions of existing machinery.

## 4.29 Collar scaling law: defect = c(beta_face) * theta — collar is log-cheap

Measured (pool of 48 polished MU vectors, min nonzero pairwise
overlap = the near-clique defect the collar must resolve):
  - theta-scan at (2.0412, pi): 0.0001 @ 0.01, 0.0002 @ 0.02,
    0.0003 @ 0.04, 0.0007 @ 0.08 — linear in theta.
  - in-face scan at fixed theta: every defect at theta=0.005 is
    EXACTLY half its theta=0.01 value across b2 in [1.94, 2.20]
    and b3 in [pi-0.2, pi+0.1] (e.g. 4.29e-4 -> 2.14e-4,
    0.81e-4 -> 0.40e-4).  So defect = c(b2,b3)*theta with smooth
    c ~ 0.005-0.075, and the IN-FACE gradient of the defect is
    itself O(theta).
Design consequence: collar tiles need theta-scaled precision only
in the theta direction (dyadic slabs, drift taxes prop. to slab
thickness); in-face resolution is theta-INDEPENDENT (set by c's
O(1)-scale variation).  Collar tile = bulk-grade coarse tile +
certified edge-deletion on the near-clique pairs (root-local
overlap lower bound  c*theta  vs taxes prop. theta — scales match
at every level).  Cost: (log2(theta_max/theta_0) ~ 2-5 levels) x
(one bulk slab), i.e. NEGLIGIBLE vs bulk.  theta_0 = tube radius
of the face-family certificate (Layer-3 walk on the theta=0
2-param group-Hadamard face — the Jaming-Matolcsi-Mora-Szollosi
locus, which the literature already excludes pointwise; we need
the tube version via K-drift).
Residual risk: c has an interior minimum (~0.005 near b2=2.10 at
b3=pi); if c -> 0 anywhere on the face that point is a deeper
stratum (codim-2) needing its own scaled chart.  Next: map c over
the face.

## 4.30 Face obstruction measured: pair-level, margin ~ 1/6 — face walk is trivial

At theta = 1e-6 (face proxy; theta=0 is the branch-chart
singularity), both the F-point (2.0412, pi) and a generic face
point (1.0, 2.0): 48 MU vectors containing exactly 8 orthonormal
bases, and the maximum mutually-unbiased clique among those bases
is 1 — NO pair of extra bases is mutually unbiased.  Since any 4th
MUB alongside {I, H} requires TWO pairwise-unbiased extra bases,
the obstruction sits at the basis-pair level with defect margin
(max entrywise ||<u,v>|^2 - 1/6|, min over the 28 pairs):
  F-point:  min 0.0917, median 0.1666 (= 1/6, the class gap of
            4.22 — the LP bridge surfaces as the face margin!)
  generic:  min 0.1649, median 0.1662
Certificate per face tile: root-list coverage (existing fat-sweep
machinery) + 8-basis enumeration + 28 certified pair-defect lower
bounds.  Margins ~ 1e-1 vs drift rates O(10) -> in-face steps and
transverse tube radius theta_0 both ~ 1e-2.  Face walk ~ 1e3-1e4
cheap tiles; dyadic collar needs ~1-2 levels (theta_max 0.02 ->
theta_0 ~ 1e-2).

PRICING, COMPLETE (all three pieces now measured):
  bulk   $30-50k CPU today, $5-15k post mechanical stack;
  collar ~ 1-2 x one bulk slab (negligible);
  face   ~ 1e3-1e4 pair-defect tiles (negligible).
The theorem's cost IS the bulk cost.  No structural unknowns
remain; remaining work is engineering (collar tiler, face walker,
certified-grade pass, driver integration) and the campaign run.

## 4.31 Collar breaking-edge law + face walk demonstrated (5/5 tiles, 5 s each)

Per-basis breaking edges (the RIGHT collar margin — one certified
non-orthogonal edge kills a 6-clique; the global min-pair c of the
4.29 map was measuring irrelevant pairs):
  max-edge overlap / theta, min over near-bases, constant in theta
  over 0.0025-0.02:
    F-point   0.419 (med 0.490)   [worst measured anywhere]
    c-map-worst (1.1,1.043)  0.74-0.85
    generic (1.0,2.0)        0.76-0.85
  So EVERY near-basis clique carries an edge with overlap
  >= 0.42*theta, ~4000x the worst global-min pair; the collar
  edge-deletion certificate has fat, perfectly linear margins.
  Deeper-collar bonus: near-basis count drops 8 -> 2 beyond
  theta ~ 0.005 at generic points, and the face unbiasedness
  defect GROWS with theta (0.166 -> 0.19 at 0.02) — margins only
  improve away from the face.

Face walk prototype (facewalk.py, sampled-grade rates RATE=0.5 =
7x the measured 0.07, PAD 1e-3): per-tile = pool + 8-basis
enumeration + max-MUB-clique sanity + 28 pair-defect lower bounds
over an hf=0.05 in-face box.  5/5 probe tiles CERTIFIED (margins
0.055-0.13), ~5 s/tile => full face = (2pi/0.1)^2 ~ 4e3 tiles
~ 5.5 CPU-HOURS.  The face piece is demonstrated end-to-end at
prototype grade; certified pass = s^2-chart rates + coverage
wiring (existing machinery).

## 4.32 Collar tile demonstrated; the pi/3 branch surface is the corner stratum

collar_tile.py, two designs tried:
  v1 (center-anchored symmetric tax): 0/4 — provably wrong shape,
  the breaking overlap c*theta HALVES toward the slab bottom while
  a symmetric tax stays put.
  v2 (bottom-anchored + FD theta-monotonicity + in-face-only
  taxes, PAD 3x): slabs [0.005,0.015] and [0.0025,0.005],
  hf 5e-3 —
    F-point:        CERTIFIED, chi 7 -> 2 (both slabs, 110 edges
                    deleted)
    generic (1,2):  CERTIFIED, chi 8 -> 5
    (1.1, 0.90):    CERTIFIED, chi 6 -> 2
    (1.1, 1.20):    CERTIFIED, chi 6 -> 2
    (1.1, pi/3):    FAILED — and diagnostically so: 1102 edges,
                    chi 42, drift median 0.59 (theta-part 0.42,
                    |S_theta| ~ 40), 132 pairs below overlap 0.01;
                    root_data2 never fails, roots are isolated but
                    ILL-CONDITIONED. Splitting hf 4x barely moves
                    it (chi 9 -> 6): the fat is intrinsic.
  Dephased-entry check: NOT a group-Hadamard locus (26 distinct
  phases). Verdict: beta3 = pi/3 is an interior ROOT-DEGENERACY
  (branch) surface — near-singular Jacobians, racing roots — the
  same fold/valley anatomy the bulk's valley-window machinery
  already handles. The collar x branch-surface corner stratum
  needs the existing valley windows composed into the collar tile
  (engineering, not new theory), mirroring the ~1-2% branch-aware
  chart class long budgeted for the bulk.
  ~5-6 s/tile at demo grade, in line with the face walk.

## 4.33 Collar census: 36/48; the ONLY failures are the beta3 = +-pi/3 surfaces

48-point demo-grade census (6 x 8 grid over the face, slab
[0.005, 0.015], hf 5e-3): 36 CERTIFIED, and ALL 12 failures sit at
beta3 in {1.043, 5.257} — the grid points nearest pi/3 and 5pi/3 —
at EVERY beta2. The degenerate set is exactly the mirror pair of
branch surfaces beta3 = +-pi/3 (beta2-independent); no other
special loci at this resolution. Adaptive closure measured: the
cell at pi/3 + 0.02 fails at hf 5e-3 but CERTIFIES at hf 2.5e-3
(chi -> 4); so the residual corner stratum is the band
|beta3 -+ pi/3| < ~0.02, to be handled by composing the bulk's
valley-window machinery into the collar tile (racing-root
territory, 4.32). Generic collar cost confirmed: ~5-6 s/tile at
demo grade, uniformly across the face.

## 4.34 The band closes: signed analytic pair rates + theta-scaled b3-thin boxes

The pi/3 band arc, in order of discovery (all measured today):
  (a) NO exact bases anywhere on the band at theta > 0 (ladder at
      (th, 1.1, pi/3), th = 0.005-0.1: 0 bases at tol 1e-7) — the
      no-triple claim is true on the walls too.
  (b) The breaking-edge law HOLDS ON THE WALL: min per-near-basis
      max-edge = (0.34-0.38) * theta, linear over th 0.005-0.04;
      off-wall max-edge ~ 0.84 * theta INDEPENDENT of distance.
      Margins never degenerate; only the certificates were failing.
  (c) Root anatomy on the wall: roots widely spaced (min 1.24) but
      column-resolved sensitivities show ONLY the wall-normal b3
      races: |S_b3| = 10.8 / theta (clean law: 2180@0.005,
      1087@0.01, 537@0.02, 258@0.04), while S_theta ~ 4.7 and
      S_b2 ~ 0.4-1.8 stay tame. root_data2's default continuation
      delta 2.5e-4 branch-jumps for ALL 48 roots (motion ~ |S| *
      delta overshoots); delta <= 1e-5 (adaptive ladder 1e-5,
      2e-6, 5e-7) recovers S for 48/48.
  (d) v3 = signed ANALYTIC pair rates from S (no FD): d<u,v>/db_l
      via phase-sensitivity differences, signed d|O|/db_l via the
      inner-product phase; theta-term DROPPED when the signed
      derivative certifies growth (bottom-anchored, the v2 trick);
      unsigned in-face taxes with b3-thin boxes hf3 ~ 0.1 *
      theta_lo (the 1/theta racing law sets the width).
  Result: wall-exact cells CERTIFY — chi 3 at (1.1, pi/3) and
  (3.0, 5pi/3) on [0.005,0.015], chi 2 on the deep slab
  [0.0025,0.005]; the marginal offset cell (pi/3 + 0.02) certifies
  at chi 2 with hf 5e-3 (blanket needed 2.5e-3); generic anchors
  certify sharper than blanket (chi 2, no deletions needed —
  certified-positive pairs never enter the adjacency).
  Band pricing: b3-strip 0.04 wide at hf3 5e-4 -> ~40 b3-tiles x
  ~630 b2-columns x 2-3 slabs ~ 6-8e4 tiles x 6 s ~ 4-5 CPU-days,
  embarrassingly parallel — negligible.
  NO valley-window composition needed. The signed-rate structure
  IS the certified-grade design (dual-AD first-order pair data +
  second-order remainders replace the sampled FD/S continuation).

## 4.35 48/48: the collar census closes completely

The two residual census cells (b3 = 1.043, b2 = 0.3 and 1.92)
were not margin failures:
  (0.3, 1.043):  closes with hf3 = 2.5e-4 (chi 5).
  (1.92, 1.043): breaking edges were fine all along (o ~ 0.0025,
      tax ~ 4e-5, theta-growth certified) and correctly deleted —
      the residual chi 6 was the GREEDY COLORING being suboptimal
      on the thinned 108-edge graph. Adding a DSATUR upper bound
      (soundness needs any proper coloring, so min(greedy, DSATUR)
      is valid) gives chi 5. Wall-exact and generic regressions
      unchanged (chi 2).
With the final configuration (signed analytic rates, adaptive
delta ladder, hf 5e-3, hf3 2.5e-4, greedy+DSATUR): the 48-point
collar census certifies 48/48 — no open cells anywhere in the
collar, band and walls included, ~6 s/tile at demo grade.

## 4.36 Seam-complete: face theta-tubes meet the collar (10/10)

facewalk.py v2 adds the theta-tube term (RATE_TH = 1.0 = the
measured theta-drift of the unbias defect ~0.01, padded 100x —
the defect is flat in theta). 10/10 face tiles CERTIFY with
tubes 0.0025 and 0.005 at hf 0.05, margins 0.050-0.126,
INCLUDING both wall points (1.1, pi/3) and (2.3, 5pi/3).
The domain decomposition is now demonstrated with OVERLAP at
every seam:
  face  theta in [0, 0.005]   (tube certificates, 10/10)
  collar theta in [0.0025, 0.02] (dyadic slabs, census 48/48,
         deepest slab [0.0025, 0.005] certified)
  bulk  theta > 0.02          (coarse tiles, 6/6 + chains)
Every piece runs on demonstrated machinery at ~5-6 s/tile; the
face-collar handoff is double-covered on [0.0025, 0.005].

## 4.37 Collar chains 4-5x; the (pi/3, pi/3) corner is the last special object

collar_chain.py: warm pools (re-polish previous roots + 500-start
top-up; count change falls back loudly to full enumeration).
Measured on [0.005, 0.015], 20 tiles per line:
  generic line (b3=2.0):  20/20 certified, 0 fallbacks,
                          1.23 s/tile (4.9x vs 6 s cold)
  wall line (b3=pi/3):    18/20, 1 fallback, 1.51 s/tile
Collar campaign re-price at 1.2-1.5 s/tile: band ~6-8e4 tiles ->
~1 CPU-day; full collar << 1 CPU-week single-core. Negligible.

The 2 wall-line failures are NOT warm-pool artifacts (identical
cold): they sit at b2 = 1.04-1.05 ~ pi/3 — the CORNER where a
b2 = pi/3 locus crosses the b3 = pi/3 wall. Probes:
  b2 = pi/3 OFF-wall (b3 = 2.0, 5.0): CERTIFIED chi 2 — the b2
      surface is not special by itself;
  corner cells (~pi/3, pi/3): FAIL box-size-independently even at
      hf 5e-4 x hf3 2.5e-4 (792-878 edges, chi 14-23, 10/48 roots
      gated at delta 5e-7) — deeper-than-1/theta degeneracy,
      codim-3 point stratum (a point in the face x the theta
      interval). Characterization running; measure-zero, bespoke
      treatment affordable whatever it is.

## 4.38 Corner law: |S_inface| = (2.2, 1.1)/theta^2, S_theta ~ 1

Corner characterization at (theta, pi/3, pi/3):
  - NOT a group-Hadamard locus (dephased dev 0.5);
  - NO exact bases at theta = 0.005-0.05 — no-triple TRUE at the
    corner as everywhere off the theta=0 face;
  - racing law: med |S_b2| = 2.17/theta^2, |S_b3| = 1.09/theta^2
    (three-digit fit over theta 0.01-0.04), while S_THETA ~ 1 —
    supertame in theta;
  - root_data2 recovers 48/48 with the adaptive delta ladder at
    theta >= 0.01.
Treatment: the corner is a self-similar codim-3 funnel walked
DOWN in theta (S_theta ~ 1 makes theta-steps O(h) cheap): per
dyadic level, a theta^2-scaled core box + a ring handing off to
wall-scaling (1/theta) boxes. Cost hinges on the basin constant
(where wall-scaling resumes) — being measured. The corner is a
POINT in the face: even a generous ring count prices in CPU-hours.

## 4.39 COMPLETE: every stratum of the domain certifies at demo grade

Basin interior follows the theta^2 law exactly: d = 0.01 certifies
at hf 1e-4, d = 0.005 at hf 2e-5, and the EXACT corner column
(pi/3, pi/3) at hf 2e-6 (chi 2, 52 edges) — the width ladder
matches |S_inface| ~ 2.2/theta^2 with the same 0.1-margin rule
used everywhere else. Basin edge measured at ~theta_hi from the
corner; ring counts O(100)/level; whole corner ~ CPU-hours.

STATUS OF THE DEMONSTRATION TIER (2026-07-31, one day):
  bulk          6/6 census + chains (amortized ~3x)      [prior]
  generic collar 48/48 census, chains 1.23 s/tile (4.9x) [today]
  walls +-pi/3  1/theta law, chi 2-3, in-census          [today]
  corner basin  theta^2 law, chi 2                        [today]
  exact corner  hf 2e-6, chi 2                            [today]
  face + tubes  10/10 incl. walls, seam double-covered    [today]
No open cells. No unpriced pieces. The three-piece architecture
of 4.28 survived contact with every special locus it met; each
special object yielded to a measured scaling law + the same
signed-rate tile. Remaining work is now genuinely mechanical:
certified-grade pass (dual-AD pair rates + 2nd-order remainders,
interval enumeration coverage), driver integration, and the bulk
campaign — the only expensive item ($5-50k depending on the
mechanical multipliers).

## 4.40 Certified signed pair rates: enclosures 0.05-1% wide, FD inside

certpair.py — the certified-grade path for the v3 primitive:
Krawczyk root enclosures (existing) -> J at center + HESS_ROW
ball pad -> dg/dbeta from dual_karlsson certified partials
(theta-ball dependence PAD-sampled, the one prototype-grade link)
-> S enclosure via midpoint-inverse residual bound -> pair-rate
interval. Measured at the breaking edges:
  generic (0.005, 1.0, 2.0):  rate (0.0802, 0.0003, 0.0484)
      +- 4.3e-5  (0.05% fat), FD reference inside;
  wall (0.005, 1.1, pi/3):    rate (0.0660, 0.00026, 0.4519)
      +- 4.3e-3  (~1% fat), FD inside — the wall's per-root
      |S| ~ 2000 does NOT poison the enclosure because the
      Krawczyk balls are ~1e-12 and kappa stays ~1.
Conclusion: the certified collar pass costs essentially nothing
in tax fatness. Remaining certified-grade links: the theta-ball
Lipschitz of dg/dbeta (replace sampled PAD with the interval
Hessian machinery) and second-order pair-rate remainders over the
box — both standard moves in the existing substrate.

## 4.41 Second-order honesty: wall beta3-curvature = 1/theta^2, hf3 ~ theta^1.5

Measured pair-overlap curvatures at the breaking edges (second
differences through re-polished roots):
  generic: curv_theta ~ 5.5, in-face ~ 0 — h^2 terms 1e-4-1e-10,
      negligible at demo widths;
  wall:    curv_theta ~ 2.2, curv_b2 ~ 0, curv_b3 ~ 4.7e4 ~
      1/theta_lo^2 — the h^2 term at the demo width hf3 = 2.5e-4
      is 2.95e-3 > the breaking overlap 1.95e-3(!). The demo
      census's first-order taxes were silently benefiting from
      the omitted remainder on the wall strip.
Honest widths: hf3 <~ sqrt(0.4 theta / (3 curv)) ~ 0.4 theta^1.5
on the wall (1e-4 at theta_lo 0.005). Verified: wall tiles at
hf3 = 1e-4 CERTIFY (chi 2, 48-64 edges) with room for the
explicit curvature charge. Cost impact: ~2.5x more beta3-tiles on
the wall strip only — pricing unchanged in substance (collar
still ~ CPU-days). Remaining certified-grade links after this:
interval Hessian for the theta-ball Lipschitz of dg/dbeta, and
wiring the curvature charge into the tile's tax line — both
existing-substrate moves.

## 4.42 Driver lesson: the corner core is 2-D theta^2-fine (~$100s parallel, not CPU-hours)

collar_driver.py walked the hardest line (b2 = pi/3 exactly,
b3 across the corner +-3e-3, slab [0.005, 0.015]) with an
auto-laddering rung schedule. Result: 200/1501 positions
certified in 72 min — every success at the finest rung, and the
ladder's bottom (hf = 1e-4, hf3 = 2e-6) was WRONG for this line:
on b2 = pi/3 the funnel needs theta^2 widths in BOTH in-face
coordinates (|S_b2| = 2.2/theta^2 is the LARGER rate; the earlier
corner success 4.39 used 2e-6 x 2e-6). Rung added.
CORRECTED corner pricing (4.38's "CPU-hours" was wrong): the 2-D
core |beta - (pi/3, pi/3)| <~ few 1e-3 at widths ~ 0.05 theta^2
is ~ 5e6-1e7 tiles ~ CPU-MONTHS single-core = ~$100-500
embarrassingly parallel at cloud CPU rates — still negligible
against the bulk, but three orders bigger than the old estimate.
Driver fixes queued: warm-start the rung ladder at the previous
position's rung (the 4-attempts-per-position tax dominated the
72 min), and treat the corner core as its own scheduled patch
rather than line-walking through it.

## 4.43 Driver validated on the wall crossing: 13/13 at 2.2 s/tile

With rung warm-start (ladder resumes one rung above the previous
success), the b2 = 1.1 wall crossing (b3 = pi/3 +- 3e-3, slab
[0.005, 0.015]) drives 13/13 CERTIFIED in 29 s — all at the
coarsest rung (5e-3, 2.5e-4), no downshifts needed. The adaptive
driver is the collar campaign's shape: warm pools + rung ladder +
loud failures + JSONL ledger. The corner core remains a separate
scheduled patch (4.42). Note the certified-grade pass will charge
the 4.41 curvature term and push wall cells one rung down — which
is exactly what the ladder is for.

## 4.44 Second-order charge wired in; dyadic slabs + theta^1.5 widths verify at all depths

collar_tile now charges t2 = curv0*span^2 + curv1*hf^2 +
curv2*hf3^2 in the signed tax line, with per-stratum constants
(CURV_GENERIC = (10,10,10); walls: slab-scaled curv3 = 4.8 /
theta_lo^2 = 4x the measured 1.2/theta^2). Verified with the
charge ON:
  generic anchor:      chi 4, hf3 2.5e-4              CERTIFIED
  wall [0.0025,0.005]: chi 4, hf3 2e-5                CERTIFIED
  wall [0.005,0.01]:   chi 4, hf3 1e-4 (both walls)   CERTIFIED
  wall [0.01,0.02]:    chi 4, hf3 2e-4                CERTIFIED
Two structural consequences, both already in the design's shape:
slabs must be DYADIC [theta, 2theta] (the span^2 charge kills the
2:3 slab [0.005,0.015] at the wall), and wall widths follow
hf3 ~ 0.2 theta^1.5. The demo tier is now honest through second
order everywhere except the corner core (its curvature constant is
unmeasured; its widths 2e-6 make any plausible charge ~1e-12-1e-6
— headroom, but measure before the campaign).

## 4.45 Corner core second-order safe; 4.42 pricing stands

Breaking-edge derivatives AT the corner (0.005, pi/3, pi/3),
o0 = 2.28e-3:
  theta: rate 0.456, curv 0.05      — supertame, growth law holds
  b2:    rate 81.5,  curv 6.0e3
  b3:    rate 43.4,  curv 7.9e6
The pair-CORRELATED rates are ~1000x below the per-root 1/theta^2
rates (8.8e4) — the signed cancellation persists at the corner.
At hf = 2e-6: first-order taxes 7.5e-4 (dominant), second-order
charges 2.4e-8 + 3.2e-5 — 70x headroom. Corner widths are
FIRST-order-limited at ~5e-6; the demonstrated 2e-6 grid and the
4.42 pricing (~$100-500 parallel) stand with second order charged.
All curvature constants in the collar are now measured; nothing
in the demo tier rests on an unmeasured second-order term.

## 4.46 HONEST CORRECTION: the corner core has a fractal tail — one genuinely open piece

8x8 corner patch at hf 2e-6 (slab [0.005, 0.01], measured-curv
charged): 41/64 at 1.61 s/tile — but the 23 failures trace CURVED
1-D sub-loci through the patch (a full b2-row, a b3-arc), and the
rung-down to 5e-7 closes only 1/5 sampled failures. Worse, on the
sub-loci root_data2's continuation itself breaks (36/48 roots
gated at one cell -> complete-graph chi 37). Interpretation: the
funnel has self-similar sub-structure — the breaking edge changes
identity across these curves and both candidates are small at the
crossings; each specialization level reveals thinner loci.
POINTWISE TILING HAS NOT BEEN SHOWN TO TERMINATE at the corner
core below ~1e-5.
Status revision of 4.39: "no open cells" holds everywhere EXCEPT
the corner core (pi/3, pi/3) at scales < ~1e-5 x theta in
(0, 0.02] — codim-3, measure-zero, but a proof needs it. Candidate
treatments, in order of plausibility:
  (i)  family-level certificates over the sub-locus curves (the
       valley/fold + s^2-chart machinery — built, unwired here);
  (ii) exact algebraic analysis at beta = (pi/3, pi/3): a highly
       special Hadamard slice, plausibly tractable by hand or by
       the Layer-3 pointwise methods that did S6/Bjorck's C;
  (iii) a deeper measured scaling law organizing the sub-loci.
This is the machine finding its own next problem — the honest
boundary of today's demonstration tier.

## 4.47 Corner treatment plan: a 12th-root arc + walked sub-curves

The corner slice H(theta, pi/3, pi/3) dephases to a 1-parameter
deformation of an EXACT 12th-roots-of-unity Hadamard: row 1 is
exactly (0, pi, pi/3, -2pi/3, pi/6, -5pi/6) and theta-INDEPENDENT;
rows 2-5 carry the theta-drift in symmetric pairs; quad-product
phase invariants sit on the pi/6 lattice. The theta -> 0 limit is
a specific 12th-root matrix (identify against the catalog — Dita
slice? — worth an hour).
Treatment (the concrete version of 4.46 option (i)+(ii)):
  1. Walk the corner COLUMN as a 1-param arc (the Layer-3
     pattern: per-theta anchors, S_theta ~ 1 is supertame,
     certified transverse tube). Measured tube radius from the
     breaking edge: r ~ 0.456 theta / (3 * 124) ~ 1.2e-3 * theta
     — at theta = 0.005 that is 6e-6, JUST short of the ~1e-5
     funnel sub-loci.
  2. The sub-loci are themselves curves (per 4.46 they are where
     the breaking edge changes identity): walk EACH as its own
     arc with its own margins — the same machinery one level in.
     If the curve census is finite per level and margins hold,
     the tail terminates; measure before building.
  3. Fallback: exact analysis of the 12th-root limit matrix
     (pointwise Layer-3 did S6 and Bjorck's C — same weight
     class).
Bounded, concrete, all on existing patterns — tomorrow's first
problem.

## 4.48 The tail decodes: THREE STRAIGHT LINES, not a fractal

16x16 map of the corner core at hf 2e-6 (slab [0.005, 0.01],
209/256 certified, 1.6 s/tile; rows = b2, cols = b3):

  ##########.#####
  ##########.##.##
  ##########.#####
  ##########.#.###
  ##########...###
  #########...####
  ................
  ##########...###
  #########..#####
  ########.#.#####
  #######.##.#####
  ######.###.#####
  #####.####.#####
  ####.#####.#####
  ###.######.#####
  ##.#######.#####

The fail set is exactly three straight lines:
  1. b2 = pi/3 - 6e-6      (theta-displaced b2-surface, full row)
  2. b3 = pi/3 + 1e-5      (theta-displaced wall, full column)
  3. b2 + b3 = 2pi/3 + 8e-6 (ANTI-DIAGONAL — a third resonance
     locus, previously unknown; displacements all ~ theta-scale)
Also: the theta = 0 anchor (pi/3, pi/3) face point has 0 exact
bases among 48 MU vectors (continua signature — the face walk
needs a special chart there; naive F6(a,b) phase-multiset fit
fails at loss 0.93, identification open).
REVISED VERDICT on 4.46: the tail is STRUCTURED, not fractal —
three smooth codim-1 curve families per level, crossings = next-
level corners (same shape recursively, displacements shrinking
with theta). Treatment: walk the three line-families with the arc
pattern (each has its own breaking margins), recurse on the
finitely many crossings; the countable tree with geometrically
shrinking measure is compactness-argument shaped. The corner
remains the one open piece, but it is now a REGULAR object with a
concrete finite-looking attack, not a pathology.

Addendum (final probe of the shift): the breaking-edge law holds
ON all three lines — min per-basis max-edge = 0.44 theta
(b2-line), 0.81 theta (b3-line), 0.63 theta (diagonal), 8
near-bases each. As at every level of this hierarchy, the margins
never degenerate; only the axis-aligned pointwise taxes do. The
line walk needs line-adapted boxes (thin transverse to each
line), and the margins to pay for them are measured and fat.

## 4.49 The three lines certify pointwise; open set shrinks to crossings + anchor

Sequence of discoveries closing 4.48's lines:
  (a) COLORING ARTIFACT: a diagonal-line cell's residual graph
      measured chi <= 4 by randomized-order coloring where greedy
      and DSATUR both said 6 (charge-free config). Randomized
      restarts (300 shuffles when >5) added to the tile's chi().
      Under the FULL charged tax the lines persist (map2: 210/256)
      — the artifact was real but not the whole story.
  (b) SLOP FLOOR: the surviving line edges sat below the fixed
      1e-3 discretization slop (polished roots are 1e-12-accurate;
      1e-3 was paranoia). Slop is now a parameter.
  (c) ALL THREE LINES CERTIFY with refined parameters, chi 2:
      diagonal:  slop 1e-5, boxes 5e-7, slab [0.005, 0.006]
      b2-line:   slop 1e-5, hf 5e-7, hf3 2e-7, slab [0.005,0.0055]
      (b3-line untested at final config; displaced-wall symmetry
      expected — confirm when wiring the line schedule.)
EVERY cell of every stratum tested today ultimately certifies at
demo grade. The remaining open set: the line-CROSSING points
(next-level corners — same recursion, thinner parameters,
expected to fall the same way), the theta = 0 corner face anchor
(MU continua — needs its own chart), and the THEOREM-side
question: a compactness/termination argument that the recursion
bottoms out in finitely many levels (the displacements shrink
with theta; margins never degenerate — the shape is there).

## 4.50 At the crossing, only the certified path survives — 48/48 enclosed

The line-crossing point (pi/3 - 6e-6, pi/3 + 1e-5): root_data2's
FD continuation gates out 48/48 roots at every delta in the
ladder — the sampled machinery is DEAD there (|S| ~ 2.7e5,
another 1/theta factor deep). But certpair.certified_S (Krawczyk
ball + analytic J + dual-AD dg/dbeta + residual bound — no FD
anywhere) encloses ALL 48 roots with ~0.3% relative width
(err ~ 6e2-9.5e2 on |S| ~ 2.7e5). The crossing tiles are a
certified-grade job by necessity, not just by preference: the
deepest stratum of the domain is reachable ONLY by the rigorous
substrate. A fitting place for the demonstration tier to end.

## 4.51 THE CROSSING CERTIFIES — every tested cell in the domain closes

collar_tile signed mode gained a certified-S fallback (certpair's
Krawczyk enclosure path) for roots where FD continuation is dead.
The line-crossing cell (pi/3 - 6e-6, pi/3 + 1e-5) — the deepest
stratum found today, where root_data2 gates out 48/48 —
CERTIFIES at the FIRST rung: chi 5, 92 edges, all 48 roots on
certified-S. Demo-grade caveat: the S enclosure err is not yet
charged into the tax (the certified pass will fold it into slop).
With this, EVERY cell tested today, at every stratum:
  bulk / generic collar / walls / basin / corner column / all
  three lines / the crossing
certifies at demo grade. The remaining open set is not
compute-shaped but paper-shaped:
  1. the theta = 0 corner face anchor (MU continua — a special
     chart or exact analysis of the 12th-root limit matrix);
  2. the termination lemma for the corner recursion (each level:
     3 lines + crossings, displacements ~ theta-scale, margins
     0.4-0.8 theta never degenerating — compactness-argument
     shaped).
The day ends with the sampled machinery dead at the bottom of the
funnel and the certified machinery — built this afternoon —
carrying the last cell. The substrate is the story.

Addendum to 4.51 (the honest last word): charging the certified-S
enclosure error into the tax (te = (err_i + err_j) * 5/6 *
(span + hf + hf3), now in the code, zero for FD-sourced roots)
FAILS the crossing (complete graph — the err*span term ~0.85
swamps everything) and the b3-line cell (42/48 on fallback,
chi 7). The 4.51 crossing certification is therefore demo-grade
with a REAL open caveat, not a formality: closing it needs
(a) error-robust growth testing (dO[0] - err > 0 gates the
span-tax drop; error must not multiply span otherwise), and
(b) tighter enclosures — err ~900 is kappa * wid with kappa =
|J0^-1| large at the crossing; the wid side (dual-AD dgdb width
+ prototype theta-ball pad) is where to push. Both are bounded
certified-substrate work. Every OTHER stratum's certification is
unaffected (FD-sourced roots carry zero err charge).

## 4.52 Overnight (1)+(2): cone-limit margin measured POSITIVE; anchor is precision-, not structure-limited

(2) TERMINATION LEMMA GROUNDWORK (PROOF_SKELETON §6 drafted —
cone compactification: rescale (s,t) = (beta - c*)/theta; the
corner region becomes a compact cylinder K; walls/lines/crossings
are FIXED subsets of K, so the "recursion" is one finite cover):
Claim C1 measured: the rescaled breaking-edge margin mu(s,t,theta)
is theta-INDEPENDENT to 3 digits along every tested ray and
radius (theta 0.01 -> 0.0025):
  s-axis 0.32-0.38, t-axis 0.38-0.45, antidiag 0.38-0.43,
  diag 0.38-0.71, generic 0.38-0.80.
Measured uniform floor mu_0 >= 0.32. The cone-limit function
exists and is positive — C1 is empirically true, and the
termination argument stands on it.

(1) ANCHOR:
  - Identification: NOT F6 and not the tested F6(a,b) points
    (permutation-invariant Haagerup multiset differs: the corner
    limit has phases at odd pi/6-multiples shifted +-0.0326 —
    F6's are pure pi/3-lattice). The anchor is its own object;
    first dephased row is exact 12th roots, full invariant is not.
  - The apparent "root collapse" below theta ~ 2e-3 on the corner
    column (LM finds nothing under 1e-10) is a DOUBLE-PRECISION
    ARTIFACT: the residual floor tracks eps_machine * cond(J)
    (cond ~ 1/theta^2: floor 8.5e-11 at 2e-3, 9e-10 at 5e-4).
    High-precision (mpmath dps=40) Newton verification running:
    if roots converge to 1e-32 residual with positive min-sv(J),
    the column has isolated ill-conditioned roots all the way
    down (and the earlier "continua" reading at the face anchor
    needs the same re-examination).

## 4.53 THE DEEP PATHOLOGY WAS FLOAT NOISE: map error = eps/theta^3 at the corner

Measured float-map error vs a 40-digit mpmath port (mp_karlsson):
  corner column: 5.1e-11 @ th=0.01, 7.6e-8 @ 1e-3, 5.4e-5 @ 1e-4,
      1.4e-2 @ 1e-5 — a clean eps_machine/theta^3 law (third-order
      catastrophic cancellation in the Moebius quotients at the
      corner).
  OFF-corner it is tiny: generic/wall/F-point <= 9e-10 even at
      th = 1e-6, <= 3.3e-8 at 1e-9.
Consequences, in order of importance:
  1. ALL slab-level demo work STANDS (th >= 0.0025: error <= 4e-9,
     orders below every margin and tax used).
  2. facewalk's TH_PROXY = 1e-6 is SAFE at generic face points
     (error <= 9e-10) — but its (pi/3, pi/3) reading ("0 bases",
     "continua") was NOISE, as were all anchor probes at
     th <= 1e-6 (the +-0.0326 Haagerup jitter, the "empty pool",
     the LM residual floors ~ eps*cond, the root_data2 gating on
     the deep column).
  3. On the TRUE map (mp Newton, dps 40): the corner column at
     th = 1e-4 has ISOLATED, WELL-CONDITIONED roots — 12 found,
     min-sv(J) = 0.05-0.09, residuals < 1e-32. No continua, no
     intrinsic ill-conditioning. The deep funnel is GENERIC-grade
     geometry evaluated through a numerically bad chart.
  4. TRUE corner limit (mp, th = 1e-12): an EXACT pi/6-lattice
     Hadamard — Haagerup phases {0, pi/6, pi/3, pi/2, 2pi/3,
     5pi/6, pi} — and NOT F6 (pure pi/3 lattice). Catalog
     comparison vs Dita/Bjorck/S6 in flight.
  5. The certified campaign needs the corner region evaluated via
     the interval map (ivkarlsson, which BOUNDS its own error) or
     a cancellation-free algebraic rewrite of the Moebius
     quotients near (pi/3, pi/3); the slab-level racing laws
     (1/theta, 1/theta^2 at th >= 0.0025) are chart-real and
     their treatments stand.

## 4.54 Item (1) RESOLVED: the corner anchor is a Fourier point with ordinary face structure

Three converging measurements on clean arithmetic:
  1. IDENTIFICATION: the true corner limit's FULL Haagerup
     multiset matches the Fourier family F6(a,b) at 27
     twelfth-root parameter points (e.g. (a,b) = (pi/6, 0),
     (0, pi/6), (pi/6, pi/6), ...). Strong-evidence grade (the
     Haagerup invariant is necessary, not sufficient; direct
     permutation equivalence is a bounded follow-up). This also
     explains the face's universal 8-basis structure: the theta=0
     face IS the Fourier locus, and the corner is its
     12th-root-parameter point.
  2. STRUCTURE (safe map, th = 3e-3 and 5e-3, float error
     ~2e-9): 48 vectors, EXACTLY 8 near-bases, max extra-MUB
     clique 1, unbias-defect min 0.1007/0.1012, median exactly
     1/6 — the corner anchor's margins are BETTER than the
     F-point's 0.0917. Face-ordinary in every measured respect.
  3. The earlier "0 bases / continua / empty pool" anchor
     readings were the eps/theta^3 float-map noise (4.53).
CONSEQUENCE: no special chart is needed at the corner anchor.
The face walk covers (pi/3, pi/3) with the STANDARD pair-defect
certificate at margin ~0.10 — only the evaluation route changes
near the corner (proxy at th >= 3e-3 with the theta-tube
covering [0, 3e-3], or ivkarlsson/mp evaluation; tube tax
RATE_TH * 3e-3 = 0.003 against margin 0.10 — trivial headroom).
With (1) resolved, the termination lemma (2) stands on: C1
measured (mu_0 >= 0.32, 4.52), the cylinder's theta = 0 face
covered by the face certificate (margin >= 0.10 at the corner
point itself), and the finite-cover conclusion of SKELETON §6.
Remaining for (2): certified-grade continuity of the rescaled
family (s^2-chart / interval map — substrate work, not geometry).

## 4.55 Identification correction and honest final state

CORRECTION to 4.54 item 1: the "27 F6(a,b) matches" used an
INVALID hand-rolled F6(a,b) (row-varying column phases — not even
Hadamard). Redone with the repo's true fourier_family(x1, x2):
  - the corner limit's full Haagerup multiset matches the TRUE
    Fourier family at 216 lattice points (all x2 at odd twelfths,
    etc.) — the invariant is constant along large family chunks
    and does not pin a point;
  - direct permutation-dephasing equivalence FAILS for the 4
    nearest candidates (F(0,1/12), F(0,1/4), F(1/12,0),
    F(1/12,1/12)) and their transposes.
Honest identification state: the anchor shares the Fourier
family's Haagerup invariant exactly; it is either an untested
Fourier parameter point or a non-Fourier Hadamard with the same
invariant. NOT LOAD-BEARING either way: the anchor's measured
no-4th-basis structure (4.54 item 2 — 8 near-bases, clique 1,
margin 0.101, median exactly 1/6) is what the face certificate
uses, and it is face-ordinary. The identification is a paper
footnote, not a proof dependency.

## 4.56 The interval map inherits eps/theta^3 at the corner — fix is bounded

iv_karlsson point-box widths at the corner: 7.0e-8 @ th=0.01,
7.1e-5 @ 1e-3, 6.8e-2 @ 1e-4, and a LOUD RuntimeError ("z2
denominator interval touches 0") at 1e-5. Diagnostics show the
z3 Moebius denominator margin collapsing ~ theta^4 (7.5e-9 ->
7.5e-13 -> 7.5e-17): at the corner both numerator and denominator
of the Moebius quotients vanish at the same order, so the naive
quotient loses eps/theta^3 in ANY fixed-precision arithmetic —
but the quotient itself is REGULAR there. Off-corner the interval
map is fine (3.9e-10 generic, 4.8e-9 wall at th = 1e-4).
Certified-campaign consequence (bounded work, two options):
  (a) evaluate corner-basin tiles in mpmath.iv at elevated dps
      (the substrate already uses mpmath.iv for transcendentals);
  (b) algebraic factorization: cancel the common theta-degeneracy
      from numerator and denominator once, exactly, and evaluate
      the regular quotient — the clean permanent fix.
Either way the failure mode is LOUD (denominator-touches-zero
raise), never silent — consistent with the machine's design
philosophy.

## 4.57 Stable Karlsson kernel: exact factorization lands (~4300x at the corner)

karlsson_stable.py implements the derivation:
  num3 = (A11 - z1 A12) * P,  den3 = (cA12 - z1 cA11) * Q,
  P/Q = 2 e^{i lam/2} [ sin(lam/2 - pi/6)
                        - sqrt3 sin^2(theta/2) sin(lam/2)
                        +- i (sqrt3/2) sin(theta) cos(lam/2+phi) ]
— each bracket term individually small and cancellation-free, and
the special loci fall out ALGEBRAICALLY: sin(lam/2 - pi/6) = 0 is
the wall lam = pi/3; cos(lam/2 + phi) = 0 is the 2phi + lam = pi
locus through the corner; sin^2(theta/2) is the even collapse.
Verified vs 40-digit reference:
  z3^2 kernel: relative error 8.9e-12 @ th=1e-4 (naive 1.5e-8),
      8.9e-10 @ 1e-6 (naive 1.3e-4), 1.1e-8 @ 1e-8 (naive O(1));
  full map at the corner: 1.2e-8 @ 1e-4 (naive 5.4e-5),
      1.2e-6 @ 1e-5 (naive 1.4e-2) — ~4300x; residual eps/theta^2
      comes from the z2 chain (its exact-delta form = documented
      completion); off-corner unchanged (~1e-15).
The same factorization applied inside ivkarlsson removes the
interval blowup of 4.56 — bounded, mechanical.
(Note for the record: a first "validation" showing uniform 5.9e-1
errors was a forgotten /sqrt6 in the test harness's mp port, not
a map bug — caught by entrywise diff.)

Addendum 4.57 (certified half): ivstable.py ports the factored
kernel to the interval substrate. Corner point-box enclosure
widths: 4.6e-10 @ th=0.01 (naive map 7.0e-8), 4.6e-6 @ 1e-4
(naive 6.8e-2), and NO RAISE down to 1e-6 (naive raises at 1e-5).
~150x tighter at every depth; residual width is the 8-ulp pad on
the pi/6 constant (house style), irrelevant at campaign box
widths. Remaining to wire: use iv_stable_z3sq inside iv_karlsson
(and the z2 chain via the same z3sq), mechanical.

## 4.58 Face walk (full domain) in flight; TWO NEW WALLS surface

facewalk_run.py: 3969 tiles (hf 0.05, tube 0.005), 7 workers.
Interim (first 200 tiles): 188 certified (margins >= 0.0591), and
the 12 failures are exactly the tiles containing b3 in
{pi/3, 2pi/3, 4pi/3, 5pi/3} — the wall set, INCLUDING two walls
(2pi/3, 4pi/3) the collar census grid never sampled (its b3 grid
missed them by ~0.2). Failure mode: at proxy theta = 1e-6 the
wall racing (|S_b3| ~ 10.8/theta ~ 1e7) breaks the enumeration
(2-6 bases found instead of 8) — same cure as the corner tiles:
safe proxy theta = 3e-3. Refinement pass over the ledger queued
for when the walk lands. By the symmetry group the full wall set
b3 = k pi/3 (k = 1, 2, 4, 5) was expected — good that the
machine surfaces it loudly on its own.

Addendum 4.57 (hybrid completes the float map): stable_karlsson
now routes z2 through an mp(dps 30) fallback when |den2| < 1e-3
(the corner basin), pending the exact next-order expansion of
N = zeta*bigden*Q - bignum*P. Full-map accuracy at the corner:
1.8e-12 @ th=1e-4 (naive 5.4e-5), 1.8e-10 @ 1e-6, 4.0e-9 @ 1e-8
— accurate at EVERY depth. Cost: 0.25 ms/call in the basin vs
0.031 ms generic (basin tiles only). The float-map story is
closed; remaining: wire iv_stable_z3sq + an interval z2 treatment
into iv_karlsson for the certified side.

## 4.59 FULL FACE EXECUTED (pass 1): 3777/3969, min margin 0.0496

facewalk_run.py + facewalk_resume.py (restart-safe via the JSONL
ledger): all 3969 tiles of the (b2, b3) face executed at hf 0.05,
theta-tube 0.005. PASS 1: 3777 certified (95.2%), min margin
0.0496; the 192 failures are exactly the four b3-wall columns
(k pi/3, k = 1, 2, 4, 5) plus shoulder tiles at the b2 = pi/3,
5pi/3 crossings (the corner rows) — the histograms are clean
delta functions on the known special set. Refine pass 1 (proxy
3e-3 alone) fixed nothing — diagnosis: find_bases' 1e-5 tolerance
can't see near-bases at a 3e-3 proxy where wall orthogonality
defects are ~0.4-0.9 * proxy. face_tile now scales the tolerance
with the proxy (30 * TH_PROXY), matching the measured wall/corner
structure (8 near-bases, margins 0.10-0.17). Refine v2 in flight.

## 4.60 FACE COMPLETE: 3969/3969 tiles certified (demo grade)

Final accounting of the full face execution (hf 0.05, theta-tube
0.005, JSONL ledgers, restart-safe):
  pass 1 (proxy 1e-6):        3777/3969, min margin 0.0496
  refine (proxy 3e-3 + proxy-scaled near-basis tolerance
          + stratum-count cap 100): 112/112 of the remaining
          tiles, min margin ~0.05 — including the four b3-walls
          (odd-k: 8 near-bases; even-k: 8-16, count fragile and
          not load-bearing), the corner rows, and 12 tiles on the
          pi/3 wall carrying the 70-near-basis F(1/6,0)-class
          Fourier stratum (richest stratum seen; pair-level
          obstruction holds there too, as JMMS requires).
TOTAL: 3969/3969. The face piece of the theorem is EXECUTED at
demo grade — the first fully executed piece of the domain
decomposition. Wall-clock ~75 min on 7 cores across passes.
Certified-pass upgrades needed per tile: dual-AD in-face rates
(replace RATE_FACE 0.5 sample-pad), interval pair defects,
enumeration coverage (fat-sweep hulls), and the stable kernel
(4.57) at the corner rows.

## 4.61 Full-circle collar chain: 628/628, 1.79 s/tile, zero fallbacks

collar_circle.py: one complete beta2 circumnavigation of the torus
at b3 = 2.0 (generic), slab [0.005, 0.01], hf 5e-3, signed tiles,
warm pools: 628/628 CERTIFIED, 0 pool fallbacks, 1.79 s/tile
amortized (1123 s wall on one core, running alongside a 7-worker
face refine). Confirms at scale: no hidden special loci in the b2
direction — the exceptional set is exactly the b3-walls and the
corners, as mapped. Collar campaign price at this rate: a full
slab (628 x 63 lines) ~ 20 CPU-hours; the 3-4 slab collar ~ 3-4
CPU-days single-core, embarrassingly parallel — consistent with
the 4.37 estimate.

## 4.62 THE BULK COLLAPSES: signed tiles certify at every depth, 1.39 s/tile chained

The v3 signed tile — invented for the collar — was never tried at
bulk theta. Tried (2026-08-01 late morning):
  ladder at (1.5349, 0.5516), slabs [th, th+0.01], hf 5e-3:
    th = 0.05, 0.2, 0.5, 1.5: chi 2, 24-50 edges     CERTIFIED
    th = 1.0 (56 roots):      chi 2, 74 edges        CERTIFIED
  chain at th-slab [0.5, 0.6] (10x-THICK slab, adaptive), b3=2.0,
    20 steps: 20/20, 0 fallbacks, 1.39 s/tile.
No valley windows, no fold machinery, no Q-tubes — margins grow
as 0.42 theta, so slab thickness and in-face widths scale with
theta (span^2 and hf^2 curvature charges bound them at
span ~ 0.1, hf ~ 0.14 theta). Estimated certificate-leg tile
count over the fundamental domain (theta in [0.02, pi/2],
adaptive widths): ~2-4e5 tiles ~ 4-8 CPU-days ~ $30-100 parallel
— THREE ORDERS below the coarse-tile campaign estimate.
CAVEAT (the honest leg): this certificate trusts enumeration at
demo grade. The certified pass needs the coverage sweep
(fat-sweep hulls -> blob-vs-root accounting) per tile; its
per-tile cost at bulk-hf boxes is being measured — coverage is
now THE cost driver of the whole theorem.

Addendum 4.62: robustness + fat tiles. The signed bulk tile
certifies at the TRUE census points (theta = 1.535 and 1.105 —
the valley-heavy anchors of the old machinery; 48 and 56 roots,
chi 2 both). Fat tiles: hf 0.02 certifies chi 2 at theta = 0.5
and 1.0; hf 0.05 certifies chi 4 at theta = 1.0 and fails only at
the theta = 1.5 equator edge (use 0.02 there). Certificate-leg
count refines to ~1e5 tiles ~ 2 CPU-days ~ $20-50 parallel.

## 4.63 Coverage-leg datum + composite bulk price: ~$100-500 all-in

The naive measurement: fat_sweep_hulls over ONE bulk signed-tile
box (hf 5e-3 at (0.5, 1.0, 2.0)) exceeded 25 CPU-minutes without
terminating (killed) — fat beta-taxes leave too many survivors.
So per-tile naive coverage is out; the viable coverage routes are
the ones already measured elsewhere:
  (a) chain-anchored coverage (coarse_chain's design: anchors at
      wmin 0.025 + D1 exclusion cache, measured ~35-40 s
      amortized per 3e-3 step -> ~60 s per 5e-3 signed tile);
  (b) the 41x GPU sweep (gputile, box-exact validated) applied to
      (a)'s anchor sweeps;
  (c) slab-level starvation sweeps (4.5's fat-tile idea) —
      unmeasured at bulk, possible further multiplier.
COMPOSITE BULK PRICE with the signed-tile redesign:
  certificate leg:  ~1e5 adaptive tiles x 1.4 s  ~ $20-50
  coverage leg:     ~1e5 x 60 s CPU / 41x GPU    ~ $50-150
  TOTAL:            ~$100-500 all-in
— a 50-100x reduction from the $25-50k coarse-tile design,
obtained by pointing the collar's signed tile at the bulk. The
whole-theorem compute (bulk + collar + face + corner patches) now
prices at ~$500-1500, i.e. the "$1M problem" of the survey has
been argued down ~three orders of magnitude by measurements and
design, not hardware.

## 4.64 External review response (5 findings, all actioned)

A code review (2026-08-01) found five issues; disposition:
  1. Coverage sampled blobs, not boxes (HIGH): FIXED at demo
     grade — _coverage_blobs and coarse_tile stage 4 now check
     EVERY member box/hull cell geometrically (cell + r inside
     ball(root, R_LOC), finer 0.025 hull cells), with the
     whisker-tail refinement sweep as the only fallback and loud
     failure otherwise. First box-wise run flushed out that the
     0.1-grid hulls were too fat for any inclusion test — the
     old cluster-sampling pass had been vacuously generous.
     The certified-grade coverage verifier is owned by the
     rigor/coverage-verifier branch (colleague).
  2. Failed anchor did not halt the chain (HIGH): FIXED —
     anchors return ok, chain_step raises on a failed state.
  3. CURV = 5.0 unproved (HIGH): marked EXPERIMENTAL in code;
     the cache shortcut is an amortization experiment until the
     interval-Hessian re-verify exists.
  4. Theta-tax dropped on a point derivative (HIGH): FIXED —
     the drop now requires derivative PERSISTENCE across the
     slab (dO_theta > curv0 * span), using the same sampled
     curvature sups charged elsewhere. Re-validated: generic,
     wall, exact corner, bulk thick-slab, census point all still
     certify (chi 2-4).
  5. Empty enumeration not fail-closed (MEDIUM): FIXED — empty
     pools and empty root sets now fail loudly in collar_tile
     and _pair_colors.
  General: [SAMPLED] grade tags added to tile verdicts and the
  README (a SAMPLED "CERTIFIED" is a validated experiment, not a
  theorem); test_smoke.py added (kernel bounds, interval
  containment, map-port agreement, fail-closed paths, one tile
  of each kind, certified-S vs FD enclosure).

Addendum 4.64 (a bug the new tests caught, in the OTHER direction):
the smoke test comparing certified_S against root_data2's FD S
failed — and arbitration by independent central differences of
the polished root showed the CERTIFIED path is right (agrees to
~1e-9, 400x inside its 1.6e-5 err bar) while root_data2's
continuation-based S carries an ~9e-5 delta-INDEPENDENT
systematic (theta-column dominated) at the tested bulk point.
Consequences: (a) certpair/4.40 enclosures stand, better than
claimed; (b) the demo tiles' sampled signed rates inherit an
~1e-4 absolute rate error — absorbed by PAD_CORR = 3x for the
rate scales in play, but the certified pass should prefer
certified_S (or central differences) over root_data2's S
everywhere; (c) the test now uses central differences as the
reference. Tracking down root_data2's systematic (its internal
continuation step?) is a small open item.

Addendum 4.64 (persistence gate follow-through): the honest
theta-tax (drop only when dO_theta > curv0 * span) fails 3 of the
48 census cells at the base rung — their earlier certificates were
leaning on the unsound point-derivative drop (the reviewer's
exact concern). All three CERTIFY on 5x-thinner slabs
([0.005, 0.006]: chi 2, 2, 4) — the persistence threshold scales
with span, so the dyadic ladder absorbs the honesty cost locally.
Collar remains 48/48 under review-hardened semantics with
adaptive slabs; census tooling should try thinner slabs as a
rung alongside thinner hf3.

## 4.65 z2-delta second layer decoded: G- = sin(theta) * g1(phi, lam)

The z2 stable form needs delta = zeta - z3sq via the numerator
N' = (t1 + t2) * G- - i * t3 * G+, where G-+ = B12^2 * bigden -+
cB11^2 * bignum. Probed at 40 digits: G- vanishes IDENTICALLY in
(phi, lam) at theta = 0 (exactly 0.0 for all offsets) and is
linear in theta with slope ~6 at the corner; G+ ~ 4.0 = O(1).
Hence G- = sin(theta) * g1(phi, lam) + O(sin^2), with g1 a closed
form obtained by differentiating at s = 0 (A11 = alpha + s mu,
A12 = beta + s nu with alpha/beta theta-only and mu/nu phi-only —
at s = 0 the phi-dependence drops, which is WHY the identity
holds). With t1 ~ (lam - pi/3), t2 ~ theta^2, t3 ~ theta *
cos(lam/2 + phi), and G- ~ theta, every product in N' carries an
explicit small factor — the z2 quotient becomes cancellation-free
one level down, same pattern as P/Q. Remaining: assemble g1
symbolically and wire delta into karlsson_stable + ivstable
(mechanical; removes the mp-hybrid and the interval z2 blowup).

Addendum 4.64 (coverage-cost verdict): the box-wise coarse-tile
validation ran 134 CPU-minutes without completing and was stopped.
Anatomy: the 0.025-grid hull-cell census produces O(1e4-1e5)
cells needing refinement at this census tile, and per-chunk
mini-sweeps (50 cells x ~45 s) put per-tile coverage at HOURS —
correct and fail-closed, but impractical as a per-tile leg.
Verdict recorded honestly: the old cluster-sampled 6/6 census
coverage claims are DOWNGRADED to sampled-grade pending sound
coverage; the box-wise semantics stand as the specification; the
practical implementation (batched streaming refinement, segment
amortization, GPU sweeps) is exactly the coverage-verifier
work (colleague's branch). The demo-tier certificate/margin
results are unaffected — coverage was always their listed gap.

## 4.66 FIRST THEOREM-GRADE SIGNED TILE (task 2 complete)

rigor_tile.py: the signed tile with every pair bound RIGOROUS —
Q-curves + TM-certified curve residuals (tmres, no sampling/PAD)
+ certified tile-local RJ_extra (rates.py) + slanted-tube
containment (tube_krawczyk with quadratic offset) + TM
inner-product overlap lower bounds (certified_overlap_lo) +
exhibited proper coloring. Result objects are
certificate_result.CertificateResult with per-dependency
Evidence; the weakest-dependency rule correctly grades the tile
SAMPLED_BOUND while the pair layer itself is RIGOROUS —
enumeration coverage (rigor/coverage-verifier) is the one
non-rigorous dependency, by design.

Measured at (0.5, 1.0, 2.0), the h-ladder (the tube-contraction
wall is |Y|-limited: Krawczyk fixpoint needs |Y| * rad_g <~
3.5e-3 with rad_g ~ h^2-scaled):
  h = 2.5e-4: tubes  0/52 (global-Lipschitz RJ was 20x fat; fixed
              to certified RJ_extra; still |Y| ~ 23 blocks)
  h = 1e-4:   tubes 41/52, 787 rigorous deletions, chi 13
  h = 5e-5:   tubes 45/52, 961, chi 9
  h = 2.5e-5: tubes 49/52, 1151, chi 5  -> CERTIFIED, ~6 s/tile
The 3 remaining tube-less roots are valley-class (large |Y|);
fail-closed as full-edge vertices they still fit chi <= 5. The
composition with valley windows (certify_tile's L4) is the
completion for larger h; the working law: rigorous-pair h ~
2.5e-5-1e-4 at bulk conditioning vs 5e-3 sampled — the rigorous
campaign economics run through chains/adaptivity exactly as the
fine-tile design always said. En route fixed: the global L_H_J
RJ bound was 20x fatter than rates.py's certified RJ_extra — the
debug trail that found it is the night's second dependency-blowup
lesson (after 4.57's).

Addendum 4.66 (validation ladder): the theorem-grade tile
certifies at 5/5 test points with adaptive h:
  theta=0.05  h=1e-5:   48/48 tubes, chi 2 (h=2.5e-5 gave chi 8)
  theta=0.2   h=2.5e-5: 48/48 tubes, chi 2
  theta=0.5   h=2.5e-5: 49/52 tubes, chi 5
  theta=1.0   h=2.5e-5: 53/56 tubes, chi 5
  census pt   h=2.5e-5: 48/48 tubes, chi 2
~6 s/tile throughout (pool-dominated). Required h tightens toward
the collar (conditioning grows as theta drops) — the adaptive-h
rung ladder is the campaign shape, as everywhere else.

## 4.67 g1 closed form VERIFIED (z2-delta layer-2 identity)

Hand-derived and verified to 1e-24 at 40 digits at arbitrary
(phi, lam):
  g1(phi, lam) = dG-/ds at s=0
    = i (sqrt3/2) [ -4 w e^{i phi} (w - z1 conj(w))
                    - 2 w^2 (e^{-i phi} - z1 e^{i phi}) ],
  w = e^{2 pi i/3}, z1 = e^{i lam}.
(Derivation collapses via conj(mu) = -nu at s = 0.) With 4.65's
structure this fully specifies the stable z2 numerator through
first order; the remaining mechanical piece is the s^2 remainder
g2 (same differentiation, c = cos(theta) now contributing) and
the symmetric treatment of num2 = cB12^2 (zeta' - z3sq) — both
routine. Until wired, the mp-hybrid (0.25 ms) covers float grade
and elevated-dps intervals cover certified corner tiles.
