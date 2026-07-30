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
