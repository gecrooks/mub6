# Roadmap: a rigorous exclusion of MU quadruples through the Karlsson family

**Target theorem.** No set of four mutually unbiased bases in C^6 contains a
Hadamard matrix from Karlsson's three-parameter family K6^(3) (equivalently,
any H2-reducible Hadamard). Combined with the existing rigorous results
(Fourier family excluded — Jaming–Matolcsi–Móra–Szöllősi–Weiner 2009 and
Matolcsi–Weiner's LP proof; S6 has no triple) and the conjectured
classification of order-6 Hadamards (S6 ∪ K6^(3) ∪ G6^(4)), this would
reduce Zauner's conjecture to the four-parameter family G6^(4).

Why this family: the McNulty–Weigert review identifies K6^(3) as precisely
where the Fourier-analytic linear program *fails* to exclude complete sets,
Gröbner computations exhaust memory, and a claimed exclusion was retracted
(erroneous lemma). It is the acknowledged soft spot, and no rigorous
quadruple statement exists for any non-affine Hadamard.

The strategy is a certified-numerics pyramid, prototyped in `certify.py`.

## Layer 0 — pointwise certificates (prototyped today)

For a fixed Hadamard H (stated for every matrix in an explicit ball around
a float H0, so exactness of H0 is not required):

1. **Certified root enumeration.** MU vectors of {I,H} are zeros of
   g_k(θ) = |⟨h_k|u(θ)⟩|² − 1/6 on the phase 5-torus (θ_0 = 0 by gauge,
   Σ_k g_k ≡ 0 for unitary H, so g_1..g_5 is a square system). Adaptive
   branch-and-bound with the exact Lipschitz bound |∂g_k/∂θ_j| ≤ 2|s_k|/6
   excludes boxes; surviving boxes cluster at roots; a Krawczyk interval-
   Newton test proves existence and uniqueness per cluster with the exact
   Hessian bounds (7/18 diagonal, 1/18 off-diagonal). Coverage is checked:
   every non-excluded box lies in a verified uniqueness box. All
   inequalities carry explicit floating-point slop (1e-11, a ~300× margin
   over worst-case IEEE error for these expression sizes); the structure is
   ball-arithmetic-ready (Arb) for publication grade.
2. **Certified graph.** Pairwise |⟨u,v⟩| intervals from the root enclosures
   (Lipschitz constant 1/6 per phase). Pairs whose interval reaches 0 are
   *possible* edges; all others are certified non-edges.
3. **Conclusion.** If the possible-edge graph has clique number < 6, no six
   MU vectors are mutually orthogonal ⟹ no orthonormal basis among MU
   vectors ⟹ **no MU triple {I,H,K}, a fortiori no quadruple containing
   {I,H}** — for every H in the stated ball.

Measured margins that make this work: root separations ~1.2 rad; smallest
nonzero overlaps 0.004–0.15; σ_min(Jacobian) 0.015–0.045; graph at generic
K3 points is a near-perfect matching (clique number 2).

Cost today: ~5M boxes, ~10 min per point in NumPy. A C/Arb implementation
with standard contractors (HC4, interval Newton on boxes rather than pure
bisection) should reach seconds per point.

## Layer 1 — parametric tiles

Replace the point H0 by a box B in the 3-dimensional K3 parameter space
(θ, φ, λ). The certifier already threads an H-ball radius (`hslop`) through
every bound; a tile certificate sets hslop = L_param · diam(B), where
L_param is an explicit Lipschitz constant of the parametrization. The
Krawczyk step then proves: for *every* parameter in B the root set consists
of N continuously varying verified roots with edges inside the fixed
possible-edge set; clique < 6 kills triples on the whole tile.

**Do NOT thread the parameter through the H-ball slop.** The naive route —
reusing the pointwise certificate with hslop = L_param·diam(B) — compounds
worst cases (|J⁻¹|·|∂g/∂H| ≈ 500 rad/rad formal root sensitivity) and
caps tiles at ~3e-5, i.e. ~10^15 tiles: dead. Measured true sensitivities
at K6(0.9, 2.1, 0.7) are two orders tamer: max |dθ*/dβ| = 7.7 rad/rad
(per-parameter: 7.7, 2.5, 2.1), pair-overlap drift ≤ 0.59/rad (median
0.08, binding pairs 0.15). A *parametric* Krawczyk — verify the root
curves θ*(β) over a parameter box using first-order sensitivities with a
certified second-order remainder — recovers the true rates, supporting
tile diameters ~1e-2 at the binding margins (0.004–0.011 / 0.6) and up to
~0.1 where a clique-robust partition certificate is used (prove no
6-clique via a ≤5-class partition whose within-class non-edges are the
fat-margin pairs — tolerates individual thin margins flipping). That
restores ~10^6–10^7 adaptive tiles for the 3-parameter family. Additional
multipliers: (a) adaptive tiling, refinement driven by certificate
failure; (b) quotienting by the family's finite equivalence group (the
analogue of the 144-triangle structure of the Fourier square); (c)
locally-constant certificates — N and the graph change only at bifurcation
strata, so neighboring tiles share almost all work.

**Measured Layer-1 law (prototype, 2026-07-25).** The working prototype
(`parametric.py`: zoned per-box slanted taxes, tube Krawczyk, partition
certificate) certifies h = 3e-4 tiles in 26 s of NumPy at the reference
point — 100x the threading radius, ~10^6x fewer tiles. Each root imposes
a ceiling h_i ~ sigma_min,i / (2 sqrt3 PAD (Hess |S_i| + g_beta)) — at the
reference tile the worst roots give 1.4e-3, the median 5e-3..2e-2 — and
just under a ceiling, axis-aligned refinement pays (L/gain)^5 per octave.
A second, tighter ceiling binds first-order tubes in practice: the
tube-boundary margin — the Krawczyk radius is pushed up by the O(h^2)
curve residual while the slant tax squeezes exclusion just outside it
(measured: h = 6e-4 fails with 42 stuck boxes at root 14's tube
boundary). SECOND-ORDER (Q-) TUBES, built 2026-07-26, eliminated that
ceiling exactly as designed — ladder: 6e-4 PASS, 1e-3 PASS, 1.4e-3 FAIL
at the predicted slope wall (1.39e-3). Measured Q-gain: ~2.4x in h, ~13x
in tiles. The slope wall itself is untouchable by curve modeling (its
constant is the e-drift of the integrand); the big levers are
fold-delegation (65,000x, above) and margin-cached continuation (~10^3x
per-tile cost). Smooth-part headroom left: direction-resolved local
taxes (~2x) and anisotropic tiles (~2x volume); rotated singular-frame
boxes near poorly conditioned roots to erase the (L/gain)^5 alignment
penalty.

**Clique-robustness is load-bearing, not optional.** A 6-point parameter
scan found binding overlap margins fluctuating from 6.3e-3 down to 1.3e-4
(K6(0.9, 2.1, 0.7) sits near an edge-creation bifurcation: one MU-vector
pair is 1.3e-4 from orthogonal). Freezing the whole graph across a tile
dies at such points; but the *conclusion* is robust — an extra edge
appearing on a near-perfect matching raises the clique number from 2 to at
most 3, nowhere near 6. Tile certificates should therefore prove "no
6-clique" via a ≤5-class partition witnessed by fat-margin non-edges only,
letting thin margins flip freely. Root conditioning also varies (σ_min
from 0.045 down to 0.005 across sampled points; the certifier's local
refinement fallback handles the resulting suspect-spread blowups).

**Family-wide ceiling statistics (2026-07-26) forced two design upgrades.**
A 24-point scan (`scan_ceilings.py`) showed the tax-slope ceiling has
quartiles 1.8e-4/5.4e-4/1.2e-3 with a tail to 1.7e-5: half the family
sits near fold strata where sigma_min -> 0, |S| ~ 1/sigma_min, and the
ceiling ~ sigma_min^2 — the tile integral across a stratum transversal
diverges. Uniform Layer-1 tiling is therefore IMPOSSIBLE in principle,
not merely slow. Upgrades:

1. **Valley windows are co-equal machinery** (built 2026-07-26,
   `fold.py`; supersedes the fold-pair hypothesis). Measurement showed
   the ill-conditioned roots are not fold pairs but racing roots in long
   shallow 1-dim trenches (phi monotone, slope ~ sigma_min, no partner).
   Certificate: singular-frame reduction (4x4 implicit-function block,
   sigma_4 healthy) + 1-dim monotone-crossing window sized so the trench
   floor clears the |s|-local beta tax; sweep collects the curved tube
   via per-root oracles; partition rows sampled along the racing dip.
   DEMONSTRATED: scan pt 8 (wall h = 1.7e-5) certified at h = 3e-4 in
   36 s — 18x beyond its Layer-1 wall. Deepest valleys at larger h need
   per-root beta-subdivision of the (cheap, 1-dim) valley certificate —
   the divergence is structural no more. The fold-PAIR normal form
   (Miranda, 0/1/2 roots) remains the design for genuine fold crossings
   when a campaign tile actually straddles one (dip enclosures with <= 2
   runs are already accepted by the code; pair bookkeeping in the
   partition is a TODO). Delegation at sigma_min cut 0.03 (12% of roots)
   lifts the tile-count integral from 4e14 to 6e9 (Result 8b).

2. **Margin-cached continuation — BUILT AND MEASURED 2026-07-26**
   (`cache.py`, Result 11): anchor stores 2.07M excluded boxes with
   excess + drift rate (95 MB float32; 260k slant-excluded boxes stored
   verbatim); chain step = 4 ms vectorized far-field re-check + patch
   mini-sweep of the failing minority + per-root stage with anchor-fixed
   guards. Measured: 12/12 chained tiles certified from one anchor
   (20.2x sweep-volume amortization), steps 3.5-10 s vs 26 s standalone;
   patch fraction grows linearly with distance (18% -> 88% over 12
   steps), setting the re-anchor cadence. In the C substrate the check
   stays O(ms) and FD sampling becomes interval evaluation — the ~10^3x
   per-tile collapse the campaign estimate assumed is real.
   ADDENDUM (deep valleys): beta-subdivision cannot fix their
   edge-clearance failure (window-edge requirement is set by the
   full-tile tax); the fix is ANISOTROPIC tiles — thin along the racing
   direction phi_beta (the stratum transversal), fat transverse (the
   trench does not move to first order) — a per-tile h-vector in the
   adaptive driver.

## Layer 2 — the failing strata are the interesting ones

Tiles where the certificate fails cluster on measure-zero sets:

- **Root-merging loci** (double roots; Krawczyk cannot separate). The
  no-6-clique conclusion does not require root *uniqueness* — switch the
  per-cluster test from Krawczyk to a Miranda/Brouwer covering certificate
  ("all roots, however many, lie in these enclosures"), which tolerates
  multiplicity. Today's F(1/6,0)-chart analysis shows exactly what these
  look like (universal degenerate spectrum 1/3, 1/3, 1/(3√2), 0, 0).
- **The Fourier subfamily F ⊂ K3** — bases exist (8–16), so no-triple is
  false there; quadruples are already excluded in the literature. Needed:
  a *stability shell* — on tiles abutting F, certify instead the
  quantitative Grassl statement: every base-pair among MU vectors has
  unbiasedness defect ≥ 1/6 − C·diam (our F6 control certifies the 1/6 wall
  at machine precision), whose Lipschitz continuation across the shell
  closes the gap to the literature result on F itself.
- **The Szöllősi subfamily X6^(2) ⊂ K3** — triples exist; this is the one
  place genuinely new mathematics is required (Layer 3).
- **Unknown strata** — any tile that fails for none of the above reasons
  would be a discovery about the family, not a failure of the program.

## Layer 3 — no fourth basis through Szöllősi triples

Each X6^(2) point carries an explicit MU triple {I, H(α), K(α)} (Zauner's
circulant-block construction, Thm 7.4–7.5 of the review). A quadruple
through X needs a vector MU to all three bases: 12 constraints on 5 phases,
massively overdetermined. Certify over the 2-parameter X-family plus 5
phases (7-dimensional parametric sweep, same machinery, faster exclusion
because 12 constraints prune boxes brutally): **no vector is MU to any
Szöllősi triple.** Our unextendibility spot-checks (0 solutions at sampled
triples, 8000 starts) and the literature's constellation searches ({5³,1}
never found) both say the margins here are healthy.

## Layer 4 — assembly and formalization

Compactness assembles: certified tiles ∪ stability shells ∪ Layer-3
certificate ∪ literature results cover all of K6^(3). Publish certificate
data (box trees, root enclosures, pair bounds) plus an independent checker.

**Rigor substrate (decided 2026-07-25): NOT a wholesale Arb port.** Scalar
arb/acb (python-flint) would cost ~300x versus the vectorized NumPy sweep
and put the campaign back at CPU-decades; and the margins (1e-4..1e-1
against double rounding 1e-16) present a trust problem, not a precision
problem. Plan instead — generator/checker separation:
1. NumPy stays the fast untrusted *searcher*; it emits a certificate file
   per tile. Algorithm design (parametric Krawczyk, partition
   certificates, Layer 3) is iterated here first; freeze the format.
2. A static rounding-error lemma (Rump-style, over the fixed expression
   DAGs of g, J, and pair overlaps) turns the 1e-11 slop constants into a
   theorem at zero runtime cost — immediate credibility for existing
   certificates. Add a scalar-Arb spot-audit of Krawczyk steps (~50/tile).
3. The trusted *checker* is a small C kernel with directed-rounding
   double intervals (~2.5 s/tile => ~a month of desktop for 10^6 tiles;
   C-level Arb would be ~200 s/tile => CPU-years). Arb (C library) is
   reserved for tile-corner enclosures of H(θ,φ,λ) entries and for
   auditing pathological clusters. Tiles landing within ~1e-12 of a
   bifurcation are split, not run at higher precision.

**Cost estimate**, revised after the working tile: the LIFO-chunked sweep
runs a full pointwise (or tile) certificate in **~2 s of NumPy** (6.9M
boxes, 447 MB peak RSS; the prototype's earlier minutes-to-hours timings
were swap-thrash of the breadth-first frontier, since eliminated). The
demonstrated tile radius is 2.9e-6 with only 7% of the margin budget
consumed — ~3e-5 tiles fit the present margins, and clique-robust
partition certificates plus parametric Krawczyk should reach the
1e-3..1e-2 sizing of the sensitivity analysis. Even at a pessimistic 1e-3
tile size, ~10^9 tiles × 2 s ≈ 70 CPU-years in Python shrinks to
CPU-weeks in C/Arb with 10-100x contractor gains — and at the
sensitivity-supported 1e-2 sizing it is ~10^6 tiles ≈ days on a desktop.
The workload is embarrassingly parallel and margin-driven. The single
biggest open design risk is the parametric Krawczyk near root-merging
strata, where sensitivities diverge; Miranda-style multiplicity-tolerant
certificates bound the refinement there.

## Beyond K3: the generic family G6^(4)

No closed-form parametrization exists, so Layer 1 must tile the Hadamard
*variety* directly: pin 4 dephased phases to get verified local charts via
interval Newton on the unitarity system (the same Krawczyk machinery, one
level up), then run the pointwise certificate per chart tile. Our numerics
make this look promising: generic G6^(4) points have *edgeless* MU-vector
graphs — the cheapest possible certificate — with zero degenerate roots.
The atlas bookkeeping is the hard part. Excluding S6 (done, both by
Lasserre in the literature and by today's independent certificate), F/F^T
(literature), K6^(3) (this program), and G6^(4) (variety atlas) would close
Zauner's conjecture modulo the classification Conjecture 7.1 — which is
itself a candidate for the same interval technology.

## Gap list to the theorem (status 2026-07-26)

Target statement: no set of four MU bases in C^6 contains an H2-reducible
Hadamard (= any member of K6^(3); by Karlsson's classification, any
Hadamard whose dephased form contains a -1). Gaps, ordered by risk:

R1. F-locus shell — DE-RISKED 2026-07-26 (`fourier_shell.py`, Result 10):
    quadruple-exclusion certified at F6 (all 16 cliques enumerated, no
    spurious ones; every clique pair certified >= 1/6 wall); wall drift
    along F(a,b) measured at ~0.05/unit-parameter => O(0.1)-thick shell
    tiles. Remaining: parametric version (same tube machinery).
R2. X-locus / Layer 3 — B-ARC COMPLETE 2026-07-27 (Result 19): the
    entire Hermitian family arc theta in [1.25, 5.03] certified
    strongly-unextendible via ~540 adaptive anchors (zero gaps; two
    Dita special points crossed — rebuild-mode anchoring required at
    special loci, where tracked partners lose their branch). Remaining:
    the same walk over Szollosi's 2-parameter X6^(2) chart, and the
    endpoint segments joining Bjorck's C.
R3. Fold-pair bookkeeping — CLOSED 2026-07-26 (Result 13): a real fold
    located by count-bisection (48 -> 52, two pairs born, orbit of 2);
    certified intra-dip self-overlap makes the single-vertex collapse
    sound for pairs; doorstep tiles certified at ordinary cost. Birth-
    STRADDLING tiles certified via phantom anchors (two-stage sweep:
    coarse divert -> far-cluster |g|-minimizer polish -> depth-thresholded
    phantom valley windows -> fine resume), demonstrated at 33 s with
    both birth sites auto-found. Root count changing inside a tile is
    routine; no open machinery questions remain in R3.
R4. Rigor substrate — CORNERSTONE BUILT 2026-07-26 (Result 14):
    directed-rounding intervals + certified Karlsson map (200/200
    containment; branch margins fat at all scan points) + dual-AD
    mean-value tile bounds that BEAT the empirical PAD*FD constants at
    campaign h (ratio 0.68 at h = 3e-4; sub-box derivative evaluation
    controls the artifact at larger h); certified g-widths 1.6e-12
    validate and replace the SLOP model. Remaining (mechanical): wire
    interval g/J into the certificate paths (curve residuals at corner
    points, tube Krawczyk data, valley phi-envelopes over t-cells);
    swap libm-faithfulness for Arb; certified 3rd-order Q-remainders.
    The flagship theorems are already upgraded on this substrate
    (Result 15): exact-point membership certified at 5.9e-14 within the
    1e-11 ball; tile radius certified at 2.0e-6; SLOP replaced by a
    derived 1e-13 static bound (100x headroom, Krawczyk needs no
    accuracy from LAPACK).
    SHARPENING (derived, unimplemented): the certified gb bound's |s|<=1
    factor can use the zone-local |s_k| <= (1/sqrt6)(1 + sqrt5 D) valid
    within l-inf distance D of a root (|grad_theta s_k|_2 <= sqrt6/6
    exactly), if the per-root taxes' validity zone is shrunk to
    D = 0.35 (beyond it the certified beta-tax excludes anyway): gb
    7.8 -> ~6.0. Same |s|-localization sharpens the Hessian row bound
    11/18 -> 0.577 (off-diagonals are exactly 1/18 — the mixed
    theta-second-derivatives of s vanish since s is a sum of
    single-phase terms). Net: coef1 down ~20%, slope walls up ~1.25x.
R5. Drivers — FIRST END-TO-END RUN 2026-07-26 (Result 16, `driver.py`):
    a 4x4 grid of margin-cached chain-lines certified a full-dimension
    block (112 tiles, 16 anchors, 17.7 min, zero failures/re-anchors) —
    theorem: no MU quadruple contains any Karlsson Hadamard from the
    certified 4.6e-3 x 3.0e-3 x 3.0e-3 box. ANISOTROPY DONE 2026-07-27
    (Result 27): h-vector plumbed through tm/rates/certify_tile, 3.0x
    volume certified at the deep-valley wall point; driver picks hv
    from center-point |u.Gb| diagnostics. Remaining: the adaptive
    campaign driver itself (per-tile hv chooser, fundamental-domain
    cover + boundary bookkeeping under the order-32 quotient,
    checkpoint/resume); the compactness root-count argument written
    out.
R6. C checker + certificate format + parallel campaign (~6e9 tiles);
    re-verification of the four pointwise theorems on the same substrate.

Out of scope for this theorem (Zauner needs more): G6^(4) via a variety
atlas, and the completeness Conjecture 7.1 (S6 u K3 u G4 = everything).

Complementary to the three published routes, and aimed exactly at their
blind spot: the Lasserre/SDP hierarchies certified S6-no-triple but scale
badly ({5³,1} out of reach); the QRAC SDP hierarchy hits a 2^(32k)-bit
memory wall; the Fourier-analytic LP provably cannot exclude K6^(3). The
certified-sweep pyramid is embarrassingly parallel, margin-driven, and its
hardest instance so far (a generic K3 point) runs on a laptop.
