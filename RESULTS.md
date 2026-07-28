# MUB-6: hunting a fourth mutually unbiased basis in dimension six

**Date:** 2026-07-24. **Verdict up front:** no counterexample (none was
realistically expected — every computation lands exactly on the known
obstructions). What this session produced: an independent working lab that
reproduces all the key published numbers, plus fresh numerical probes of the
*unexplored* corner (generic / non-affine Hadamards), which came back strongly
supporting Zauner's conjecture that only 3 MUBs exist in d = 6.

## The problem

Zauner's conjecture (1999): the maximum number of pairwise mutually unbiased
bases in C^6 is 3 (vs the d+1 = 7 upper bound achieved in prime-power
dimensions). A counterexample = 4 MUBs {I, H1, H2, H3} where each Hi is a
complex Hadamard matrix and each Hi†Hj is Hadamard too. Weiner's no-gap
result plus known constructions leave max#MUB(6) ∈ {3, 4, 5, 7}.

State of the art (McNulty–Weigert review, arXiv:2410.23997, accepted in
Quantum 2026): rigorously excluded from any quadruple: the entire Fourier
family F(a,b) and transpose (Jaming et al.; Matolcsi–Weiner LP), Björck's C,
Diţă D(0), Butson BH(6,12), Schmidt-rank ≤ 2, three product columns; S6 has
no triple. NOT rigorously excluded: non-affine families B, M, X (Szöllősi),
K (Karlsson), and the generic four-parameter family G6^(4) (a published
exclusion was retracted — erroneous lemma; Gröbner computations exhaust
memory there). The Fourier-analytic LP also fails to exclude complete sets
built from K6^(3).

## What we ran (code in this directory)

`mub.py` — constructors for all order-6 Hadamard families with explicit forms
(F(x1,x2), F^T, Diţă D(x), Björck C, Tao S6, Beauchamp–Nicoara B(θ)), plus:
MU-vector enumeration via multistart Levenberg–Marquardt on phase variables,
orthonormal-basis (6-clique) search, unbiased-pair testing, Tadej–Życzkowski
defect, and a global 4-MUB least-squares search in the dephased-phase
parametrization (unbiasedness to I exact by construction).

All constructors verified Hadamard to 1e-9; defects match theory (S6: 0,
family members: 4).

## Result 1 — Grassl's obstruction, reproduced in 3.6 s

Vectors MU to both I and F6 (= biunimodular sequences for the 6-point DFT):
**exactly 48** (machine precision 8e-15), forming **16 orthonormal bases**,
with **zero mutually unbiased base-pairs** — so no quadruple contains
{I, F6}. Matches Grassl 2004 (MAGMA computer algebra) in every particular.
The nearest any two of the 16 bases get to unbiasedness is a defect of
exactly 1/6.

## Result 2 — the global floor, independently reproduced to 7 digits

Multistart LM over {I, B1, B2, B3} (85 phase parameters after gauge fixing):

| setting | starts | best cost | D̄₄ of best | exact sets found |
|---|---|---|---|---|
| d = 5 control | 8 | 3.4e-31 | 1.0000000 | yes (half of starts) |
| **d = 6** | 200 | 2.562461e-2 | **0.9982917** | none |
| d = 6, F6 pinned | 100 | 8.946667e-2 | 0.9940356 | none (rigorously impossible) |
| d = 7 control | 400 | 5.77e-2 | 0.9967956 | none in 400 starts |

The d = 6 best value **D̄₄ = 0.9982917 equals the published "four most
distant bases" record** (Butterley–Hall 2007; Raynal–Lü–Englert 2011; also
the see-saw SDP / Monte Carlo / Bell-inequality searches) to all printed
digits, found by 48/200 starts. Identity check: our least-squares cost and
their distance measure obey cost = 15·(1 − D̄₄); 0.0256246/15 = 0.00170831 =
1 − 0.9982917 exactly. Same optimum, two coordinate systems.

The d = 6 landscape is *quantized*: 200 random starts collapse onto a small
discrete menu of local floors (0.02562 ×48, 0.06403 ×92, 0.08078 ×15,
0.12440 ×14, …) — a rigid landscape with no soft directions, very unlike a
near-miss scenario. Pinning F6 pushes the floor up threefold, as the
rigorous exclusion demands.

Structural identification: all three Hadamards of our best constellation
carry the transposed-Fourier fingerprint (three dephased columns containing
−1; 24–30/36 entries near 6th roots), matching Raynal–Lü–Englert's analytic
identification of the "four most distant bases" — we reproduce the value
*and* the structure.

Honesty note, and what actually carries the evidence. In d = 7 four MU
bases *provably exist* (eight do), yet 400 random starts found zero exact
quadruples: the d ≥ 7 landscape is glassy and multistart search is nearly
blind there. So "search failed in d = 6" is, by itself, weak. What is
strong: (i) the d = 6 landscape is *not* glassy — 24% of independent starts
funnel to one identical quantized value D̄₄ = 0.9982917, (ii) that value is
exactly what five independent published methods find (LM, steepest ascent,
see-saw SDP, Monte Carlo, Bell-inequality optimization), and (iii) the
rigorous family exclusions above. A hidden exact quadruple would have to
hide in a landscape that everywhere else advertises its optima loudly.

## Result 3 — the unexplored corner: generic and non-affine Hadamards

Counting vectors MU to {I, H} and their structure across families
(3000-start enumeration per point):

| H | MU vectors | bases (→ triples) | MU base-pairs (→ quadruples) |
|---|---|---|---|
| F(0,0) = F6 | 48 | 16 | 0 |
| F(a,b) generic | 48 | 8 | 0 |
| F(1/6,0) (this chart; root-merging point) | 42 distinct (=48 w/ mult.) | 7 | 0 |
| Diţă D(0) | 120 | 10 | 0 |
| D(0.05) | 72 | 4 | 0 |
| D(1/8) | 48 | 0 | 0 |
| Björck C | 54 distinct (see Result 4) | 1 | 0 |
| Tao S6 | 90 | 0 | 0 |
| B(θ), 4 values | 72, 72, 60, 60 | **1 each** | 0 |
| random Hadamards ×11 (G6^(4) territory) | 48–54 | **0** | 0 |

All literature values reproduced (120/72/48 along Diţă, 56 for C, 90 for S6,
no bases for S6). New territory:

- **Generic random Hadamards** (defect 4, not H2-reducible — i.e. generic
  points where no rigorous result exists): ~48 MU vectors whose
  orthogonality graph is **edgeless** — largest orthogonal set = 1. Not a
  single orthogonal *pair* among the MU vectors, versus 16 full bases for
  F6. Extension structure does not fade gradually away from the affine
  families — it disappears entirely. (One random draw converged onto the
  isolated Tao matrix S6 — defect 0, 90 vectors — a nice self-check.)
- **Hermitian family B(θ)**: exactly one basis among its MU vectors at each
  sampled θ — the triple exists (consistent with B ⊂ X6^(2), Szöllősi) and
  is unique; no second basis, so nothing to pair into a quadruple. The
  MU-vector counts (72/60 with a transition between θ ≈ 1.80 and 2.07)
  appear not to be in the literature.

Everything is consistent with Conjecture 8.1 (McNulty–Weigert): {I,H}
extends to a triple only for H in F, F^T, or X6^(2) — and with Zauner: no
quadruple anywhere, with huge margins in every direction probed.

**Karlsson family K6^(3) — the LP-resistant weak spot, probed directly.**
The review singles out K6^(3) as the family the Fourier-analytic linear
program cannot exclude from complete sets. Sampling it via its defining
H2-reducibility (all nine 2×2 blocks Hadamard — 18 extra residuals in the
random-Hadamard solver; no closed form needed): 8 samples, 48–60 MU vectors
each, **zero MU base-pairs**; the only samples admitting any bases (48/8,
min defect ≈ 0.10 — the Fourier-generic signature) turn out to *be*
Fourier-family points (three dephased rows containing −1; 24/36 entries at
exact 6th roots — F(a,b)'s fingerprint). Genuine generic K3 points: no
bases at all. Conjecture 8.1 survives its first direct numerical probe
inside the LP-resistant family. (Notable solver sociology: least-squares is
attracted to high-symmetry strata — it found Fourier points inside K3 and
Tao's S6 from generic starts.)

## Result 4 — double roots and a numerical trap (methodological)

At special symmetric points the MU-vector polynomial system has *isolated
degenerate roots*: solutions where the 5-parameter Jacobian drops to rank 3,
with a universal singular-value pattern (1/3, 1/3, 1/(3√2), 0, 0). Verified
by predictor–corrector walks that snap back (not continua). Consequences:

- **Björck C**: 57 raw solutions collapse to **54 distinct** (52 simple + 2
  degenerate). The literature's "56" is consistent with a count *with
  multiplicity* (52 + 2×2). After polishing to machine floor, exactly **one
  orthonormal basis** appears among the 54 — the (unique) MU triple through
  {I, C}, as B ⊂ X6^(2) extendability predicts. Sloppy tolerance had hidden
  it: near a degenerate root, residual 1e-10 still means phase error ~1e-5,
  which silently deletes orthogonality-graph edges at test tolerance 1e-8.
- **F(1/6,0) in the Bengtsson chart**: same pathology — 42 distinct
  vectors after 6 pairwise root mergers (= 48 with multiplicity; a
  248,832-seed lattice enumeration confirms 42 is complete), 7 bases. A
  6×4 scan of the fundamental triangle shows every generic point at
  exactly 48/8/0 (matching the literature) with two special points: F(0,0)
  (16 bases, Grassl) and this root-merging point. The literature's
  "70 bases at (1/6,0)" refers to the McNulty–Weigert chart of the family
  (their Eq. 7.15), a different parametrization from the Bengtsson
  F(x1,x2) implemented here — their special point sits at an exceptional
  parameter (codimension ≥ 1) our coarse grid does not hit. At that point
  too the literature reports zero quadruples.
- **Generic random Hadamards have no degenerate roots at all** (all 48
  simple, well-conditioned) — the generic-corner conclusions above are
  unaffected. The recurring count of 48 at every generic point (Fourier or
  random) looks like the universal generic root count of the system.

Moral for anyone hunting numerically: the interesting special points are
exactly where double-precision multistart quietly undercounts, and the
failure mode (missing vectors AND missing graph edges) biases *against*
finding extensions — a counterexample hunter must polish to machine floor
and match tolerances to conditioning, or they will discard their needle.

## Result 6 — certified (rigorous) exclusions: from numerics to theorems

`certify.py` implements a computer-assisted proof pipeline (see
`PROOF_ROADMAP.md` for the full program): adaptive branch-and-bound over
the phase 5-torus with exact Lipschitz bounds excludes regions containing
no MU vectors; Krawczyk interval-Newton verifies each surviving cluster
holds exactly one root; coverage is checked; certified pairwise-overlap
intervals bound the orthogonality graph; clique number < 6 then *proves*
no MU triple {I,H,K} exists — for every matrix H in an explicit ball
around the target (so float representation is immaterial). All
inequalities carry explicit slop (1e-11, ~300x worst-case IEEE error);
structure is Arb-ready.

Certified results (each ~5M boxes, ~10-15 min in NumPy):

- **F6 control**: exactly 48 MU vectors verified with full coverage;
  300-edge graph, min non-edge margin 0.1547; 6-cliques correctly found
  (bases exist — control passes in both directions). A certified
  re-derivation of Grassl's enumeration.
- **S6 (Tao)**: certified — exactly 90 MU vectors, **edgeless** graph
  (min non-edge margin 0.1545 ≈ sin(π/10)/2, a golden-ratio value), clique
  number 1. **Theorem: no MU triple contains {I, S6}** — caveat-free,
  since S6's exact entries (third roots of unity) lie within trivial float
  slop of the target, well inside the certified ball. Independently
  reproduces the known Lasserre-SDP certification by an unrelated method.
- **K3 seed 3003 (LM-sampled generic Karlsson)**: certified — exactly 48
  MU vectors; possible-edge graph is a perfect matching (24 edges, min
  non-edge margin 0.0112); clique number 2. **Theorem: no MU triple
  {I,H,K} for any H in the 1e-11 ball** (existence of an exact K3 member
  inside this ball is numerically overwhelming but not itself certified;
  the exact-map run below closes that gap).
- **K6(0.9, 2.1, 0.7) (exact Karlsson map, arXiv:1003.4177 transcribed and
  verified 200/200)**: certified — exactly 56 MU vectors (two
  ill-conditioned roots, σ_min down to 0.0104, verified via the local
  refinement fallback); possible-edge graph a perfect matching (28 edges);
  clique number 2; the thinnest non-edge margin, 1.3e-4 (this point sits
  near an edge-creation bifurcation), still certified positive since root
  enclosures are ~1e-9. **Theorem: no MU triple — a fortiori no MU
  quadruple — contains the exact Karlsson Hadamard K6(0.9, 2.1, 0.7)**
  (exact membership via Karlsson's construction: the transcribed map's
  float error ≲ 1e-13 lies well inside the 1e-11 certified ball). To our
  knowledge the first rigorous exclusion at a generic point of the
  non-affine (LP-resistant) territory.

- **Layer-1 tile at K6(5.9785, 4.0075, 1.6328)** (the well-conditioned
  point from a 6-point margin scan): certified with a fattened ball
  hslop = 3e-6 — 48 roots, perfect matching, clique 2, margin degraded
  only 0.0063 → 0.0059. **Theorem: no MU triple through {I,H} for every H
  in the 3e-6 ball, i.e. a parameter tile of radius ≥ 2.89e-6 in all
  three Karlsson parameters.** The Layer-1 mechanism of PROOF_ROADMAP.md,
  demonstrated end-to-end; the margin budget says ~10x fatter tiles
  already fit.

Engineering notes (they change the economics): the O(n²) suspect
clustering was replaced by grid-hash union-find (20 min → 1 s at 61k
suspects); the breadth-first sweep was rewritten as LIFO-chunked (100k
boxes per evaluation), which eliminated both the OOM kills *and* the
runtime — the earlier 12-min-to-2.5-h sweeps were swap-thrash, not
compute: **a full certified sweep is ~2 s of NumPy** (6.9M boxes,
peak RSS 447 MB). A pointwise certificate is therefore seconds, not
minutes — the full-family estimate in PROOF_ROADMAP.md drops accordingly.
Sweeps checkpoint suspects (`*_suspects.npz`) so verification reruns skip
the sweep.

## Result 7 — parametric-Krawczyk prototype (Layer 1) and its measured law

`parametric.py` prototypes tile certificates: "no MU triple through
{I, H(b)} for ALL b in a parameter box of half-width h." Design that
survived five iterations of failure analysis: per root, a slanted tube
theta*(b) = theta0 + S(b-b0) +- rho with S = -J^{-1} dg/db (kills the
first-order residual); a single ZONED sweep where each box's exclusion tax
is the minimum of the plain threading bound and any nearby root's slanted
bound (reach-dependent, per box); handoff to tube-Krawczyk boxes at the
contraction radius; first-order pair-overlap models; and a <=5-coloring
partition certificate (clique <= colors < 6). Map-side constants are
sampled+padded (EMPIRICAL — the rigorous version swaps in interval
enclosures of the Karlsson map, one Arb evaluation per tile).

Certified (prototype-grade) at K6(5.9785, 4.0075, 1.6328):
**h = 3e-4 in 26 s** (5.0M boxes; 48 tubes, worst radius 5.8e-4; 24
conflicts, 2 colors) — a 100x radius gain over naive threading (2.9e-6),
i.e. ~10^6x fewer tiles in 3 parameters.

The measured Layer-1 law: each root imposes TWO ceilings, and the binding
one is the tighter:

1. Tax-slope ceiling h_i ~ sigma_min,i / (2 sqrt3 PAD (Hess|S_i| + g_b)) —
   where the slant tax grows with distance as fast as |g| rises along the
   root's worst direction (worst roots here: 1.4e-3; median 5e-3..2e-2).
   Just underneath it, axis-aligned boxes pay (L/gain)^5 per octave
   (h = 1e-3 is inside the ceiling but costs ~10^7 boxes per bad root).
2. Tube-boundary ceiling — the Krawczyk radius is pushed UP by the O(h^2)
   curve residual (R_K >~ |Y| rad_g) while the slant tax squeezes the
   exclusion margin just outside it. At h = 6e-4 root 14 fails with
   exactly 42 stuck wmin-boxes at the tube boundary (measured residual
   margin ~3e-5): the all-in first-order ceiling at this point is ~5e-4.

Design lessons, each found by a distinct failure mode: (1) breadth-first
frontiers swap-thrash (fixed: LIFO chunks); (2) any nonzero parameter tax
makes the un-excludable band 5-dimensionally thick — never refine it to
wmin (fixed: inner-region collection); (3) spherical guard balls leak
along the worst singular direction (fixed: per-box reach-dependent slanted
tax — the zoned sweep); (4) uniform per-level cascade taxes overtax the
outer annulus (fixed: single-level zoned pass); (5) the near-root
Lipschitz constant must be the local 2|s|/6, not the global cap (a 2.4x
that fifth-powers); (6) the tube-boundary margin — not the tax slope — is
what actually binds first-order tubes.

**Second-order (Q-) tubes: built and confirmed (2026-07-26).**
`root_data2` tracks the curve Hessian Q by FD root-following; the
on-curve residual is sampled along the quadratic curve (h^2 -> h^3); the
zoned-sweep taxes take the sampled curve residual in place of the
analytic h^2 term. Ladder at the reference point:

    h = 3e-4   PASS  29 s   (regression)
    h = 6e-4   PASS  35 s   (failed with first-order tubes)
    h = 1e-3   PASS  93 s   (hopeless with first-order tubes)
    h = 1.4e-3 FAIL         (box blowout AT the predicted slope wall 1.39e-3)
    h = 2e-3   FAIL

Both ceilings behave exactly as derived: Q-tubes eliminated the
tube-boundary ceiling (~5e-4 -> wall) and the tax-slope wall stands.
Honest revision of the earlier "~10x ceiling, ~1000x fewer tiles" hope:
the measured Q-tube gain is ~2.4x in h (~13x in tile count) — the slope
wall, h ~ sigma_min/(2 sqrt3 PAD (Hess|S| + g_b)), is untouchable by
higher-order curve modeling since its constant comes from the e-drift of
the integrand, not curve error. The real levers are elsewhere (Result 8:
fold-delegation, 65,000x; margin-cached continuation, ~10^3x on per-tile
cost). Remaining smooth-part headroom: direction-resolved local taxes
(local |s| in place of the global Hessian row bound, ~2x) and anisotropic
tiles matched to per-direction sensitivities (~2x volume).

## Result 8 — family-wide ceiling statistics: folds are load-bearing

`scan_ceilings.py`, 24 random Karlsson points: per-root tax-slope ceilings
h_i = sigma_min,i / (2 sqrt3 PAD (Hess |S_i| + g_b)), whose min bounds the
Layer-1 tile size at that point.

    h_slope quartiles across the family:  1.8e-4 / 5.4e-4 / 1.2e-3
    range: 1.7e-5 .. 4.0e-3;  points below 5e-4: 12/24
    per-point MEDIAN-root ceiling: 3.9e-3 .. 8.7e-3 (uniformly healthy)

The reference point (1.4e-3) sits in the top quartile — it was lucky, not
typical. The bad tail is structural: points with h_slope ~ 1e-5..1e-4
have worst-root |S| = 26..79 and sigma_min = 0.003..0.007 — they are near
FOLD strata (root pairs created/annihilated as beta varies; consistent
with the observed 48/52/56/60 vector-count plateaus across the family).
Since S = -J^{-1} dg/dbeta ~ 1/sigma_min, the ceiling scales like
sigma_min^2 near a stratum, so the tile integral over a transversal
DIVERGES: uniform-in-kind Layer-1 tiles cannot cross fold strata at any
budget — campaign integral N ~ 4e14 (Q-tubes) / 1e16 (first-order) tiles.

Conclusion: fold-aware certificates (Miranda/normal-form tubes treating a
near-singular root PAIR as one object — the pair's discriminant varies
smoothly even when each root's Jacobian degenerates) are not Layer-2
polish; they carry ~half the parameter volume. The graph conclusion is
fold-robust: merging roots approach each other (mutual overlap -> 1, not
0), so a fold-tube enters the partition certificate as one color class
holding up to two near-identical vertices. With ill-conditioned roots
delegated to fold-tubes, the tile size is set by the well-conditioned
roots — per-point median-root ceilings of 4e-3..9e-3 — restoring a finite,
tractable campaign.

**8b — the fold-delegation trade curve** (same 24 points, 1280 roots;
delegate roots with sigma_min < cut to fold-tubes, tile ceiling from the
rest):

    cut 0.00:  folded   0/1280   h_med 5.4e-4   N ~ 4.0e14 tiles
    cut 0.01:  folded  19/1280   h_med 7.3e-4   N ~ 2.6e11
    cut 0.02:  folded  63/1280   h_med 1.2e-3   N ~ 2.7e10
    cut 0.03:  folded 150/1280   h_med 1.9e-3   N ~ 6.1e9

Delegating the worst 12% of roots wins a factor 65,000 in tile count.
At ~6e9 tiles: Python 2 s/tile is out; a C interval kernel (~10 ms) is
~months-on-a-cluster; margin-cached continuation (below) at ~0.1-1 ms
amortized brings it to WEEKS ON A DESKTOP. The remaining structural cost
lever, sketched in PROOF_ROADMAP.md: cache the sweep's per-box exclusion
margins once per region and re-certify tiles against stored margins
(m_box > L_beta * step), re-sweeping only when a chain of tiles exhausts
its floor — turning the 5M-box sweep from a per-tile cost into an
amortized one.

## Result 9 — valley tubes break the stratum wall (the divergence is gone)

Building the "fold tubes" of Result 8 produced a discovery and a working
certificate (`fold.py`):

**Discovery.** The ill-conditioned roots are (at every point examined)
NOT fold pairs: the bifurcation scalar phi along the singular direction is
MONOTONE through the root with slope ~ sigma_min, and no partner root
exists anywhere along the window. The structure is a long shallow VALLEY:
a 1-dim trench in |g| of depth < the exclusion taxes over >1 radian of
arc, along which the root races as beta varies (|S| ~ phi_beta/sigma).
Layer-1 tiles died not from folds but from trying to pay 5-dim sweep
taxes along an un-excludable 1-dim trench.

**Certificate** (valley window): rotate to the singular frame; solve the
4 well-conditioned directions by implicit function (residuals ~ 1e-16,
sigma_4 ~ 0.2-0.33 uniformly healthy); on the remaining scalar, certify
phi monotone with a single sign change over a window [-T, T] sized so the
trench floor at the ends clears the local beta-tax; endpoint signs clear
the tile's beta-drift => exactly one root in the window for every b in
the tile. The zoned sweep collects exactly the curved valley TUBE via a
per-root oracle (frame + y-polyline + tube radius; the thin-to-fat shell
is root-free by the sigma_4 bound); overlap rows for the partition come
from vectors sampled along the racing dip.

**Enabling fix — |s|-local beta taxes.** The global threading tax
(|s| <= 1) can never be cleared by a trench floor. But near MU points
|s_k| ~ 1/sqrt6, and the correct local bound |dg_k| <= 2|s_k| ||dh_k||
|dbeta| is 2.4x smaller — implemented per-box per-component in the zoned
sweep (benefits all exclusion, not just valleys).

**Demo at scan pt 8** (worst point: sigma_min = 0.0028, |S| = 78, Layer-1
wall h = 1.68e-5): **TILE CERTIFIED at h = 3e-4 in 36 s** — 18x beyond
its Layer-1 wall — with 5 valley windows (T = 0.35..1.05, all monotone)
+ 47 Q-tubes, 27 conflicts, 2 colors. At h = 1e-3 only the deepest
valley (sigma 0.0028) fails: its trench floor cannot clear the h-scaled
tax by T = 1.6. Clearly-scoped fix: per-root beta-subdivision of the
valley certificate (a ~3 s 1-dim object; 8 subtiles ~ 24 s for that one
root), plus longer windows. With valleys handled, the binding ceilings
revert to the well-conditioned roots (Result 8b: h ~ 1.6-2.4e-3 typical),
and the stratum divergence that made uniform tiling impossible is gone.

## Result 10 — both research risks de-risked; the 1/6 wall is universal

The two items where the K3 program could still have died (PROOF_ROADMAP
gap list R1, R2) both fell in one session, with maximal margins:

**R2 — Layer 3 (X-triple unextendability), certified.** `layer3.py`
constructs the exact triple {I, B(1.6), K} (Hermitian family point, B in
X6^(2); 120 MU vectors, unique basis K, triple exact to 2.5e-16) and
sweeps the phase torus against BOTH Hadamards' constraint stacks (10
components): **908k boxes, 0.4 s, zero suspects** — certified: no vector
is MU to the triple, for any matrices in 1e-9 balls. The 10-constraint
overdetermination makes Layer 3 cheaper than a pointwise certificate.
Measured unextendability margin over the 114 non-basis MU vectors:
**min 0.1664 = 1/6** — maximal.

**R1 — Fourier shell (quadruple exclusion where triples exist),
certified at F6.** `fourier_shell.py` enumerates ALL 6-cliques of the
certified possible-edge graph: exactly 16 — Grassl's bases, no spurious
cliques — and certifies for every clique pair a lower bound on the max
cross-deviation from 1/6: **0.166667** (the 1/6 wall, now a certificate
rather than a numerical observation). Drift along the Fourier family:
0.166667 -> 0.166177 over eps = 0.01 — rate ~ 0.05/unit-parameter, so
the shell certificate tolerates O(0.1)-thick tiles; the wall is
effectively indestructible near F.

**Observation worth a theorem hunt: the same 1/6 = 1/d margin appears in
all three obstructions** — Grassl's F6 wall, the X-triple's extension
gap, and the base-pair deviations. The deviation spectrum at the
symmetric points is quantized in units of 1/6, which is exactly the
integrality structure of the Matolcsi-Weiner Fourier-analytic linear
constraints — suggesting the campaign's measured walls and the
LP-integrality route are two views of one mechanism, and an exact-proof
bridge may exist.

## Result 11 — margin-cached continuation works: 20x sweep amortization

`cache.py` implements the amortization lever sketched in the roadmap.
The anchor sweep stores every box excluded by the |s|-local beta tax
with its exclusion EXCESS and per-unit-distance drift RATE (2.07M boxes,
95 MB float32; boxes excluded by frame-dependent root-slant taxes —
260k — are stored verbatim for per-step re-sweeping). A chain step to a
nearby tile center then re-certifies the entire cached far field with
ONE vectorized comparison E > R * dist — measured at 4 MILLISECONDS —
patches the failing minority plus the slant boxes with a mini-sweep, and
re-runs the per-root machinery (tubes, partition) at the new center with
anchor-fixed guards padded by root motion.

Measured chain at the reference point, h = 3e-4, steps of 1.6h along
theta: **12/12 steps certified from one anchor** — a certified sausage
of 13 overlapping tiles, theta-half-length 6.06e-3 (20.2x one tile's
sweep volume) — with per-step cost 3.5-10.1 s vs 26 s standalone: the
step cost is dominated by the patch sweep, whose size grows linearly
with distance (18% of cache at step 1, 88% at step 12), mapping the
natural re-anchor cadence at roughly the point where patch ~ full sweep.
Per-root FD sampling (~1-2 s/step) becomes interval evaluation in the C
substrate; the vectorized check is the part that stays O(ms). In the
campaign this collapses per-tile sweep cost by ~10^3 for chained
traversal, exactly as the roadmap estimated.

**Deep-valley anisotropy (design decision recorded).** The h=1e-3
failure mode at the deepest valleys is edge-clearance: the trench floor
must beat a tax ~ h while the frame caps the window length. Beta-
subdivision does NOT fix this (the window-edge requirement is set by the
full-tile sweep tax); the correct fix is ANISOTROPIC tiles thin along
the 1-dim racing direction phi_beta (the stratum transversal) and fat in
the two transverse directions (which do not move the trench to first
order) — restoring full-size tile volume near strata. Slots into the
adaptive driver as a per-tile h-vector; plumbing change, no new theory.

## Result 12 — parametric Layer 3: certified theta-intervals in seconds

`layer3_param.py` extends the X-triple unextendability certificate to
parameter INTERVALS: sampled drift rates |dB/dtheta| = 0.408 and
|dK/dtheta| = 0.389 (K tracked by warm-started column polishing; triple
defect stays 1.9e-16), ball radius PAD*(rates)*delta threaded through
the 10-component sweep. Delta-ladder outcome: **theta in [1.594, 1.606]
certified in 1 s** (3.5M boxes, zero suspects) — no vector is MU to any
triple {I, B(theta), K(theta)} in the interval. Drowning threshold at
sweep slop ~ 0.096 vs the 1/6 margin (the binding neighborhoods are the
120 MU vectors of {I,B}); the |s|-local tax refinement (as in the zoned
sweep) would widen the interval ~2x. At 0.012 per 1-second anchor, the
full B-circle is ~524 anchors (~minutes); the 2-parameter X-family at
this rate is ~10^5 anchors (~days in Python, hours in C) — Layer 3 is
campaign-ready.

## Result 13 — a real fold, found and certified over (R3 closed for pairs)

`fold_hunt.py` bisected the MU-vector count along a path between scan
points (48 <-> 60, with plateaus 48/52/56/60 crossed en route — fold
sheets are plentiful): a birth event pinned at path parameter
s = 0.29822, where the count jumps 48 -> 52 — **two root pairs born
simultaneously** (a symmetry orbit). Just past it: pair separations
0.039/0.055, member sigma_min 0.0016-0.0025, and pair-vector overlaps
0.9998 — near-parallel, as fold theory demands.

The clique-soundness insight that makes pairs cheap: all vectors in a
merged fold-dip are mutually non-orthogonal (certified via the dip's
span-sample self-overlap), so ANY clique uses at most one of them — the
single-collapsed-vertex bookkeeping with worst-case overlap rows is
sound for pairs with NO new enclosure logic. Implemented as a certified
self-overlap floor (reject if < 0.05; measured ~0.99).

**Certified: tiles at the fold's doorstep** (52 roots incl. both
near-born pairs) at h = 1e-4 AND h = 3e-4, ~35 s each — six valley
windows (the four pair members flag non-monotone, i.e. the windows see
the genuine dip structure with the partner inside) + 46 Q-tubes,
2-colored. Fold-pair tiles cost the same as ordinary tiles.

**Birth-straddling tiles: CERTIFIED (R3 fully closed).** Phantom
anchors implemented as a two-stage sweep: stage A runs coarse
(wmin 0.02) and diverts unresolved boxes; far-from-root clusters are
LM-polished to their local |g|-minimizers, and those with certified-
small depth (below the tax scale) become phantom-anchored valley
windows — their dips certify wherever 0, 1, or 2 roots may live across
the tile, with no center root required; shallow pockets and near-guard
shells go to stage B, which resumes exactly on the diverted boxes at
full resolution with the augmented oracles. Measured at a tile centered
on the 48-root side with the birth 2.9e-4 away (inside the tile):
stage A diverted 40k coarse boxes (62 clusters, 13 far); BOTH birth
sites found (|g|min = 1.14e-4 each — the two pairs of the orbit) and
phantom-windowed (T = 0.70, non-monotone dips); 11 shallow pockets
correctly triaged to stage B; **certified in 33 s** — the same cost as
an ordinary tile. Root count changing inside a tile is now routine.

## Result 14 — the rigor substrate stands up (R4 cornerstone)

`interval.py` (directed-rounding real intervals + complex rectangles;
libm faithful-rounding assumed at <= 2 ulp, padded 4 — the one
documented gap to swap for Arb/crlibm), `ivkarlsson.py` (certified
Karlsson map over parameter boxes with explicit branch obligations),
`dual.py` (forward AD over the interval types -> certified first
partials and mean-value tile bounds). Findings:

1. **The interval map is correct**: 200/200 float evaluations contained
   in point-width enclosures (width 1.9e-11).
2. **The SLOP model is validated by computation**: certified interval
   g-evaluations have width <= 1.6e-12 against the assumed 1e-11 —
   the "static rounding-error lemma" is now a measurement, and in the
   production checker the interval evaluation simply replaces the slop
   model outright.
3. **Branch tracking is a non-issue at campaign scale**: zero
   cut/denominator failures across all 24 scan points at h = 3e-4;
   min |den|^2 = 0.255, min cut margin 0.025.
4. **Naive interval tile enclosures are unusable (32-39x dependency
   blowup) — and the mean-value form via dual AD fixes it**: at
   h = 3e-4 the certified bound is 7.3e-4 vs the prototype's empirical
   1.08e-3 — ratio 0.68, i.e. the CERTIFIED constant is TIGHTER than
   the PAD*FD guess it replaces (PAD = 2 confirmed conservative). The
   derivative enclosure's own artifact re-enters linearly (ratio 2.4 at
   h = 1e-3, blown at 3e-3); the standard cure — evaluating the
   derivative over 2^3-subdivided boxes, at microseconds per map eval —
   halves the artifact per level and keeps certified taxes at
   prototype magnitude through campaign tile sizes.

With this, every EMPIRICAL constant in the certificates has a
demonstrated certified replacement path: map values and Lipschitz taxes
(this result), g/J evaluations (interval, widths ~1e-12), curve
residuals (interval evaluation at corner points), and the elementary-
function faithfulness assumption isolated as the single remaining
substitution (Arb) for publication grade.

## Result 15 — the flagship theorems, upgraded to certified

**Exact-point theorem, now self-contained.** The interval map certifies
that the EXACT Karlsson Hadamard K6(0.9, 2.1, 0.7) lies within
**5.9e-14** of the float matrix the Result-6 certificate swept (branch
obligations verified: |den|^2 >= 3.35, cut margin >= 0.032) — well
inside its certified 1e-11 ball. The "float error ~1e-13 (estimated)"
caveat is gone: no MU triple contains the exact K6(0.9, 2.1, 0.7),
with every constant in the chain certified (modulo the libm-
faithfulness assumption isolated in interval.py).

**Tile theorem, certified radius.** For the 3e-6-ball tile at
K6(5.9785, 4.0075, 1.6328): anchor point ball 1.9e-13, certified
mean-value rate sum 1.444 (dual AD over the box) => **certified
parameter radius r = 2.0e-6** in all three Karlsson parameters
(dev = 1.444 r + 1.9e-13 <= 3e-6). The FD-claimed 2.89e-6 takes a 31%
haircut and becomes rigorous.

**Static rounding-error lemma (replacing SLOP by derivation).** Under
IEEE-754 doubles (eps = 2^-53) with faithfully-rounded libm sin/cos
(<= 2 ulp): the sweep's computed g~_k satisfies |g~_k - g_k| <= 1e-13
for all inputs. Sketch: the u-entries carry <= 3.5 eps; the 6-term
complex dot has error <= gamma_14 * sum|H||u| + propagation <= 22 eps
= 4.9e-15 (sum|H_jk||u_j| = 1); the modulus-square step brings the
total below 2e-14; 5x rounding-up gives 1e-13. The margin arithmetic
adds <= 1e-14. Hence SLOP = 1e-11 over-covers by >= 100x. The same DAG
shape bounds the analytic J entries and the Krawczyk residual products
(<= 1e-13 each); the approximate inverse Y needs NO accuracy bound —
the Krawczyk operator is valid for arbitrary Y, so LAPACK stays a
black box. Measured cross-check: certified interval g-widths 1.6e-12
(Result 14), consistent with the derived bound's headroom.

## Result 16 — the pipeline runs end-to-end: a certified 3-D block

`driver.py` ran the first mini-campaign: 16 margin-cached chain-lines
(4x4 grid in (phi, lam), spacing 1.8h) x 7 tiles each along theta
(step 1.6h), h = 4e-4. Outcome: **112 tiles, 16 anchors, 17.7 minutes,
zero failures, zero mid-line re-anchors** (96/96 chain steps OK).

    THEOREM (prototype-grade). For every beta in
      theta [5.978103, 5.982743] x phi [4.007135, 4.010095]
                                 x lam [1.632365, 1.635325]
    no MU triple {I, H(beta), K} exists — hence no MU quadruple
    contains any Karlsson Hadamard H(beta) from this box.

Extents 4.6e-3 x 3.0e-3 x 3.0e-3, volume 4.1e-8 — 79x one tile, the
first certified region of full dimension. Measured throughput 9.5
s/tile all-in (anchors amortized). Extrapolation at this uniform h:
~1800 CPU-years for the family in Python — which is precisely the gap
the remaining engineering closes: adaptive h (~50-100x fewer tiles),
the C kernel (~100-1000x per tile), anisotropy near strata; the product
lands at the weeks-on-a-desktop scale the roadmap projects. Every stage
of that projection is now backed by a measured run.

## Result 17 — certified taxes wired into the live certificate

`rates.py` derives every sweep-tax constant rigorously from the dual-AD
map enclosures (derivations in the docstring): per-COLUMN beta taxes
BR_k = (2/sqrt6) h sum_j ||dH_{.k}/db_j||_1 * |s_k|-local, the mixed-
derivative gb bound for the slant-tax slope, the Krawczyk J-drift, and
the |s| box-drift — no sampling, no PAD. Measured at the reference tile
(h = 3e-4): **the certified beta/far taxes are 2.7x TIGHTER than the
prototype's padded FD constants** (ratio 0.37 — the honest constants
beat the guessed ones); only the gb bound is fatter (7.8 vs sampled
0.3-0.7, being worst-case over all components; sharpenable with
box-local |s|).

Wired into `certify_tile(use_certified=True)`:
- **h = 3e-4 CERTIFIED with fully certified taxes in 25 s** — faster
  than the empirical run (sweep 3.8M boxes vs 5.0M). The fatter gb
  re-routed four moderate roots into valley windows automatically; the
  machinery composed correctly around the new constants (44 tubes + 4
  valleys, 2-colored, no phantoms needed).
- h = 1e-3 fails at a valley y-shell margin: the certified ceiling at
  this point sits ~25% below the empirical wall, as the honest gb
  dictates. Remaining EMPIRICAL after this wiring: curve-residual
  corner sampling, root sensitivities Sn/Q, overlap gradients — each
  with a designed certified replacement (Result 14 / roadmap R4).

## Result 18 — the 3-D block re-certified on certified constants, faster

The full Result-16 campaign re-ran with `use_certified=True` end to end
(rates threaded through anchors, chain steps, and the cache): **112
tiles, 16 anchors, 10.6 minutes — 40% faster than the empirical run —
96/96 chain steps, zero re-anchors, zero failures**, same certified
box (4.6e-3 x 3.0e-3 x 3.0e-3, 79x one tile).

Why rigor made it FASTER: the certified per-column beta taxes are 2.7x
tighter than the padded FD constants, shrinking sweeps and patches; and
they exclude everything the frame-dependent root-slant taxes used to
handle — the "+slant" re-sweep category dropped to ZERO boxes, so chain
steps fell to 2.1-3.5 s. Soundness across the chain span is kept by
sup-ing the per-column unit drift rates over a train of h-sized
sub-boxes along the chain line (tight constants per sub-box, coverage
by the max — `chain_certified_rates`).

The campaign's flagship demonstration now runs on certified tax
constants throughout. Remaining EMPIRICAL in the loop: curve-residual
corner sampling, Sn/Q sensitivities, overlap gradients (roadmap R4
tail), plus the libm-faithfulness assumption under everything.

## Result 19 — the first family-complete certificate: the whole B-arc

`layer3_family.py` + `layer3_patch.py`: an adaptive anchor walk (delta
grows 1.3x on success, halves on drowning, gap-and-skip at the floor;
K warm-tracked between anchors and the triple re-verified at each)
certified, across ~540 anchors and ~50 min in three runs:

    THEOREM (prototype-grade). For every theta in [1.25, 5.03] — the
    Hermitian (Beauchamp-Nicoara) arc, minus the endpoint segments that
    degenerate to Bjorck's C (excluded pointwise, Result 5) — no vector
    is mutually unbiased to the triple {I, B(theta), K(theta)}: the
    ENTIRE Hermitian family of MU triples is strongly unextendible.

Zero gaps after patching. Discoveries made by the walk itself:

1. **The arc crosses two Dita-equivalent special points** — theta = pi
   and theta = 3pi/2 (y = -1 and y = -i in the family chart, exactly
   the special values Bengtsson et al. list) — where the MU-vector
   count spikes (108-120 vs 56-72 generic) and certification margins
   tighten (delta ground to ~1e-3 through pi).
2. **K-tracking silently loses its branch at a special point** (at Dita
   the partner-basis structure reorganizes: 10 bases live there). The
   first pass through 3pi/2 produced a spurious rate blowup (29.7) and
   a 0.15-rad honest-gap trail from a broken tracked K; rebuild-mode
   anchors (fresh multistart triple per anchor) crossed the same zone
   at rates 0.9-1.4 with delta = 2e-3, zero gaps, in 6 min. CAMPAIGN
   LESSON: continuation must re-anchor by rebuild at special loci —
   tracked objects can lose their branch exactly where the structure is
   richest, and the certificate correctly refuses rather than lies.

With Result 12 (theta-interval mechanics) this completes Layer 3 for
one of the three known infinite triple families; the remaining X6^(2)
territory is the same walk over Szollosi's two-parameter chart.

## Result 20 — the B-arc closed end to end (Bjorck to Bjorck)

`endpoints.py` completed Result 19's arc:

- **Bjorck's C certified**: 54 MU vectors, unique partner basis K
  (triple defect 1.4e-9), sweep 842k boxes / 0.3 s / zero suspects at
  ball radius 5e-9 — which contains the exact triple. Result 5's
  multistart observation is now a certificate.
- **Both endpoint segments walked with zero gaps** (rebuild-mode, 26
  anchors each, ~4.5 min): coverage ground to theta = 1.1965 and
  5.0858, with drift rates rising to ~15 at the last anchors — the
  chart's sqrt-branch divergence, handled by delta ~ 6.6e-4.

Combined statement: for every theta in [1.1965, 5.0858] — plus the
endpoint Bjorck matrices themselves — no vector is MU to the triple
{I, B(theta), K(theta)}. The two residual slivers at the branch points
(total width < 1.8e-3, where the chart's rates diverge) are the honest
boundary term; closing them rigorously wants the standard cure of
reparametrizing by the branch variable (theta = theta_min + s^2), noted
for the substrate pass.

## Result 21 — Layer 3 over the Dita arc (multi-triple anchors)

`layer3_dita.py`: the same rebuild-mode walk over the Dita family's
fundamental domain x in [-1/8, 1/8] — the second explicit slice of
X6^(2), and the first with MULTIPLE triples per point: **118 anchors,
up to 4 bases each (every basis certified per anchor), zero gaps,
32.3 min**:

    THEOREM (prototype-grade). For every x in the Dita arc, every MU
    triple {I, D(x), K} is strongly unextendible.

The walk crossed the maximally rich D(0) region (120 MU vectors)
without special handling — the multi-basis loop simply certified each
of the 4 surviving bases (the 10-basis structure lives only exactly at
the special point, as with the Fourier family's 16 -> 8). With Results
19/20 (the B-arc, Bjorck to Bjorck), both explicit 1-parameter slices
of the Szollosi set are now family-complete; the remaining Layer-3
territory is the generic 2-parameter X6^(2) interior, which needs the
Szollosi-chart transcription (eq 7.14 + Zauner's circulant triple
construction) and then this same walker in two dimensions.

## Result 22 — the Szollosi chart, and Layer 3 into the generic X-interior

`szollosi.py` transcribes eq 7.14 with the cubic-root rule ((x,y) roots
of f_alpha, (u,v) of f_{-alpha}, alpha in the deltoid intersection D):
**validated 40/40 sampled alpha at unitarity defect 5e-15**, one
canonical root-pair choice throughout; rows 3-5 of the transcription are
identically orthogonal to the flat row, a structural self-check. The
last family chart of the program is implemented.

`layer3_x.py` then took Layer 3 into the generic X-interior:
- **8/8 spot certificates** at random interior alpha (60 MU vectors,
  unique triple each, all certified unextendible);
- **a chain line across D certified: 241 anchors, 33.7 min**, with 2
  gap anchors at t ~ -0.2596 where even 40k multistarts find only
  54/59 of the 60 vectors — a root-birth on the line collapsing the
  K-columns' basins — **closed by continuation**: the neighbor's K
  tracked across the gap verifies to 2e-16 and certifies (zero
  suspects).

Failure-mode symmetry worth recording: on the B-arc, TRACKING lost its
branch at a special point and rebuild-mode saved coverage; on the
X-line, REBUILD (enumeration) went blind at a root-birth and tracking
saved it. The campaign walker must run both modes and cross-verify —
each mode's blind spot is the other's strength.

Full-D coverage at the measured anchor economics (~5.6e-3 spacing,
8.4 s/anchor): ~1.1e5 anchors ~ 11 days of Python, ~3 h with the C
kernel — a campaign line-item, now costed. With Results 19-22, Layer 3
has touched every known triple territory: Fourier (literature),
Hermitian arc (complete), Dita arc (complete), generic X (spot + line +
costed cover).

## Result 23 — the full X-family covered, on cloud GPUs

The program's first cloud campaign. `gpusweep.py` ports the certified
triple sweep to CuPy (array-module-generic; box-for-box identical
verdicts to the NumPy path, and the static rounding-error lemma applies
unchanged since GPU FP64 is IEEE-compliant). Measured on a Modal A100:
**0.14-0.24 s per sweep vs 1.7-8.7 s CPU — 12-36x**, ~33M boxes/s.
`modal_app.py` then mapped the X-cover line walker across Modal
containers (10-GPU account cap; walker on container CPU, sweeps on
GPU): **115 lines dispatched, all lines inside the deltoid region
walked in ~35 min wall-clock** — the same job the laptop had projected
at 21 hours. Rough spend: ~$20 of A100 time.

Combined with the 10 local lines:

    THEOREM (prototype-grade, modulo slivers). Across the ENTIRE
    Szollosi region D — 95 chain lines, 7,412 certified anchors,
    line spacing 2.5e-2, vertical extent |Im alpha| < ~1.18 — every
    anchored triple {I, X(alpha), K(alpha)} is strongly unextendible.
    27 lines carry isolated gap slivers (~1-2 each, ~2e-3 wide in
    Re alpha; total measure < ~0.1% of D), concentrated symmetrically
    near the axis band and the sub-family crossings — root-birth loci,
    patchable by the Result-22 continuation pass.

With this, Layer 3 has certified anchors across every known triple
territory in dimension six: the Fourier family (literature), the
Hermitian arc (complete), the Dita arc (complete), and now the full
two-parameter Szollosi region (complete modulo slivers). The
Layer-3 hole in the eventual K3 theorem is reduced to sliver-patching
plus the anchor-interval-to-full-cover bookkeeping.

GPU lesson for the campaign: 40 lines of CuPy bought 12-36x on the
dominant kernel with zero rigor loss; the walker's CPU parts now
dominate anchor cost, and batching those (many anchors' polishes per
GPU call) is the next order of magnitude.

## Result 24 — the X-cover closed gap-free; Layer 3 complete

The Modal patch pass (27 gap lines re-walked with delta_min 2e-4, 20k
rebuilds, 5-anchor cadence; ~$8) collapsed the 27 sliver-bearing lines
to THREE residual points; local continuation anchors (11 per sliver at
delta 1.3e-4, tracked-K) then closed all three — 33/33 anchors, zero
residual. **The Szollosi region D is covered by certified anchors with
no gaps.**

The three stubborn loci were identified before closing: alpha = 0 is a
FOURIER-FAMILY crossing (fingerprint 36/36 sixth-roots with the -1-rich
F structure — the F and X families intersect at the region's center),
and alpha ~ +-(0.874 + 0.65i) are rational special points of the same
kind (14/36 sixth roots). The walker had stumbled on enumeration
artifacts at the crossings; the certificates themselves pass cleanly.

**Layer-3 status: COMPLETE across every known triple territory in
dimension six.** Fourier family — literature (Thm 7.7). Hermitian arc —
Results 19-20 (Bjorck to Bjorck). Dita arc — Result 21. Szollosi region
— Results 22-24, gap-free. Remaining bookkeeping for the eventual
theorem: the anchor certificates cover overlapping parameter intervals
by construction; promoting interval-chains to a stated continuum
theorem is assembly, not computation. Total cloud spend for the entire
X-campaign: ~$29.

## Result 25 — Taylor-model curve residuals: the last big EMPIRICAL
## category certified

`tm.py` (degree-2 Taylor models in the three tile variables with
rigorous remainders — coefficient roundoff Gamma-folded, elementary
series with explicit tails, Moebius/csqrt via bounded expansions) and
`tmres.py` (the Karlsson map over TMs, shared per tile, plus the
g-along-the-curve composition) replace `curve_residual`'s corner
sampling + PAD with a certified bound:

- **Soundness simplification discovered en route: S and Q need no
  certification at all** — they merely define the tube curve; rigor
  lives entirely in the certified residual along the defined curve. The
  FD-status of the sensitivities stops being a rigor gap.
- Certified R(h) scales exactly as h^3 (measured 37x from h = 3e-4 to
  1e-3, matching (10/3)^3), runs in 0.01 s/root — FASTER than the
  sampling it replaces — and lands 4+ orders below the binding taxes
  (2e-8 vs 1e-3 at h = 3e-4), so certificates lose nothing despite the
  bound's 15-400x headroom over sampled values.
- Wired into `certify_tile(use_certified=True)`; the reference tile
  re-certifies in 24 s with TM residuals and no PAD on rad_g.

**And the sole survivor fell in the same session**: `certified_overlap_lo`
bounds every Q-tube pair's overlap over the tile by TM inner products of
the u-curves (0.15 s for all ~1100 pairs), replacing the FD overlap
gradients + PAD; fold/valley rows keep their span-sampled bounds.
Regressions: the reference tile (24 s) and the fold-doorstep tile
(31 s, valleys + TM path coexisting) both re-certify.

**Status: for Q-tube content, the certified tile path is free of
sampled constants end-to-end** — taxes (dual-AD), curve residuals (TM),
pair overlaps (TM), map values (certified enclosures), rounding
(derived lemma, interval-validated). Remaining assumptions: libm
faithfulness (Arb swap), and the fold/valley windows' 1-dim sampled
envelopes (same TM treatment available, not yet built).

## Result 26 — the symmetry quotient measured: divide the campaign by 32

`quotient.py` discovered the Karlsson parametrization's equivalence
group numerically (Haagerup-invariant screening + full-multiset
confirmation at multiple generic points):

- **Period pi in every axis**: theta -> theta+pi, phi -> phi+pi,
  lam -> lam+pi are all Hadamard equivalences (2^3 = 8x);
- the central inversion (-theta, -phi, -lam) and the reflection
  (theta, pi-phi, pi-lam), closing with the translations into a
  point-group of 4;
- **every H(beta) is equivalent to its own conjugate AND transpose**
  (so the anti-linear quadruple-statement symmetries are already
  internal — structurally notable: the family is self-dual under both).

Orbit size at generic points: **32** — the fundamental domain is 1/32
of the naive parameter box, and every campaign tile-count divides
accordingly: ~6e9 tiles -> **~1.9e8**, which at the GPU-measured sweep
economics projects to **~10-30 GPU-days ~ days at the 10-GPU cap ~
$500-1500** for the main K3 campaign — down from CPU-years at the
session's start. (Fundamental-domain boundary bookkeeping — half-open
faces under the group action — goes with the adaptive driver.)

## Result 27 — anisotropic tiles: 3x volume at the deep-valley wall

The h-vector is now plumbed through the whole certified stack: `tm.py`
Taylor models carry per-variable half-widths (per-axis series
remainders), `rates.py` weights every certified rate per-direction
(BR_k = (2/sqrt6) sum_j c1[j,k] hv_j; RJ and s_drift likewise), and
`certify_tile` accepts an hv vector on the certified path — slope terms
use hv-weighted row sums sum_k |S_ik| hv_k, q_offset is
0.5 sum |Q| hv_k hv_l, curve-residual and overlap TMs get the true box,
valley corners broadcast. All formulas reduce *exactly* to the scalar
ones at isotropic hv: the reference tile reproduces its prior
certificate bit-for-bit-comparable (48 roots, 3.84M boxes, 24 s).

Demonstration at the deep-valley scan point beta = (1.386501, 1.678339,
4.019655) (52 roots, 12 near-singular, racing axes *mixed* across
valleys — beta1 for some, beta2/beta0 for others):

- isotropic ceiling under certified taxes: **h = 3e-4 certifies;
  6e-4 fails** (valley 20 y-shell margin -1.8e-3);
- the cheap direction is the weakly-coupled beta0 (|u.Gb|_0 = 0.06 vs
  0.23 on beta1 for the binding valley 20);
- **hv = (9e-4, 3e-4, 3e-4) CERTIFIES: 3.0x the isotropic ceiling
  volume**; (1e-3, 3e-4, 3e-4) fails by only -2.3e-4 on the same
  margin, and trading beta1 down does not rescue larger beta0 —
  pinning the binding term as the *global* far tax
  sum_j c1[j] hv_j inflating rho_y, not the corner spread.

Campaign consequence: deep-valley strata (the tile-count-dominant
regions) get ~3x fewer tiles from shape alone; the identified next
lever, if ever needed, is a trench-local per-column valley tax in
place of the global max. The adaptive driver should pick hv from the
center-point diagnostics (|u.Gb| axis profile of the flagged roots),
which cost one SVD sweep per tile — negligible.

## Result 28 — zoned tile sweep on GPU: box-exact, 41x

`gputile.py` ports `parametric.zoned_sweep` (the certify_tile inner
engine: zones, per-root slant taxes, guard/oracle collection,
valley-tube membership, chain cache emission) to the xp-generic
NumPy/CuPy pattern of `gpusweep.py`. Validation is EXACT, not
statistical: replaying the reference tile's captured stage-A inputs,
xp=numpy reproduces the CPU sweep's (stuck, D0, total) to the bit
(3595 / 3,813,469 boxes, D0 diff 0.0), and the A100 CuPy path
returns the *same exact counts* — the FP64-IEEE rigor-identity
claim observed, not assumed, on the full tile machinery including
oracle interpolation.

A100 timings (Modal, reference tile stage A): container CPU 52.9 s,
GPU cold 91 s (JIT, once per worker), **GPU warm 1.28 s = 41x** —
at the tile pipeline's inherited 20k chunk; larger chunks raise it
further. Ops note: Modal calls from this Mac get their in-flight
input cancelled when the local client's connection drops (~3 min at
night — the first "failures" were also container OOM at default
memory, fixed with memory=16384); campaign dispatches must use
`.spawn()` + app-log/ledger collection, never a live `.remote()`.

## Where a proof (not a counterexample) might come from

Per the review's closing strategy list: Lasserre/SDP hierarchies (memory
wall ~2^32k bits at level k for the QRAC route; symmetry reduction helps —
Gribling–Polak got d+2 exclusion for d ≤ 8), the Jaming-style rigorous
discretization generalized beyond the Fourier family, positive-definite
functions à la Matolcsi–Weiner, or completing the classification of order-6
Hadamards (Conjecture 7.1: S6 ∪ K6^(3) ∪ G6^(4) is everything) and then
excluding family by family. The acknowledged weak spot: the LP bound cannot
yet kill complete sets built from K6^(3).

## Result 5 — strong unextendibility of triples (capstone)

For two explicitly constructed MU triples — the Grassl-type {I, F6, B0} and
the Björck triple {I, C, K} (the latter recovered only after the Result-4
polishing fix) — searches for even a *single* vector unbiased to all three
bases find **zero solutions** (8000 starts each). This reproduces, on
concrete triples, the constellation-level finding ({5^3,1} never found)
that MU triples in d = 6 are strongly unextendible: the gap to a fourth
basis is not "six missing vectors" but *one* missing vector, six times over.

## Files

- `mub.py` — library; `grassl.py` — Result 1; `search4.py` — Result 2
  (logs `run_d6.log`, `run_d7_long.log`, `run_d6_fixF.log`, `best4_*.npz`);
  `families.py` — Result 3 (`run_families.log`, `families_results.npy`);
  `karlsson.py` — K6^(3) probe + exact Karlsson map (arXiv:1003.4177);
  `degenerate.py` — Result 4 machinery; `certify.py` — Result 6 pipeline
  (logs `cert_*.log`; targets `K3_seed3003.npy`, `K6_0.9_2.1_0.7.npy`);
  `PROOF_ROADMAP.md` — the full-family program;
  `f160_vectors.npy` — lattice-enumerated vectors at F(1/6,0).
- Key sources: McNulty–Weigert review arXiv:2410.23997 (accepted Quantum,
  2026-03); explicit matrices from Bengtsson–Bruzda–Ericsson–Larsson–Tadej–
  Życzkowski quant-ph/0610161; Grassl quant-ph/0406175; Raynal–Lü–Englert
  arXiv:1103.1025; Brierley–Weigert arXiv:0901.4051.
