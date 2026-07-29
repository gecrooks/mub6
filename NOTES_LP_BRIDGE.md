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

## 5. What would falsify the bridge

- R-B LPs with optimum strictly below the certified gap (relaxation
  too weak; the walls are genuinely nonlinear).
- A certified point with defect wall NOT at a 1/6 multiple (none
  found in ~30 results, but the survey's h<3e-4 territory has not
  been wall-measured).
