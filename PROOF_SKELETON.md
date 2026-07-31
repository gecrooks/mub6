# Proof skeleton: no MU quadruple contains a Karlsson Hadamard

R10 deliverable. The complete logical chain of the computer-assisted
theorem, each lemma tagged with its code counterpart and rigor status
(CERT = certified, no sampling/PAD; SLOP = float + rounding lemma;
PROTO = sampled + PAD, marked for upgrade). This document is the map
a referee (or a future session) reads first.

## 0. Statement

**Theorem (target).** No four mutually unbiased bases of C^6 contain,
after the standard equivalences, a complex Hadamard matrix from the
Karlsson family K6^(3) (equivalently: any H2-reducible Hadamard).

Reduction (standard, literature): a MU quadruple {I, H, ...} may be
normalized so H is dephased; if H is H2-reducible it lies in the
3-parameter Karlsson family beta = (theta, phi, lam) in [0, 2pi)^3.
A fourth-basis vector must be unbiased to {I, H}, i.e. a MU vector:
u on the phase 5-torus with g_k(u) := |<u, h_k>|^2 - 1/6 = 0 for all
k. A quadruple needs 6 mutually orthogonal MU vectors twice over; it
suffices to show no set of 6 mutually orthogonal MU vectors exists —
and in fact the certificates prove the stronger clique bound
omega(G_MU-orthogonality) < 6 for every beta.

## 1. The certificate for one tile

Fix a tile T = {beta0 + db : |db_j| <= hv_j}. The certificate is a
finite object whose validity implies: for every beta in T, the MU
vectors of H(beta) admit no 6-clique of mutual orthogonality.

Layers, in dependency order:

L1 (map enclosure). Certified interval/dual-AD enclosures of H(beta)
    and dH/dbeta over T; branch obligations (csqrt cut clearance,
    denominator margins) are checked, not assumed.
    [ivkarlsson.py, dual.py, interval.py — CERT; transcendentals via
    mpmath.iv (Result 30), sqrt IEEE; no libm.]

L2 (rates). From L1: per-column beta-tax rates BR_k, unit-drift BU_k,
    gb-rate, Jacobian drift RJ, s-drift — sup-bounds over T.
    [rates.py — CERT (Result 25/27).]

L3 (global sweep). LIFO-chunked branch-and-bound over the 5-torus.
    A box is discarded only if (a) some |g_k| at its center exceeds
    the certified theta-Lipschitz spread + beta tax (+SLOP), or (b)
    it lies inside a declared collection zone (root guard, Q-tube
    oracle, valley oracle). Anything else splits; boxes below wmin
    are STUCK and fail the certificate.
    [parametric.zoned_sweep / gputile.py — SLOP for the pointwise
    evaluations (the one remaining libm scope), CERT for all
    constants; box-exact CPU/GPU agreement (Result 28).]

L4 (local certificates). Every collection zone carries its own
    exclusion-or-structure proof:
    - Q-tubes: theta*(db) = th0 + S db + Q[db,db]/2 with TM-certified
      residual bound rad_g and Krawczyk containment: every root in
      the guard zone lies in the tube. [tm.py/tmres.py,
      certify_root_tube — CERT.]
    - Valley/fold windows: ball-coverage dichotomy floors (local
      frames, adaptive 4/1 vs 3/2 split) certify that roots in the
      collected window live only in the dip enclosures; certified
      tube radius where the growth quadratic solves. [fold.py —
      CERT for 19/22 validation valleys; near-cusp tube handoff
      PROTO (open R8b: the sampled shell drops sqrt5 / l2->linf
      handoff factors).]

L5 (clique bound). Pairwise overlap lower bounds: TM inner products
    for Q-tube pairs (CERT), interval inner products over certified
    dip boxes for valley rows (CERT where consistent, PROTO
    fallback otherwise), drift bounds for the rest. Pairs with
    lo > 0 are certified non-orthogonal; the conflict graph
    (lo <= 0) is greedily <= 5-colored. Same-color classes are
    pairwise non-orthogonal; a valley dip contributes at most ONE
    clique vertex (certified intra-dip self-overlap > 0). Hence any
    orthogonality clique has size <= #colors < 6. [color_conflicts,
    certified_dip_rows — CERT/PROTO as above.]

## 2. The compactness / root-count argument (why enumeration is
##    never trusted)

The center enumeration (find_mu_vectors, multistart) SUGGESTS the
root list; the certificate never relies on its completeness:

**Lemma (cover-completeness).** If the L3 sweep terminates with zero
stuck boxes, then every point of the 5-torus is in a discarded box
or a collection zone. Discarded boxes contain no MU vector of any
beta in T (L3 test + rates). Collection zones contain MU vectors
only inside their certified structures (tubes/dips) by L4. Hence the
COMPLETE MU-vector set of every H(beta), beta in T, is contained in
the union of tubes and dip enclosures indexed by the enumerated
roots — compactness of the torus is used only through termination of
the box subdivision (finite tree, wmin floor).

Corollary: a root missed by enumeration is impossible-or-visible —
either its neighborhood is excluded (so it does not exist), or the
sweep leaves stuck boxes / an uncollected region and the certificate
FAILS loudly. Root birth inside a tile (count changes between tile
centers) is covered the same way: phantom anchors (two-stage sweep)
attach windows to un-excludable blobs with no center root.

**Lemma (clique transfer).** For beta in T, any 6 mutually
orthogonal MU vectors map to 6 clique vertices in the certificate's
structures (each tube one vertex; each dip at most one vertex by
certified self-overlap). The <= 5-coloring of the conflict graph
contradicts a 6-clique. QED (tile).

## 3. Assembling tiles into the theorem

- Fundamental domain: F = [0, pi/2] x [0, pi/2] x [0, pi] suffices
  by the order-32 equivalence group (Result 26; Haagerup-invariant
  screened, multiset-confirmed — the group action is exact algebra).
  Closed tiles may overhang F (over-cover is sound).
- Coverage: the campaign ledger must witness a finite closed cover
  of F by certified tiles (chain tiles + adaptive standalone tiles).
  [campaign.py/dispatch.py — machinery validated; campaign not run.]
- Layer-3 families (triples containing specific Hadamard arcs) are
  separate, completed certificates: B-arc end-to-end incl. branch
  points (Results 20/31), Dita arc (21), Szollosi region gap-free
  (23/24), S6 and pointwise anchors (6..13). These are corollaries
  of the same sweep semantics at fixed Hadamards.

## 4. Assumption ledger (current)

1. Sweep pointwise floats: CLOSED (Result 32) — in-repo certified
   trig kernel (E_TRIG = 3e-13, proven by exact coefficient algebra)
   in every sweep; ~10 calls/test <= 3e-12 vs SLOP = 1e-11.
2. Near-cusp valley tube handoff: CLOSED (Result 33) — thin oracle
   (actual root territory), pure-phi ball floors standalone, shell
   boxes self-exclude or stick loudly (cover-completeness). 22/22
   validation valleys certified, zero fallbacks.
3. Layer-3 family rates (B/K arc walks): FD + PAD. NARROWED
   (bnrates.py): certified dB/dtheta is DONE and tight (dual-AD B-N
   map, no dependency blowup); what remains is K-column motion for
   ILL-CONDITIONED columns (||J^-1|| ~ 5e2) — the 1-parameter
   frame-split continuation (strong-block Krawczyk contracts at
   q ~ 0.02; weak direction needs signed interval dg with
   cancellation + a 1-dim crossing argument).
4. tau-coverage: CLOSED (ball-coverage lemma, Result 29-2).
5. libm in certified constants: CLOSED (Result 30, mpmath.iv).

Everything else in the certificate chain is certified arithmetic.

## 5. The two-regime architecture (2026-07-30/31 program)

The uniform-tile campaign of section 3 is superseded by a measured
two-regime decomposition (full diary: NOTES_LP_BRIDGE.md 4.5-4.26):

**Bulk (>= ~99% of the fundamental domain).** Coarse tiles at
h ~ 1.5e-3-3e-3 prove the Layer-1 claim (no MU triple through
{I, H(beta)} — strictly stronger than the theorem needs, and TRUE
off the triple-carrying strata). Certificate composition per tile:
coarse sweep with first-order taxes (Result 34) -> S-sloped
correlated tubes (pair coloring, typically 2-5 colors) -> batched
valley windows with one-level beta-recursion -> wild roots as
colored vertices with measured localization -> w^4 whisker-tail
refinement fallback -> blob-coverage against the root list
(enumeration is never trusted; missed roots surface as uncovered
blobs). Validated 6/6 census points; chained along beta with a
first-order exclusion cache (measured amortization ~3x, re-anchor
cadence ~5; window sub-caching and GPU sweeps are the remaining
mechanical multipliers).

**Collar (< ~1%: neighborhoods of the triple-carrying strata,
concentrated at the theta=0 branch face = the group-Hadamard
locus).** Here triples EXIST, the Layer-1 claim is false, and the
statement switches to the quadruple level: no fourth basis extends
the strata's triples. This is the Layer-3 family machinery
(complete, gap-free along the B-arc/Dita/Szollosi slices)
THICKENED transversally via the s^2 branch chart and the certified
K-drift substrate (bnrates). Status: ingredients built and
measured; the thickened certificates are the remaining new build.

**Seam.** The regimes self-delimit: bulk tiles fail loudly exactly
where triples appear, so the collar's extent is discovered by the
machine, not assumed.

Assumption ledger deltas vs section 4: the coarse composition
reintroduces PROTOTYPE-grade elements (batched-Newton anchors are
definitional and residual-checked; measured wild localizations and
the refinement fallback are sampled-grade pending their certified
pass). The fine-tile chain (sections 1-4) remains assumption-free.
