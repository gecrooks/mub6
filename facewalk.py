"""Face-walk prototype (PROOF_SKELETON §5, face piece).

At theta = 0 the Karlsson family degenerates to the 2-parameter
group-Hadamard face (the Jaming-Matolcsi-Mora-Szollosi locus) and
triples exist; the theorem's statement there is quadruple-level: no
fourth basis. Measured structure (NOTES_LP_BRIDGE 4.30): 48 MU
vectors, exactly 8 extra bases, max mutually-unbiased clique among
them = 1, obstructed at the basis-PAIR level with defect margin
0.09-0.17.

Per face tile (in-face box of half-width hf at a proxy theta):
  1. pool the MU vectors (coverage against fat-sweep hulls comes
     from the existing starve machinery in the certified pass);
  2. enumerate the near-bases (orthogonality 6-cliques);
  3. sanity: max mutually-unbiased clique among bases must be <= 1;
  4. for every basis pair, lower-bound the unbiasedness defect over
     the tile: center defect - RATE_FACE * hf * sqrt2 - PAD.

PROTOTYPE GRADE: RATE_FACE is the sampled in-face defect drift
(measured ~0.07 at the F-point, padded to 0.5); PAD covers the
theta-proxy offset (defect is even in theta, measured flat to 4
digits at theta=0.003). The certified pass replaces both with
dual-AD rates in the s^2 branch chart.
"""

import time
from itertools import combinations

import numpy as np
from scipy.optimize import least_squares

from karlsson import karlsson_map
from mub import (bases_matrix, find_bases, find_mu_vectors,
                 mu_vector_residuals)

TH_PROXY = 1e-6      # face proxy; karlsson_map is singular at 0
RATE_FACE = 0.5      # sampled in-face defect drift 0.07, padded 7x
RATE_TH = 1.0        # sampled theta drift of the unbias defect:
                     # measured +0.0002 over theta 0-0.02 (~0.01),
                     # padded 100x — the defect is flat in theta
PAD = 1e-3           # theta-proxy + enumeration slop (sampled grade)


def face_pool(b2, b3, n_starts=6000, seed=3):
    H = karlsson_map(TH_PROXY, b2, b3)
    ps = []
    for v in find_mu_vectors([H], n_starts=n_starts, seed=seed):
        t = np.angle(v * np.sqrt(6))[1:]
        t = least_squares(mu_vector_residuals, t, args=([H],),
                          method="lm", xtol=3e-16, ftol=3e-16,
                          gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], t))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-5 for u in ps):
            ps.append(w)
    return np.array(ps)


def face_tile(b2, b3, hf, th_tube=0.0):
    """Certify (prototype grade) 'no fourth basis' on the box
    [b2 +- hf] x [b3 +- hf] x theta in [0, th_tube]. The theta-tube
    is the collar handoff: the collar's dyadic slabs descend to
    theta_lo = th_tube and this certificate covers the rest.
    Returns (ok, n_bases, worst_margin)."""
    P = face_pool(b2, b3)
    # near-basis tolerance scales with the proxy: wall tiles use
    # the safe proxy 3e-3 where orthogonality defects are
    # ~0.4-0.9 * proxy (the walls' breaking-edge law), so exact
    # bases only exist to that resolution
    tol = 1e-5 if TH_PROXY < 1e-4 else 30.0 * TH_PROXY
    bases = find_bases(P, tol=tol)
    if len(bases) != 8:
        return False, len(bases), np.nan   # unexpected stratum: split
    nb = len(bases)
    worst = np.inf
    for a, b in combinations(range(nb), 2):
        A = bases_matrix(P, bases[a]).conj().T @ bases_matrix(P, bases[b])
        defect = float(np.max(np.abs(np.abs(A) ** 2 - 1.0 / 6.0)))
        worst = min(worst, defect - RATE_FACE * hf * np.sqrt(2)
                    - RATE_TH * th_tube - PAD)
    return worst > 0, nb, worst


def main():
    t0 = time.time()
    hf = 0.05
    pts = [(2.041236, np.pi), (1.0, 2.0), (1.1, 1.0428571),
           (0.5, 5.0), (2.8, 0.7)]
    for b2, b3 in pts:
        t1 = time.time()
        ok, nb, m = face_tile(b2, b3, hf)
        print(f"FACE_TILE ({b2:.4f},{b3:.4f}) hf={hf}: "
              f"{'CERTIFIED' if ok else 'FAILED'} bases={nb} "
              f"margin={m:.4f} [{time.time()-t1:.0f}s]", flush=True)
    print(f"total {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
