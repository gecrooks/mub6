"""Multi-locus rigidity census (LP-bridge program, note section 4.7).

Maps the clique/orthogonality structure of the fourth-basis candidate
pool across the fundamental domain: (a) a coarse grid over the
theta ~ 0 branch face (where the Fourier locus lives) to find every
rigidity hot spot; (b) rays off the face at the known F-point in
several directions. Multistart enumeration grade — the certified
census comes later if the picture holds.
"""

import sys
import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from karlsson import karlsson_map
from mub import find_mu_vectors, mu_vector_residuals

warnings.filterwarnings("ignore")


def pool_of(H, n=4000, seed=7):
    pool = find_mu_vectors([H], n_starts=n, seed=seed)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([H],),
                           method="lm", xtol=3e-16, ftol=3e-16,
                           gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-6 for u in pol):
            pol.append(w)
    return np.array(pol)


def cliq(adj):
    best = [0]
    order = np.argsort(-adj.sum(axis=1))

    def bb(cand, cur):
        nonlocal best
        if len(cur) + len(cand) <= len(best):
            return
        if not cand:
            if len(cur) > len(best):
                best = list(cur)
            return
        v = cand[0]
        bb([u for u in cand[1:] if adj[v, u]], cur + [v])
        bb(cand[1:], cur)

    bb(list(order), [])
    return len(best)


def census(beta):
    t0 = time.time()
    P = pool_of(karlsson_map(*beta))
    G = np.abs(P.conj() @ P.T) ** 2
    iu = np.triu_indices(len(P), 1)
    v = G[iu]
    adj = (G < 1e-6) & ~np.eye(len(P), dtype=bool)
    print(f"CENSUS {beta[0]:.4f} {beta[1]:.4f} {beta[2]:.4f} "
          f"n={len(P)} orth={int((v < 1e-6).sum())} "
          f"clique={cliq(adj)} near16={int(((v > 0.9/6) & (v < 1.1/6)).sum())}"
          f" [{time.time()-t0:.0f}s]", flush=True)


def main():
    # (a) the branch face: theta = 0.002 fixed, coarse (phi, lam) grid
    for phi in np.linspace(0.15, 1.45, 6):
        for lam in np.linspace(0.2, 3.0, 8):
            census((0.002, float(phi), float(lam)))
    # (b) rays off the known F-point in +phi, -phi, +lam, diag
    bF = np.array([0.0, 2.041236, np.pi])
    for d in ([0, 1, 0], [0, -1, 0], [0, 0, -1], [1, 1, 1] / np.sqrt(3)):
        d = np.asarray(d, float)
        for t in (0.003, 0.03, 0.3):
            b = bF + t * d
            b[0] = max(b[0], 0.002)
            census(tuple(b))


if __name__ == "__main__":
    main()
