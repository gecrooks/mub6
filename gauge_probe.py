"""Gauge optimization for the graded LP (note 4.13): find diagonal
row phases (and discrete pair swaps) minimizing the off-graded mass
of the Karlsson unbiasedness kernels h_k h_k^H."""
import numpy as np
from scipy.optimize import minimize
from karlsson import karlsson_map

H = karlsson_map(5.978503016422594, 4.007534549834652, 1.6327649325136653)
Ms = [np.outer(H[:, k], H[:, k].conj()) for k in range(6)]


def graded_proj(M):
    P = np.zeros_like(M)
    for m in range(3):
        for mp in range(3):
            B = M[2*m:2*m+2, 2*mp:2*mp+2]
            d0 = 0.5 * (B[0, 0] + B[1, 1])
            d1 = 0.5 * (B[0, 1] + B[1, 0])
            P[2*m:2*m+2, 2*mp:2*mp+2] = [[d0, d1], [d1, d0]]
    return P


def offmass(phi, swaps):
    d = np.exp(1j * np.concatenate(([0.0], phi)))
    perm = np.arange(6)
    for m in range(3):
        if swaps & (1 << m):
            perm[2*m], perm[2*m+1] = perm[2*m+1], perm[2*m]
    tot = 0.0
    for M in Ms:
        Mg = (d[:, None] * M * d.conj()[None, :])[np.ix_(perm, perm)]
        tot += np.linalg.norm(Mg - graded_proj(Mg))**2
    return tot / sum(np.linalg.norm(M)**2 for M in Ms)


best = (1e9, None, None)
rng = np.random.default_rng(1)
for swaps in range(8):
    for _ in range(24):
        r = minimize(offmass, rng.uniform(0, 2*np.pi, 5), args=(swaps,),
                     method="Nelder-Mead",
                     options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        if r.fun < best[0]:
            best = (r.fun, swaps, r.x)
print(f"min off-graded mass fraction: {best[0]:.6e} (swaps={best[1]})")
print("phases/pi:", np.round(best[2] / np.pi, 4))
