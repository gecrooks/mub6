import numpy as np
from scipy.optimize import minimize
from karlsson import karlsson_map


def graded_proj(M):
    P = np.zeros_like(M)
    for m in range(3):
        for mp in range(3):
            B = M[2*m:2*m+2, 2*mp:2*mp+2]
            d0 = 0.5 * (B[0, 0] + B[1, 1])
            d1 = 0.5 * (B[0, 1] + B[1, 0])
            P[2*m:2*m+2, 2*mp:2*mp+2] = [[d0, d1], [d1, d0]]
    return P


def u2(p):
    a, b, c, d = p
    ph = np.exp(1j * a)
    return ph * np.array([[np.cos(b) * np.exp(1j * c),
                           np.sin(b) * np.exp(1j * d)],
                          [-np.sin(b) * np.exp(-1j * d),
                           np.cos(b) * np.exp(-1j * c)]])


def offmass_u2(p, Ms, den):
    V = np.zeros((6, 6), complex)
    for m in range(3):
        V[2*m:2*m+2, 2*m:2*m+2] = u2(p[4*m:4*m+4])
    tot = 0.0
    for M in Ms:
        Mg = V @ M @ V.conj().T
        tot += np.linalg.norm(Mg - graded_proj(Mg))**2
    return tot / den


def probe(beta, tag, n_starts=30):
    H = karlsson_map(*beta)
    Ms = [np.outer(H[:, k], H[:, k].conj()) for k in range(6)]
    den = sum(np.linalg.norm(M)**2 for M in Ms)
    rng = np.random.default_rng(2)
    best = 1e9
    for _ in range(n_starts):
        r = minimize(offmass_u2, rng.uniform(0, 2*np.pi, 12),
                     args=(Ms, den), method="Nelder-Mead",
                     options=dict(xatol=1e-10, fatol=1e-12,
                                  maxiter=8000))
        best = min(best, r.fun)
    print(f"{tag}: U(2)^3 min off-graded mass = {best:.6e}", flush=True)


probe((5.978503016422594, 4.007534549834652, 1.6327649325136653), "bulk ref")
probe((0.01, 2.041236, np.pi), "near F-face th=0.01")
probe((0.1, 2.041236, np.pi), "off-face th=0.1")
