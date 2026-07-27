"""Probe the Karlsson family K6^(3) — the family the Fourier-analytic LP
cannot exclude from complete sets (McNulty-Weigert Sec. 10.2).

Sampling trick: K6^(3) = H2-reducible Hadamards (Thm 7.2: all nine 2x2
blocks are themselves Hadamard), so instead of Karlsson's closed form we add
9 block-orthogonality residuals to the random-Hadamard least squares.

Findings (2026-07-24): 8 samples, 48-60 MU vectors each, zero MU base-pairs.
Only samples that are secretly Fourier-family points (fingerprint: three
dephased rows contain -1, 24/36 entries at exact 6th roots) admit bases at
all (48 vectors / 8 bases). Generic K3 points: no bases. Consistent with
Conjecture 8.1 and Zauner.
"""

import warnings

import numpy as np
from scipy.optimize import least_squares

from mub import dephase
from families import probe

warnings.filterwarnings("ignore")


F2 = np.array([[1, 1], [1, -1]], dtype=complex)


def karlsson_map(theta, phi, lam):
    """Exact K6^(3) member, Karlsson arXiv:1003.4177 (Prop 4, Cor 6, Eq 4.2).

    Lambda(theta,phi) self-adjoint unitary; A = F2(-e/2 + i sqrt(3)/2 Lambda),
    B = -F2 - A; z1 = e^{i lam}; z3^2, z4^2, z2^2 via the unit-circle-
    preserving Moebius maps of z1^2 given by the unimodularity relations.
    Verified: 200/200 random parameters give exact Hadamards to 1e-15.
    """
    L = np.array([[np.cos(theta), np.exp(1j * phi) * np.sin(theta)],
                  [np.exp(-1j * phi) * np.sin(theta), -np.cos(theta)]],
                 dtype=complex)
    A = F2 @ (-0.5 * np.eye(2) + 1j * (np.sqrt(3) / 2) * L)
    B = -F2 - A
    A11, A12 = A[0, 0], A[0, 1]
    B11, B12 = B[0, 0], B[0, 1]
    z1 = np.exp(1j * lam)
    w = z1 ** 2
    z3sq = (A11 ** 2 - w * A12 ** 2) / (np.conj(A12) ** 2 - w * np.conj(A11) ** 2)
    z4sq = (B11 ** 2 - w * B12 ** 2) / (np.conj(B12) ** 2 - w * np.conj(B11) ** 2)
    z3, z4 = np.sqrt(z3sq), np.sqrt(z4sq)
    z2 = np.sqrt((B11 ** 2 - z3sq * np.conj(B12) ** 2)
                 / (B12 ** 2 - z3sq * np.conj(B11) ** 2))
    Z1 = np.array([[1, 1], [z1, -z1]], dtype=complex)
    Z2 = np.array([[1, 1], [z2, -z2]], dtype=complex)
    Z3 = np.array([[1, z3], [1, -z3]], dtype=complex)
    Z4 = np.array([[1, z4], [1, -z4]], dtype=complex)
    H = np.zeros((6, 6), dtype=complex)
    H[0:2, 0:2] = F2
    H[0:2, 2:4] = Z1
    H[0:2, 4:6] = Z2
    H[2:4, 0:2] = Z3
    H[4:6, 0:2] = Z4
    H[2:4, 2:4] = 0.5 * Z3 @ A @ Z1
    H[2:4, 4:6] = 0.5 * Z3 @ B @ Z2
    H[4:6, 2:4] = 0.5 * Z4 @ B @ Z1
    H[4:6, 4:6] = 0.5 * Z4 @ A @ Z2
    return H / np.sqrt(6)


def h2_residuals(phases, d=6):
    E = np.exp(1j * phases.reshape(d - 1, d - 1))
    H = np.ones((d, d), dtype=complex)
    H[1:, 1:] = E
    G = H.conj().T @ H / d - np.eye(d)
    iu = np.triu_indices(d, 1)
    res = [G[iu].real, G[iu].imag]
    for bj in range(3):
        for bk in range(3):
            B = H[2 * bj:2 * bj + 2, 2 * bk:2 * bk + 2]
            z = B[0, 0] * np.conj(B[1, 0]) + B[0, 1] * np.conj(B[1, 1])
            res.append(np.array([z.real, z.imag]))
    return np.concatenate(res)


def random_k3(seed):
    local = np.random.default_rng(seed)
    for _ in range(200):
        p0 = local.uniform(0, 2 * np.pi, size=25)
        r = least_squares(h2_residuals, p0, method="lm",
                          xtol=3e-16, ftol=3e-16, gtol=3e-16)
        if r.cost < 1e-22:
            E = np.exp(1j * r.x.reshape(5, 5))
            H = np.ones((6, 6), dtype=complex)
            H[1:, 1:] = E
            return H / np.sqrt(6)
    return None


def fourier_fingerprint(H):
    """(rows with -1, cols with -1, #entries at exact 6th roots) of dephased H."""
    D = dephase(H) * np.sqrt(6)
    ang = np.angle(D) / (2 * np.pi) * 6
    n6 = int(np.sum(np.abs(ang - np.round(ang)) < 1e-6))
    rows = sum(1 for r in range(6) if np.min(np.abs(D[r, :] + 1)) < 1e-7)
    cols = sum(1 for c in range(6) if np.min(np.abs(D[:, c] + 1)) < 1e-7)
    return rows, cols, n6


if __name__ == "__main__":
    found = 0
    for s in range(40):
        H = random_k3(3000 + s)
        if H is None:
            continue
        found += 1
        r, c, n6 = fourier_fingerprint(H)
        print(f"  fingerprint: rows-with--1={r} cols-with--1={c} sixth-roots={n6}/36")
        probe(f"K3 sample #{found} (seed {3000+s})", H, n_starts=4000, seed=s)
        if found >= 8:
            break
