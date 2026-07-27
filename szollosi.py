"""Szollosi's two-parameter family X6^(2) (review eq 7.14): the last
Layer-3 territory. (x, y) are roots of f_alpha(z) = z^3 - alpha z^2 +
conj(alpha) z - 1 and (u, v) roots of f_{-alpha}, alpha in the deltoid
intersection D (both cubics all-unimodular; the cubics are
self-inversive, so root products are 1 and pair-sums are conjugate
sums automatically).

Transcription of the dephased form (validated numerically in main):
    [1,  1,    1,     1,       1,     1   ]
    [1,  x²y,  xy²,   xy/uv,   uxy,   vxy ]
    [1,  x/y,  x²y,   x/u,     x/v,   uvx ]
    [1,  uvx,  uxy,   -1,      -uxy,  -uvx]
    [1,  x/u,  vxy,   -x/u,    -1,    -vxy]
    [1,  x/v,  xy/uv, -xy/uv,  -x/v,  -1  ]
Rows 3-5 are identically orthogonal to row 0; the (x,y,u,v) constraints
from the cubic supply the rest.
"""

import itertools
import warnings

import numpy as np

warnings.filterwarnings("ignore")


def _roots(alpha, sign=+1):
    a = sign * alpha
    r = np.roots([1.0, -a, np.conj(a), -1.0])
    return r


def in_D(alpha, tol=1e-9):
    return (np.max(np.abs(np.abs(_roots(alpha, +1)) - 1)) < tol and
            np.max(np.abs(np.abs(_roots(alpha, -1)) - 1)) < tol)


def x_matrix(x, y, u, v):
    r0 = [1, 1, 1, 1, 1, 1]
    r1 = [1, x * x * y, x * y * y, x * y / (u * v), u * x * y, v * x * y]
    r2 = [1, x / y, x * x * y, x / u, x / v, u * v * x]
    r3 = [1, u * v * x, u * x * y, -1, -u * x * y, -u * v * x]
    r4 = [1, x / u, v * x * y, -x / u, -1, -v * x * y]
    r5 = [1, x / v, x * y / (u * v), -x * y / (u * v), -x / v, -1]
    return np.array([r0, r1, r2, r3, r4, r5], dtype=complex) / np.sqrt(6)


def szollosi_map(alpha, choice=None):
    """Build X6(alpha). If choice is None, search root-pair choices for
    one that yields a Hadamard; else use the given ((i,j),(k,l))."""
    rx = _roots(alpha, +1)
    ru = _roots(alpha, -1)
    tries = ([choice] if choice is not None else
             [(p, q) for p in itertools.permutations(range(3), 2)
              for q in itertools.permutations(range(3), 2)])
    best = None
    for (p, q) in tries:
        H = x_matrix(rx[p[0]], rx[p[1]], ru[q[0]], ru[q[1]])
        d = np.max(np.abs(H.conj().T @ H - np.eye(6)))
        if best is None or d < best[0]:
            best = (d, (p, q), H)
        if d < 1e-9:
            return H, (p, q), d
    return best[2], best[1], best[0]


def main():
    rng = np.random.default_rng(3)
    print("=== Szollosi chart validation ===")
    # sample D
    samples = []
    while len(samples) < 40:
        a = rng.uniform(-1.5, 1.5) + 1j * rng.uniform(-1.5, 1.5)
        if in_D(a):
            samples.append(a)
    n_ok = 0
    worst = 0.0
    choice_counts = {}
    for a in samples:
        H, ch, d = szollosi_map(a)
        if d < 1e-9:
            n_ok += 1
            choice_counts[ch] = choice_counts.get(ch, 0) + 1
        worst = max(worst, d)
    print(f"Hadamard at {n_ok}/40 sampled alpha in D "
          f"(worst unitarity defect {worst:.2e})")
    top = sorted(choice_counts.items(), key=lambda kv: -kv[1])[:3]
    print(f"root-pair choices that worked: {top}")


if __name__ == "__main__":
    main()
