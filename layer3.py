"""Layer 3 prototype: certified unextendability of an X-family MU triple.

Statement certified: for the concrete triple {I, B, K} (B in the
Hermitian family B(theta) subset X6^(2), K its unique MU partner basis),
NO unit vector unbiased to the identity basis is simultaneously unbiased
to the columns of B' and K', for ANY B', K' in explicit balls around the
float matrices. Since exact X-triples exist within the balls (Zauner's
circulant construction; float residuals ~1e-13), this certifies strong
unextendability of the exact triple.

Method: chunked LIFO branch-and-bound over the phase 5-torus, exclusion
per component over BOTH Hadamards' constraint stacks (10 components),
with the local-|s| Lipschitz bound and ball slop. Expect full exclusion
(zero suspects) -- the margin is the min over the torus of the largest
excludable component.
"""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from certify import SLOP, L_H_G, _uvec
from mub import (beauchamp_nicoara, bases_matrix, find_bases,
                 find_mu_vectors, mu_vector_residuals, unbiasedness_defect)

warnings.filterwarnings("ignore")


def build_triple(theta=1.6, seed=11, n_starts=20000):
    """Construct {I, B(theta), K}: enumerate MU vectors of {I, B}, polish,
    find the (unique) orthonormal basis K among them."""
    B = beauchamp_nicoara(theta)
    pool = find_mu_vectors([B], n_starts=n_starts, seed=seed)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([B],), method="lm",
                          xtol=3e-16, ftol=3e-16, gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-4 for u in pol):
            pol.append(w)
    pol = np.array(pol)
    bases = find_bases(pol, tol=1e-6)
    if len(bases) != 1:
        raise RuntimeError(f"expected unique basis, found {len(bases)}")
    K = bases_matrix(pol, bases[0])
    return B, K, pol, bases[0]


def margin_stats(B, K, pol, basis_idx):
    """Unextendability margin: for each vector MU to {I,B}, its worst
    unbiasedness defect to K (excluding K's own members)."""
    defs = []
    for i, v in enumerate(pol):
        if i in basis_idx:
            continue
        P = np.abs(K.conj().T @ v) ** 2
        defs.append(float(np.max(np.abs(P - 1 / 6))))
    return np.array(defs)


def certified_triple_sweep(B, K, hslop=1e-9, wmin=2e-3, chunk=100_000,
                           max_boxes=4e8, verbose=True):
    """Exclude the whole phase torus for the 10-component system.

    Returns (n_suspects, min_excl_margin, boxes). Zero suspects certifies
    the statement; min_excl_margin is the weakest certified inequality.
    """
    stacks = [B.conj(), K.conj()]
    stack_C = [np.full((1, 5), np.pi)]
    stack_W = [np.full((1, 5), np.pi)]
    n_sus = 0
    total = 0
    min_margin = np.inf
    t0 = time.time()
    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        total += len(C)
        if total > max_boxes:
            raise RuntimeError(f"exceeded {max_boxes:g} boxes")
        u = _uvec(C)
        sw = W.sum(axis=1)
        best = np.full(len(C), -np.inf)
        for Hc in stacks:
            s = u @ Hc
            g = np.abs(s) ** 2 - 1.0 / 6.0
            L = 2.0 * np.minimum(np.abs(s) + sw[:, None] / 6.0, 1.0) / 6.0
            marg = np.abs(g) - L * sw[:, None] - (SLOP + L_H_G * hslop)
            best = np.maximum(best, marg.max(axis=1))
        excl = best > 0
        min_margin = min(min_margin,
                         float(best[excl].min()) if excl.any() else np.inf)
        keep = ~excl
        Ck, Wk = C[keep], W[keep]
        small = Wk.max(axis=1) <= wmin
        n_sus += int(small.sum())
        Cb, Wb = Ck[~small], Wk[~small]
        if len(Cb):
            j = np.argmax(Wb, axis=1)
            rows = np.arange(len(Cb))
            Wn = Wb.copy()
            Wn[rows, j] /= 2.0
            Cl, Cr = Cb.copy(), Cb.copy()
            Cl[rows, j] -= Wn[rows, j]
            Cr[rows, j] += Wn[rows, j]
            stack_C.append(np.vstack([Cl, Cr]))
            stack_W.append(np.vstack([Wn, Wn]))
        if verbose and total % 20_000_000 < chunk <= total:
            print(f"    ...{total} boxes, {n_sus} suspects "
                  f"[{time.time()-t0:.0f} s]", flush=True)
    return n_sus, min_margin, total


def main():
    theta = 1.6
    print(f"=== Layer 3 prototype: triple {{I, B({theta}), K}} ===")
    B, K, pol, bidx = build_triple(theta)
    print(f"MU vectors of {{I,B}}: {len(pol)}; unique basis found")
    print(f"triple defect |B vs K|: {unbiasedness_defect(B, K):.2e}")
    # polish K to a nearby exact unitary for tighter balls
    U_, _, Vt_ = np.linalg.svd(K)
    K = U_ @ Vt_
    print(f"after unitary polish:   {unbiasedness_defect(B, K):.2e}")

    defs = margin_stats(B, K, pol, bidx)
    print(f"unextendability margin over {len(defs)} non-basis MU vectors: "
          f"min {defs.min():.4f}  median {np.median(defs):.4f}")

    print("certified sweep over the phase torus (10 components)...")
    t0 = time.time()
    n_sus, min_margin, total = certified_triple_sweep(B, K)
    dt = time.time() - t0
    print(f"  boxes {total}, suspects {n_sus}, min margin {min_margin:.5f} "
          f"[{dt:.0f} s]")
    if n_sus == 0:
        print("  ==> THEOREM (prototype-grade): no vector is MU to the "
              "triple {I, B', K'} for any B', K' in the 1e-9 balls -- "
              "the X-family triple is strongly unextendible.")


if __name__ == "__main__":
    main()
