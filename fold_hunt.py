"""Locate a genuine fold (root birth/death) along a path in the Karlsson
family by bisecting the MU-vector count, then examine the near-fold pair
structure for the R3 (fold-pair) certificate."""

import warnings

import numpy as np

from certify import _g_and_J, _torus_delta
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import polish_root

warnings.filterwarnings("ignore")

# scan pt 2 (48 roots) and pt 0 (60 roots), regenerated from the seed
rng = np.random.default_rng(20260726)
PTS = [tuple(rng.uniform(0.25, 2 * np.pi - 0.25, 3)) for _ in range(3)]
B48, B60 = np.array(PTS[2]), np.array(PTS[0])


def count_at(s, n_starts=6000, seed=0):
    beta = tuple((1 - s) * B48 + s * B60)
    H = karlsson_map(*beta)
    vecs = find_mu_vectors([H], n_starts=n_starts, seed=seed)
    return len(vecs), beta, vecs


def main():
    print("=== fold hunt along pt2(48) -> pt0(60) ===")
    lo, hi = 0.0, 1.0
    n_lo = count_at(lo)[0]
    n_hi = count_at(hi)[0]
    print(f"s=0: {n_lo} vecs   s=1: {n_hi} vecs")
    # first, a coarse scan to find one transition bracket
    grid = np.linspace(0, 1, 9)
    counts = []
    for s in grid:
        n, beta, _ = count_at(s)
        counts.append(n)
        print(f"  s={s:.3f}: {n} vecs", flush=True)
    k = next(i for i in range(8) if counts[i] != counts[i + 1])
    lo, hi = grid[k], grid[k + 1]
    n_lo, n_hi = counts[k], counts[k + 1]
    print(f"bracket [{lo:.3f}, {hi:.3f}]: {n_lo} -> {n_hi}")

    for _ in range(9):
        mid = 0.5 * (lo + hi)
        n_mid = count_at(mid, n_starts=9000)[0]
        print(f"  s={mid:.5f}: {n_mid}", flush=True)
        if n_mid == n_lo:
            lo = mid
        else:
            hi = mid
    print(f"fold near s = {0.5*(lo+hi):.5f}")

    # pair structure on the MANY side, just past the fold
    side = hi if n_hi > n_lo else lo
    n, beta, vecs = count_at(side, n_starts=12000)
    H = karlsson_map(*beta)
    roots = [polish_root(H, np.angle(v * np.sqrt(6))[1:]) for v in vecs]
    R = np.array(roots)
    m = len(R)
    print(f"many side (s={side:.5f}): {n} vecs; nearest pairs and sigmas:")
    D = np.full((m, m), np.inf)
    for a in range(m):
        for b in range(a + 1, m):
            D[a, b] = np.max(np.abs(_torus_delta(R[a], R[b])))
    pairs = np.dstack(np.unravel_index(np.argsort(D, axis=None), D.shape))[0]
    shown = 0
    for a, b in pairs:
        if D[a, b] > 0.6 or shown >= 4:
            break
        sa = np.linalg.svd(_g_and_J(H, R[a])[1], compute_uv=False)[-1]
        sb = np.linalg.svd(_g_and_J(H, R[b])[1], compute_uv=False)[-1]
        # overlap of the pair's vectors
        ua = np.exp(1j * np.concatenate(([0.0], R[a]))) / np.sqrt(6)
        ub = np.exp(1j * np.concatenate(([0.0], R[b]))) / np.sqrt(6)
        ov = abs(np.vdot(ua, ub))
        print(f"  pair ({a},{b}): sep {D[a,b]:.4f}  sigmas {sa:.4f}/{sb:.4f}"
              f"  |<u,v>| = {ov:.4f}")
        shown += 1
    np.save("fold_point.npy", dict(beta=beta, s=side), allow_pickle=True)


if __name__ == "__main__":
    main()
