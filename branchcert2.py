"""Fat-Q-tube branch certificate (corrected design, note 4.15)."""
import time
import warnings

import numpy as np

from certify import SLOP, _g_and_J
from karlsson import karlsson_map
from parametric import polish_root, root_data2
from starve import fat_sweep
from tmres import certified_curve_residual, tm_karlsson

warnings.filterwarnings("ignore")


def union_find_branches(C, pitch=0.05):
    """Components of survivor boxes via face-adjacency on a snap grid."""
    keys = np.round(C / pitch).astype(np.int64)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    n = len(uniq)
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    lut = {tuple(k): i for i, k in enumerate(uniq)}
    for ax in range(5):
        for d in (-1, 1):
            nb = uniq.copy()
            nb[:, ax] += d
            for i, k in enumerate(map(tuple, nb)):
                j = lut.get(k)
                if j is not None:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        parent[ri] = rj
    roots = np.array([find(i) for i in range(n)])
    _, comp = np.unique(roots, return_inverse=True)
    return comp[inv], int(comp.max()) + 1


def certify_point(beta, h=3e-3, smax=40.0, sig_min=0.05):
    t0 = time.time()
    hv = np.array([h, h, h])
    C, W, _tot = fat_sweep(beta, h)
    comp, ncomp = union_find_branches(C)
    H0 = karlsson_map(*beta)
    Htm = tm_karlsson(beta, hv)
    tubes, wild = [], 0
    for b in range(ncomp):
        m = comp == b
        cen = C[m].mean(axis=0)
        th = polish_root(H0, cen)
        _, J = _g_and_J(H0, th)
        sig = np.linalg.svd(J, compute_uv=False)[-1]
        if sig < sig_min:
            wild += 1
            continue
        try:
            S, Q, defect = root_data2(beta, th)
        except Exception:
            wild += 1
            continue
        if np.abs(S).max() > smax or defect > 1e-6:
            wild += 1
            continue
        try:
            resid = certified_curve_residual(beta, hv, th, S, Q, Htm=Htm)
        except Exception:
            wild += 1
            continue
        loc = np.abs(S) @ hv + resid / max(sig, 1e-3) + 1e-3
        tubes.append((th, S, float(loc.sum())))
    n = len(tubes)
    cols = 0
    if n:
        cen = np.array([t[0] for t in tubes])
        Ss = np.array([t[1] for t in tubes])
        locs = np.array([t[2] for t in tubes])
        inv6 = 1.0 / np.sqrt(6.0)
        u = np.empty((n, 6), complex)
        u[:, 0] = inv6
        u[:, 1:] = np.exp(1j * cen) * inv6
        O = np.abs(u.conj() @ u.T)
        # pair drift: sum_la |S_i - S_j| h / 6 + residual locs
        Sd = np.abs(Ss[:, None, :, :] - Ss[None, :, :, :]).sum(axis=(2, 3))
        lo = O - Sd * h / 6.0 - (locs[:, None] + locs[None, :]) / 60.0 - SLOP
        adj = (lo <= 0) & ~np.eye(n, dtype=bool)
        order = np.argsort(-adj.sum(axis=1))
        color = -np.ones(n, dtype=int)
        for v in order:
            used = set(color[adj[v]]) - {-1}
            cx = 0
            while cx in used:
                cx += 1
            color[v] = cx
        cols = int(color.max()) + 1
    bound = cols + wild
    print(f"TUBE beta=({beta[0]:.3f},{beta[1]:.3f},{beta[2]:.3f}): "
          f"{'CERTIFIED' if bound <= 5 else 'FAILED'} — branches {ncomp}, "
          f"tubes {n}, wild {wild}, colors {cols}, bound {bound} "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    certify_point((1.5349059690692832, 0.5515593746413914,
                   2.9470193319949907))
