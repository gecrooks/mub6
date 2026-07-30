"""beta-correlated blob certificate, stage 1 (note 4.12): branch
decomposition + hull-covering parametric Krawczyk (class R) + same-
beta pair bounds + coloring, at bulk points, h = 3e-3.

R-branch semantics: the Krawczyk box COVERS the branch hull, so the
verified statement is "exactly one root in this branch for every
beta in the tile", and its enclosure box bounds that root's position
uniformly in beta. Pair bounds then compare enclosure boxes — valid
at every fixed beta simultaneously. Branches failing (deep sigma,
fat hulls, no contraction) are WILD; the tile's clique bound is
colors(R-graph) + n_wild (crude Wild handling = stage 2)."""
import time
import warnings

import numpy as np

from certify import SLOP, cluster_suspects, krawczyk_verify
from karlsson import karlsson_map
from parametric import polish_root, _g_and_J
from rates import certified_rates
from starve import fat_sweep

warnings.filterwarnings("ignore")


def branches(C, W, link=0.06):
    return cluster_suspects(C, W, link=link)


def certify_point(beta, h=3e-3):
    t0 = time.time()
    C, W, total = fat_sweep(beta, h)
    cl = branches(C, W)
    H0 = karlsson_map(*beta)
    r = certified_rates(beta, (h, h, h))
    # certified max entry drift of H over the beta box
    hs = float(sum(np.abs(r["dH0"][j]).max() * h for j in range(3))
               + r["WD"].max() * np.sqrt(6.0) * h) + 1e-9
    R_list, wild = [], 0
    for cc, rr, _ in cl:
        th = polish_root(H0, np.asarray(cc))
        _, J = _g_and_J(H0, th)
        sig = np.linalg.svd(J, compute_uv=False)[-1]
        hull = np.maximum(np.asarray(rr) + 0.02, 0.06)
        if sig < 0.08:
            wild += 1
            continue
        ok, c, r_enc, _ = krawczyk_verify(H0, th, hull, hslop=hs)
        if not ok:
            wild += 1
            continue
        R_list.append((c, r_enc))
    n = len(R_list)
    cols = 0
    if n:
        cen = np.array([c for c, _ in R_list])
        enc = np.array([e for _, e in R_list])
        inv6 = 1.0 / np.sqrt(6.0)
        u = np.empty((n, 6), complex)
        u[:, 0] = inv6
        u[:, 1:] = np.exp(1j * cen) * inv6
        O = np.abs(u.conj() @ u.T)
        r1 = enc.sum(axis=1)
        lo = O - (r1[:, None] + r1[None, :]) / 6.0 - SLOP
        adj = (lo <= 0) & ~np.eye(n, dtype=bool)
        order = np.argsort(-adj.sum(axis=1))
        color = -np.ones(n, dtype=int)
        for v in order:
            used = set(color[adj[v]]) - {-1}
            cix = 0
            while cix in used:
                cix += 1
            color[v] = cix
        cols = int(color.max()) + 1
    bound = cols + wild
    print(f"BRANCH beta=({beta[0]:.3f},{beta[1]:.3f},{beta[2]:.3f}): "
          f"{'CERTIFIED' if bound <= 5 else 'FAILED'} — "
          f"{len(cl)} branches, R={n}, wild={wild}, colors_R={cols}, "
          f"clique bound {bound}, hs={hs:.1e}, "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    for beta in [(1.5349059690692832, 0.5515593746413914, 2.9470193319949907),
                 (1.1045934240780566, 0.6254576854126334, 0.8763425369356228),
                 (0.6429646757586674, 0.57148158011095, 2.607748463725592)]:
        certify_point(beta)
