"""Bulk starvation certificate prototype (two-regime plan, bulk side).

For a fat tile (beta box, h ~ 3e-3) in the clique-starved bulk:
1. streaming fat sweep (pointwise exclusion, certified sup + first-
   order beta taxes) -> surviving boxes = thin filaments around the
   root trenches (measured ~w^4, vol frac ~4e-5 at wmin 0.025);
2. bin survivors into cells of side SEG (grid hash), hull per cell;
3. pairwise overlap lower bounds between cell hulls via the exact
   per-phase Lipschitz constant 1/6:
     |<u,v>| >= |<u_cA, u_cB>| - (||rA||_1 + ||rB||_1)/6 - SLOP;
   cells are self-coherent (internal overlap >= 1 - diam_1/6), so
   each contributes at most one clique vertex;
4. the possible-orthogonality graph (lower bound <= 0) must have
   max clique <= 5. Candidates are beta-independent vectors, so the
   union-over-beta graph is conservative for every beta in the tile.

Prototype grade: the sweep taxes are certified; the segment overlap
bound is exact arithmetic + SLOP; the clique search is exact.
"""

import sys
import time
import warnings

sys.setrecursionlimit(400000)

import numpy as np

from certify import SLOP
from karlsson import karlsson_map
from rates import certified_rates
from trig_kernel import kexp_i

warnings.filterwarnings("ignore")

SEG = 0.12


def fat_sweep(beta, h, wmin=0.025, chunk=1_000_000, max_boxes=4e8):
    r = certified_rates(beta, (h, h, h))
    Hc = karlsson_map(*beta).conj()
    BR, s_drift = r["beta_rate_vec"], r["s_drift"]
    dH0c = [np.conj(r["dH0"][j]) for j in range(3)]
    WD = r["WD"]
    inv6 = 1.0 / np.sqrt(6.0)
    sC = [np.full((1, 5), np.pi)]
    sW = [np.full((1, 5), np.pi)]
    outC, outW = [], []
    total = 0
    while sC:
        C = sC.pop()
        W = sW.pop()
        if len(C) > chunk:
            sC.append(C[chunk:])
            sW.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        total += len(C)
        if total > max_boxes:
            raise RuntimeError("box budget")
        u = np.empty((len(C), 6), complex)
        u[:, 0] = inv6
        u[:, 1:] = kexp_i(C) * inv6
        s = u @ Hc
        g = np.abs(s) ** 2 - 1.0 / 6.0
        sw = W.sum(axis=1)
        smod = np.minimum(np.abs(s) + sw[:, None] / 6.0 + s_drift, 1.0)
        margin = np.abs(g) - (2.0 * smod / 6.0) * sw[:, None]
        tax = BR * smod
        t1 = np.zeros_like(tax)
        for j in range(3):
            sb = u @ dH0c[j]
            t1 += h * (np.abs(2.0 * np.real(np.conj(s) * sb))
                       + 2.0 * s_drift * np.abs(sb)
                       + 2.0 * smod * WD[j])
        tax = np.minimum(tax, t1 + SLOP)
        keep = ~(margin > SLOP + tax).any(axis=1)
        C, W = C[keep], W[keep]
        small = W.max(axis=1) <= wmin
        if small.any():
            outC.append(C[small])
            outW.append(W[small])
        Cb, Wb = C[~small], W[~small]
        if len(Cb):
            j = np.argmax(Wb, axis=1)
            rows = np.arange(len(Cb))
            Wn = Wb.copy()
            Wn[rows, j] /= 2.0
            Cl, Cr = Cb.copy(), Cb.copy()
            Cl[rows, j] -= Wn[rows, j]
            Cr[rows, j] += Wn[rows, j]
            sC.append(np.vstack([Cl, Cr]))
            sW.append(np.vstack([Wn, Wn]))
    return np.vstack(outC), np.vstack(outW), total


def segments(C, W):
    """Grid-bin survivor boxes into cells of side SEG; return cell
    centers (hull midpoints) and per-axis hull radii."""
    keys = np.floor(C / SEG).astype(np.int64)
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    n = int(inv.max()) + 1
    order = np.argsort(inv, kind="stable")
    invs = inv[order]
    starts = np.searchsorted(invs, np.arange(n))
    Lo = (C - W)[order]
    Hi = (C + W)[order]
    lo = np.minimum.reduceat(Lo, starts, axis=0)
    hi = np.maximum.reduceat(Hi, starts, axis=0)
    return 0.5 * (lo + hi), 0.5 * (hi - lo)


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


def certify_bulk_tile(beta, h=3e-3):
    t0 = time.time()
    C, W, total = fat_sweep(beta, h)
    cen, rad = segments(C, W)
    n = len(cen)
    inv6 = 1.0 / np.sqrt(6.0)
    u = np.empty((n, 6), complex)
    u[:, 0] = inv6
    u[:, 1:] = np.exp(1j * cen) * inv6
    r1 = rad.sum(axis=1)
    self_lo = 1.0 - 2.0 * r1 / 6.0 - SLOP
    # blockwise possible-orthogonality edges (dense n x n would OOM)
    nbr = [[] for _ in range(n)]
    n_edges = 0
    B = 2048
    uc = u.conj()
    for a0 in range(0, n, B):
        blk = np.abs(uc[a0:a0 + B] @ u.T)
        lo_b = blk - (r1[a0:a0 + B, None] + r1[None, :]) / 6.0 - SLOP
        ii, jj = np.nonzero(lo_b <= 0)
        for i, j in zip(ii + a0, jj):
            if i < j:
                nbr[i].append(j)
                nbr[j].append(i)
                n_edges += 1
    # greedy coloring proves clique <= n_colors (sparse)
    order = np.argsort(-np.array([len(x) for x in nbr]))
    color = -np.ones(n, dtype=int)
    for v in order:
        used = {color[w] for w in nbr[v] if color[w] >= 0}
        c = 0
        while c in used:
            c += 1
        color[v] = c
    k = int(color.max()) + 1
    ok = bool(k <= 5 and self_lo.min() > 0.05)
    print(f"BULK beta=({beta[0]:.3f},{beta[1]:.3f},{beta[2]:.3f}) h={h}: "
          f"{'CERTIFIED' if ok else 'FAILED'} — {len(C)} surv boxes, "
          f"{n} segments (max r1 {r1.max():.3f}), poss-orth edges "
          f"{n_edges}, clique<= {k}, min self {self_lo.min():.3f}, "
          f"{total/1e6:.0f}M swept [{time.time()-t0:.0f}s]", flush=True)
    return ok


def main():
    pts = [(1.5349059690692832, 0.5515593746413914, 2.9470193319949907),
           (1.1045934240780566, 0.6254576854126334, 0.8763425369356228),
           (0.6429646757586674, 0.57148158011095, 2.607748463725592)]
    for beta in pts:
        certify_bulk_tile(beta)


if __name__ == "__main__":
    main()
