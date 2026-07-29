"""Feasibility probe for the FAT-TILE starvation certificate
(NOTES_LP_BRIDGE.md two-regime plan, bulk regime).

At a bulk beta, sweep the 5-torus over a fat beta-box (h ~ 0.01-0.03)
using only pointwise exclusion with certified beta-taxes (sup-tax
min'd with the first-order signed-gradient tax of Result 34) — no
tubes, no valleys, no guards. The un-excluded remainder is where MU
vectors of ANY beta in the box can live. Measures: survivor volume,
blob count/diameters (grid-hash clustering), and the minimum
pairwise overlap between blob representatives — the quantities that
decide whether blob-level non-orthogonality (clique starvation over
fat tiles) is certifiable.
"""

import sys
import time
import warnings

import numpy as np

from certify import SLOP, cluster_suspects
from karlsson import karlsson_map
from rates import certified_rates
from trig_kernel import kexp_i

warnings.filterwarnings("ignore")


def fat_sweep(beta, h, wmin=0.1, chunk=200_000, max_boxes=1.2e8):
    r = certified_rates(beta, (h, h, h))
    H0 = karlsson_map(*beta)
    Hc = H0.conj()
    BR = r["beta_rate_vec"]
    s_drift = r["s_drift"]
    dH0c = [np.conj(r["dH0"][j]) for j in range(3)]
    WD = r["WD"]
    inv6 = 1.0 / np.sqrt(6.0)
    stack_C = [np.full((1, 5), np.pi)]
    stack_W = [np.full((1, 5), np.pi)]
    surv_C, surv_W = [], []
    total = 0
    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
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
        L = 2.0 * smod / 6.0
        margin = np.abs(g) - L * sw[:, None]
        tax = BR * smod
        t1 = np.zeros_like(tax)
        for j in range(3):
            sb = u @ dH0c[j]
            d0g = 2.0 * np.real(np.conj(s) * sb)
            t1 += h * (np.abs(d0g) + 2.0 * s_drift * np.abs(sb)
                       + 2.0 * smod * WD[j])
        tax = np.minimum(tax, t1 + SLOP)
        excl = (margin > SLOP + tax).any(axis=1)
        C, W = C[~excl], W[~excl]
        small = W.max(axis=1) <= wmin
        if small.any():
            surv_C.append(C[small])
            surv_W.append(W[small])
        Cb, Wb = C[~small], W[~small]
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
    return (np.vstack(surv_C) if surv_C else np.zeros((0, 5)),
            np.vstack(surv_W) if surv_W else np.zeros((0, 5)), total)


def probe(beta, h):
    t0 = time.time()
    C, W, total = fat_sweep(beta, h)
    vol = float((2 * W).prod(axis=1).sum()) if len(C) else 0.0
    frac = vol / (2 * np.pi) ** 5
    if len(C) == 0:
        print(f"STARVE beta={beta} h={h}: EMPTY survivor set "
              f"({total} boxes) — vacuously starved", flush=True)
        return
    clusters = cluster_suspects(C, W, link=0.12)
    reps = np.array([cc for cc, _r, _i in clusters])
    diam = [float(2 * _r.max()) for _c, _r, _i in clusters]
    inv6 = 1.0 / np.sqrt(6.0)
    u = np.empty((len(reps), 6), complex)
    u[:, 0] = inv6
    u[:, 1:] = np.exp(1j * reps) * inv6
    O = np.abs(u.conj() @ u.T)
    iu = np.triu_indices(len(reps), 1)
    print(f"STARVE beta=({beta[0]:.3f},{beta[1]:.3f},{beta[2]:.3f}) "
          f"h={h}: {len(C)} boxes, vol frac {frac:.2e}, "
          f"{len(clusters)} blobs, max diam {max(diam):.3f}, "
          f"min rep-overlap {O[iu].min() if len(iu[0]) else 1:.3f}, "
          f"{total} swept [{time.time()-t0:.0f}s]", flush=True)


def main():
    pts = [(1.5349059690692832, 0.5515593746413914, 2.9470193319949907),
           (1.1045934240780566, 0.6254576854126334, 0.8763425369356228)]
    for beta in pts:
        for h in (0.003, 0.01):
            probe(beta, h)


if __name__ == "__main__":
    main()
