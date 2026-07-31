"""coarse_certify_tile v1: the assembled fat-tile certificate.

Composition of tonight's measured parts, at h ~ 3e-3:
  1. enumerate + polish roots (multistart; chain mode will warm-track);
  2. per root: batched valley window (110 ms); failures recurse ONE
     octant level (8 sub-windows, measured sufficient for 11/12 at
     the worst census point); residual failures -> wild list;
  3. S-fat tube pair bounds + greedy coloring (correlated, colors ~2);
  4. coarse sweep (first-order taxes, wmin 5e-3): stuck blobs must
     polish into the certified root list (coverage semantics), else
     the tile FAILS loudly.
Verdict: CERTIFIED iff colors + wild <= 5 and coverage holds.
Prototype-grade composition; per-stage timings reported.
"""
import itertools
import time
import warnings

import numpy as np

from certify import SLOP, cluster_suspects
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import polish_root, root_data2, _g_and_J
from fold import valley_certificate
from rates import certified_rates
from starve import fat_sweep

warnings.filterwarnings("ignore")


def _window_ok(beta, th, h):
    hv = np.array([h] * 3)
    cr = certified_rates(beta, hv)
    vt = cr["far_tax"] * (1.0 / np.sqrt(6.0) + 0.06)
    Hs = karlsson_map(*beta)
    th_s = polish_root(Hs, th)
    try:
        c = valley_certificate(beta, th_s, hv, vt, cert_rates=cr,
                               corners_mode="lite")
        cc = c.get("cert")
        return (cc is not None and cc.get("consistent")
                and min(c["self_mins"]) > 0.05)
    except RuntimeError:
        return False


def coarse_certify_tile(beta, h=3e-3):
    t0 = time.time()
    H0 = karlsson_map(*beta)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:])
             for v in vecs]
    t_enum = time.time() - t0
    # stage 2: windows with one-level recursion
    t0w = time.time()
    slow = rec = 0
    wild = []
    for th in roots:
        if _window_ok(beta, th, h):
            slow += 1
            continue
        offs = [-h / 2, h / 2]
        all_ok = all(
            _window_ok((beta[0] + dx, beta[1] + dy, beta[2] + dz),
                       th, h / 2)
            for dx, dy, dz in itertools.product(offs, offs, offs))
        if all_ok:
            rec += 1
        else:
            wild.append(th)
    t_win = time.time() - t0w
    # stage 3: S-tube pair coloring
    t0p = time.time()
    wild_set = {tuple(np.round(w_, 6)) for w_ in wild}
    tubes = []
    for th in roots:
        is_wild = tuple(np.round(th, 6)) in wild_set
        S = None
        if not is_wild:
            try:
                S, Q, defect = root_data2(beta, th)
                if defect > 1e-6 or np.abs(S).max() >= 60:
                    S = None
            except Exception:
                S = None
        tubes.append((th, S))          # S=None -> wild vertex,
                                        # blanket localization
    n = len(tubes)
    cols = 0
    if n:
        cen = np.array([t[0] for t in tubes])
        WILD_LOC = 0.75          # blanket theta-l1 localization for
                                 # wild vertices (their blob extent)
        drift = np.array([
            (np.abs(t[1]).sum() * h if t[1] is not None else WILD_LOC)
            for t in tubes])
        inv6 = 1.0 / np.sqrt(6.0)
        u = np.empty((n, 6), complex)
        u[:, 0] = inv6
        u[:, 1:] = np.exp(1j * cen) * inv6
        O = np.abs(u.conj() @ u.T)
        lo = O - (drift[:, None] + drift[None, :]) / 6.0 \
            - 5e-3 / 6.0 - SLOP
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
    t_pair = time.time() - t0p
    # stage 4: coverage — coarse sweep stuck blobs must polish into roots
    t0s = time.time()
    C, W, swept = fat_sweep(beta, h, wmin=0.025)
    cl = cluster_suspects(C, W, link=0.1)
    uncovered = 0
    R = np.array(roots)
    for cc_, rr_, _ in cl:
        th = polish_root(H0, np.asarray(cc_))
        d = np.abs((th - R + np.pi) % (2 * np.pi) - np.pi).max(axis=1)
        if d.min() > 0.05:
            uncovered += 1
    t_swp = time.time() - t0s
    bound = cols                 # wilds are IN the coloring now
    ok = bound <= 5 and uncovered == 0
    print(f"COARSE_TILE ({beta[0]:.3f},{beta[1]:.3f},{beta[2]:.3f}) "
          f"h={h}: {'CERTIFIED' if ok else 'FAILED'} — roots "
          f"{len(roots)}, slow {slow}, recursed {rec}, wild "
          f"{len(wild)}, colors {cols}, bound {bound}, uncovered-blobs "
          f"{uncovered} | enum {t_enum:.0f}s win {t_win:.0f}s pair "
          f"{t_pair:.1f}s sweep {t_swp:.0f}s total "
          f"{time.time()-t0:.0f}s", flush=True)
    return ok


if __name__ == "__main__":
    for beta in [(1.33, 1.23, 1.42),
                 (0.51, 0.52, 2.93)]:
        coarse_certify_tile(beta)
