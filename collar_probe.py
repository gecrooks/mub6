"""Collar validation: coarse tiles near the theta=0 branch face with
anisotropic hv (theta thin, phi/lam fat). Measures whether the
assembled coarse machinery certifies in the rigidity collar, and at
what anisotropy."""
import itertools
import time
import warnings

import numpy as np

from certify import SLOP
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import polish_root, root_data2
from fold import valley_certificate
from rates import certified_rates

warnings.filterwarnings("ignore")


def _win_ok(beta, th, hv, cr=None):
    if cr is None:
        cr = certified_rates(beta, hv)
    vt = cr["far_tax"] * (1.0 / np.sqrt(6.0) + 0.06)
    th_s = polish_root(karlsson_map(*beta), th)
    try:
        c = valley_certificate(beta, th_s, np.asarray(hv), vt,
                               cert_rates=cr, corners_mode="lite")
        cc = c.get("cert")
        return (cc is not None and cc.get("consistent")
                and min(c["self_mins"]) > 0.05)
    except RuntimeError:
        return False


def collar_tile(beta, hv):
    t0 = time.time()
    hv = np.asarray(hv, float)
    H0 = karlsson_map(*beta)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:])
             for v in vecs]
    cr = certified_rates(beta, hv)
    slow = rec = 0
    wild, locs = [], {}
    for th in roots:
        if _win_ok(beta, th, hv, cr):
            slow += 1
            continue
        all_ok = True
        spread = 0.0
        for sx, sy, sz in itertools.product((-1, 1), repeat=3):
            sb = (beta[0] + sx * hv[0] / 2, beta[1] + sy * hv[1] / 2,
                  beta[2] + sz * hv[2] / 2)
            th_s = polish_root(karlsson_map(*sb), th)
            spread = max(spread, float(np.abs(
                (th_s - th + np.pi) % (2 * np.pi) - np.pi).sum()))
            if all_ok and not _win_ok(sb, th, hv / 2):
                all_ok = False
        if all_ok:
            rec += 1
        else:
            wild.append(th)
            locs[tuple(np.round(th, 6))] = 2.0 * spread + 0.15
    # pair coloring
    n = len(roots)
    drift = np.zeros(n)
    for i, th in enumerate(roots):
        key = tuple(np.round(th, 6))
        if key in locs:
            drift[i] = locs[key]
        else:
            try:
                S, _q, d = root_data2(beta, th)
                drift[i] = (np.abs(S) @ hv).sum() if d < 1e-6 else 0.75
            except Exception:
                drift[i] = 0.75
    inv6 = 1.0 / np.sqrt(6.0)
    u = np.empty((n, 6), complex)
    u[:, 0] = inv6
    u[:, 1:] = np.exp(1j * np.array(roots)) * inv6
    O = np.abs(u.conj() @ u.T)
    lo = O - (drift[:, None] + drift[None, :]) / 6.0 - SLOP
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
    print(f"COLLAR ({beta[0]:.4f},{beta[1]:.3f},{beta[2]:.3f}) "
          f"hv=({hv[0]:g},{hv[1]:g},{hv[2]:g}): "
          f"{'OK' if cols <= 5 else 'FAIL'} roots {len(roots)} "
          f"slow {slow} rec {rec} wild {len(wild)} colors {cols} "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    collar_tile((0.02, 2.0412, np.pi - 0.01), (2e-4, 3e-3, 3e-3))
    collar_tile((0.008, 2.0412, np.pi - 0.01), (1e-4, 2e-3, 2e-3))
