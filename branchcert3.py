"""Fat-Q-tube certificate, roots-first architecture (4.15 corrected)."""
import time
import warnings

import numpy as np

from certify import SLOP, _g_and_J
from karlsson import karlsson_map
from parametric import polish_root, root_data2
from starve import fat_sweep
from tmres import certified_curve_residual, tm_karlsson

warnings.filterwarnings("ignore")


def certify_point(beta, h=3e-3, stride=None, smax=40.0, sig_min=0.05):
    t0 = time.time()
    hv = np.array([h, h, h])
    C, W, _ = fat_sweep(beta, h)
    print(f"  sweep: {len(C)} survivors [{time.time()-t0:.0f}s]", flush=True)
    H0 = karlsson_map(*beta)
    stride = stride or max(1, len(C) // 4000)
    roots, seen = [], []
    for c in C[::stride]:
        th = polish_root(H0, c)
        key = tuple(np.round(th / 2e-3).astype(int))
        if key in seen:
            continue
        seen.append(key)
        roots.append(th)
    # dedup finer (torus-aware quick pass)
    ded = []
    for th in roots:
        if not any(np.max(np.abs((th - o + np.pi) % (2*np.pi) - np.pi))
                   < 5e-3 for o in ded):
            ded.append(th)
    roots = ded
    print(f"  {len(roots)} roots from {len(C[::stride])} polishes "
          f"[{time.time()-t0:.0f}s]", flush=True)
    Htm = tm_karlsson(beta, hv)
    tubes, wild = [], 0
    for th in roots:
        _, J = _g_and_J(H0, th)
        sig = np.linalg.svd(J, compute_uv=False)[-1]
        try:
            S, Q, defect = root_data2(beta, th)
            resid = certified_curve_residual(beta, hv, th, S, Q, Htm=Htm)
        except Exception:
            wild += 1
            continue
        if sig < sig_min or np.abs(S).max() > smax or defect > 1e-6:
            wild += 1
            continue
        loc = np.abs(S) @ hv + 5e-3 + resid
        tubes.append((th, S, loc))
    print(f"  tubes {len(tubes)}, wild-roots {wild} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    # coverage: every survivor box must lie in some tube hull or count
    # as wild volume
    assigned = np.zeros(len(C), dtype=bool)
    for th, S, loc in tubes:
        d = np.abs((C - th + np.pi) % (2*np.pi) - np.pi)
        assigned |= (d <= loc + W + 0.02).all(axis=1)
    n_un = int((~assigned).sum())
    print(f"  coverage: {n_un} unassigned boxes "
          f"({100*n_un/len(C):.2f}%) [{time.time()-t0:.0f}s]", flush=True)
    n = len(tubes)
    cen = np.array([t[0] for t in tubes])
    Ss = np.array([t[1] for t in tubes])
    locs = np.array([t[2].sum() for t in tubes])
    inv6 = 1.0 / np.sqrt(6.0)
    u = np.empty((n, 6), complex)
    u[:, 0] = inv6
    u[:, 1:] = np.exp(1j * cen) * inv6
    O = np.abs(u.conj() @ u.T)
    Sd = np.abs(Ss[:, None] - Ss[None, :]).sum(axis=(2, 3))
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
    print(f"TUBE beta=({beta[0]:.3f},..): "
          f"{'CERTIFIED' if bound <= 5 and n_un == 0 else 'PARTIAL'} — "
          f"tubes {n}, wild {wild}, colors {cols}, bound {bound}, "
          f"uncovered {n_un} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    certify_point((1.5349059690692832, 0.5515593746413914,
                   2.9470193319949907))
