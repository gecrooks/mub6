"""Collar-tile demonstrator (PROOF_SKELETON §5, collar piece).

In the open collar (0 < theta < ~0.02) triples do not exist
(NOTES_LP_BRIDGE 4.28) but the bulk tile's pair coloring fails:
each near-basis inherited from the face is a 6-clique of
near-orthogonal pairs (overlaps ~ c_edge * theta, c_edge >= 0.42).
The collar fix is EDGE DELETION: a pair whose overlap is
certifiably positive over the beta-box is not an orthogonality
edge, and deleting the certifiable edges must bring the chromatic
number back to <= 5.

Key structure (4.31): the breaking overlap c * theta is MONOTONE
INCREASING in theta (c constant to 3 digits over an 8x theta
range), so its minimum over a theta-slab sits on the BOTTOM face.
The certificate therefore anchors at theta_lo and pays taxes only
for in-face drift (measured correlated rate ~ 4e-3, i.e. theta *
dc/dbeta — two orders below the anchored overlap), plus an FD
monotonicity check across the slab. A center-anchored symmetric
tax provably cannot work here: the overlap halves toward the
bottom while the tax stays put — that was the first (failed)
version of this file.

Demonstrator per slab [theta_lo, theta_hi] x in-face box (b2, b3)
+- hf, all sampled grade (FD rates x PAD; the certified pass
derives signed rates from the S-tube first-order data):
  1. pool the roots at the BOTTOM center (theta_lo, b2, b3);
  2. adjacency: possibly-orthogonal pairs under per-root phase
     drift (root_data2 S when sane, else 2x measured FD spread);
  3. in-face correlated rates from 4 in-face FD probes;
  4. theta-monotonicity: overlap must increase at theta mid + top;
  5. delete edges with overlap_lo - PAD * rate_if * 2 hf > 0 that
     pass the monotonicity check; greedy-color the remainder;
     CERTIFIED (demo grade) iff chi <= 5.
"""

import time

import numpy as np
from scipy.optimize import least_squares

from karlsson import karlsson_map
from mub import find_mu_vectors, mu_vector_residuals
from parametric import root_data2

PAD_CORR = 3.0    # pad on FD in-face correlated rates
MONO_SLOP = 0.2   # overlap must grow by >20% of itself across slab


def _polish(H, th):
    return least_squares(mu_vector_residuals, th, args=([H],),
                         method="lm", xtol=3e-16, ftol=3e-16,
                         gtol=3e-16).x


def _pool_phases(beta, n_starts=6000, seed=3):
    H = karlsson_map(*beta)
    ph = []
    for v in find_mu_vectors([H], n_starts=n_starts, seed=seed):
        t = _polish(H, np.angle(v * np.sqrt(6))[1:])
        if not any(np.max(np.abs(np.exp(1j * t) - np.exp(1j * u)))
                   < 1e-5 for u in ph):
            ph.append(t)
    return np.array(ph)


def _uvecs(ph):
    u = np.empty((len(ph), 6), complex)
    u[:, 0] = 1.0
    u[:, 1:] = np.exp(1j * ph)
    return u / np.sqrt(6.0)


def _repolish_pool(beta, ph0):
    H = karlsson_map(*beta)
    return np.array([_polish(H, t) for t in ph0])


def collar_tile(theta_lo, theta_hi, b2, b3, hf):
    t0 = time.time()
    b_lo = (theta_lo, b2, b3)
    ph0 = _pool_phases(b_lo)
    n = len(ph0)
    U0 = _uvecs(ph0)
    O0 = np.abs(U0.conj() @ U0.T)

    # per-root localization for the adjacency
    h_dir = np.array([theta_hi - theta_lo, hf, hf])
    probes, spreads = [], np.zeros(n)
    probe_betas = [(theta_lo, b2 + s * hf, b3) for s in (1, -1)] + \
                  [(theta_lo, b2, b3 + s * hf) for s in (1, -1)]
    for bp in probe_betas:
        php = _repolish_pool(bp, ph0)
        d = np.abs(np.exp(1j * php) - np.exp(1j * ph0)).max(axis=1)
        spreads = np.maximum(spreads, d)
        probes.append(np.abs(_uvecs(php).conj()
                             @ _uvecs(php).T))
    drift = np.empty(n)
    for i in range(n):
        try:
            S, _q, defect = root_data2(b_lo, ph0[i])
            if defect < 1e-6 and np.abs(S).max() < 60:
                drift[i] = float(np.abs(S) @ h_dir).real \
                    if np.ndim(np.abs(S) @ h_dir) == 0 \
                    else float((np.abs(S) @ h_dir).sum())
                continue
        except Exception:
            pass
        drift[i] = 2.0 * spreads[i] + 0.05      # measured, padded
    lo = O0 - (drift[:, None] + drift[None, :]) / 6.0 - 1e-3 / 6.0
    adj = (lo <= 0) & ~np.eye(n, dtype=bool)
    n_edges = int(adj.sum() // 2)

    # in-face correlated overlap rates (4 probes at slab bottom)
    rate_if = np.zeros((n, n))
    for Op in probes:
        rate_if = np.maximum(rate_if, np.abs(Op - O0) / hf)

    # theta-monotonicity across the slab (mid + top)
    grow = np.ones((n, n), dtype=bool)
    for th in (0.5 * (theta_lo + theta_hi), theta_hi):
        php = _repolish_pool((th, b2, b3), ph0)
        Op = np.abs(_uvecs(php).conj() @ _uvecs(php).T)
        grow &= Op >= (1.0 + MONO_SLOP) * O0

    dele = adj & grow & (O0 - PAD_CORR * rate_if * 2.0 * hf > 0)
    adj2 = adj & ~(dele | dele.T)
    n_del = int((dele | dele.T).sum() // 2)

    def chi(a):
        order = np.argsort(-a.sum(axis=1))
        col = -np.ones(n, dtype=int)
        for v in order:
            used = set(col[a[v]]) - {-1}
            c = 0
            while c in used:
                c += 1
            col[v] = c
        return int(col.max()) + 1

    c_before, c_after = chi(adj), chi(adj2)
    ok = c_after <= 5
    print(f"COLLAR_TILE [{theta_lo:g},{theta_hi:g}] "
          f"({b2:.4f},{b3:.4f}) hf={hf:g}: "
          f"{'CERTIFIED' if ok else 'FAILED'} — roots {n}, edges "
          f"{n_edges}, deleted {n_del}, chi {c_before} -> {c_after} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return ok


def main():
    for (b2, b3) in ((2.041236, np.pi), (1.0, 2.0), (1.1, 1.0428571)):
        collar_tile(0.005, 0.015, b2, b3, hf=5e-3)
    collar_tile(0.0025, 0.005, 2.041236, np.pi, hf=5e-3)


if __name__ == "__main__":
    main()
