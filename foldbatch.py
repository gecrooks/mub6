"""Batched certified valley floors: vectorizes the per-anchor loop of
fold.certified_valley_floors (anchor evals, local re-polish Newton,
5x5 SVDs, floor arithmetic) with numpy stacked linalg. Same
mathematics, same constants; per-window cost target ~5-10x down."""
import numpy as np

from fold import SLOPV, SQ5
from karlsson import karlsson_map

INV6 = 1.0 / np.sqrt(6.0)


def _g_J_batch(Hc, TH):
    """g[1:6] and 5x5 J for a batch of phase points TH (n,5)."""
    n = len(TH)
    u = np.empty((n, 6), complex)
    u[:, 0] = INV6
    u[:, 1:] = np.exp(1j * TH) * INV6
    s = u @ Hc                                   # (n, 6)
    g = np.abs(s) ** 2 - 1.0 / 6.0
    ds = 1j * u[:, 1:, None] * Hc[None, 1:, :]   # (n,5,6) dds_k/dth_j
    J = 2.0 * np.real(np.conj(s)[:, None, :] * ds)  # (n,5,6) J[j,k]
    return u, s, g[:, 1:6], np.swapaxes(J[:, :, 1:6], 1, 2)  # J (n,5,5)


def floors_batch(beta, hv, th0, w, Wc, tgrid, Yc, cert_rates, rho_tube):
    H0 = karlsson_map(*beta)
    Hc = H0.conj()
    hvv = np.broadcast_to(np.asarray(hv, float), (3,))
    hmag = cert_rates["hmag"]
    c1 = cert_rates["c1"]
    s_beta = (c1 * hvv[:, None]).sum(axis=0) * INV6
    n_t = len(tgrid)
    r_t = 0.5 * (tgrid[1] - tgrid[0])
    TH = th0[None, :] + tgrid[:, None] * w[None, :] + Yc @ Wc.T
    # local re-polish: 3 batched Newton steps in the local frame
    for _ in range(3):
        _u, _s, gv, J = _g_J_batch(Hc, TH)
        U, sv, Vt = np.linalg.svd(J)
        P = U[:, :, :4]                          # (n,5,4)
        Wl = np.swapaxes(Vt[:, :4, :], 1, 2)     # (n,5,4)
        F = np.einsum("nij,ni->nj", P, gv)       # (n,4)
        Jy = np.einsum("nij,nik,nkl->njl", P, J, Wl)
        step = np.linalg.solve(Jy, F[:, :, None])[:, :, 0]
        TH = TH - np.einsum("nij,nj->ni", Wl, step)
    u, s_c, gv, J = _g_J_batch(Hc, TH)
    U, sv, Vt = np.linalg.svd(J)
    # adaptive split per anchor
    m3 = sv[:, 3] < 0.05
    drift = np.zeros(n_t)
    drift[1:] = np.linalg.norm(np.diff(TH, axis=0), axis=1)
    R_ball = r_t + rho_tube + 0.5 * drift
    s_hat = np.abs(s_c)
    s_enc = np.minimum(s_hat + hmag * INV6 * SQ5 * R_ball[:, None]
                       + s_beta, 1.0)
    Dk = 2.0 * (hmag ** 2 / 6.0 + hmag * INV6 * s_enc)
    O = 2.0 * hmag ** 2 / 6.0
    qk = 0.5 * (Dk + 4.0 * O)[:, 1:6]            # (n,5)
    rate = (2.0 * INV6) * s_enc[:, None, :] * c1[None, :, :]
    bcol = (rate * hvv[None, :, None]).sum(axis=1)[:, 1:6]
    if "dH0" in cert_rates:
        bfo = np.zeros((n_t, 6))
        uc = u.conj()
        for j in range(3):
            sbj = uc @ cert_rates["dH0"][j]
            d0g = 2.0 * np.real(np.conj(s_c) * sbj)
            bfo += hvv[j] * (np.abs(d0g)
                             + 2.0 * s_beta.max() * np.abs(sbj)
                             + 2.0 * s_enc * cert_rates["WD"][j])
        bcol = np.minimum(bcol, bfo[:, 1:6] + SLOPV)
    floors = np.full(n_t, np.nan)
    f_ok = np.zeros(n_t, dtype=bool)
    r_stars = np.full(n_t, np.nan)
    for k in range(n_t):                          # cheap scalar tail
        m = 3 if m3[k] else 4
        u_k = U[k][:, m:5]
        P_k = U[k][:, :m]
        phi_hat = float(np.linalg.norm(u_k.T @ gv[k]))
        Fres = float(np.linalg.norm(P_k.T @ gv[k]))
        sv4, sv5 = sv[k, m - 1], sv[k, m]
        q_phi = float(np.linalg.norm(np.abs(u_k).T @ qk[k]))
        q_F = float(np.linalg.norm(qk[k]))
        cb_F = float(np.linalg.norm(bcol[k]))
        cb_phi = float(np.linalg.norm(np.abs(u_k).T @ bcol[k]))
        Rb = R_ball[k]
        r_star = 1.1 * (Fres + cb_F + q_F * Rb ** 2 + SLOPV + 1e-4) \
            / max(sv4, 1e-6)
        f_ok[k] = r_star <= Rb
        r_stars[k] = min(r_star, Rb)
        floors[k] = (phi_hat - sv5 * Rb - q_phi * Rb ** 2
                     - cb_phi - SLOPV)
    low = floors <= 0
    runs, k = [], 0
    while k < n_t:
        if low[k]:
            j = k
            while j + 1 < n_t and low[j + 1]:
                j += 1
            runs.append((max(0, k - 1), min(n_t - 1, j + 1)))
            k = j + 1
        else:
            k += 1
    merged = []
    for a, b in runs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    edge_ok = (not low[0]) and (not low[-1])
    return dict(floors=floors, f_ok=f_ok, r_stars=r_stars, dips=merged,
                edge_ok=edge_ok,
                edge_margin=float(min(floors[0], floors[-1])),
                r_t=r_t, anchors=TH, all_f_ok=bool(f_ok.all()))


def profile_batch(beta, th0, w, Wc, tgrid, corners, n_newton=10):
    """Batched replacement for valley_certificate's profile(): for
    each beta-corner, solve the transverse block at all t stations
    with vectorized Newton (init 0), return phi (n_c, n_t), Y
    (n_c, n_t, 4), res_max. Advisory in floors-first mode; anchors
    are definitional and residuals are measured."""
    n_t = len(tgrid)
    phi = np.zeros((len(corners), n_t))
    Y = np.zeros((len(corners), n_t, 4))
    res_max = 0.0
    from fold import fold_frame
    for ci, db in enumerate(corners):
        Hb = karlsson_map(beta[0] + db[0], beta[1] + db[1],
                          beta[2] + db[2])
        Hc = Hb.conj()
        TH = th0[None, :] + tgrid[:, None] * w[None, :]
        Yl = np.zeros((n_t, 4))
        from certify import _g_and_J
        _, J0 = _g_and_J(Hb, th0)
        U0, _sv, Vt0 = np.linalg.svd(J0)
        P = U0[:, :4]
        u5 = U0[:, 4]
        # block-warm-started Newton: stations ordered outward from
        # t=0, chunks of 16, each chunk warm-started from its inward
        # neighbor (branch tracking with batched speed)
        order = np.argsort(np.abs(tgrid), kind="stable")
        prev_y = {}
        B = 16
        for b0 in range(0, n_t, B):
            idx = order[b0:b0 + B]
            for i in idx:
                # nearest inward already-solved station
                cand = [j for j in (i - 1, i + 1) if j in prev_y]
                if cand:
                    Yl[i] = prev_y[cand[0]]
            for _ in range(n_newton):
                _u, _s, gv, J = _g_J_batch(Hc, TH[idx] + Yl[idx] @ Wc.T)
                F = gv @ P
                Jy = np.einsum("ik,nkl,lj->nij", P.T, J, Wc)
                Yl[idx] = Yl[idx] - np.linalg.solve(
                    Jy, F[:, :, None])[:, :, 0]
            for i in idx:
                prev_y[int(i)] = Yl[i]
        _u, _s, gv, J = _g_J_batch(Hc, TH + Yl @ Wc.T)
        F = gv @ P
        res_max = max(res_max, float(np.abs(F).max()))
        phi[ci] = gv @ u5
        Y[ci] = Yl
    return phi, Y, res_max
