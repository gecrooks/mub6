"""PROTOTYPE parametric Krawczyk: certify "no MU triple through {I, H(b)}"
for ALL b in a parameter tile of half-width h around a Karlsson point.

Architecture (Layer 1 of PROOF_ROADMAP.md):
  1. Pointwise pipeline at the tile center b0: roots + Krawczyk data.
  2. Per root: sensitivity S = -J^{-1} dg/db kills the first-order residual
     along the predicted curve theta*(b) ~ theta0 + S (b - b0); the root is
     enclosed in a slanted TUBE around that curve for all b in the tile.
  3. Far field: plain parameter-threading exclusion (tax ~ L_G*L_H*|db|).
  4. Near field per root: cascade of slanted sweeps whose beta-tax shrinks
     with the region radius D (drift of [J S + dg/db] over D), until the
     remainder fits a tube-Krawczyk uniqueness box.
  5. Pairs: first-order overlap models o(b) ~ o0 + G.db plus tube drift;
     partition certificate: greedy <=5-coloring whose same-color pairs are
     certified non-orthogonal over the tile => clique < 6 => no triple.

PROTOTYPE STATUS: constants marked EMPIRICAL are estimated by sampling the
Karlsson map (finite differences / tile-corner sampling), padded by PAD;
the rigorous version replaces them with certified interval enclosures of
H(beta) and its derivatives (one Arb evaluation per tile). Everything else
(Lipschitz/Hessian bounds on g in theta, sweep logic, Krawczyk algebra) is
inherited from the certified pointwise pipeline in certify.py.
"""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from certify import SLOP, L_H_G, _g_and_J, _torus_delta, _uvec, sweep
from karlsson import karlsson_map
from mub import find_mu_vectors, mu_vector_residuals

warnings.filterwarnings("ignore")

PAD = 2.0                    # padding on all EMPIRICAL sampled constants
HESS_ROW_TH = 11.0 / 18.0    # rigorous row bound on d2g/dtheta2 (certify.py)
SQ3 = np.sqrt(3.0)
SEP_GUARD = 0.55             # max allowed root-neighborhood radius


# ---------------------------------------------------------------------------
# map-side quantities (EMPIRICAL: finite differences of the Karlsson map)
# ---------------------------------------------------------------------------

def g_at(H, th):
    v = _uvec(th[None, :])[0]
    s = v @ H.conj()
    return (np.abs(s) ** 2 - 1.0 / 6.0)[1:6]


def map_lipschitz(beta, delta=1e-6):
    """max |dH/dbeta_k| entrywise (FD). EMPIRICAL."""
    H0 = karlsson_map(*beta)
    L = 0.0
    for k in range(3):
        b = list(beta)
        b[k] += delta
        L = max(L, np.max(np.abs(karlsson_map(*b) - H0)) / delta)
    return L


def dg_dbeta(beta, th, delta=1e-6):
    """5x3 matrix dg/dbeta at fixed theta (FD through the map). EMPIRICAL."""
    out = np.zeros((5, 3))
    g0 = g_at(karlsson_map(*beta), th)
    for k in range(3):
        b = list(beta)
        b[k] += delta
        out[:, k] = (g_at(karlsson_map(*b), th) - g0) / delta
    return out


def polish_root(H, th):
    return least_squares(mu_vector_residuals, th, args=([H],), method="lm",
                         xtol=3e-16, ftol=3e-16, gtol=3e-16).x


# ---------------------------------------------------------------------------
# per-root parametric data
# ---------------------------------------------------------------------------

def root_data(beta, th0):
    """Sensitivity S (5x3) and empirical first-order defect at the center."""
    H0 = karlsson_map(*beta)
    _, J = _g_and_J(H0, th0)
    Gb = dg_dbeta(beta, th0)
    S = -np.linalg.solve(J, Gb)
    defect = float(np.max(np.abs(J @ S + Gb)))       # ~ FD noise
    return S, defect


def root_data2(beta, th0, delta=2.5e-4):
    """S (5x3) plus curve Hessian Q (5x3x3) via FD tracking of the polished
    root, and the center defect. EMPIRICAL. Tracking guard: each stencil
    polish must stay within 0.1 of its linear prediction (same branch)."""
    S, defect = root_data(beta, th0)

    def track(db):
        H = karlsson_map(beta[0] + db[0], beta[1] + db[1], beta[2] + db[2])
        th = polish_root(H, th0 + S @ db)
        if np.max(np.abs(th - (th0 + S @ db))) > 0.1:
            raise RuntimeError("root tracking jumped branch")
        return th

    Q = np.zeros((5, 3, 3))
    thp, thm = {}, {}
    for k in range(3):
        e = np.zeros(3)
        e[k] = delta
        thp[k] = track(e)
        thm[k] = track(-e)
        Q[:, k, k] = (thp[k] + thm[k] - 2 * th0) / delta ** 2
    for k in range(3):
        for l in range(k + 1, 3):
            e = np.zeros(3)
            e[k] = delta
            e[l] = delta
            Qkl = (track(e) - thp[k] - thp[l] + th0) / delta ** 2
            Q[:, k, l] = Qkl
            Q[:, l, k] = Qkl
    return S, Q, defect


def curve_residual(beta, th0, S, Q, h, quadratic=True):
    """max |g| along the (linear or quadratic) predicted root curve over the
    26 corner/face points of the tile. This directly bounds the on-curve
    part of the slant tax: O(h^2) for the linear curve, O(h^3) with Q.
    EMPIRICAL (sampled, padded at use site)."""
    worst = 0.0
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if sx == sy == sz == 0:
                    continue
                db = h * np.array([sx, sy, sz], dtype=float)
                H = karlsson_map(beta[0] + db[0], beta[1] + db[1],
                                 beta[2] + db[2])
                th = th0 + S @ db
                if quadratic:
                    th = th + 0.5 * np.einsum("ikl,k,l->i", Q, db, db)
                worst = max(worst, float(np.max(np.abs(g_at(H, th)))))
    return worst


def q_offset(Q, h):
    """Per-coordinate bound on the quadratic curve's deviation from the
    linear frame: 0.5 sum_kl |Q_ikl| h_k h_l (h scalar or (3,))."""
    hv = np.broadcast_to(np.asarray(h, float), (3,))
    return 0.5 * np.einsum("ikl,k,l->i", np.abs(Q), hv, hv)


def sampled_tube_residual(beta, th0, S, h):
    """max over tile corners of |g(theta0 + S db, b0 + db)| — the
    second-order residual along the predicted curve. EMPIRICAL."""
    worst = 0.0
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                db = h * np.array([sx, sy, sz], dtype=float)
                H = karlsson_map(beta[0] + db[0], beta[1] + db[1],
                                 beta[2] + db[2])
                worst = max(worst, float(np.max(np.abs(g_at(H, th0 + S @ db)))))
    return worst


def sampled_J_drift(beta, th0, S, h, n=4, seed=1):
    """max sampled beta-drift of the Jacobian along the slanted curve.
    (theta-drift within the tube is bounded rigorously at use site.)
    EMPIRICAL."""
    rng = np.random.default_rng(seed)
    _, J0 = _g_and_J(karlsson_map(*beta), th0)
    worst = 0.0
    for _ in range(n):
        db = h * rng.choice([-1.0, 1.0], size=3)
        H = karlsson_map(beta[0] + db[0], beta[1] + db[1], beta[2] + db[2])
        _, J = _g_and_J(H, th0 + S @ db)
        worst = max(worst, float(np.max(np.abs(J - J0))))
    return worst


def sampled_gb_drift(beta, roots, D=0.3, n=5, seed=2):
    """max drift rate of dg/dbeta over |e| <= D around sampled roots,
    per unit |e|. EMPIRICAL."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for th0 in roots[:3]:
        G0 = dg_dbeta(beta, th0)
        for _ in range(n):
            e = rng.normal(size=5)
            e = D * e / np.max(np.abs(e))
            worst = max(worst, np.max(np.abs(dg_dbeta(beta, th0 + e) - G0)) / D)
    return worst


# ---------------------------------------------------------------------------
# slanted cascade + tube Krawczyk
# ---------------------------------------------------------------------------

def slant_tax(D, h, S_norm, gb_rate, defect):
    """Uniform bound on |g(theta0 + S db + e, b) - g(theta0 + e, b0)| for
    |e|_inf <= D, |db|_inf <= h: mean-value along beta with the integrand
    J.S + dg/db vanishing at the center; its drift over the region is
    (Hess_theta*reach*|S| + gb_rate*reach + defect), reach = D + |S db|."""
    reach = D + S_norm * SQ3 * h
    rate = HESS_ROW_TH * reach * S_norm + gb_rate * reach + defect
    return PAD * rate * SQ3 * h


def zoned_sweep(H0, roots, coef0, coef1, guards, far_tax,
                zone_R=0.85, wmin=1e-4, chunk=20_000, max_boxes=6e7,
                init_C=None, init_W=None, s_drift=0.0, oracles=None,
                beta_rate=None, cache=None, beta_unit=0.0,
                stuck_out=None):
    """Global sweep with per-box tax = min(plain threading, slanted bound of
    any root within zone_R): tax_i(box) = coef0_i + coef1_i * reach_i.
    Boxes componentwise inside a root's guard box are collected for that
    root's cascade. Returns (stuck, per-root interior reach D0 (n x 5),
    boxes_processed).

    Exclusion soundness: a box is excluded only if some |g_k(c)| exceeds
    its theta-Lipschitz spread (local |s| widened by s_drift for the tile's
    beta-variation) plus a valid bound on the beta-variation of g — plain
    L_G*L_H*|db| always valid; the slanted bound of root i (coef0_i =
    on-curve sampled residual + defect term, coef1_i = integrand drift
    rate) is valid with the slant applied, and the guard handoff keeps the
    slanted frames consistent with the near-zone pass (same tax, same
    frame).
    """
    R = np.array(roots)
    n = len(R)
    Hc = H0.conj()
    if init_C is None:
        stack_C = [np.full((1, 5), np.pi)]
        stack_W = [np.full((1, 5), np.pi)]
    else:
        stack_C = [np.asarray(init_C, float).reshape(-1, 5)]
        stack_W = [np.asarray(init_W, float).reshape(-1, 5)]
    D0 = np.zeros((n, 5))
    stuck = 0
    total = 0
    coef0 = np.asarray(coef0, float)
    coef1 = np.asarray(coef1, float)
    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        m = len(C)
        total += m
        if total > max_boxes:
            raise RuntimeError(f"zoned sweep exceeded {max_boxes:g} boxes")
        u = _uvec(C)
        s = u @ Hc
        g = np.abs(s) ** 2 - 1.0 / 6.0
        sw = W.sum(axis=1)
        # distances to all roots (m x n), Linf in theta
        delta = np.abs(_torus_delta(C[:, None, :], R[None, :, :]))
        dist = delta.max(axis=2)
        # per-box tax: plain vs best slanted root bound within zone
        reach = dist + sw[:, None]
        taxes = coef0[None, :] + coef1[None, :] * reach
        taxes = np.where(dist <= zone_R, taxes, np.inf)
        tax_box = np.minimum(far_tax, taxes.min(axis=1))
        smod_w = np.abs(s) + sw[:, None] / 6.0 + s_drift
        L = 2.0 * np.minimum(smod_w, 1.0) / 6.0
        margin = np.abs(g) - L * sw[:, None]
        if beta_rate is not None:
            # per-component beta tax: |dg_k| <= 2|s_k| ||dh_k|| |dbeta|,
            # with |s_k| widened over the box and tile
            beta_tax = beta_rate * np.minimum(smod_w, 1.0)
            tax_arr = np.minimum(beta_tax, tax_box[:, None])
        else:
            beta_tax = None
            tax_arr = tax_box[:, None]
        excl = (margin > SLOP + tax_arr).any(axis=1)
        if cache is not None and excl.any():
            # cacheable = excluded by the beta tax alone (drift semantics
            # valid for other tile centers); the rest (root-slant/far
            # exclusions) is stored for per-step re-sweeping
            if beta_tax is not None:
                exb = margin - SLOP - beta_tax
                kb = np.argmax(exb, axis=1)
                rows_i = np.arange(len(C))
                eb = exb[rows_i, kb]
                cacheable = excl & (eb > 0)
            else:
                cacheable = np.zeros(len(C), dtype=bool)
            if cacheable.any():
                sm = np.minimum(smod_w[rows_i, kb], 1.0)
                bu = (np.asarray(beta_unit)[kb[cacheable]]
                      if np.ndim(beta_unit) else beta_unit)
                cache["C"].append(C[cacheable].astype(np.float32))
                cache["W"].append(W[cacheable].astype(np.float32))
                cache["E"].append(eb[cacheable].astype(np.float32))
                cache["R"].append((bu * sm[cacheable]).astype(np.float32))
            rest = excl & ~cacheable
            if rest.any():
                cache["SC"].append(C[rest].astype(np.float32))
                cache["SW"].append(W[rest].astype(np.float32))
        keep = ~excl
        C, W, delta = C[keep], W[keep], delta[keep]
        if len(C):
            # guard collection: componentwise inside root i's guard box;
            # oracle roots use the fine valley-tube membership test instead
            inside_hull = (delta + W[:, None, :]
                           <= guards[None, :, :]).all(axis=2)
            collected = np.zeros(len(C), dtype=bool)
            for i in range(n):
                sel = inside_hull[:, i]
                if not sel.any():
                    continue
                if oracles is not None and oracles.get(i) is not None:
                    o = oracles[i]
                    d = _torus_delta(C[sel], o["th0"])
                    tb = d @ o["w"]
                    yb = d @ o["Wc"]
                    tW = W[sel] @ np.abs(o["w"])
                    yW = W[sel] @ np.abs(o["Wc"])
                    yc = np.stack([np.interp(tb, o["tgrid"], o["Yc"][:, j])
                                   for j in range(4)], axis=1)
                    fine = ((np.abs(tb) + tW <= o["T"]) &
                            ((np.abs(yb - yc) + yW <= o["rho_y"]).all(axis=1)))
                    collected[np.where(sel)[0][fine]] = True
                else:
                    collected |= sel
                    D0[i] = np.maximum(
                        D0[i], (delta[sel, i, :] + W[sel]).max(axis=0))
            if collected.any():
                C, W = C[~collected], W[~collected]
        if len(C):
            small = W.max(axis=1) <= wmin
            stuck += int(small.sum())
            if stuck_out is not None and small.any():
                stuck_out.append((C[small].copy(), W[small].copy()))
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
    return stuck, D0, total


def tube_krawczyk(H0, th0, R, rad_g, RJ_extra, iters=5):
    """Krawczyk in the slanted frame: unique root of g(.,b) within
    theta0 + S db +- R (scalar or per-coordinate vector) for every b in the
    tile. rad_g bounds the residual along the predicted curve; RJ_extra the
    slanted beta-drift of J."""
    r = np.full(5, float(R)) if np.ndim(R) == 0 else np.asarray(R, float).copy()
    g, J = _g_and_J(H0, th0)
    try:
        Y = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        return False, r
    rad_g_tot = rad_g + float(np.max(np.abs(g))) + SLOP
    ok = False
    for _ in range(iters):
        RJ_col = r / 3.0 + r.sum() / 18.0 + RJ_extra + SLOP
        M_mid = np.eye(5) - Y @ J
        M_rad = np.abs(Y).sum(axis=1)[:, None] * RJ_col[None, :]
        K_off = -Y @ g
        K_rad = (np.abs(Y) @ np.full(5, rad_g_tot)
                 + (np.abs(M_mid) + M_rad) @ r + SLOP)
        tot = np.abs(K_off) + K_rad
        if np.all(tot < r - SLOP):
            ok = True
            r = np.minimum(tot + SLOP, r)
        else:
            break
    return ok, r


def certify_root_tube(H0, th0, coef0_i, coef1_i, rad_g, RJ_extra, qoff, D0,
                      s_drift=0.0, beta_rate=None):
    """Near-zone pass for one root: per-box-tax sweep from the guard box
    down to a tube-Krawczyk box (widened per-coordinate by the quadratic
    curve offset qoff). Stuck wmin boxes => h too big at this root."""
    # largest tube radius at which Krawczyk contracts (pre-test candidates)
    R_K = None
    for cand in (0.06, 0.04, 0.02, 0.01, 5e-3, 2.5e-3, 1.2e-3):
        ok, rho = tube_krawczyk(H0, th0, cand + qoff, rad_g=rad_g,
                                RJ_extra=RJ_extra)
        if ok:
            R_K = cand
            break
    if R_K is None:
        return False, None, "no contracting tube radius"
    D = np.minimum(np.asarray(D0, float), SEP_GUARD)
    try:
        stuck, _D0_out, nb = zoned_sweep(
            H0, [th0], np.array([coef0_i]), np.array([coef1_i]),
            (R_K + qoff)[None, :], far_tax=np.inf,
            init_C=th0[None, :], init_W=D[None, :], max_boxes=5e7,
            wmin=3e-5, s_drift=s_drift, beta_rate=beta_rate)
    except RuntimeError as e:
        return False, None, str(e)
    if stuck:
        return False, None, f"near zone: {stuck} stuck boxes (D={D.max():.3f})"
    ok, rho = tube_krawczyk(H0, th0, R_K + qoff, rad_g=rad_g,
                            RJ_extra=RJ_extra)
    if not ok:
        return False, None, "tube contraction lost"
    return True, rho, dict(R_K=R_K, boxes=nb)


# ---------------------------------------------------------------------------
# pair overlap models and partition certificate
# ---------------------------------------------------------------------------

def overlap_gradients(beta, roots, delta=2e-5):
    """G[a,b,:] = d|<u_a,u_b>|/dbeta by tracking polished roots. EMPIRICAL."""
    U0 = _uvec(np.array(roots))
    O0 = np.abs(U0.conj() @ U0.T)
    n = len(roots)
    G = np.zeros((n, n, 3))
    for k in range(3):
        b = list(beta)
        b[k] += delta
        Hk = karlsson_map(*b)
        Uk = _uvec(np.array([polish_root(Hk, th) for th in roots]))
        G[:, :, k] = (np.abs(Uk.conj() @ Uk.T) - O0) / delta
    return O0, G


def color_conflicts(lo, n_colors=5, verbose=True):
    """Greedy <=5-coloring of the conflict graph {lo <= 0}. Same-color
    pairs are certified non-orthogonal over the tile => clique < 6."""
    n = lo.shape[0]
    conflict = [set() for _ in range(n)]
    n_conf = 0
    for a in range(n):
        for b in range(a + 1, n):
            if lo[a, b] <= 0:
                conflict[a].add(b)
                conflict[b].add(a)
                n_conf += 1
    colors = [-1] * n
    for v in sorted(range(n), key=lambda x: -len(conflict[x])):
        used = {colors[u] for u in conflict[v] if colors[u] >= 0}
        c = next((cc for cc in range(n_colors) if cc not in used), None)
        if c is None:
            return False, n_conf, None
        colors[v] = c
    if verbose:
        k = len(set(colors))
        print(f"    pairs: {n_conf} conflicts; colored with {k} colors "
              f"(clique <= {k} < 6)", flush=True)
    return True, n_conf, colors


def partition_certificate(O0, G, tubes, h, n_colors=5, verbose=True):
    """Certified pairwise lower bounds + coloring (no-fold path)."""
    n = O0.shape[0]
    lo = np.empty((n, n))
    for a in range(n):
        for b in range(n):
            drift = PAD * float(np.abs(G[a, b]).sum()) * h
            tube = (tubes[a].sum() + tubes[b].sum()) / 6.0
            lo[a, b] = O0[a, b] - drift - tube - SLOP
    ok, n_conf, colors = color_conflicts(lo, n_colors, verbose)
    return ok, n_conf, colors


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def certify_tile(beta, h, verbose=True, fold_cut=0.03, use_certified=False):
    """Attempt the full tile certificate at half-width h. Roots whose own
    slope ceiling is within 2.5x of the requested h (or sigma_min <
    fold_cut outright) are certified by valley windows (fold.py) instead
    of slanted Q-tubes. use_certified=True swaps the EMPIRICAL FD/sampled
    tax constants for the certified dual-AD rates of rates.py (remaining
    EMPIRICAL after the swap: curve residual corner-sampling, root
    sensitivities Sn/Q, overlap gradients)."""
    from fold import fold_overlap_rows, valley_certificate

    t0 = time.time()
    hv = np.broadcast_to(np.asarray(h, float), (3,)).copy()
    if float(np.ptp(hv)) > 0 and not use_certified:
        raise ValueError("anisotropic tiles require use_certified=True")
    h = float(np.linalg.norm(hv) / np.sqrt(3.0))   # l2-equivalent scalar
    L_map = map_lipschitz(beta)
    L_H = PAD * L_map
    cert_rates = None
    Htm = None
    if use_certified:
        from rates import certified_rates
        from tmres import tm_karlsson
        cert_rates = certified_rates(beta, hv)
        Htm = tm_karlsson(beta, hv)
        far_tax = cert_rates["far_tax"]
        beta_rate = cert_rates["beta_rate_vec"]
        valley_tax = cert_rates["far_tax"] * (1.0 / np.sqrt(6.0) + 0.06)
    else:
        far_tax = L_H_G * L_H * SQ3 * h
        # per-component beta-tax rate: PAD * 2 * ||dh_k||_2 * |dbeta|_2,
        # multiplied by the box-local |s_k| bound inside the sweep
        beta_rate = PAD * 2.0 * np.sqrt(6.0) * L_map * SQ3 * h
        # trench-local tax for valleys: |s_k| ~ 1/sqrt6 near MU points
        valley_tax = beta_rate * (1.0 / np.sqrt(6.0) + 0.06)

    H0 = karlsson_map(*beta)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:]) for v in vecs]
    if verbose:
        print(f"  tile h={h:g}: {len(roots)} roots at center", flush=True)

    # per-root parametric data: S, Q, taxes, guards -- or fold windows
    if cert_rates is not None:
        gb_rate = cert_rates["gb"]
        s_drift = cert_rates["s_drift"]
    else:
        gb_rate = sampled_gb_drift(beta, roots)
        s_drift = 2.5 * L_H * SQ3 * h
    n = len(roots)
    per_root = [None] * n
    fold_certs = {}
    coef0 = np.zeros(n)
    coef1 = np.zeros(n)
    guards = np.zeros((n, 5))
    oracles = {}
    for i, th0 in enumerate(roots):
        _, J = _g_and_J(H0, th0)
        _U, sv, Vt = np.linalg.svd(J)
        Gb = dg_dbeta(beta, th0)
        Sn_est = float(np.max(
            np.abs(np.linalg.lstsq(J, Gb, rcond=None)[0]) @ hv)) / h
        h_ceiling = sv[-1] / (PAD * SQ3 * (HESS_ROW_TH * Sn_est + gb_rate))
        if sv[-1] < fold_cut or h_ceiling < 2.5 * h:
            try:
                cert = valley_certificate(beta, th0, hv, valley_tax,
                                          cert_rates=cert_rates)
            except RuntimeError as e:
                return dict(ok=False, h=h, seconds=time.time() - t0,
                            reason=f"valley {i}: {e}")
            if min(cert["self_mins"]) <= 0.05:
                return dict(ok=False, h=h, seconds=time.time() - t0,
                            reason=f"valley {i}: dip self-overlap "
                                   f"{min(cert['self_mins']):.3f} too small "
                                   f"(clique collapse unsound)")
            fold_certs[i] = cert
            oracles[i] = cert["oracle"]
            guards[i] = cert["guard"]
            coef0[i] = np.inf
            coef1[i] = np.inf
            if verbose:
                print(f"    root {i}: VALLEY T={cert['T']:.2f} "
                      f"sig={sv[-1]:.4f} edge={cert['edge_margin']:.4f} "
                      f"rho_y={cert['rho_y']:.3f} "
                      f"mono={cert['monotone']}", flush=True)
                if cert.get("cert") is not None:
                    cc = cert["cert"]
                    print(f"      R7 floors: edge={cc['edge_margin']:.4f} "
                          f"dips={cc['dips']} consistent="
                          f"{cc['consistent']}", flush=True)
            continue
        S, Q, defect = root_data2(beta, th0)
        Sn = float(np.max(np.sum(np.abs(S), axis=1)))
        Sn_hv = float(np.max(np.abs(S) @ hv))     # = Sn*h when isotropic
        qoff = q_offset(Q, hv)
        if cert_rates is not None:
            from tmres import certified_curve_residual
            Rcurve = certified_curve_residual(beta, hv, th0, S, Q, Htm=Htm)
            rad_g = Rcurve + defect * SQ3 * h     # certified: no PAD
            coef0_res = Rcurve
            RJx = cert_rates["RJ_extra"]
        else:
            Rcurve = curve_residual(beta, th0, S, Q, h, quadratic=True)
            rad_g = PAD * Rcurve + defect * SQ3 * h
            coef0_res = PAD * Rcurve
            RJx = PAD * sampled_J_drift(beta, th0, S, h)
        coef1[i] = PAD * SQ3 * (HESS_ROW_TH * Sn_hv + gb_rate * h)
        coef0[i] = coef0_res + PAD * defect * SQ3 * h + coef1[i] * qoff.max()
        per_root[i] = dict(S=S, Q=Q, defect=defect, Sn=Sn, qoff=qoff,
                           rad_g=rad_g, RJ_extra=RJx)
        tax_est = coef0[i] + coef1[i] * 0.6 + 4e-4
        gj = 1.8 * np.abs(Vt.T) @ (tax_est / sv)
        guards[i] = np.clip(gj + 0.01, 0.03, 0.5)

    # STAGE A: coarse sweep (wmin 0.02). Unresolved coarse boxes are
    # either near-guard shells (stage B refines them) or un-excludable
    # birth blobs with no center root -> phantom anchors.
    stuck_boxes = []
    try:
        stuck, D0, nboxes = zoned_sweep(H0, roots, coef0, coef1, guards,
                                        far_tax, s_drift=s_drift,
                                        oracles=oracles, beta_rate=beta_rate,
                                        stuck_out=stuck_boxes, wmin=0.02)
    except RuntimeError as e:
        return dict(ok=False, h=h, seconds=time.time() - t0,
                    reason=f"zoned sweep A: {e}")
    phantoms = []
    if stuck:
        from certify import cluster_suspects
        pC = np.vstack([c for c, _ in stuck_boxes])
        pW = np.vstack([w for _, w in stuck_boxes])
        clusters = cluster_suspects(pC, pW, link=0.05)
        R_arr = np.array(roots)
        far_clusters = [cc for cc, _cr, _idx in clusters
                        if np.min(np.max(np.abs(
                            _torus_delta(R_arr, cc)), axis=1)) > 0.12]
        if verbose:
            print(f"    stage A: {stuck} coarse boxes, {len(clusters)} "
                  f"clusters ({len(far_clusters)} far from roots)",
                  flush=True)
        for cc in far_clusters:
            r = least_squares(mu_vector_residuals, cc, args=([H0],),
                              method="lm", xtol=3e-16, ftol=3e-16,
                              gtol=3e-16)
            th_p = r.x
            gmin = float(np.sqrt(2.0 * r.cost))
            if gmin > 3.0 * valley_tax + 1e-3:
                continue          # shallow pocket: stage B refines it away
            if np.min(np.max(np.abs(_torus_delta(R_arr, th_p)),
                             axis=1)) < 0.05:
                continue          # converged onto a real root's zone
            if any(np.max(np.abs(_torus_delta(np.array(p), th_p))) < 0.05
                   for p in phantoms):
                continue
            try:
                cert = valley_certificate(beta, th_p, hv, valley_tax,
                                          cert_rates=cert_rates)
            except RuntimeError as e:
                return dict(ok=False, h=h, seconds=time.time() - t0,
                            reason=f"phantom at |g|min "
                                   f"{np.sqrt(2*r.cost):.2e}: {e}")
            if min(cert["self_mins"]) <= 0.05:
                return dict(ok=False, h=h, seconds=time.time() - t0,
                            reason="phantom dip self-overlap too small")
            phantoms.append(th_p)
            i_ph = len(roots)
            roots.append(th_p)
            fold_certs[i_ph] = cert
            oracles[i_ph] = cert["oracle"]
            guards = np.vstack([guards, cert["guard"][None, :]])
            coef0 = np.append(coef0, np.inf)
            coef1 = np.append(coef1, np.inf)
            if verbose:
                print(f"    phantom {i_ph}: |g|min {np.sqrt(2*r.cost):.2e} "
                      f"T={cert['T']:.2f} dips={cert['n_enclosures']} "
                      f"mono={cert['monotone']}", flush=True)
        # STAGE B: refine exactly the diverted boxes at full resolution
        # with the (possibly) augmented oracles
        try:
            stuck, D0b, nb2 = zoned_sweep(
                H0, roots, coef0, coef1, guards, far_tax, s_drift=s_drift,
                oracles=oracles, beta_rate=beta_rate,
                init_C=pC, init_W=pW)
            nboxes += nb2
        except RuntimeError as e:
            return dict(ok=False, h=h, seconds=time.time() - t0,
                        reason=f"zoned sweep B: {e}")
        if stuck:
            return dict(ok=False, h=h, seconds=time.time() - t0,
                        reason=f"stage B: {stuck} stuck boxes")
        D0 = np.maximum(np.vstack([D0, np.zeros((len(phantoms), 5))]),
                        D0b)
        per_root.extend([None] * len(phantoms))
        n = len(roots)
    if verbose:
        print(f"    zoned sweep: {nboxes} boxes, max guard "
              f"{guards.max():.3f}, far tax {far_tax:.4f}, "
              f"folds {len(fold_certs)} (phantoms {len(phantoms)})",
              flush=True)

    rho_arr = np.zeros((n, 5))
    n_tubes = 0
    for i, th0 in enumerate(roots):
        if i in fold_certs:
            continue
        pr = per_root[i]
        ok, rho, info = certify_root_tube(
            H0, th0, coef0[i], coef1[i], pr["rad_g"], pr["RJ_extra"],
            pr["qoff"], D0=D0[i] + 0.01, s_drift=s_drift,
            beta_rate=beta_rate)
        if not ok:
            if verbose:
                print(f"    root {i}: TUBE FAILED ({info})", flush=True)
            return dict(ok=False, h=h, seconds=time.time() - t0,
                        reason=f"tube {i}: {info}")
        rho_arr[i] = rho
        n_tubes += 1
    if verbose:
        print(f"    {n_tubes} tubes + {len(fold_certs)} fold windows "
              f"certified", flush=True)

    # partition certificate: Q-tube rows analytic, fold rows span-sampled
    lo = np.empty((n, n))
    if cert_rates is not None:
        # fully certified pair bounds via TM inner products of u-curves
        from tmres import certified_overlap_lo, u_curve_tms
        curves = [None] * n
        for i in range(n):
            if i not in fold_certs and per_root[i] is not None:
                pr = per_root[i]
                curves[i] = u_curve_tms(hv, roots[i], pr["S"], pr["Q"])
        lo_tm = certified_overlap_lo(hv, curves, rho_arr, slop=SLOP)
        O0, G = overlap_gradients(beta, roots)   # fills non-TM pairs
        for a in range(n):
            for b in range(n):
                if a != b and not np.isnan(lo_tm[a, b]):
                    lo[a, b] = lo_tm[a, b]
                else:
                    drift = PAD * float(np.abs(G[a, b]) @ hv)
                    tube = (rho_arr[a].sum() + rho_arr[b].sum()) / 6.0
                    lo[a, b] = O0[a, b] - drift - tube - SLOP
    else:
        O0, G = overlap_gradients(beta, roots)
        for a in range(n):
            for b in range(n):
                drift = PAD * float(np.abs(G[a, b]).sum()) * h
                tube = (rho_arr[a].sum() + rho_arr[b].sum()) / 6.0
                lo[a, b] = O0[a, b] - drift - tube - SLOP
    if fold_certs:
        from fold import certified_dip_rows
        centers = _uvec(np.array(roots))
        for i, cert in fold_certs.items():
            rows = fold_overlap_rows(cert, centers)
            row = rows.min(axis=0) - rho_arr.sum(axis=1) / 6.0 - SLOP
            cc = cert.get("cert")
            if cert_rates is not None:
                # R7 ENFORCING: the certified floors must stand and the
                # interval rows/self-overlaps replace the sampled ones.
                if cc is None or not cc.get("consistent"):
                    return dict(ok=False, h=h, seconds=time.time() - t0,
                                reason=f"valley {i}: certified floors "
                                       f"inconsistent (R7)")
                crows, cselfs = certified_dip_rows(cc, centers)
                if min(cselfs) <= 0.05:
                    return dict(ok=False, h=h, seconds=time.time() - t0,
                                reason=f"valley {i}: certified self-"
                                       f"overlap {min(cselfs):.3f}")
                row = np.min(np.vstack(crows), axis=0) \
                    - rho_arr.sum(axis=1) / 6.0 - SLOP
                if verbose:
                    print(f"    R7 rows[{i}]: certified min "
                          f"{row.min():.3f} self {min(cselfs):.3f}",
                          flush=True)
            lo[i, :] = row
            lo[:, i] = row
            lo[i, i] = 1.0
    ok, n_conf, _ = color_conflicts(lo)
    dt = time.time() - t0
    if verbose:
        print(f"  ==> {'TILE CERTIFIED (prototype)' if ok else 'FAILED at coloring'}"
              f" at h={h:g}  [{dt:.0f} s]", flush=True)
    return dict(ok=ok, h=h, n_roots=n, n_folds=len(fold_certs),
                n_conflicts=n_conf, seconds=dt)


def main():
    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    print("=== parametric tile prototype at K6(5.9785, 4.0075, 1.6328) ===")
    print("(EMPIRICAL constants: FD map derivatives, corner-sampled residuals"
          f" and drifts, padded x{PAD}; rigorous version swaps in interval"
          " enclosures of the map.)", flush=True)
    results = []
    for h in (3e-4, 1e-3, 3e-3, 1e-2):
        res = certify_tile(beta, h)
        results.append(res)
        if not res.get("ok"):
            print(f"  tile h={h:g} FAILED: {res.get('reason', 'coloring')}",
                  flush=True)
            break
    good = [r for r in results if r.get("ok")]
    if good:
        hmax = max(r["h"] for r in good)
        t_tile = good[-1]["seconds"]
        n_tiles = (2 * np.pi / (2 * hmax)) ** 3
        print(f"\nmax certifying half-width: h = {hmax:g}  "
              f"(naive-threading equivalent was 2.9e-6)")
        print(f"naive full-family tile count: {n_tiles:.2e}; "
              f"cost/tile {t_tile:.0f} s => "
              f"{n_tiles * t_tile / 86400:.0f} CPU-days in Python "
              f"(before symmetry quotient, adaptivity, C kernel)")


if __name__ == "__main__":
    main()
