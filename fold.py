"""PROTOTYPE fold/ridge tubes: certify ill-conditioned roots over a
parameter tile via singular-frame reduction, replacing the (divergent)
slanted-tube machinery near fold strata.

Frame: at a flagged root theta0 with Jacobian J = U S V^T, let u = U[:,-1]
(near-left-null), w = V[-1] (near-right-null / ridge direction), Wc = the
complementary right frame (5x4), P = U[:, :4] (well-conditioned left
frame). Coordinates theta = theta0 + t w + Wc y.

Split g = 0 into:
  regular block  F(t, y, b) = P^T g = 0   (4x4, conditioning sigma_4 ~ 0.2)
     -> y = Y(t, b) by least-squares polish (implicit function; prototype
        certifies via sampled residuals, conditioning, and PAD)
  bifurcation scalar  phi(t, b) = u^T g(theta0 + t w + Wc Y(t,b), b)
     -> 1-dim: sample phi on a t-grid x beta-corners; root enclosures =
        grid cells where the sampled envelope dips below the certified
        floor; exclusion elsewhere; at most 2 enclosures accepted.

The certificate covers the WINDOW {|t| <= T} x {|y|_inf <= y_max} for all
b in the tile; its guard box (handed to the zoned sweep for geometric
collection) is the axis-aligned hull of the window, widened so that
outside the guard the plain far-field tax excludes: needs the certified
|phi| edge margin at |t| = T to exceed the far tax, and the y-shell
growth sigma_4 * dy to exceed it within the y-widening.

EMPIRICAL constants (sampled, PAD-padded) as everywhere in the prototype;
the rigorous version replaces samples with interval enclosures.
"""

import warnings

import numpy as np
from scipy.optimize import least_squares

from certify import _g_and_J, _uvec
from karlsson import karlsson_map
from parametric import PAD, g_at

warnings.filterwarnings("ignore")


def fold_frame(H0, th0):
    _, J = _g_and_J(H0, th0)
    U, sv, Vt = np.linalg.svd(J)
    return dict(u=U[:, -1], P=U[:, :4], w=Vt[-1], Wc=Vt[:4].T, sv=sv)


SLOPV = 1e-9        # float-eval + SVD-backward-error allowance per cell
SQ5 = np.sqrt(5.0)


def certified_valley_floors(beta, hv, th0, w, Wc, tgrid, Yc, cert_rates,
                            rho_tube):
    """Certified per-cell dichotomy floors for a valley window (R7).

    For each t-cell around anchor theta_hat_k = th0 + t_k w + Wc Yc[k]
    (anchors are DEFINITIONAL; only evaluations at them are certified):
    local SVD J_k = U S V^T, u_k = U[:,-1], P_k = U[:,:4]; decompose
    dtheta = tau w_k + Wc_k eta. With certified constants (hmag, c1
    from `rates.certified_rates`) the two sides are, over the cell
    (|tau| <= r_t) x (||eta||_2 <= rho) x (|dbeta_j| <= hv_j):

      quad(R)   = sqrt5 * (hmag^2/6 + hmag*(2/sqrt6 + hmag)/2...) — per
                  component 0.5*(diag+offdiag) row form; projected l2.
      Cb        = sum_j hv_j * ||(2/sqrt6) c1[j,:]||_2
      phi-side  |u.g| >= |phi_hat| - sv5*R2 - QP*R2^2 - Cb - SLOPV,
                  R2^2 = r_t^2 + r*^2       (valid for ||eta|| <= r*)
      F-side    ||P^T g|| >= sv4*rho - ||P^T gamma_hat|| - QP*(r_t^2 +
                  rho^2) - Cb - SLOPV       (valid at ||eta|| = rho;
                  concave in rho -> check endpoints r* and rho_tube)

    A cell is EXCLUDED if phi-floor > 0 and F-side > 0 on [r*, rho_tube];
    a DIP cell has F-side > 0 but phi-floor <= 0: its enclosure is
    cell x {||eta|| <= r*}. Returns floors, dip runs, edge clearance."""
    from karlsson import karlsson_map as kmap
    H0 = kmap(*beta)
    hvv = np.broadcast_to(np.asarray(hv, float), (3,))
    hmag = cert_rates["hmag"]
    c1 = cert_rates["c1"]
    inv6 = 1.0 / np.sqrt(6.0)
    s_beta = (c1 * hvv[:, None]).sum(axis=0) * inv6   # |s| drift from beta
    n_t = len(tgrid)
    r_t = 0.5 * (tgrid[1] - tgrid[0])
    floors = np.full(n_t, np.nan)
    f_ok = np.zeros(n_t, dtype=bool)
    r_stars = np.full(n_t, np.nan)
    anchors = np.zeros((n_t, 5))
    for k in range(n_t):
        th_k = th0 + tgrid[k] * w + Wc @ Yc[k]
        # local re-polish: the global-frame polyline leaves a residual
        # in the LOCAL regular frame (trench value bleeds into F at the
        # window ends); one transverse solve in the local frame kills it.
        _gam0, J0 = _g_and_J(H0, th_k)
        _U0, _sv0, Vt0 = np.linalg.svd(J0)
        Wc_loc, P_loc = Vt0[:4].T, _U0[:, :4]
        y_loc, _r = _solve_y(H0, th_k, Vt0[-1], Wc_loc, P_loc, 0.0,
                             np.zeros(4))
        th_k = th_k + Wc_loc @ y_loc
        anchors[k] = th_k
        gam, J = _g_and_J(H0, th_k)
        U, sv, _Vt = np.linalg.svd(J)
        # adaptive split: near-cusp cells (sv4 also small) use a 3/2
        # split — strong 3-dim regular block, 2-dim bifurcation vector
        # with l2 floor (roots need BOTH components zero; dip logic
        # ports, monotonicity is not used by the dichotomy)
        m = 3 if sv[3] < 0.1 else 4
        u_k = U[:, m:5]                       # (5, 5-m) bifurcation block
        P_k = U[:, :m]
        phi_hat = float(np.linalg.norm(u_k.T @ gam))
        Fres = float(np.linalg.norm(P_k.T @ gam))
        sv4, sv5 = sv[m - 1], sv[m]           # cond / linear-loss rates
        uvec = np.concatenate(([1.0], np.exp(1j * th_k))) * inv6
        s_hat = np.abs(uvec.conj() @ H0)          # (6,)

        def losses(rho):
            # |s_k| enclosure over cell x rho-ball x beta box
            R2 = float(np.sqrt(r_t ** 2 + rho ** 2))
            s_enc = np.minimum(
                s_hat + hmag * inv6 * SQ5 * R2 + s_beta, 1.0)
            # per-COMPONENT certified Hessian quad coefs, |s|-local:
            # q_k <= 0.5 (D_k + 4 O) R2^2; project by |u| (phi side)
            # and l2 (F side) instead of sqrt5-uniform
            Dk = 2.0 * (hmag ** 2 / 6.0 + hmag * inv6 * s_enc)   # (6,)
            O = 2.0 * hmag ** 2 / 6.0
            qk = 0.5 * (Dk + 4.0 * O)[1:6]                       # (5,)
            q_phi = float(np.linalg.norm(np.abs(u_k).T @ qk))
            q_F = float(np.linalg.norm(qk))
            # per-component beta losses, |s|-local; g-components 1..5
            rate = (2.0 * inv6) * s_enc[None, :] * c1        # (3, 6)
            bcol = (rate * hvv[:, None]).sum(axis=0)[1:6]    # (5,)
            return R2, q_phi, q_F, bcol

        # split radius: minimal r* clearing the F-side, small headroom;
        # one fixed-point pass (s_enc grows mildly with rho)
        r_star = 0.0
        for _ in range(2):
            _R2, _qp, qF, bcol = losses(r_star)
            cb_F = float(np.linalg.norm(bcol))
            r_star = 1.1 * (Fres + cb_F + qF * (r_t ** 2 + r_star ** 2)
                            + SLOPV + 1e-4) / max(sv4, 1e-6)
        ok = True
        for rho in (r_star, rho_tube):
            if rho < r_star - 1e-15:
                continue
            _R2r, _qpr, qFr, bcolr = losses(rho)
            m = sv4 * rho - Fres - qFr * (r_t ** 2 + rho ** 2) \
                - float(np.linalg.norm(bcolr)) - SLOPV
            if m <= 0:
                ok = False
        f_ok[k] = ok
        r_stars[k] = r_star
        R2s, qps, _qFs, bcols = losses(r_star)
        cb_phi = float(np.linalg.norm(np.abs(u_k).T @ bcols))
        floors[k] = (phi_hat - sv5 * R2s - qps * R2s ** 2
                     - cb_phi - SLOPV)
    # dip runs from certified floors
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
    edge_ok = (not low[0]) and (not low[-1]) and f_ok[0] and f_ok[-1]
    return dict(floors=floors, f_ok=f_ok, r_stars=r_stars, dips=merged,
                edge_ok=edge_ok, edge_margin=float(min(floors[0],
                                                       floors[-1])),
                r_t=r_t, anchors=anchors,
                all_f_ok=bool(f_ok.all()))


def _solve_y(H, th0, w, Wc, P, t, y0):
    def res(y):
        return P.T @ g_at(H, th0 + t * w + Wc @ y)

    r = least_squares(res, y0, method="lm", xtol=3e-16, ftol=3e-16,
                      gtol=3e-16)
    return r.x, float(np.max(np.abs(r.fun)))


def fold_certificate(beta, th0, h, T, n_t=41, verbose=False):
    """Certify the window {|t|<=T} x y-tube around a flagged root, for all
    b in the tile of half-width h. Returns dict or raises on failure."""
    H0 = karlsson_map(*beta)
    fr = fold_frame(H0, th0)
    u, P, w, Wc, sv = fr["u"], fr["P"], fr["w"], fr["Wc"], fr["sv"]
    sigma4 = sv[3]

    tgrid = np.linspace(-T, T, n_t)
    corners = [np.zeros(3)] + [h * np.array(c, float)
                               for c in [(1, 1, 1), (1, 1, -1), (1, -1, 1),
                                         (-1, 1, 1), (1, -1, -1),
                                         (-1, 1, -1), (-1, -1, 1),
                                         (-1, -1, -1)]]
    phi = np.zeros((len(corners), n_t))
    ymax = 0.0
    res_max = 0.0
    for ci, db in enumerate(corners):
        Hb = karlsson_map(beta[0] + db[0], beta[1] + db[1], beta[2] + db[2])
        y = np.zeros(4)
        # walk outward from t=0 in both directions, warm-starting y
        order = np.argsort(np.abs(tgrid), kind="stable")
        ys = {}
        for k in order:
            y0 = ys.get(k - 1 if tgrid[k] > 0 else k + 1, y)
            yk, rk = _solve_y(Hb, th0, w, Wc, P, tgrid[k], y0)
            ys[k] = yk
            res_max = max(res_max, rk)
            ymax = max(ymax, float(np.max(np.abs(yk))))
            phi[ci, k] = float(
                u @ g_at(Hb, th0 + tgrid[k] * w + Wc @ yk))
        y = ys[order[0]]

    # certified floor: sampling pad = PAD * (grid variation + beta spread)
    dt = tgrid[1] - tgrid[0]
    grid_var = np.max(np.abs(np.diff(phi, axis=1)))
    beta_spread = np.max(phi.max(axis=0) - phi.min(axis=0))
    floor = PAD * (0.5 * grid_var + 0.5 * beta_spread + res_max + 1e-9)

    env = np.min(np.abs(phi), axis=0)          # envelope over beta samples
    low = env < floor
    # enclosures: contiguous low runs, padded one cell
    runs = []
    k = 0
    while k < n_t:
        if low[k]:
            j = k
            while j + 1 < n_t and low[j + 1]:
                j += 1
            runs.append((max(0, k - 1), min(n_t - 1, j + 1)))
            k = j + 1
        else:
            k += 1
    # merge overlapping padded runs
    merged = []
    for a, b in runs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    if len(merged) > 2:
        raise RuntimeError(f"fold window has {len(merged)} phi-dips (>2)")
    if merged and (merged[0][0] == 0 or merged[-1][1] == n_t - 1):
        raise RuntimeError("phi-dip touches window edge (T too small)")

    edge_margin = float(min(env[0], env[-1]))
    # vectors representing the enclosure(s), overlap bounds handled upstream
    reps = []
    for a, b in merged:
        k = (a + b) // 2
        yk, _ = _solve_y(H0, th0, w, Wc, P, tgrid[k], np.zeros(4))
        reps.append(th0 + tgrid[k] * w + Wc @ yk)
    # sampled vectors across each enclosure + tile corners, for overlap rows
    span_vecs = []
    for a, b in merged:
        vs = []
        for k in range(a, b + 1):
            for ci, db in enumerate(corners[:5]):
                Hb = karlsson_map(beta[0] + db[0], beta[1] + db[1],
                                  beta[2] + db[2])
                yk, _ = _solve_y(Hb, th0, w, Wc, P, tgrid[k], np.zeros(4))
                vs.append(_uvec((th0 + tgrid[k] * w + Wc @ yk)[None, :])[0])
        span_vecs.append(np.array(vs))

    guard = np.abs(w) * T + np.abs(Wc) @ np.full(4, ymax) + 0.02
    return dict(n_enclosures=len(merged), reps=reps, span_vecs=span_vecs,
                guard=guard, edge_margin=edge_margin, floor=floor,
                sigma4=sigma4, ymax=ymax, res_max=res_max,
                enclosure_cells=merged, tgrid=tgrid)


def fold_overlap_rows(cert, other_vecs):
    """Certified-lower-bound overlaps between each fold enclosure's swept
    vectors and the other roots' center vectors. Returns (n_enc, n_other)
    of min |<span, other>| minus PAD spread."""
    rows = []
    for vs in cert["span_vecs"]:
        O = np.abs(vs.conj() @ np.asarray(other_vecs).T)
        lo = O.min(axis=0) - (PAD - 1.0) * (O.max(axis=0) - O.min(axis=0))
        rows.append(lo)
    return np.array(rows)


def _rho_y_estimate(fr, far_tax):
    """Conservative (large) rho_tube estimate for the in-loop certified
    F-side check; the final floors re-run with the true rho_y."""
    sig4 = fr["sv"][3]
    return max(0.05, PAD * (far_tax + 0.01) / max(0.6 * sig4, 0.02) + 0.02)


def valley_certificate(beta, th0, h, far_tax, n_t=65, T_cap=1.6,
                       verbose=False, cert_rates=None):
    """Certify a shallow-valley (or fold) root over the tile.

    Searches a window half-length T such that the trench floor |phi| at
    both window ends exceeds the sweep's exclusion needs (PAD-safe above
    far_tax), then certifies on a t-grid x beta-corner sample:
      - either phi is monotone with a single sign change (racing root), or
      - phi has <= 2 certified dips (fold pair) not touching the edges;
      - the y-block solve stays tight (residuals ~ machine, |y| bounded).
    Exports a tube ORACLE (frame + Y polyline + rho_y) so the zoned sweep
    can collect exactly the valley-tube region. EMPIRICAL constants,
    PAD-padded, as everywhere in the prototype."""
    H0 = karlsson_map(*beta)
    fr = fold_frame(H0, th0)
    u, P, w, Wc = fr["u"], fr["P"], fr["w"], fr["Wc"]

    corners = [np.zeros(3)] + [h * np.array(c, float)
                               for c in [(1, 1, 1), (1, 1, -1), (1, -1, 1),
                                         (-1, 1, 1), (1, -1, -1),
                                         (-1, 1, -1), (-1, -1, 1),
                                         (-1, -1, -1)]]

    def profile(T):
        tgrid = np.linspace(-T, T, n_t)
        phi = np.zeros((len(corners), n_t))
        Y = np.zeros((len(corners), n_t, 4))
        res_max = 0.0
        for ci, db in enumerate(corners):
            Hb = karlsson_map(beta[0] + db[0], beta[1] + db[1],
                              beta[2] + db[2])
            order = np.argsort(np.abs(tgrid), kind="stable")
            ys = {}
            for k in order:
                prev = k - 1 if tgrid[k] > 0 else k + 1
                y0 = ys.get(prev, np.zeros(4))
                yk, rk = _solve_y(Hb, th0, w, Wc, P, tgrid[k], y0)
                ys[k] = yk
                res_max = max(res_max, rk)
                Y[ci, k] = yk
                phi[ci, k] = float(
                    u @ g_at(Hb, th0 + tgrid[k] * w + Wc @ yk))
        return tgrid, phi, Y, res_max

    if cert_rates is not None:
        T_cap = max(T_cap, 2.0)   # certified edges get room to escalate
    T = 0.35
    accepted = None
    cert_T = None
    sampled_ok = None            # first window the SAMPLED test accepts
    while T <= T_cap:
        tgrid, phi, Y, res_max = profile(T)
        env = np.min(np.abs(phi), axis=0)
        # far_tax arrives pre-padded (trench-local rate); modest headroom
        need = 1.2 * far_tax + PAD * res_max + 1e-4
        if env[0] > need and env[-1] > need:
            if sampled_ok is None:
                sampled_ok = (tgrid, phi, Y, res_max)
            if cert_rates is not None:
                # R7 (report-only): escalate T until the PAD-free floors
                # also clear; if the cap runs out, fall back to the
                # sampled window with consistent=False recorded.
                cT = certified_valley_floors(
                    beta, np.broadcast_to(np.asarray(h, float), (3,)),
                    th0, w, Wc, tgrid, Y[0], cert_rates,
                    _rho_y_estimate(fr, far_tax))
                if not (cT["edge_ok"] and cT["all_f_ok"]
                        and all(a > 0 and b < len(tgrid) - 1
                                for a, b in cT["dips"])):
                    cert_T = None
                    T += 0.35
                    continue
                cert_T = cT
            accepted = (tgrid, phi, Y, res_max)
            break
        T += 0.35
    if accepted is None:
        accepted = sampled_ok    # certified escalation exhausted the cap
    if accepted is None:
        raise RuntimeError(f"valley floor never clears far tax by T={T_cap}")
    tgrid, phi, Y, res_max = accepted

    dt = tgrid[1] - tgrid[0]
    grid_var = np.max(np.abs(np.diff(phi, axis=1)))
    beta_spread = np.max(phi.max(axis=0) - phi.min(axis=0))
    floor = PAD * (0.5 * grid_var + 0.5 * beta_spread + res_max + 1e-9)

    # enclosure(s): cells where the envelope dips below the floor
    env = np.min(np.abs(phi), axis=0)
    low = env < floor
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
    if not merged:
        raise RuntimeError("no phi-dip found (lost the root?)")
    if len(merged) > 2:
        raise RuntimeError(f"{len(merged)} phi-dips (>2)")
    if merged[0][0] == 0 or merged[-1][1] == n_t - 1:
        raise RuntimeError("phi-dip touches window edge")

    # acceptance: monotone single-crossing OR clean dips
    mono = np.all(np.diff(phi, axis=1) > 0) or np.all(np.diff(phi, axis=1) < 0)
    if not mono and len(merged) == 1:
        # single dip without monotonicity: still fine (near-tangency);
        # enclosure logic already isolates it
        pass

    # y-tube: center polyline + radius over beta-samples. The collected
    # (fat) tube must extend to where the sigma4-growth of the regular
    # block beats the far tax, so no un-excludable, un-collected shell
    # remains; the thin-to-fat shell is root-free by the sigma4 bound.
    Yc = Y[0]                                  # center-beta polyline
    rho_thin = float(np.max(np.abs(Y - Yc[None]))) + 0.01
    sig4 = fr["sv"][3]
    rho_y = max(rho_thin,
                PAD * (far_tax + floor) / max(0.6 * sig4, 0.02) + 0.02)
    shell_margin = 0.6 * sig4 * rho_y - 0.5 * (11.0 / 18.0) * rho_y ** 2 \
        - PAD * res_max
    if shell_margin <= 0:
        raise RuntimeError(f"y-shell margin {shell_margin:.2e} <= 0")
    ymax = float(np.max(np.abs(Y)))

    # spanned vectors over each dip for overlap rows; also the certified
    # intra-dip self-overlap (all vectors in a dip mutually non-orthogonal
    # => any clique uses at most ONE vector from the dip => single-vertex
    # collapse with worst-case rows is sound even for fold PAIRS)
    span_vecs = []
    self_mins = []
    for a, b in merged:
        vs = []
        for k in range(a, b + 1):
            for ci in range(min(5, len(corners))):
                th = th0 + tgrid[k] * w + Wc @ Y[ci, k]
                vs.append(_uvec(th[None, :])[0])
        vs = np.array(vs)
        span_vecs.append(vs)
        O = np.abs(vs.conj() @ vs.T)
        iu = np.triu_indices(len(vs), 1)
        smin = float(O[iu].min()) if len(iu[0]) else 1.0
        spread = float(O[iu].max() - O[iu].min()) if len(iu[0]) else 0.0
        self_mins.append(smin - (PAD - 1.0) * spread)

    guard_hull = np.abs(w) * tgrid[-1] + np.abs(Wc) @ np.full(4, ymax + rho_y) \
        + 0.02
    oracle = dict(th0=th0.copy(), w=w.copy(), Wc=Wc.copy(),
                  tgrid=tgrid.copy(), Yc=Yc.copy(), rho_y=rho_y,
                  T=float(tgrid[-1]))
    cert = None
    if cert_rates is not None:
        # R7 certified dichotomy floors (report-only for now): same
        # window, PAD-free. Consistency: certified dips must lie within
        # one cell of the sampled enclosures.
        cert = certified_valley_floors(beta, np.broadcast_to(
            np.asarray(h, float), (3,)), th0, w, Wc, tgrid, Yc,
            cert_rates, rho_y)
        # consistency: certified analysis stands alone (edge + F-side
        # everywhere); each sampled enclosure must lie INSIDE a certified
        # dip (else the samples found a crossing the floors missed). A
        # certified dip may CONTAIN several sampled dips (fold pair whose
        # floors merge) — the <=2-dip cap and the downstream self-overlap
        # collapse still apply to the certified boxes.
        cert["consistent"] = (
            cert["edge_ok"] and cert["all_f_ok"]
            and len(cert["dips"]) <= 2
            and all(any(ca <= sa and sb <= cb
                        for (ca, cb) in cert["dips"])
                    for (sa, sb) in merged))
    return dict(n_enclosures=len(merged), span_vecs=span_vecs,
                self_mins=self_mins,
                guard=guard_hull, oracle=oracle,
                edge_margin=float(min(env[0], env[-1])), floor=floor,
                sigma4=fr["sv"][3], ymax=ymax, rho_y=rho_y,
                res_max=res_max, monotone=bool(mono), T=float(tgrid[-1]),
                cert=cert)
