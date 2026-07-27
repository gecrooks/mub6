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


def valley_certificate(beta, th0, h, far_tax, n_t=65, T_cap=1.6,
                       verbose=False):
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

    T = 0.35
    accepted = None
    while T <= T_cap:
        tgrid, phi, Y, res_max = profile(T)
        env = np.min(np.abs(phi), axis=0)
        # far_tax arrives pre-padded (trench-local rate); modest headroom
        need = 1.2 * far_tax + PAD * res_max + 1e-4
        if env[0] > need and env[-1] > need:
            accepted = (tgrid, phi, Y, res_max)
            break
        T += 0.35
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
    return dict(n_enclosures=len(merged), span_vecs=span_vecs,
                self_mins=self_mins,
                guard=guard_hull, oracle=oracle,
                edge_margin=float(min(env[0], env[-1])), floor=floor,
                sigma4=fr["sv"][3], ymax=ymax, rho_y=rho_y,
                res_max=res_max, monotone=bool(mono), T=float(tgrid[-1]))
