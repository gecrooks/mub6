"""Coarse-mode chain: amortize coarse_certify_tile along a beta line.

Anchor: zoned_sweep at coarse params (wmin 5e-3, fo taxes, inf slant
coefs, zero guards -> pure tax exclusion) with the E/R exclusion
cache; windows classified and cached (class + recursion level +
wild_loc); stuck blobs recorded.
Step: (1) warm-polish cached roots (coverage: stuck blobs re-polish
into the root list, births appear as new blobs -> loud); (2) cache
re-verify E > R*dist, failures + anchor-stuck re-swept from their
boxes only; (3) windows re-run using cached class (known-fast skip
the doomed full-tile attempt); (4) pair coloring (cheap).
Measures the amortization factor = anchor_time / step_time.
"""
import itertools
import time
import warnings

import numpy as np

import parametric
from certify import cluster_suspects
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import polish_root, root_data2, zoned_sweep, SLOP
from fold import valley_certificate
from rates import certified_rates

warnings.filterwarnings("ignore")


def _win_ok(beta, th, h, cr=None):
    hv = np.array([h] * 3)
    if cr is None:
        cr = certified_rates(beta, hv)
    vt = cr["far_tax"] * (1.0 / np.sqrt(6.0) + 0.06)
    th_s = polish_root(karlsson_map(*beta), th)
    try:
        c = valley_certificate(beta, th_s, hv, vt, cert_rates=cr,
                               corners_mode="lite")
        cc = c.get("cert")
        return (cc is not None and cc.get("consistent")
                and min(c["self_mins"]) > 0.05)
    except RuntimeError:
        return False


def _classify(beta, th, h, cr):
    if _win_ok(beta, th, h, cr):
        return "slow", 0.0
    offs = [-h / 2, h / 2]
    spread = 0.0
    all_ok = True
    for dx, dy, dz in itertools.product(offs, offs, offs):
        sb = (beta[0] + dx, beta[1] + dy, beta[2] + dz)
        th_s = polish_root(karlsson_map(*sb), th)
        spread = max(spread, float(np.abs(
            (th_s - th + np.pi) % (2 * np.pi) - np.pi).sum()))
        if all_ok and not _win_ok(sb, th, h / 2):
            all_ok = False
    return ("rec" if all_ok else "wild"), 2.0 * spread + 0.15


def _pair_colors(beta, roots, classes, locs, h):
    n = len(roots)
    cen = np.array(roots)
    drift = np.zeros(n)
    for i, th in enumerate(roots):
        if classes[i] == "wild":
            drift[i] = locs[i]
        else:
            try:
                S, _q, d = root_data2(beta, th)
                drift[i] = (np.abs(S).sum() * h if d < 1e-6
                            else locs[i] or 0.75)
            except Exception:
                drift[i] = 0.75
    inv6 = 1.0 / np.sqrt(6.0)
    u = np.empty((n, 6), complex)
    u[:, 0] = inv6
    u[:, 1:] = np.exp(1j * cen) * inv6
    O = np.abs(u.conj() @ u.T)
    lo = O - (drift[:, None] + drift[None, :]) / 6.0 - 5e-3 / 6.0 - SLOP
    adj = (lo <= 0) & ~np.eye(n, dtype=bool)
    order = np.argsort(-adj.sum(axis=1))
    color = -np.ones(n, dtype=int)
    for v in order:
        used = set(color[adj[v]]) - {-1}
        cx = 0
        while cx in used:
            cx += 1
        color[v] = cx
    return int(color.max()) + 1


def anchor(beta, h):
    t0 = time.time()
    H0 = karlsson_map(*beta)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:])
             for v in vecs]
    cr = certified_rates(beta, (h, h, h))
    classes, locs = [], []
    for th in roots:
        cl, lc = _classify(beta, th, h, cr)
        classes.append(cl)
        locs.append(lc)
    cols = _pair_colors(beta, roots, classes, locs, h)
    n = len(roots)
    cache = dict(C=[], W=[], E=[], R=[], SC=[], SW=[], D1=[])
    stuck = []
    fo = dict(dH0c=[np.conj(cr["dH0"][j]) for j in range(3)],
              WD=cr["WD"], hv=np.array([h] * 3), s_drift=cr["s_drift"])
    st, _D0, nb = zoned_sweep(
        H0, roots, np.full(n, np.inf), np.full(n, np.inf),
        np.zeros((n, 5)), cr["far_tax"], s_drift=cr["s_drift"],
        beta_rate=cr["beta_rate_vec"], beta_unit=cr["beta_unit_vec"],
        cache=cache, stuck_out=stuck, wmin=0.025, fo=fo,
        max_boxes=4e8)
    C = np.vstack(cache["C"]) if cache["C"] else np.zeros((0, 5))
    W = np.vstack(cache["W"]) if cache["W"] else np.zeros((0, 5))
    E = np.concatenate(cache["E"]) if cache["E"] else np.zeros(0)
    R = np.concatenate(cache["R"]) if cache["R"] else np.zeros(0)
    D1 = np.vstack(cache["D1"]) if cache["D1"] else np.zeros((0, 3))
    # coverage at anchor (blobs persist along the chain)
    blobs = _blobs(stuck)
    unc = _coverage_blobs(H0, roots, blobs, beta=beta, h=h, cr=cr)
    ok = cols <= 5 and unc == 0
    print(f"  ANCHOR ({beta[0]:.3f},..): {'OK' if ok else 'FAIL'} "
          f"roots {n} colors {cols} stuck {st} cached {len(C)} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return dict(beta=np.array(beta), h=h, roots=roots,
                classes=classes, locs=locs, C=C, W=W, E=E, R=R,
                D1=D1, blobs=blobs, t_anchor=time.time() - t0)


def _blobs(stuck_bag):
    if not stuck_bag:
        return []
    C = np.vstack([c for c, w in stuck_bag])
    W = np.vstack([w for c, w in stuck_bag])
    cl = cluster_suspects(C, W, link=0.1)
    out = []
    rng = np.random.default_rng(3)
    for cc_, _r, idx in cl:
        m = np.asarray(idx)
        picks = C[m[rng.integers(0, len(m), size=min(5, len(m)))]]             if len(m) else np.asarray(cc_)[None, :]
        out.append((np.asarray(cc_), picks, C[m], W[m]))
    return out


def _refine_out(beta, h, boxes_C, boxes_W, cr):
    """Mini-sweep of a blob's members at a finer floor: True iff the
    whole blob EXCLUDES (whisker-parked w^4 tail, not root land)."""
    n = 1
    fo = dict(dH0c=[np.conj(cr["dH0"][j]) for j in range(3)],
              WD=cr["WD"], hv=np.array([h] * 3),
              s_drift=cr["s_drift"])
    stuck = []
    H = karlsson_map(*beta)
    st, _d, _nb = zoned_sweep(
        H, [np.zeros(5)], np.full(1, np.inf), np.full(1, np.inf),
        np.zeros((1, 5)), cr["far_tax"], s_drift=cr["s_drift"],
        beta_rate=cr["beta_rate_vec"], init_C=boxes_C, init_W=boxes_W,
        stuck_out=stuck, wmin=1.5e-3, fo=fo, max_boxes=2e6)
    return st == 0


def _coverage_blobs(H0, roots, blobs, beta=None, h=None, cr=None):
    R = np.array(roots)
    unc = 0
    for cc_, picks, bC, bW in blobs:
        covered = False
        for p in [cc_] + list(picks):
            th = polish_root(H0, np.asarray(p))
            d = np.abs((th - R + np.pi) % (2 * np.pi)
                       - np.pi).max(axis=1)
            if d.min() <= 0.05:
                covered = True
                break
        if not covered and cr is not None:
            covered = _refine_out(beta, h, bC, bW, cr)
        if not covered:
            unc += 1
    return unc


def chain_step(state, beta_new):
    t0 = time.time()
    h = state["h"]
    H1 = karlsson_map(*beta_new)
    roots1 = [polish_root(H1, th) for th in state["roots"]]
    cr = certified_rates(beta_new, (h, h, h))
    # windows with cached classes: known-fast skip the full attempt
    slow = rec = 0
    wild = []
    locs1 = []
    classes1 = []
    for th, cl0, lc0 in zip(roots1, state["classes"], state["locs"]):
        if cl0 == "slow" and _win_ok(beta_new, th, h, cr):
            classes1.append("slow")
            locs1.append(0.0)
            slow += 1
            continue
        if cl0 == "rec":
            # known fast mover: skip the doomed full-tile attempt
            offs = [-h / 2, h / 2]
            ok8 = all(_win_ok((beta_new[0] + dx, beta_new[1] + dy,
                               beta_new[2] + dz), th, h / 2)
                      for dx in offs for dy in offs for dz in offs)
            if ok8:
                classes1.append("rec")
                locs1.append(lc0)
                rec += 1
                continue
        if cl0 == "wild":
            # cached wild: re-measure spread cheaply (polish only)
            spread = 0.0
            offs = [-h / 2, h / 2]
            for dx in offs:
                for dy in offs:
                    for dz in offs:
                        sb = (beta_new[0] + dx, beta_new[1] + dy,
                              beta_new[2] + dz)
                        th_s = polish_root(karlsson_map(*sb), th)
                        spread = max(spread, float(np.abs(
                            (th_s - th + np.pi) % (2 * np.pi)
                            - np.pi).sum()))
            classes1.append("wild")
            locs1.append(2.0 * spread + 0.15)
            wild.append(th)
            continue
        cl, lc = _classify(beta_new, th, h, cr)
        classes1.append(cl)
        locs1.append(lc)
        if cl == "rec":
            rec += 1
        elif cl == "wild":
            wild.append(th)
    cols = _pair_colors(beta_new, roots1, classes1, locs1, h)
    # sweep re-verify from cache
    db = np.asarray(beta_new) - state["beta"]
    dist = float(np.linalg.norm(db))
    # first-order cache re-verify: E + D1.db - curv*dist^2, sup fallback
    CURV = 5.0
    fo_ok = (state["E"] + state["D1"] @ db.astype(np.float32)
             - CURV * dist * dist > 1e-4)
    okm = fo_ok | (state["E"] > state["R"] * dist * 1.05 + 1e-9)
    n_fail = int((~okm).sum())
    iC = state["C"][~okm]
    iW = state["W"][~okm]
    n = len(roots1)
    fo = dict(dH0c=[np.conj(cr["dH0"][j]) for j in range(3)],
              WD=cr["WD"], hv=np.array([h] * 3), s_drift=cr["s_drift"])
    stuck = []
    st, _D0, nb = zoned_sweep(
        H1, roots1, np.full(n, np.inf), np.full(n, np.inf),
        np.zeros((n, 5)), cr["far_tax"], s_drift=cr["s_drift"],
        beta_rate=cr["beta_rate_vec"], init_C=iC, init_W=iW,
        stuck_out=stuck, wmin=0.025, fo=fo, max_boxes=1e8)
    # persistent blobs: re-associate reps; NEW stuck (from cache
    # failures) would join as new blobs — count them loudly
    unc = _coverage_blobs(H1, roots1, state["blobs"],
                          beta=beta_new, h=h, cr=cr)
    if stuck:
        unc += _coverage_blobs(H1, roots1, _blobs(stuck),
                               beta=beta_new, h=h, cr=cr)
    ok = cols <= 5 and unc == 0
    dt = time.time() - t0
    print(f"  STEP -> ({beta_new[0]:.4f},..): {'OK' if ok else 'FAIL'} "
          f"colors {cols} cachefail {n_fail} resweep {nb} stuck {st} "
          f"unc {unc} [{dt:.0f}s]", flush=True)
    return ok, dt


if __name__ == "__main__":
    beta0 = (1.5349059690692832, 0.5515593746413914, 2.9470193319949907)
    h = 3e-3
    st = anchor(beta0, h)
    times = []
    for k in range(1, 7):
        bn = (beta0[0] + k * 1.6 * h, beta0[1], beta0[2])
        ok, dt = chain_step(st, bn)
        times.append(dt)
        if not ok:
            break
    if times:
        print(f"CHAIN: anchor {st['t_anchor']:.0f}s, steps mean "
              f"{np.mean(times):.0f}s -> amortization "
              f"{st['t_anchor']/np.mean(times):.1f}x", flush=True)
