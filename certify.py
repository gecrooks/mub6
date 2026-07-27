"""Certified enumeration of MU vectors and rigorous non-extendability.

Given a (near-)Hadamard matrix H0, this module produces a computer-assisted
proof of the statement:

    For EVERY matrix H with |H - H0|_entrywise <= HSLOP, the set of unit
    vectors mutually unbiased to both the identity basis and the columns
    of H consists of exactly N vectors, each localized in a verified
    enclosure; the orthogonality graph on these vectors is a subgraph of
    an explicit "possible-edge" graph. If that graph has clique number < 6,
    no orthonormal basis of MU vectors exists, hence no MU triple {I,H,K}
    and a fortiori no MU quadruple containing {I,H}.

Rigor model: pure IEEE-754 double arithmetic with explicit worst-case slop
constants added to every certified inequality. Expression graphs here are
short (<200 flops on O(1) quantities => absolute error < 3e-14); we carry
SLOP = 1e-11, a ~300x safety margin. A publication-grade version would swap
the slop model for ball arithmetic (Arb); the structure is Arb-ready.

Mathematical bounds used (u_j = exp(i th_j)/sqrt(6), th_0 = 0):
  s_k(th) = sum_j conj(H_jk) u_j,   g_k = |s_k|^2 - 1/6,
  |ds_k/dth_j| = 1/6 exactly,  |s_k| <= 1,
  |dg_k/dth_j| <= 2 |s_k| / 6   (local Lipschitz; global <= 1/3),
  |d2g_k/dth_j dth_l| <= 7/18 (diag), 1/18 (offdiag); row sum <= 11/18,
  |dg_k/dH_jk| <= 2 |s_k| |u_j| <= 2/sqrt(6) < 0.82  (6 entries per k),
  identity sum_k g_k = |u|^2 - 1 = 0 for exactly unitary H
    => the square subsystem g_1..g_5 = 0 suffices for exact Hadamard H.
"""

import numpy as np

SQ6 = np.sqrt(6.0)
SLOP = 1e-11            # float-evaluation slack per certified inequality
L_H_G = 5.0             # |dg/dH| * (6 entries): 6*0.82 => use 5.0
L_H_J = 6.0             # crude bound for Jacobian entry sensitivity to H
HESS_ROW = 11.0 / 18.0  # row-sum bound on Hessian of g_k


def _uvec(C):
    """Phase matrix (m,5) -> unit vectors (m,6) with first entry 1/sqrt6."""
    m = C.shape[0]
    u = np.empty((m, 6), dtype=complex)
    u[:, 0] = 1.0 / SQ6
    u[:, 1:] = np.exp(1j * C) / SQ6
    return u


def sweep(H, hslop=1e-11, wmin=5e-4, verbose=True, chunk=100_000,
          init_C=None, init_W=None, extra_slop=0.0,
          inner_C=None, inner_R=None, max_boxes=None):
    """Adaptive branch-and-bound over the phase 5-torus (or given boxes).

    LIFO-chunked (depth-first-ish) processing keeps peak memory bounded by
    ~depth*chunk boxes regardless of total workload. Returns
    (suspect_C, suspect_W, stats): every point of the start region is
    rigorously either inside some excluded box (no MU vector there, for any
    H' in the hslop-ball, up to the caller-supplied extra_slop allowance),
    inside a returned suspect box, or — when inner balls are given — inside
    a collected interior box. extra_slop is an additive exclusion-threshold
    term used by parametric certificates to absorb residual drift over a
    parameter tile.

    inner_C (k x 5) with inner_R (scalar): boxes entirely inside any inner
    ball (Linf) are collected WITHOUT further refinement and returned in
    stats['int_C'], stats['int_W'] — this is what keeps annulus sweeps from
    tiling thick un-excludable regions at wmin resolution. max_boxes: abort
    with RuntimeError beyond this workload (fail fast instead of OOM).
    """
    Hc = H.conj()
    if init_C is None:
        stack_C = [np.full((1, 5), np.pi)]
        stack_W = [np.full((1, 5), np.pi)]
    else:
        stack_C = [np.asarray(init_C, dtype=float).reshape(-1, 5)]
        stack_W = [np.asarray(init_W, dtype=float).reshape(-1, 5)]
    inner = None
    if inner_C is not None:
        inner = np.asarray(inner_C, dtype=float).reshape(-1, 5)
    int_C, int_W = [], []
    sus_C, sus_W = [], []
    n_excl = 0
    total_boxes = 0
    n_sus = 0
    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:                 # keep working set bounded
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        total_boxes += len(C)
        u = _uvec(C)
        s = u @ Hc                        # s[m,k] = sum_j u_j conj(H_jk)
        g = np.abs(s) ** 2 - 1.0 / 6.0
        sw = W.sum(axis=1)
        smod = np.abs(s)
        L = 2.0 * np.minimum(smod + sw[:, None] / 6.0 + L_H_G * hslop, 1.0) / 6.0
        margin = np.abs(g) - L * sw[:, None]
        excl = (margin > SLOP + L_H_G * hslop + extra_slop).any(axis=1)
        n_excl += int(excl.sum())
        keep = ~excl
        Ck, Wk = C[keep], W[keep]
        if inner is not None and len(Ck):
            # box fully inside some inner region -> collect, don't refine.
            # inner_R: scalar (Linf ball) or per-region (k,5) box radii.
            Rarr = np.asarray(inner_R, dtype=float)
            is_int = np.zeros(len(Ck), dtype=bool)
            for idx_p, p in enumerate(inner):
                d = np.abs(_torus_delta(Ck, p))
                r = Rarr if Rarr.ndim == 0 else Rarr[idx_p]
                is_int |= (d + Wk <= r).all(axis=1)
            if is_int.any():
                int_C.append(Ck[is_int])
                int_W.append(Wk[is_int])
                Ck, Wk = Ck[~is_int], Wk[~is_int]
        if max_boxes is not None and total_boxes > max_boxes:
            raise RuntimeError(f"sweep exceeded max_boxes={max_boxes}")
        small = Wk.max(axis=1) <= wmin
        if small.any():
            sus_C.append(Ck[small])
            sus_W.append(Wk[small])
            n_sus += int(small.sum())
        Cb, Wb = Ck[~small], Wk[~small]
        if len(Cb):
            j = np.argmax(Wb, axis=1)
            rows = np.arange(len(Cb))
            Wn = Wb.copy()
            Wn[rows, j] = Wb[rows, j] / 2.0
            Cl = Cb.copy()
            Cr = Cb.copy()
            Cl[rows, j] -= Wn[rows, j]
            Cr[rows, j] += Wn[rows, j]
            stack_C.append(np.vstack([Cl, Cr]))
            stack_W.append(np.vstack([Wn, Wn]))
        if verbose and total_boxes % 5_000_000 < chunk <= total_boxes:
            import resource
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 2**20
            print(f"    ...{total_boxes} boxes, {n_sus} suspects, "
                  f"stack {sum(len(x) for x in stack_C)}, peakRSS {rss}MB",
                  flush=True)
    sus_C = np.vstack(sus_C) if sus_C else np.zeros((0, 5))
    sus_W = np.vstack(sus_W) if sus_W else np.zeros((0, 5))
    stats = dict(excluded=n_excl, suspects=len(sus_C),
                 boxes_processed=total_boxes)
    if inner is not None:
        stats["int_C"] = np.vstack(int_C) if int_C else np.zeros((0, 5))
        stats["int_W"] = np.vstack(int_W) if int_W else np.zeros((0, 5))
    if verbose:
        print(f"  sweep: {total_boxes} boxes, "
              f"{n_excl} excluded, {len(sus_C)} suspect", flush=True)
    return sus_C, sus_W, stats


def _torus_delta(a, b):
    d = (a - b + np.pi) % (2 * np.pi) - np.pi
    return d


def cluster_suspects(sus_C, sus_W, link=0.02):
    """Grid-hash single-linkage clustering of suspect boxes on the torus.

    Two suspects land in one cluster if their grid cells (size=link) are
    equal or adjacent (torus-wrapped) — an over-approximation of the
    distance-link rule, safe because root separations are >> link.
    """
    import itertools as _it

    n = len(sus_C)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    ncells = max(1, int(round(2 * np.pi / link)))
    keys = np.floor(sus_C / link).astype(int) % ncells
    buckets = {}
    for i, k in enumerate(map(tuple, keys)):
        if k in buckets:
            union(i, buckets[k])
        else:
            buckets[k] = i
    offsets = [o for o in _it.product((-1, 0, 1), repeat=5) if o != (0,) * 5]
    for k, rep in buckets.items():
        for o in offsets:
            k2 = tuple((k[j] + o[j]) % ncells for j in range(5))
            if k2 in buckets:
                union(rep, buckets[k2])
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    clusters = []
    for idx in groups.values():
        base = sus_C[idx[0]]
        pts = np.array([base + _torus_delta(sus_C[i], base) for i in idx])
        wds = sus_W[idx]
        lo = (pts - wds).min(axis=0)
        hi = (pts + wds).max(axis=0)
        clusters.append(((lo + hi) / 2.0, (hi - lo) / 2.0, idx))
    return clusters


def _g_and_J(H, th):
    """g_1..g_5 and its 5x5 Jacobian at phase point th (analytic)."""
    Hc = H.conj()
    u = _uvec(th[None, :])[0]
    s = u @ Hc                                  # (6,)
    g = np.abs(s) ** 2 - 1.0 / 6.0
    ds = 1j * u[1:, None] * Hc[1:, :]           # ds[j-1,k] = dds_k/dth_j
    J = 2.0 * np.real(np.conj(s)[None, :] * ds).T   # J[k,j-1]
    return g[1:6], J[1:6, :]


def krawczyk_verify(H, center, radius, hslop=1e-11, iters=6):
    """Krawczyk test on the box (polished center) +- radius (vector).

    On success returns (True, c, r_encl, R_box): for every H' in the
    hslop-ball there is exactly one MU-vector phase point in the uniqueness
    box c +- R_box, and it lies within the tighter enclosure c +- r_encl.
    """
    c = center.copy()
    R = radius.copy()
    # float Newton polish of the center first (not part of the certificate)
    for _ in range(50):
        g, J = _g_and_J(H, c)
        try:
            step = np.linalg.solve(J, g)
        except np.linalg.LinAlgError:
            return False, c, R, R
        c = c - step
        if np.max(np.abs(step)) < 1e-14:
            break
    r = R.copy()
    ok = False
    for _ in range(iters):
        g, J = _g_and_J(H, c)
        rad_g = SLOP + L_H_G * hslop
        try:
            Y = np.linalg.inv(J)
        except np.linalg.LinAlgError:
            return False, c, r, R
        # entrywise Jacobian remainder over the box:
        # |J_kj(x)-J_kj(c)| <= (7/18) r_j + (1/18) sum_{l!=j} r_l
        RJ_col = r / 3.0 + r.sum() / 18.0 + L_H_J * hslop + SLOP
        M_mid = np.eye(5) - Y @ J
        M_rad = np.abs(Y).sum(axis=1)[:, None] * RJ_col[None, :]
        K_mid_off = -Y @ g
        K_rad = (np.abs(Y) @ np.full(5, rad_g)
                 + (np.abs(M_mid) + M_rad) @ r + SLOP)
        # K box = c + K_mid_off +- K_rad must sit strictly inside c +- r
        if np.all(np.abs(K_mid_off) + K_rad < r - SLOP):
            ok = True
            r = np.minimum(np.abs(K_mid_off) + K_rad + SLOP, r)
        else:
            break
    return ok, c, r, R


def _local_refine(H, boxes_C, boxes_W, hslop, wmin):
    """Re-run the branch-and-bound restricted to the given boxes."""
    Hc = H.conj()
    C, W = boxes_C.copy(), boxes_W.copy()
    out_C, out_W = [], []
    while len(C):
        u = _uvec(C)
        s = u @ Hc
        g = np.abs(s) ** 2 - 1.0 / 6.0
        sw = W.sum(axis=1)
        smod = np.abs(s)
        L = 2.0 * np.minimum(smod + sw[:, None] / 6.0 + L_H_G * hslop, 1.0) / 6.0
        margin = np.abs(g) - L * sw[:, None]
        excl = (margin > SLOP + L_H_G * hslop).any(axis=1)
        keep = ~excl
        Ck, Wk = C[keep], W[keep]
        small = Wk.max(axis=1) <= wmin
        if small.any():
            out_C.append(Ck[small])
            out_W.append(Wk[small])
        Cb, Wb = Ck[~small], Wk[~small]
        if len(Cb) == 0:
            break
        j = np.argmax(Wb, axis=1)
        rows = np.arange(len(Cb))
        Wn = Wb.copy()
        Wn[rows, j] = Wb[rows, j] / 2.0
        Cl, Cr = Cb.copy(), Cb.copy()
        Cl[rows, j] -= Wn[rows, j]
        Cr[rows, j] += Wn[rows, j]
        C = np.vstack([Cl, Cr])
        W = np.vstack([Wn, Wn])
    if out_C:
        return np.vstack(out_C), np.vstack(out_W)
    return np.zeros((0, 5)), np.zeros((0, 5))


def verify_all(H, sus_C, sus_W, hslop=1e-11, verbose=True):
    """Cluster suspects, Krawczyk each cluster, check coverage.

    On Krawczyk failure the cluster is locally re-swept 16x deeper (its
    suspect spread shrinks proportionally) and retried with a tighter box.
    Returns (enclosures, ok) where enclosures is a list of (c, r).
    """
    clusters = cluster_suspects(sus_C, sus_W)
    enclosures = []   # (c, r_encl, R_box)
    covered = np.zeros(len(sus_C), dtype=bool)
    for cc, cr, idx in clusters:
        R0 = cr + 8e-3
        ok, c, r_enc, Rbox = krawczyk_verify(H, cc, R0, hslop=hslop)
        if not ok:  # retry with a smaller uniqueness box
            ok, c, r_enc, Rbox = krawczyk_verify(H, cc, cr + 2e-3, hslop=hslop)
        refined_ok = False
        if not ok:  # local refinement: shrink the suspect spread 16x
            if verbose:
                print(f"  refining cluster at {np.round(cc, 3)} "
                      f"(spread {np.max(cr):.4f})")
            wloc = max(np.max(sus_W[idx]) / 16.0, 1e-6)
            rc_C, rc_W = _local_refine(H, sus_C[idx], sus_W[idx], hslop, wloc)
            if len(rc_C):
                # union bounding box of the refined (much tighter) suspects
                base = rc_C[0]
                pts = base + _torus_delta(rc_C, base)
                lo = (pts - rc_W).min(axis=0)
                hi = (pts + rc_W).max(axis=0)
                cc2, cr2 = (lo + hi) / 2.0, (hi - lo) / 2.0
                ok, c, r_enc, Rbox = krawczyk_verify(H, cc2, cr2 + 1.5e-3,
                                                     hslop=hslop)
                if ok:
                    # the refined sweep rigorously excluded the rest of the
                    # original boxes; coverage holds if every refined suspect
                    # lies inside the verified uniqueness box
                    inside = np.all(np.abs(_torus_delta(pts, c)) + rc_W
                                    <= Rbox + SLOP)
                    refined_ok = bool(inside)
        if not ok:
            if verbose:
                print(f"  KRAWCZYK FAILED for cluster at {np.round(cc, 3)}")
            return [], False
        if refined_ok:
            covered[list(idx)] = True
        else:
            # coverage: every suspect box of the cluster inside uniqueness box
            for i in idx:
                d = np.abs(_torus_delta(sus_C[i], c)) + sus_W[i]
                if np.all(d <= Rbox + SLOP):
                    covered[i] = True
        # merge if this root was already verified (overlapping clusters)
        for (c2, _r2, _R2) in enclosures:
            if np.max(np.abs(_torus_delta(c, c2))) < 1e-6:
                break
        else:
            enclosures.append((c, r_enc, Rbox))
    if not covered.all():
        if verbose:
            print(f"  COVERAGE GAP: {int((~covered).sum())} suspect boxes uncovered")
        return enclosures, False
    if verbose:
        print(f"  verified {len(enclosures)} unique roots; all suspects covered")
    return [(c, r) for c, r, _ in enclosures], True


def certified_graph(enclosures, verbose=True):
    """Certified pairwise overlap intervals between verified MU vectors.

    Edge (possible orthogonality) iff the interval around |<u,v>| reaches 0.
    Returns (possible_edges, min_nonedge_margin).
    """
    n = len(enclosures)
    cs = np.array([c for c, _ in enclosures])
    rs = np.array([r for _, r in enclosures])
    U = _uvec(cs)
    Z = np.abs(U.conj() @ U.T)
    edges = []
    min_margin = np.inf
    for a in range(n):
        for b in range(a + 1, n):
            drift = (rs[a].sum() + rs[b].sum()) / 6.0 + SLOP
            lo = Z[a, b] - drift
            if lo <= 0:
                edges.append((a, b))
            else:
                min_margin = min(min_margin, lo)
    if verbose:
        print(f"  possible-edge graph: {n} vertices, {len(edges)} edges "
              f"(min certified non-edge margin {min_margin:.4f})")
    return edges, min_margin


def clique_number(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    best = 0

    def extend(cl, cand):
        nonlocal best
        best = max(best, len(cl))
        if len(cl) + len(cand) <= best:
            return
        for v in sorted(cand):
            extend(cl + [v], cand & adj[v] & set(range(v + 1, n)))

    extend([], set(range(n)))
    return best


def certify_no_triple(name, H, hslop=1e-11, wmin=1e-4, save_suspects=None,
                      chunk=100_000):
    """Full pipeline. Returns dict with the certificate summary."""
    print(f"=== certifying {name} (H-ball radius {hslop:g}) ===")
    sus_C, sus_W, stats = sweep(H, hslop=hslop, wmin=wmin, chunk=chunk)
    if save_suspects:
        np.savez(save_suspects, sus_C=sus_C, sus_W=sus_W, H=H)
    enclosures, ok = verify_all(H, sus_C, sus_W, hslop=hslop)
    if not ok:
        print("  CERTIFICATION FAILED (refine wmin and retry)")
        return dict(name=name, ok=False, **stats)
    edges, margin = certified_graph(enclosures)
    omega = clique_number(len(enclosures), edges)
    print(f"  clique number of possible-edge graph: {omega}")
    if omega < 6:
        print(f"  ==> THEOREM: no MU triple {{I,H,K}} exists for any H in the "
              f"{hslop:g}-ball around {name}.")
    else:
        print(f"  (cliques of size {omega} present -- bases may exist; "
              f"no-triple NOT certified, as expected for this H)")
    return dict(name=name, ok=True, n_roots=len(enclosures),
                n_edges=len(edges), clique=omega, margin=margin, **stats)
