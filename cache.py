"""Margin-cached continuation: amortize the global sweep across a chain
of parameter tiles.

Anchor: one zoned sweep at beta_a stores, for every box excluded by the
|s|-local beta tax, the exclusion EXCESS (margin beyond the anchor tile's
own tax) and the box's beta-drift RATE (per unit l2 parameter distance).
Boxes excluded by root-slant taxes (frame-dependent) are stored verbatim.

Chain step to a tile centered beta_1 (same half-width h):
  1. vectorized check: cached box stays excluded iff E > R * dist(beta_1,
     anchor) * pad. (The cached excess already covers the anchor tile's
     own +-h box; moving the center by dist adds at most R*dist drift.)
  2. failing cached boxes + all slant boxes -> mini zoned sweep at beta_1
     (fresh taxes, guards anchored, oracles/roots re-polished).
  3. per-root stage at beta_1: Q-tubes / valley windows with D0 = anchor
     guard + root-motion pad; partition certificate.
Guards stay FIXED at the anchor (with slack); the chain re-anchors when
root motion approaches the slack or the patch fraction grows.

PROTOTYPE rigor model as everywhere (EMPIRICAL sampled constants, PAD).
"""

import time
import warnings

import numpy as np

from certify import SLOP, _g_and_J, _uvec
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import (HESS_ROW_TH, PAD, SQ3, certify_root_tube,
                        color_conflicts, curve_residual, dg_dbeta, g_at,
                        map_lipschitz, overlap_gradients, polish_root,
                        q_offset, root_data2, sampled_J_drift,
                        sampled_gb_drift, zoned_sweep)
from certify import L_H_G

warnings.filterwarnings("ignore")


def _per_root_stage(beta, h, roots, guards, coef0, coef1, per_root, D0pad,
                    s_drift, beta_rate, verbose=False):
    """Tubes + partition at center beta, guards given. Returns (ok, info)."""
    H0 = karlsson_map(*beta)
    n = len(roots)
    rho_arr = np.zeros((n, 5))
    for i, th0 in enumerate(roots):
        pr = per_root[i]
        ok, rho, info = certify_root_tube(
            H0, th0, coef0[i], coef1[i], pr["rad_g"], pr["RJ_extra"],
            pr["qoff"], D0=D0pad[i], s_drift=s_drift, beta_rate=beta_rate)
        if not ok:
            return False, f"tube {i}: {info}", None
        rho_arr[i] = rho
    O0, G = overlap_gradients(beta, roots)
    lo = np.empty((n, n))
    for a in range(n):
        for b in range(n):
            drift = PAD * float(np.abs(G[a, b]).sum()) * h
            tube = (rho_arr[a].sum() + rho_arr[b].sum()) / 6.0
            lo[a, b] = O0[a, b] - drift - tube - SLOP
    ok, n_conf, _ = color_conflicts(lo, verbose=verbose)
    return ok, f"{n_conf} conflicts", rho_arr


def _root_data_all(beta, h, roots, gb_rate):
    """Q-tube data for every root (chain demo assumes no valley roots)."""
    H0 = karlsson_map(*beta)
    n = len(roots)
    per_root = [None] * n
    coef0 = np.zeros(n)
    coef1 = np.zeros(n)
    sig_min = np.zeros(n)
    for i, th0 in enumerate(roots):
        S, Q, defect = root_data2(beta, th0)
        Sn = float(np.max(np.sum(np.abs(S), axis=1)))
        Rcurve = curve_residual(beta, th0, S, Q, h, quadratic=True)
        qoff = q_offset(Q, h)
        Jdrift = sampled_J_drift(beta, th0, S, h)
        coef1[i] = PAD * SQ3 * h * (HESS_ROW_TH * Sn + gb_rate)
        coef0[i] = PAD * Rcurve + PAD * defect * SQ3 * h \
            + coef1[i] * qoff.max()
        per_root[i] = dict(S=S, Q=Q, defect=defect, Sn=Sn, qoff=qoff,
                           rad_g=PAD * Rcurve + defect * SQ3 * h,
                           RJ_extra=PAD * Jdrift)
        _, J = _g_and_J(H0, th0)
        sig_min[i] = np.linalg.svd(J, compute_uv=False)[-1]
    return per_root, coef0, coef1, sig_min


def anchored_tile(beta_a, h, verbose=True):
    """Full certified tile at the anchor, exclusion cache collected."""
    t0 = time.time()
    L_map = map_lipschitz(beta_a)
    beta_rate = PAD * 2.0 * np.sqrt(6.0) * L_map * SQ3 * h
    beta_unit = PAD * 2.0 * np.sqrt(6.0) * L_map
    far_tax = L_H_G * PAD * L_map * SQ3 * h
    s_drift = 2.5 * PAD * L_map * SQ3 * h

    H0 = karlsson_map(*beta_a)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:]) for v in vecs]
    gb_rate = sampled_gb_drift(beta_a, roots)
    per_root, coef0, coef1, sig_min = _root_data_all(beta_a, h, roots,
                                                     gb_rate)
    if np.min(sig_min) < 0.012:
        raise RuntimeError("valley root present -- chain demo expects "
                           "Q-tube-only points")
    guards = np.zeros((len(roots), 5))
    for i, th0 in enumerate(roots):
        _, J = _g_and_J(H0, th0)
        _U, sv, Vt = np.linalg.svd(J)
        tax_est = coef0[i] + coef1[i] * 0.6 + 4e-4
        gj = 1.8 * np.abs(Vt.T) @ (tax_est / sv)
        guards[i] = np.clip(gj + 0.02, 0.04, 0.5)

    cache = dict(C=[], W=[], E=[], R=[], SC=[], SW=[])
    stuck, D0, nboxes = zoned_sweep(
        H0, roots, coef0, coef1, guards, far_tax, s_drift=s_drift,
        beta_rate=beta_rate, cache=cache, beta_unit=beta_unit)
    if stuck:
        raise RuntimeError(f"anchor sweep stuck boxes: {stuck}")
    C = np.vstack(cache["C"])
    W = np.vstack(cache["W"])
    E = np.concatenate(cache["E"])
    R = np.concatenate(cache["R"])
    SC = np.vstack(cache["SC"]) if cache["SC"] else np.zeros((0, 5),
                                                             np.float32)
    SW = np.vstack(cache["SW"]) if cache["SW"] else np.zeros((0, 5),
                                                             np.float32)
    D0pad = D0 + 0.01
    ok, msg, _ = _per_root_stage(beta_a, h, roots, guards, coef0, coef1,
                                 per_root, D0pad, s_drift, beta_rate)
    if not ok:
        raise RuntimeError(f"anchor per-root stage failed: {msg}")
    dt = time.time() - t0
    if verbose:
        mb = (C.nbytes + W.nbytes + E.nbytes + R.nbytes) / 2 ** 20
        print(f"  anchor: {nboxes} boxes, cached {len(C)} "
              f"(+{len(SC)} slant) [{mb:.0f} MB], {msg}, {dt:.0f} s",
              flush=True)
    return dict(beta_a=np.array(beta_a), h=h, roots=roots, guards=guards,
                D0pad=D0pad, C=C, W=W, E=E, R=R, SC=SC, SW=SW,
                gb_rate=gb_rate, L_map=L_map, anchor_seconds=dt)


def chain_step(state, beta_new, verbose=True):
    """Certify the tile at beta_new from the anchored cache."""
    t0 = time.time()
    h = state["h"]
    L_map = state["L_map"]
    beta_rate = PAD * 2.0 * np.sqrt(6.0) * L_map * SQ3 * h
    far_tax = L_H_G * PAD * L_map * SQ3 * h
    s_drift = 2.5 * PAD * L_map * SQ3 * h

    dist = float(np.linalg.norm(np.array(beta_new) - state["beta_a"]))
    ok_mask = state["E"] > state["R"] * dist * 1.05 + 1e-9
    n_fail = int((~ok_mask).sum())
    t_check = time.time() - t0

    # re-polish roots at beta_new (warm start from anchor roots)
    H1 = karlsson_map(*beta_new)
    roots1 = [polish_root(H1, th) for th in state["roots"]]
    motion = max(np.max(np.abs(np.asarray(r1) - np.asarray(r0)))
                 for r1, r0 in zip(roots1, state["roots"]))

    per_root, coef0, coef1, sig_min = _root_data_all(
        beta_new, h, roots1, state["gb_rate"])
    if np.min(sig_min) < 0.012:
        return dict(ok=False, reason="valley root emerged on chain")

    # patch sweep: failing cached boxes + all slant boxes, fresh taxes
    patch_C = np.vstack([state["C"][~ok_mask], state["SC"]]).astype(float)
    patch_W = np.vstack([state["W"][~ok_mask], state["SW"]]).astype(float)
    stuck, _D0, nboxes = zoned_sweep(
        H1, roots1, coef0, coef1, state["guards"], far_tax,
        s_drift=s_drift, beta_rate=beta_rate,
        init_C=patch_C, init_W=patch_W, max_boxes=3e7)
    if stuck:
        return dict(ok=False, reason=f"patch sweep stuck: {stuck}")
    t_patch = time.time() - t0 - t_check

    # per-root stage with anchor guards + motion pad
    D0pad = state["D0pad"] + motion + 0.005
    ok, msg, _ = _per_root_stage(beta_new, h, roots1, state["guards"],
                                 coef0, coef1, per_root, D0pad, s_drift,
                                 beta_rate)
    dt = time.time() - t0
    if verbose:
        print(f"    step dist={dist:.2e}: check {t_check*1e3:.0f} ms, "
              f"fail {n_fail} (+{len(state['SC'])} slant, "
              f"{nboxes} patch boxes, {t_patch:.1f} s), motion {motion:.4f},"
              f" {'OK' if ok else 'FAIL'} ({msg}) [{dt:.1f} s]", flush=True)
    return dict(ok=ok, reason=msg, seconds=dt, n_fail=n_fail,
                patch_boxes=nboxes, motion=motion, dist=dist)


def main():
    beta_a = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    h = 3e-4
    print(f"=== margin-cached chain at h={h:g}, anchor K6{beta_a} ===")
    state = anchored_tile(beta_a, h)

    step = 1.6 * h            # tile centers overlap slightly
    n_ok = 0
    for k in range(1, 13):
        beta_new = (beta_a[0] + k * step, beta_a[1], beta_a[2])
        res = chain_step(state, beta_new)
        if not res.get("ok"):
            print(f"    chain stopped at step {k}: {res.get('reason')}")
            break
        n_ok += 1
    span = n_ok * step + h
    print(f"\nchain: {n_ok} steps certified from one anchor; certified "
          f"theta-interval half-length ~ {span:.2e} "
          f"(vs single tile {h:g}) -- amortization factor ~ "
          f"{span/h:.1f}x on sweep volume")


if __name__ == "__main__":
    main()
