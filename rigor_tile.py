"""Theorem-grade signed tile: every ingredient except enumeration
coverage at RIGOROUS grade (PROOF_SKELETON §7 upgrade path).

The chain, per pair of roots over the beta-box B = beta0 +- hv:
  1. hslop: certified sup entry-drift of H over B from the box
     interval dual map (mean-value, no sampling).
  2. Parametric Krawczyk tubes: krawczyk_verify with this hslop
     gives, for EVERY beta in B, a unique root in c +- r_encl.
  3. Interval J over (tube x box): center J plus the rigorous
     Hessian-row remainder (11/18 row bound) + map-drift terms —
     the same algebra krawczyk_verify uses internally.
  4. Interval dg/dbeta over (tube x box): u over the tube (phases
     in c +- r via iv_exp_i), H and dH/dbeta from the box dual.
  5. Interval S = -J^{-1} dg/dbeta via midpoint solve + residual
     norm bound.
  6. Interval pair derivative dip_l = sum_k conj(u_i) u_j i
     (S_j - S_i)[k, l] / 6 — the CORRELATED rate; its magnitude is
     a certified sup of |d<u_i, u_j>/dbeta_l| over B (mean-value:
     NO curvature constants anywhere).
  7. Overlap enclosure over B: center value +- sum_l mag_l h_l
     (plus tube slop on the center evaluation).

Edge semantics (rigorous): overlap_lo > 0  => the pair is
certifiably non-orthogonal over ALL of B (edge deleted);
otherwise the pair is a conservative adjacency edge. chi <= 5 of
the adjacency graph over the enumerated roots then bounds the
orthogonality clique among ENUMERATED roots — enumeration
completeness (coverage) remains the SAMPLED_BOUND dependency,
owned by rigor/coverage-verifier.

Returns certificate_result.CertificateResult with per-dependency
evidence.
"""

import time

import numpy as np

from certificate_result import (CertificateGrade, CertificateResult,
                                Evidence)
from certify import SLOP, L_H_G, L_H_J, _g_and_J, krawczyk_verify
from collar_tile import _pool_phases
from dual import dual_karlsson
from interval import CIV, IV, iv_exp_i
from karlsson import karlsson_map

HESS_ROW = 11.0 / 18.0


def _box_dual(beta0, hv):
    return dual_karlsson(IV(beta0[0] - hv[0], beta0[0] + hv[0]),
                         IV(beta0[1] - hv[1], beta0[1] + hv[1]),
                         IV(beta0[2] - hv[2], beta0[2] + hv[2]))


def _civ_mag(z):
    a2 = z.abs2()
    return float(np.sqrt(max(a2.hi, 0.0)))


def _hslop(Hd, hv):
    """Certified sup entry drift of H over the box (mean value on
    the box-interval partials)."""
    m = 0.0
    for i in range(6):
        for k in range(6):
            s = sum(_civ_mag(Hd[i][k].d[j]) * hv[j]
                    for j in range(3))
            m = max(m, s)
    return m


def rigorous_signed_tile(theta_lo, theta_hi, b2, b3, hf, hf3=None,
                         pool=None, n_starts=6000):
    """Theorem-grade signed tile: pair layer on the CERT-blessed
    L4/L5 path — Q-curves with TM-certified residuals
    (tmres.certified_curve_residual), slanted-tube containment
    (tube_krawczyk with quadratic offset, as certify_root_tube),
    and TM inner-product overlap lower bounds
    (tmres.certified_overlap_lo). Every dependency RIGOROUS except
    enumeration coverage (rigor/coverage-verifier's lane) and the
    definitional float S, Q (containment is the certified claim).
    The beta-box is the cube |db|_inf <= h with
    h = max(span/2, hf, hf3) (over-cover is sound)."""
    import parametric as _pm
    from tmres import (certified_curve_residual, certified_overlap_lo,
                       tm_karlsson, u_curve_tms)
    t0 = time.time()
    hf3 = hf if hf3 is None else hf3
    span = theta_hi - theta_lo
    h = float(max(0.5 * span, hf, hf3))
    beta0 = np.array([0.5 * (theta_lo + theta_hi), b2, b3])
    hv = np.full(3, h)
    ph0 = (_pool_phases((theta_lo, b2, b3), n_starts=n_starts)
           if pool is None else np.asarray(pool))
    n = len(ph0)
    deps = [Evidence("map+rates+TM-residuals",
                     CertificateGrade.RIGOROUS,
                     "tm_karlsson TMs; certified_curve_residual "
                     "(no sampling, no PAD)"),
            Evidence("enumeration-coverage",
                     CertificateGrade.SAMPLED_BOUND,
                     "multistart pool; sound coverage owned by "
                     "rigor/coverage-verifier")]
    if n == 0:
        return CertificateResult(
            False, CertificateGrade.SAMPLED_BOUND, tuple(deps),
            reason="empty enumeration (fail-closed)")
    H0 = karlsson_map(*beta0)
    from rates import certified_rates
    cr = certified_rates(beta0, (h, h, h))
    RJx = cr["RJ_extra"]
    Htm = tm_karlsson(beta0, h)
    curves = []
    rho_arr = []
    n_tube = n_fail = 0
    for th in ph0:
        try:
            S, Q, defect = _pm.root_data2(beta0,
                                          np.asarray(th, float))
            Rcurve = certified_curve_residual(beta0, np.full(3, h),
                                              th, S, Q, Htm=Htm)
            rad_g = Rcurve + defect * np.sqrt(3.0) * h
            qoff = _pm.q_offset(Q, np.full(3, h))
            ok = False
            for cand in (0.02, 0.01, 5e-3, 2.5e-3, 1.2e-3, 6e-4):
                ok, rho = _pm.tube_krawczyk(
                    H0, np.asarray(th, float), cand + qoff,
                    rad_g=rad_g, RJ_extra=RJx)
                if ok:
                    break
            if not ok:
                raise RuntimeError("no contracting tube")
            curves.append(u_curve_tms(h, np.asarray(th, float),
                                      S, Q))
            rho_arr.append(rho)
            n_tube += 1
        except RuntimeError:
            curves.append(None)
            rho_arr.append(np.full(5, np.nan))
            n_fail += 1
    lo = certified_overlap_lo(h, curves, rho_arr)
    adj = np.zeros((n, n), dtype=bool)
    n_rig = 0
    for i in range(n):
        for j in range(i + 1, n):
            v = lo[i, j]
            if np.isfinite(v) and v > 0.0:
                n_rig += 1
            else:
                adj[i, j] = adj[j, i] = True
    deps.append(Evidence(
        "pair-overlap-lo", CertificateGrade.RIGOROUS,
        f"{n_rig} pairs deleted by TM inner-product bounds; "
        f"{n_fail} tube-less roots fail-closed to edges"))
    best = n
    rng = np.random.default_rng(1)
    for t in range(400):
        order = (np.argsort(-adj.sum(axis=1)) if t == 0
                 else rng.permutation(n))
        col = -np.ones(n, dtype=int)
        for v in order:
            used = set(col[adj[v]]) - {-1}
            cx = 0
            while cx in used:
                cx += 1
            col[v] = cx
        best = min(best, int(col.max()) + 1)
        if best <= 5:
            break
    deps.append(Evidence("coloring", CertificateGrade.RIGOROUS,
                         f"proper {best}-coloring exhibited"))
    ok = best <= 5
    res = CertificateResult(
        ok, CertificateGrade.SAMPLED_BOUND, tuple(deps),
        reason=f"chi<={best} over {n} enumerated roots",
        metadata=dict(n_roots=n, chi=best, n_rig_deleted=n_rig,
                      n_tube=n_tube, n_tube_fail=n_fail,
                      h=h, dt=time.time() - t0))
    print(f"RIGOR_TILE h={h:g} ({b2:.4f},{b3:.4f}) "
          f"th_lo={theta_lo:g}: {'CERTIFIED' if ok else 'FAILED'}"
          f"[{res.grade.name}] — roots {n}, chi {best}, tubes "
          f"{n_tube}, tube-fails {n_fail}, rig-deletions {n_rig} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return res


if __name__ == "__main__":
    rigorous_signed_tile(0.5, 0.5005, 1.0, 2.0,
                         hf=2.5e-4, hf3=2.5e-4)
