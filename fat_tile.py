"""Coverage-only fat tile (PROOF_SKELETON §8): the ball-vertex
lemma makes the pair layer nearly free; the coarse coverage sweep
is the sole real cost. Prototype measures the end-to-end cost of
the design that could take the campaign to ~10^2-10^3 GPU-hours.

Per tile at box half-width h:
  1. pool the roots (enumeration; coverage makes it sound);
  2. certified coverage: fat_sweep_hulls with certified taxes
     (mv rates feed the same constants) — every surviving hull
     cell must fit inside ball(root, R_LOC) (box-wise, loud);
  3. pairs: exact interval phasor sums with per-coordinate width
     2 R_LOC + BU_max h (certified blanket drift);
  4. chi <= 5 of the ball graph.
Returns CertificateResult; coverage evidence is RIGOROUS-track
(certified taxes + box-wise inclusion) pending the
coverage-verifier's resume/frontier accounting.
"""

import time

import numpy as np

from certificate_result import (CertificateGrade, CertificateResult,
                                Evidence)
from collar_tile import _pool_phases
from interval import IV, iv_cos, iv_sin


def _pair_lo(ti, tj, w):
    """Certified lower bound on |<u_i, u_j>| when each phase
    coordinate may move by +-w around the listed roots: interval
    phasor sum (ball-vertex budget)."""
    re_acc = IV(1.0)
    im_acc = IV(0.0)
    for k in range(5):
        d0 = float(((tj[k] - ti[k]) + np.pi) % (2 * np.pi) - np.pi)
        ph = IV(d0 - 2 * w, d0 + 2 * w)
        re_acc = re_acc + iv_cos(ph)
        im_acc = im_acc + iv_sin(ph)
    a2 = re_acc * re_acc + im_acc * im_acc
    return max(a2.lo, 0.0) ** 0.5 / 6.0


def fat_tile(beta, h, R_LOC=0.06, n_starts=4000):
    import rates as _rates
    from starve import fat_sweep_hulls
    t0 = time.time()
    beta = tuple(float(b) for b in beta)
    ph0 = _pool_phases(beta, n_starts=n_starts)
    n = len(ph0)
    deps = [Evidence("ball-vertex-lemma", CertificateGrade.RIGOROUS,
                     "SKELETON §8: r <= 0.88 >> R_LOC"),
            Evidence("enumeration", CertificateGrade.SAMPLED_BOUND,
                     "multistart pool; soundness via coverage")]
    if n == 0:
        return CertificateResult(False,
                                 CertificateGrade.SAMPLED_BOUND,
                                 tuple(deps),
                                 reason="empty pool (fail-closed)")
    # certified blanket drift: per-unit l2 drift BU (mv rates)
    _rates.MV_DEFAULT = True
    cr = _rates.certified_rates(beta, (h, h, h))
    BU = float(np.max(cr["beta_unit_vec"]))
    drift = BU * h * 3.0          # l1 over the box, conservative
    w = R_LOC + drift             # per-coordinate phase budget
    t_rates = time.time() - t0
    # coverage sweep (certified taxes inside fat_sweep_hulls)
    t1 = time.time()
    cen_h, rad_h, _cnt, swept = fat_sweep_hulls(beta, h,
                                                wmin=0.025,
                                                cell=0.025)
    C_h = np.atleast_2d(cen_h)
    r_h = np.atleast_2d(rad_h)
    uncovered = 0
    if len(C_h) and n:
        R = np.asarray(ph0)
        d = np.abs((C_h[:, None, :] - R[None, :, :] + np.pi)
                   % (2 * np.pi) - np.pi).max(axis=2)
        inside = (d.min(axis=1) + r_h.max(axis=1)) <= R_LOC
        uncovered = int((~inside).sum())
    t_cov = time.time() - t1
    deps.append(Evidence(
        "coverage", CertificateGrade.SAMPLED_BOUND,
        f"{len(C_h)} hull cells, {uncovered} uncovered "
        f"(box-wise, R_LOC={R_LOC}); certified taxes; "
        f"frontier accounting = coverage-verifier"))
    # pair layer: interval phasor bounds
    t2 = time.time()
    adj = np.zeros((n, n), dtype=bool)
    n_del = 0
    for i in range(n):
        for j in range(i + 1, n):
            if _pair_lo(ph0[i], ph0[j], w) > 0.0:
                n_del += 1
            else:
                adj[i, j] = adj[j, i] = True
    deps.append(Evidence("pairs", CertificateGrade.RIGOROUS,
                         f"{n_del} pairs non-orthogonal by "
                         f"interval phasor sum, budget w={w:.3f}"))
    best = n
    rng = np.random.default_rng(1)
    for t in range(400):
        order = (np.argsort(-adj.sum(axis=1)) if t == 0
                 else rng.permutation(n))
        col = -np.ones(n, dtype=int)
        for v in order:
            used = set(col[adj[v]]) - {-1}
            c = 0
            while c in used:
                c += 1
            col[v] = c
        best = min(best, int(col.max()) + 1)
        if best <= 5:
            break
    ok = best <= 5 and uncovered == 0
    res = CertificateResult(
        ok, CertificateGrade.SAMPLED_BOUND, tuple(deps),
        reason=f"chi<={best}, uncovered={uncovered}",
        metadata=dict(n_roots=n, chi=best, drift=drift, w=w,
                      swept=int(swept), t_rates=t_rates,
                      t_cov=t_cov, t_pairs=time.time() - t2,
                      dt=time.time() - t0))
    print(f"FAT_TILE {beta} h={h:g}: "
          f"{'CERTIFIED' if ok else 'FAILED'}[SAMPLED] — roots "
          f"{n}, chi {best}, uncovered {uncovered}, drift "
          f"{drift:.3f}, sweep {swept/1e6:.1f}M boxes "
          f"(cov {t_cov:.0f}s, total {time.time()-t0:.0f}s)",
          flush=True)
    return res


if __name__ == "__main__":
    fat_tile((0.5, 1.0, 2.0), 0.01)
