"""Certified tax rates for the tile certificates, derived from dual-AD
enclosures of the Karlsson map — replacing the EMPIRICAL FD/sampled
constants (map_lipschitz, sampled_gb_drift, sampled_J_drift).

Notation: s_k(theta, b) = sum_i conj(H_ik(b)) u_i(theta), |u_i| = 1/sqrt6,
g_k = |s_k|^2 - 1/6. With certified per-entry derivative magnitudes
e1[j][i][k] >= max over the tile box of |dH_ik/db_j| (dual AD, mean-value
valid pointwise over the box) and column l1-norms c1[j][k] = sum_i e1:

  |ds_k/db_j|      <= (1/sqrt6) c1[j][k]
  |dg_k/db_j|      <= 2 |s_k| (1/sqrt6) c1[j][k]
  beta tax (|db_j| <= h each):
      tax_k <= BR_k * min(|s_k|_widened, 1),  BR_k = (2/sqrt6) h sum_j c1[j][k]

  |ds_k/dth_l| = |H_lk|/sqrt6 <= hmag/sqrt6   (hmag >= max |H| entry)
  |d2 s_k/dth_l db_j| <= e1[j][l][k]/sqrt6
  theta-gradient of dg/db (row sum over l, |s|<=1):
      sum_l |d2 g_k/dth_l db_j| <= (2/sqrt6) c1[j][k] (hmag sqrt6/6 + 1)
  => certified gb rate  GB_k = sum_j h-free row sums (used as the slant
     tax's D-slope addition, in place of sampled_gb_drift)

  J_kl = 2 Re(conj(s_k) ds_k/dth_l):
  |dJ_kl/db_j| <= 2 [ (1/sqrt6) c1[j][k] hmag/sqrt6 + e1[j][l][k]/sqrt6 ]
  => RJ_cert = max_kl sum_j |dJ_kl/db_j| h   (replaces PAD*sampled_J_drift)

All quantities are upper bounds valid for EVERY b in the tile box; no
sampling, no PAD.
"""

import warnings

import numpy as np

from dual import cdual_mag, dual_karlsson
from interval import IV

warnings.filterwarnings("ignore")

SQ6 = np.sqrt(6.0)


MV_DEFAULT = False   # module switch: mean-value dual map for
                     # ALL certified_rates calls (4.76)


def certified_rates(beta, h, mv=None):
    """h: scalar or per-direction (3,) half-widths (anisotropic
    tiles). mv=True uses the mean-value dual map (4.76) — required
    near the branch surfaces (walls/corner), tighter everywhere.
    mv=None defers to the module switch MV_DEFAULT."""
    hv = np.broadcast_to(np.asarray(h, float), (3,)).copy()
    if mv is None:
        mv = MV_DEFAULT
    if mv:
        from dual import dual_karlsson_mv
        Hd = dual_karlsson_mv(beta, hv)
    else:
        Hd = dual_karlsson(IV(beta[0] - hv[0], beta[0] + hv[0]),
                           IV(beta[1] - hv[1], beta[1] + hv[1]),
                           IV(beta[2] - hv[2], beta[2] + hv[2]))
    e1 = np.zeros((3, 6, 6))
    hmag = 0.0
    for i in range(6):
        for k in range(6):
            hmag = max(hmag, Hd[i][k].v.mag().hi)
            for j in range(3):
                e1[j, i, k] = cdual_mag(Hd[i][k].d[j])
    c1 = e1.sum(axis=1)                          # c1[j, k]
    BR = (2.0 / SQ6) * (c1 * hv[:, None]).sum(axis=0)   # per-column tax
    # per-unit l2-distance drift rate per column (Cauchy-Schwarz over j)
    BU = (2.0 / SQ6) * np.sqrt((c1 ** 2).sum(axis=0))
    GB = (2.0 / SQ6) * c1.sum(axis=0) * (hmag * SQ6 / 6.0 + 1.0)
    dJj = np.zeros((3, 6, 6))
    for k in range(6):
        for l in range(6):
            for j in range(3):
                dJj[j, k, l] = 2.0 * ((1.0 / SQ6) * c1[j, k] * hmag / SQ6
                                      + e1[j, l, k] / SQ6)
    RJ = float((dJj * hv[:, None, None]).sum(axis=0).max())
    s_drift = float((1.0 / SQ6) * (c1 * hv[:, None]).sum(axis=0).max())
    # first-order sweep-tax data: center beta-gradient of H (point
    # evaluation) + certified deviation of the gradient over the tile
    # (enclosure minus point, componentwise; the curvature term needs
    # no second derivatives — the dual enclosure IS the bound)
    Hp = dual_karlsson(IV(beta[0]), IV(beta[1]), IV(beta[2]))
    dH0 = np.zeros((3, 6, 6), complex)
    WD = np.zeros((3, 6))
    for j in range(3):
        for i in range(6):
            for k in range(6):
                pe = Hp[i][k].d[j]
                pc = complex(0.5 * (pe.re.lo + pe.re.hi),
                             0.5 * (pe.im.lo + pe.im.hi))
                dH0[j, i, k] = pc
                te = Hd[i][k].d[j]
                dev = (max(te.re.hi - pc.real, pc.real - te.re.lo)
                       + max(te.im.hi - pc.imag, pc.imag - te.im.lo))
                WD[j, k] += dev / SQ6
    return dict(beta_rate_vec=BR, beta_unit_vec=BU, far_tax=float(BR.max()),
                gb=float(GB.max()), RJ_extra=RJ, hmag=hmag,
                c1=c1, s_drift=s_drift, dH0=dH0, WD=WD)


def chain_certified_rates(beta, h, span, axis=0):
    """Sup of the per-column unit drift rates BU over a chain segment
    [beta_axis, beta_axis + span], via a train of h-sized sub-boxes
    (each keeps tight h-scale constants; their max covers the span)."""
    n_sub = max(1, int(np.ceil((span + 2 * h) / (2 * h))))
    BU = np.zeros(6)
    for i in range(n_sub):
        c = list(beta)
        c[axis] = beta[axis] + (2 * i) * h
        r = certified_rates(tuple(c), h)
        BU = np.maximum(BU, r["beta_unit_vec"])
    return BU


def main():
    from parametric import map_lipschitz, PAD, SQ3
    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    h = 3e-4
    r = certified_rates(beta, h)
    L = map_lipschitz(beta)
    proto_beta_rate = PAD * 2.0 * SQ6 * L * SQ3 * h
    proto_far = 5.0 * PAD * L * SQ3 * h
    print(f"=== certified vs prototype rates at reference, h={h:g} ===")
    print(f"beta tax (at |s|=1): certified {r['far_tax']:.2e}  "
          f"prototype {proto_beta_rate:.2e}  ratio "
          f"{r['far_tax']/proto_beta_rate:.2f}")
    print(f"far tax:            certified {r['far_tax']:.2e}  "
          f"prototype {proto_far:.2e}  ratio {r['far_tax']/proto_far:.2f}")
    print(f"gb rate:            certified {r['gb']:.2f}  "
          f"(sampled+PAD was ~0.6-1.3)")
    print(f"RJ drift:           certified {r['RJ_extra']:.2e}")
    print(f"s_drift:            certified {r['s_drift']:.2e}")
    print(f"hmag = {r['hmag']:.6f} (1/sqrt6 = {1/SQ6:.6f})")


if __name__ == "__main__":
    main()
