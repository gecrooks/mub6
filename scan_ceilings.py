"""Ceiling statistics across the Karlsson family.

For ~24 random K6(3) points, measure each root's tax-slope ceiling

    h_i = sigma_min,i / (2 sqrt3 PAD (HESS |S_i| + g_beta)),

whose minimum over roots bounds the first-order parametric tile size at
that point. Outputs the distribution (is the reference point's ~1.4e-3
worst-root wall typical?) and the implied adaptive campaign integral
N ~ V * E[1/(2h)^3].
"""

import warnings

import numpy as np

from certify import _g_and_J, _torus_delta
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import (HESS_ROW_TH, PAD, SQ3, polish_root, root_data,
                        sampled_gb_drift)

warnings.filterwarnings("ignore")


def scan_point(beta, seed):
    H = karlsson_map(*beta)
    vecs = find_mu_vectors([H], n_starts=3000, seed=seed)
    roots = [polish_root(H, np.angle(v * np.sqrt(6))[1:]) for v in vecs]
    gb = sampled_gb_drift(beta, roots)
    ceils, sigs, Sns = [], [], []
    for th in roots:
        S, _defect = root_data(beta, th)
        Sn = float(np.max(np.sum(np.abs(S), axis=1)))
        _, J = _g_and_J(H, th)
        sig = float(np.linalg.svd(J, compute_uv=False)[-1])
        sigs.append(sig)
        Sns.append(Sn)
        ceils.append(sig / (PAD * SQ3 * (HESS_ROW_TH * Sn + gb)))
    R = np.array(roots)
    n = len(R)
    sep = min(np.max(np.abs(_torus_delta(R[a], R[b])))
              for a in range(n) for b in range(a + 1, n))
    return dict(n_roots=n, gb=gb, min_sep=sep,
                min_sig=min(sigs), max_Sn=max(Sns),
                h_slope=min(ceils), h_median_root=float(np.median(ceils)))


def main():
    rng = np.random.default_rng(20260726)
    rows = []
    for t in range(24):
        beta = tuple(rng.uniform(0.25, 2 * np.pi - 0.25, 3))
        try:
            r = scan_point(beta, seed=t)
        except Exception as e:
            print(f"pt {t}: FAILED ({type(e).__name__}: {e})", flush=True)
            continue
        rows.append(r)
        print(f"pt {t}: roots={r['n_roots']:2d} minsep={r['min_sep']:.2f} "
              f"minsig={r['min_sig']:.4f} maxSn={r['max_Sn']:5.1f} "
              f"gb={r['gb']:.2f} h_slope={r['h_slope']:.2e} "
              f"h_med={r['h_median_root']:.2e}", flush=True)

    hs = np.array([r["h_slope"] for r in rows])
    print(f"\n=== {len(rows)} points ===")
    print(f"h_slope quartiles: {np.percentile(hs, 25):.2e} / "
          f"{np.percentile(hs, 50):.2e} / {np.percentile(hs, 75):.2e}  "
          f"min {hs.min():.2e}  max {hs.max():.2e}")
    print(f"points with h_slope < 5e-4: {(hs < 5e-4).sum()}/{len(hs)}")
    V = (2 * np.pi) ** 3
    for label, hb in [("Q-tubes (h=h_slope)", hs),
                      ("first-order (h=h_slope/3)", hs / 3)]:
        for cap in (2e-3, 5e-3):
            hbc = np.minimum(hb, cap)
            N = V * np.mean(1.0 / (2 * hbc) ** 3)
            print(f"campaign integral, {label}, cap {cap:g}: "
                  f"N ~ {N:.2e} tiles")


if __name__ == "__main__":
    main()
