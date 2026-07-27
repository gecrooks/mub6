"""Gap patch for the B-arc Layer-3 walk: the zone around theta = 3pi/2
(a Dita-equivalent special point) breaks K-tracking (the partner basis
bifurcates — Dita's matrix carries 10 bases) and has diverging drift
rates. Patch mode: REBUILD the triple fresh at every anchor (no
tracking), certify EVERY basis found (statement: no vector is MU to ANY
triple {I, B(theta), K} at these theta), delta floor 2e-4."""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from layer3 import certified_triple_sweep
from layer3_param import track_K
from mub import (bases_matrix, beauchamp_nicoara, find_bases,
                 find_mu_vectors, mu_vector_residuals, unbiasedness_defect)
from parametric import PAD

warnings.filterwarnings("ignore")


def triples_at(theta, n_starts=8000, seed=0):
    B = beauchamp_nicoara(theta)
    pool = find_mu_vectors([B], n_starts=n_starts, seed=seed)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([B],), method="lm",
                          xtol=3e-16, ftol=3e-16, gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-4 for u in pol):
            pol.append(w)
    pol = np.array(pol)
    bases = find_bases(pol, tol=1e-6)
    Ks = []
    for b in bases:
        K = bases_matrix(pol, b)
        U_, _, Vt_ = np.linalg.svd(K)
        Ks.append(U_ @ Vt_)
    return B, Ks, len(pol)


def main(th_lo=4.75, th_hi=4.92, delta0=1.5e-3, delta_min=2e-4):
    t0 = time.time()
    print(f"=== Layer-3 gap patch (rebuild mode) over [{th_lo}, {th_hi}] ===")
    th = th_lo + delta0
    delta = delta0
    covered_to = th_lo
    anchors = 0
    gaps = []
    eps = 5e-5
    while covered_to < th_hi:
        B, Ks, nvec = triples_at(th, seed=anchors)
        if not Ks:
            gaps.append((th, "no basis found"))
            covered_to = th + delta_min
            th = covered_to + delta_min
            continue
        ok_all = True
        rate_max = 0.0
        for K in Ks:
            dtri = unbiasedness_defect(B, K)
            if dtri > 1e-9:
                ok_all = False
                gaps.append((th, f"triple defect {dtri:.1e}"))
                break
            B1 = beauchamp_nicoara(th + eps)
            K1 = track_K(B1, K)
            rate = (np.max(np.abs(B1 - B)) + np.max(np.abs(K1 - K))) / eps
            rate_max = max(rate_max, rate)
            got = False
            d = delta
            while d >= delta_min:
                try:
                    n_sus, _m, _nb = certified_triple_sweep(
                        B, K, hslop=PAD * rate * d, wmin=2e-3,
                        max_boxes=1.2e7, verbose=False)
                except RuntimeError:
                    d *= 0.5
                    continue
                if n_sus == 0:
                    got = True
                    break
                d *= 0.5
            if not got:
                ok_all = False
                gaps.append((th, f"uncertifiable (rate {rate:.1f})"))
                break
            delta = d
        if ok_all:
            anchors += 1
            covered_to = th + delta
            th = covered_to + 0.95 * delta
            delta = min(delta * 1.3, 2e-3)
            print(f"  theta {covered_to:.4f}  nvec={nvec} nbases={len(Ks)} "
                  f"rate={rate_max:.1f} delta={delta:.2e} "
                  f"[{(time.time()-t0)/60:.1f} min]", flush=True)
        else:
            covered_to = th + delta_min
            th = covered_to + delta_min
            delta = max(delta, 4 * delta_min)
    print(f"\nPATCH DONE: {anchors} anchors, {len(gaps)} gaps, "
          f"{(time.time()-t0)/60:.1f} min")
    for g_th, why in gaps:
        print(f"  GAP theta ~ {g_th:.4f}: {why}")


if __name__ == "__main__":
    main()
