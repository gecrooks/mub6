"""Parametric Layer 3: certified unextendability of the X-family triple
{I, B(theta), K(theta)} over a theta-INTERVAL.

The pointwise sweep (layer3.py) showed margin 1/6. Here the sweep carries
an extra additive slop = (rate_B + rate_K) * delta, where the rates are
sampled (PAD-padded) drift rates of the two Hadamards along theta:
B directly from the closed form, K by tracking its polished columns.
A certificate at slop_theta covers every theta in [t0-delta, t0+delta].
"""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from certify import L_H_G
from layer3 import build_triple, certified_triple_sweep
from mub import beauchamp_nicoara, mu_vector_residuals, unbiasedness_defect
from parametric import PAD

warnings.filterwarnings("ignore")


def track_K(B1, K0):
    """Re-polish K's columns as MU vectors of the new B."""
    cols = []
    for j in range(6):
        v = K0[:, j]
        ph0 = np.angle(v / v[0])[1:]
        r = least_squares(mu_vector_residuals, ph0, args=([B1],),
                          method="lm", xtol=3e-16, ftol=3e-16, gtol=3e-16)
        cols.append(np.exp(1j * np.concatenate(([0.0], r.x))) / np.sqrt(6))
    K1 = np.array(cols).T
    # restore the global column phases lost by dephasing: align to K0
    for j in range(6):
        z = np.vdot(K1[:, j], K0[:, j])
        K1[:, j] *= z / abs(z)
    return K1


def main():
    t0c = 1.6
    print(f"=== parametric Layer 3 around theta = {t0c} ===")
    B0, K0, pol, bidx = build_triple(t0c)
    U_, _, Vt_ = np.linalg.svd(K0)
    K0 = U_ @ Vt_

    # sampled drift rates along theta
    eps = 1e-4
    B1 = beauchamp_nicoara(t0c + eps)
    K1 = track_K(B1, K0)
    rate_B = float(np.max(np.abs(B1 - B0))) / eps
    rate_K = float(np.max(np.abs(K1 - K0))) / eps
    print(f"sampled drift rates: |dB/dth| {rate_B:.3f}  |dK/dth| {rate_K:.3f}")
    print(f"triple defect at th+eps: {unbiasedness_defect(B1, K1):.2e}")

    for delta in (0.02, 0.012, 0.006, 0.003):
        hslop = PAD * (rate_B + rate_K) * delta
        print(f"delta = {delta}: ball {hslop:.4f}, "
              f"sweep slop {L_H_G * hslop:.4f} (vs 1/6 margin)", flush=True)
        t0 = time.time()
        try:
            n_sus, _mm, total = certified_triple_sweep(
                B0, K0, hslop=hslop, wmin=2e-3, max_boxes=4e7, verbose=False)
        except RuntimeError as e:
            print(f"  drowned ({e}) [{time.time()-t0:.0f} s]", flush=True)
            continue
        dt = time.time() - t0
        print(f"  boxes {total}, suspects {n_sus} [{dt:.0f} s]", flush=True)
        if n_sus == 0:
            print(f"==> THEOREM (prototype-grade): for every theta in "
                  f"[{t0c-delta}, {t0c+delta}], no vector is MU to the "
                  f"X-family triple {{I, B(theta), K(theta)}} -- "
                  f"one anchor covers a theta-interval of length "
                  f"{2*delta:g}; the whole family needs ~{np.pi/delta:.0f} "
                  f"anchors at this rate.")
            break


if __name__ == "__main__":
    main()
