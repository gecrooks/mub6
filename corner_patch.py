"""Corner-patch scheduler: 2-D snake over the theta^2 core around
(pi/3, pi/3) at the corner rung (NOTES_LP_BRIDGE 4.42, 4.45).

The corner core is first-order-limited at widths ~5e-6 with
measured curvature headroom 70x (4.45). The patch tiles an n x n
grid at hf = hf3 = 2e-6 (steps 4e-6), snaking so adjacent tiles
share warm pools. The measured s/tile here is the constant in the
corner campaign price (~5e6-1e7 tiles, embarrassingly parallel).
"""

import time

import numpy as np

from collar_chain import _warm_pool
from collar_tile import _pool_phases, collar_tile

CURV_CORNER = (10.0, 2.4e4, 3.2e7)   # 4x the measured 4.45 sups


def patch(theta_lo, theta_hi, b2_c, b3_c, n=8, hf=2e-6):
    t00 = time.time()
    n_cert = n_tiles = 0
    ph_prev = None
    for i in range(n):
        b2 = b2_c + (i - n / 2 + 0.5) * 2 * hf
        js = range(n) if i % 2 == 0 else range(n - 1, -1, -1)
        for j in js:
            b3 = b3_c + (j - n / 2 + 0.5) * 2 * hf
            if ph_prev is None:
                ph = _pool_phases((theta_lo, b2, b3))
            else:
                ph, _n, _f = _warm_pool((theta_lo, b2, b3), ph_prev)
            ok = collar_tile(theta_lo, theta_hi, b2, b3, hf,
                             adjacency="signed", hf3=hf, pool=ph,
                             curv=CURV_CORNER)
            n_tiles += 1
            n_cert += bool(ok)
            ph_prev = ph
    dt = time.time() - t00
    print(f"CORNER_PATCH {n}x{n} at ({b2_c:.5f},{b3_c:.5f}) slab "
          f"[{theta_lo:g},{theta_hi:g}]: {n_cert}/{n_tiles} "
          f"certified, {dt/n_tiles:.2f} s/tile [{dt:.0f}s]",
          flush=True)
    return n_cert, n_tiles


def main():
    p3 = float(np.pi / 3)
    patch(0.005, 0.01, p3, p3, n=8)


if __name__ == "__main__":
    main()
