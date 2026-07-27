"""Mini-campaign driver: certify a 3-dimensional parameter BLOCK by
margin-cached chains — anchored tiles chained along theta, a grid of
chain-lines across (phi, lam), re-anchoring on failure or cache
exhaustion. The certified region is the union rectangle of the tile
grid (overlapping by construction: chain step 1.6h < 2h, line spacing
1.8h < 2h)."""

import time
import warnings

import numpy as np

from cache import anchored_tile, chain_step

warnings.filterwarnings("ignore")


def certify_block(beta0, h=4e-4, n_steps=6, n_phi=4, n_lam=4,
                  step_fac=1.6, line_fac=1.8, use_certified=False):
    step = step_fac * h
    spacing = line_fac * h
    t0 = time.time()
    tiles = 0
    anchors = 0
    for ip in range(n_phi):
        for il in range(n_lam):
            base = (beta0[0], beta0[1] + ip * spacing,
                    beta0[2] + il * spacing)
            print(f"  line (phi+{ip}, lam+{il}):", flush=True)
            state = anchored_tile(base, h, verbose=True,
                                  use_certified=use_certified)
            anchors += 1
            tiles += 1
            k = 1
            while k <= n_steps:
                target = (base[0] + k * step, base[1], base[2])
                res = chain_step(state, target, verbose=True)
                if res.get("ok"):
                    tiles += 1
                    k += 1
                    continue
                print(f"    re-anchoring at step {k} "
                      f"({res.get('reason')})", flush=True)
                state = anchored_tile(target, h, verbose=True,
                                      use_certified=use_certified)
                anchors += 1
                tiles += 1
                k += 1
    dt = time.time() - t0
    ext_t = n_steps * step + 2 * h
    ext_p = (n_phi - 1) * spacing + 2 * h
    ext_l = (n_lam - 1) * spacing + 2 * h
    print(f"\nBLOCK CERTIFIED (prototype): {tiles} tiles, {anchors} anchors,"
          f" {dt/60:.1f} min")
    print(f"certified box: theta [{beta0[0]-h:.6f}, "
          f"{beta0[0]-h+ext_t:.6f}] x phi [{beta0[1]-h:.6f}, "
          f"{beta0[1]-h+ext_p:.6f}] x lam [{beta0[2]-h:.6f}, "
          f"{beta0[2]-h+ext_l:.6f}]")
    print(f"extents {ext_t:.2e} x {ext_p:.2e} x {ext_l:.2e} "
          f"= volume {ext_t*ext_p*ext_l:.2e} "
          f"({ext_t*ext_p*ext_l/(2*h)**3:.0f}x one tile)")
    return dict(tiles=tiles, anchors=anchors, minutes=dt / 60)


if __name__ == "__main__":
    import sys
    cert = "--certified" in sys.argv
    beta0 = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    print(f"=== mini-campaign block at K6{np.round(beta0, 4)}, h=4e-4"
          f"{' [CERTIFIED TAXES]' if cert else ''} ===")
    certify_block(beta0, use_certified=cert)
