"""Adaptive collar driver: walks a beta3 line across strata,
auto-laddering box widths (NOTES_LP_BRIDGE 4.28-4.41).

Each position tries rungs coarse-to-fine; the first certifying
rung sets the step (2 * hf3 of that rung). Failures ladder down;
a position that exhausts the ladder is recorded FAILED and stepped
past at the finest width (loud, resumable). Warm pools chain
between adjacent tiles; a root-count change falls back to a cold
pool automatically (collar_chain machinery).

The demo walk is the hardest line in the domain: b2 = pi/3
EXACTLY, crossing the (pi/3, pi/3) corner in beta3. The driver
must downshift into the theta^2 funnel and upshift out of it
without guidance.

Rungs are (hf, hf3): generic/wall-honest (4.41), funnel, core.
"""

import json
import time

import numpy as np

from collar_chain import _warm_pool
from collar_tile import _pool_phases, collar_tile

RUNGS = [(5e-3, 2.5e-4), (5e-3, 1e-4), (1e-3, 2e-5), (1e-4, 2e-6),
         (2e-6, 2e-6)]   # last rung: the 2-D theta^2 corner core
                         # (BOTH in-face widths race as 1/theta^2
                         # on the b2 = pi/3 line — 4.42)


def drive_line(theta_lo, theta_hi, b2, b3_start, b3_end,
               ledger_path=None):
    t00 = time.time()
    b3 = b3_start
    ph_prev = None
    ledger = []
    n_tiles = n_cert = 0
    while b3 < b3_end:
        if ph_prev is None:
            ph = _pool_phases((theta_lo, b2, b3))
        else:
            ph, _n, _fell = _warm_pool((theta_lo, b2, b3), ph_prev)
        ok_rung = None
        for hf, hf3 in RUNGS:
            if collar_tile(theta_lo, theta_hi, b2, b3, hf,
                           adjacency="signed", hf3=hf3, pool=ph):
                ok_rung = (hf, hf3)
                break
        n_tiles += 1
        if ok_rung is None:
            step = 2.0 * RUNGS[-1][1]
            ledger.append((b3, step, None))
        else:
            n_cert += 1
            step = 2.0 * ok_rung[1]
            ledger.append((b3, step, ok_rung))
        b3 += step
        ph_prev = ph
    dt = time.time() - t00
    print(f"DRIVE b2={b2:.4f} b3 [{b3_start:.4f},{b3_end:.4f}] "
          f"slab [{theta_lo:g},{theta_hi:g}]: {n_cert}/{n_tiles} "
          f"tiles certified, {dt:.0f}s "
          f"({dt/max(n_tiles,1):.1f} s/tile)", flush=True)
    if ledger_path:
        with open(ledger_path, "w") as f:
            for b3v, step, rung in ledger:
                f.write(json.dumps({"b3": b3v, "step": step,
                                    "rung": rung}) + "\n")
    return n_cert, n_tiles, ledger


def main():
    p3 = float(np.pi / 3)
    # the hardest line: b2 = pi/3 exactly, crossing the corner
    drive_line(0.005, 0.015, p3, p3 - 3e-3, p3 + 3e-3,
               ledger_path="collar_drive_corner.jsonl")


if __name__ == "__main__":
    main()
