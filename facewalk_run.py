"""Full face walk — the first fully executed piece of the domain
decomposition (PROOF_SKELETON §5 face piece, demo grade).

Tiles the full (b2, b3) face on a step-0.1 grid (hf = 0.05,
overlapping closed boxes), theta-tube 0.005 to overlap the collar's
deepest slab. Runs in a worker pool; every tile's verdict goes to a
JSONL ledger. Special handling: the corner-adjacent tiles evaluate
at the safe proxy theta = 3e-3 (float map error there ~2e-9;
NOTES 4.53) — the standard proxy 1e-6 is fine everywhere else
(off-corner map error <= 9e-10).

Failures are recorded, not fatal: unexpected basis counts mark
SPLIT tiles for the refinement pass (the corner tile is expected
to need the safe proxy; anything else failing is news).
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HF = 0.05
TUBE = 0.005


def tile_job(args):
    b2, b3 = args
    import numpy as np
    from facewalk import face_tile, TH_PROXY
    import facewalk
    p3 = np.pi / 3
    # corner-adjacent tiles: safe proxy (map noise at 1e-6 there)
    near_corner = min(abs(b2 - p3), abs(b2 - 5 * p3)) < 0.12 and \
        min(abs(b3 - p3), abs(b3 - 5 * p3)) < 0.12
    if near_corner:
        facewalk.TH_PROXY = 3e-3
    else:
        facewalk.TH_PROXY = 1e-6
    t0 = time.time()
    try:
        ok, nb, m = face_tile(b2, b3, HF, th_tube=TUBE)
    except Exception as e:
        return {"b2": b2, "b3": b3, "ok": False, "err": str(e)[:60],
                "dt": time.time() - t0}
    return {"b2": b2, "b3": b3, "ok": bool(ok), "bases": int(nb),
            "margin": float(m) if m == m else None,
            "proxy": facewalk.TH_PROXY, "dt": time.time() - t0}


def main():
    grid = [(round(float(b2), 4), round(float(b3), 4))
            for b2 in np.arange(0.05, 2 * np.pi, 0.1)
            for b3 in np.arange(0.05, 2 * np.pi, 0.1)]
    print(f"face walk: {len(grid)} tiles, hf {HF}, tube {TUBE}",
          flush=True)
    t0 = time.time()
    n_ok = n_bad = 0
    with open("facewalk_ledger.jsonl", "w") as led:
        with ProcessPoolExecutor(max_workers=7) as ex:
            for i, r in enumerate(ex.map(tile_job, grid,
                                         chunksize=8)):
                led.write(json.dumps(r) + "\n")
                n_ok += r["ok"]
                n_bad += not r["ok"]
                if (i + 1) % 200 == 0:
                    led.flush()
                    print(f"  [{i+1}/{len(grid)}] ok {n_ok} "
                          f"bad {n_bad} "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"FACE WALK DONE: {n_ok}/{len(grid)} certified, "
          f"{n_bad} failed, {time.time()-t0:.0f}s total",
          flush=True)


if __name__ == "__main__":
    main()
