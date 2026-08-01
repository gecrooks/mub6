"""Resume facewalk_run from its ledger (spawn-safe driver file)."""

import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from facewalk_run import tile_job


def main():
    done = set()
    for line in open("facewalk_ledger.jsonl"):
        d = json.loads(line)
        done.add((d["b2"], d["b3"]))
    grid = [(round(float(b2), 4), round(float(b3), 4))
            for b2 in np.arange(0.05, 2 * np.pi, 0.1)
            for b3 in np.arange(0.05, 2 * np.pi, 0.1)]
    todo = [g for g in grid if g not in done]
    print(f"resume: {len(done)} done, {len(todo)} to go",
          flush=True)
    n_ok = n_bad = 0
    t0 = time.time()
    with open("facewalk_ledger.jsonl", "a") as led:
        with ProcessPoolExecutor(max_workers=7) as ex:
            for r in ex.map(tile_job, todo, chunksize=4):
                led.write(json.dumps(r) + "\n")
                n_ok += r["ok"]
                n_bad += not r["ok"]
    print(f"RESUME DONE: +{n_ok} ok +{n_bad} bad "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
