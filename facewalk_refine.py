"""Refinement pass over the face-walk ledger: rerun failed tiles
at the safe proxy theta = 3e-3 (wall/corner tiles break the
1e-6-proxy enumeration via the 1/theta wall racing — NOTES 4.58).
Appends results to facewalk_ledger3.jsonl."""

import json
import time
from concurrent.futures import ProcessPoolExecutor

HF = 0.05
TUBE = 0.005


def tile_job(args):
    b2, b3 = args
    import facewalk
    from facewalk import face_tile
    facewalk.TH_PROXY = 3e-3
    t0 = time.time()
    try:
        ok, nb, m = face_tile(b2, b3, HF, th_tube=TUBE)
    except Exception as e:
        return {"b2": b2, "b3": b3, "ok": False,
                "err": str(e)[:60], "dt": time.time() - t0}
    return {"b2": b2, "b3": b3, "ok": bool(ok), "bases": int(nb),
            "margin": float(m) if m == m else None,
            "proxy": 3e-3, "dt": time.time() - t0}


def main():
    fails = []
    for line in open("facewalk_ledger2_pass1.jsonl"):
        d = json.loads(line)
        if not d["ok"]:
            fails.append((d["b2"], d["b3"]))
    print(f"refining {len(fails)} failed tiles at proxy 3e-3",
          flush=True)
    n_ok = n_bad = 0
    t0 = time.time()
    with open("facewalk_ledger3.jsonl", "w") as led:
        with ProcessPoolExecutor(max_workers=7) as ex:
            for r in ex.map(tile_job, fails, chunksize=4):
                led.write(json.dumps(r) + "\n")
                n_ok += r["ok"]
                n_bad += not r["ok"]
    print(f"REFINE DONE: {n_ok}/{len(fails)} now certified, "
          f"{n_bad} still failing [{time.time()-t0:.0f}s]",
          flush=True)
    if n_bad:
        for line in open("facewalk_ledger3.jsonl"):
            d = json.loads(line)
            if not d["ok"]:
                print("  ", {k: d[k] for k in d
                             if k in ("b2", "b3", "bases", "err")},
                      flush=True)


if __name__ == "__main__":
    main()
