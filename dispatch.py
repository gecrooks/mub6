"""Multi-line campaign dispatcher: spawn `campaign_line` Modal calls up
to the GPU cap, collect ledgers into the local campaign ledger, resume
across crashes and disconnects.

State file (JSON): {line_key: {"spec": {...}, "fc": id|null,
"done": bool, "covered": float|null}}. Spawned call ids are persisted
BEFORE polling, so a killed dispatcher never loses a paid-for result —
rerun and it re-attaches to in-flight calls by id. The local client
being flaky is fine: spawned calls run server-side to completion.

Usage:
  dispatch.py --h 3e-4 --phi0 X --lam0 Y --n-lines 8 [--th-lo A --th-hi B]
              [--cap 8] [--state campaign_state.json]
              [--ledger campaign_ledger.jsonl]

The (phi, lam) grid matches campaign.py: lines at lam0 + k*1.8h, fixed
phi0 (v1 slab). Respect the account cap: never more than --cap (<= 10)
calls in flight.
"""

import argparse
import json
import os
import time

import modal

PI = 3.141592653589793


def load_state(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_state(path, st):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(st, f, indent=1)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=3e-4)
    ap.add_argument("--phi0", type=float, required=True)
    ap.add_argument("--lam0", type=float, required=True)
    ap.add_argument("--th-lo", type=float, default=0.0)
    ap.add_argument("--th-hi", type=float, default=PI / 2)
    ap.add_argument("--n-lines", type=int, default=8)
    ap.add_argument("--cap", type=int, default=8)
    ap.add_argument("--state", default="campaign_state.json")
    ap.add_argument("--ledger", default="campaign_ledger.jsonl")
    args = ap.parse_args()
    cap = min(args.cap, 10)                      # hard account cap

    fn = modal.Function.from_name("mub6-xcover", "campaign_line")
    st = load_state(args.state)
    spacing = 1.8 * args.h
    for k in range(args.n_lines):
        key = f"{args.phi0:.9f},{args.lam0 + k * spacing:.9f}"
        if key not in st:
            st[key] = dict(spec=dict(phi=args.phi0,
                                     lam=args.lam0 + k * spacing,
                                     th_lo=args.th_lo, th_hi=args.th_hi,
                                     h=args.h),
                           fc=None, done=False, covered=None)
    save_state(args.state, st)

    ledger = open(args.ledger, "a")
    t0 = time.time()
    while True:
        in_flight = [k for k, v in st.items() if v["fc"] and not v["done"]]
        pending = [k for k, v in st.items() if not v["fc"] and not v["done"]]
        # collect finished calls
        for k in list(in_flight):
            try:
                fc = modal.FunctionCall.from_id(st[k]["fc"])
                r = fc.get(timeout=0)
            except TimeoutError:
                continue
            except Exception as e:                # call crashed: respawn
                print(f"line {k}: call failed ({type(e).__name__}: "
                      f"{str(e)[:120]}) — will respawn", flush=True)
                st[k]["fc"] = None
                save_state(args.state, st)
                continue
            for rec in r.get("ledger", "").splitlines():
                ledger.write(rec + "\n")
            ledger.flush()
            st[k]["done"] = True
            st[k]["covered"] = r["covered"]
            full = r["covered"] >= st[k]["spec"]["th_hi"]
            save_state(args.state, st)
            print(f"line {k}: {'COMPLETE' if full else 'STALLED'} "
                  f"covered={r['covered']:.6f} n={r['n']}", flush=True)
        # top up
        in_flight = [k for k, v in st.items() if v["fc"] and not v["done"]]
        pending = [k for k, v in st.items() if not v["fc"] and not v["done"]]
        for k in pending[:max(0, cap - len(in_flight))]:
            call = fn.spawn(st[k]["spec"])
            st[k]["fc"] = call.object_id
            save_state(args.state, st)
            print(f"line {k}: spawned {call.object_id}", flush=True)
        if all(v["done"] for v in st.values()):
            break
        time.sleep(20)
    n_done = sum(v["done"] for v in st.values())
    stalled = [k for k, v in st.items()
               if v["done"] and v["covered"] < v["spec"]["th_hi"]]
    print(f"\n{n_done} lines done ({len(stalled)} stalled: {stalled}) "
          f"in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
