"""Full-D X-cover driver: 2-D Layer-3 mini-campaign over the Szollosi
deltoid region. Lines at fixed c = Im(alpha), walked in t = Re(alpha)
with the dual-mode walker (rebuild at line start and every 8th anchor,
tracked-K with per-anchor defect verification between — Result 22's
lesson); |s|-local ball slop (delta cruise ~1.5e-2). Lines are walked
center-out so partial runs leave a contiguous certified band. Per-line
results appended to xcover_lines.log."""

import sys
import time
import warnings

import numpy as np

from layer3 import certified_triple_sweep
from layer3_param import track_K
from layer3_x import triples_at, CHOICE
from mub import unbiasedness_defect
from parametric import PAD
from szollosi import in_D, szollosi_map

warnings.filterwarnings("ignore")

DELTA0 = 8e-3
DELTA_MIN = 1e-3
DELTA_MAX = 2e-2
REBUILD_EVERY = 8


def walk_line(c, delta_line=1.5e-2, log=None, delta_min=DELTA_MIN,
              n_starts=9000, rebuild_every=REBUILD_EVERY):
    # find D-extent at this c by scanning
    ts = np.arange(-1.7, 1.7, 5e-3)
    inside = [t for t in ts if in_D(t + 1j * c)]
    if not inside:
        return dict(c=c, anchors=0, gaps=0, note="line outside D")
    t_lo, t_hi = min(inside) + 2e-3, max(inside) - 2e-3
    t0c = time.time()
    t = t_lo + DELTA0
    delta = DELTA0
    covered = t_lo
    anchors = gaps = since_rebuild = 0
    gap_ts = []
    K = None
    while covered < t_hi:
        a = t + 1j * c
        if not in_D(a):
            covered = t + delta
            t = covered + delta
            continue
        need_rebuild = (K is None) or (since_rebuild >= rebuild_every)
        ok = False
        try:
            H, _ch, dH = szollosi_map(a, choice=CHOICE)
            if dH > 1e-9:
                raise RuntimeError("chart defect")
            if not need_rebuild:
                K = track_K(H, K)
                if unbiasedness_defect(H, K) > 1e-9:
                    need_rebuild = True
            if need_rebuild:
                _H, Ks, _nv = triples_at(a, n_starts=n_starts, seed=anchors)
                if Ks:
                    K = Ks[0]
                    since_rebuild = 0
                elif K is not None:
                    K = track_K(H, K)      # enumeration blind: continue
                else:
                    raise RuntimeError("no basis, no tracked K")
            if unbiasedness_defect(H, K) > 1e-9:
                raise RuntimeError("K invalid")
            eps = 5e-5
            H1, _c1, _d1 = szollosi_map(a + eps, choice=CHOICE)
            K1 = track_K(H1, K)
            rate = (np.max(np.abs(H1 - H)) + np.max(np.abs(K1 - K))) / eps
            d = delta
            while d >= delta_min:
                try:
                    n_sus, _m, _nb = certified_triple_sweep(
                        H, K, hslop=PAD * rate * d, wmin=2e-3,
                        max_boxes=2.5e7, verbose=False)
                except RuntimeError:
                    d *= 0.5
                    continue
                if n_sus == 0:
                    ok = True
                    break
                d *= 0.5
            delta = max(d, delta_min)
        except RuntimeError:
            ok = False
        if ok:
            anchors += 1
            since_rebuild += 1
            covered = t + delta
            t = covered + 0.95 * delta
            delta = min(delta * 1.4, DELTA_MAX)
        else:
            gaps += 1
            gap_ts.append(round(float(t), 5))
            K = None                      # force rebuild next anchor
            covered = t + delta_min
            t = covered + delta_min
            delta = max(delta, 4 * delta_min)
    res = dict(c=c, t_lo=t_lo, t_hi=t_hi, anchors=anchors, gaps=gaps,
               gap_ts=gap_ts,
               minutes=(time.time() - t0c) / 60)
    if log:
        with open(log, "a") as f:
            f.write(f"c={c:+.4f} t=[{t_lo:+.3f},{t_hi:+.3f}] "
                    f"anchors={anchors} gaps={gaps} "
                    f"min={res['minutes']:.1f}\n")
    return res


def main(cs):
    t0 = time.time()
    tot_a = tot_g = 0
    for i, c in enumerate(cs):
        r = walk_line(c, log="xcover_lines.log")
        tot_a += r["anchors"]
        tot_g += r.get("gaps", 0)
        print(f"[{i+1}/{len(cs)}] c={c:+.4f}: anchors={r['anchors']} "
              f"gaps={r.get('gaps', 0)} ({r.get('minutes', 0):.1f} min) "
              f"[total {(time.time()-t0)/60:.0f} min]", flush=True)
    print(f"DONE: {tot_a} anchors, {tot_g} gaps, "
          f"{(time.time()-t0)/60:.0f} min", flush=True)


if __name__ == "__main__":
    half = sys.argv[1] if len(sys.argv) > 1 else "pos"
    spacing = 2.5e-2
    if half == "pos":
        cs = list(np.arange(0.0, 1.56, spacing))
    else:
        cs = list(np.arange(-spacing, -1.56, -spacing))
    main(cs)
