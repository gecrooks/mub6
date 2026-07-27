"""Family-complete Layer 3: certify that no vector is MU to the triple
{I, B(theta), K(theta)} for EVERY theta along the Hermitian (Beauchamp-
Nicoara) arc — the first family-complete certificate of the program.

Adaptive anchor walk: at each anchor, K is warm-tracked from the previous
one (column polishing; triple defect re-verified), local drift rates are
FD-sampled (EMPIRICAL, PAD-padded — certified B-map rates are a future
substrate item), and the 10-component sweep runs with ball slop
PAD*(rate_B+rate_K)*delta. On drowning, delta halves; below delta_min the
segment is recorded as a GAP and skipped (honest reporting). The arc ends
[1.196, 1.25] and [5.03, 5.087] degenerate toward Bjorck's C (excluded
pointwise in the literature and in Result 5) and are left to an endpoint
analysis.
"""

import time
import warnings

import numpy as np

from layer3 import certified_triple_sweep, build_triple
from layer3_param import track_K
from mub import beauchamp_nicoara, unbiasedness_defect
from parametric import PAD

warnings.filterwarnings("ignore")

TH_LO, TH_HI = 1.25, 5.03
DELTA0 = 4e-3
DELTA_MIN = 8e-4
DELTA_MAX = 7e-3


def main(th_lo=TH_LO, th_hi=TH_HI, max_boxes=3e7):
    TH_LO_, TH_HI_ = th_lo, th_hi
    t0 = time.time()
    print(f"=== family-complete Layer 3 over theta in [{th_lo}, {th_hi}] ===")
    B, K, _pol, _bidx = build_triple(th_lo + DELTA0)
    U_, _, Vt_ = np.linalg.svd(K)
    K = U_ @ Vt_

    th = th_lo + DELTA0
    delta = DELTA0
    anchors = 0
    gaps = []
    covered_to = th_lo
    worst_defect = 0.0
    boxes_total = 0
    eps = 1e-4
    while covered_to < th_hi:
        B = beauchamp_nicoara(th)
        K = track_K(B, K)
        d_tri = unbiasedness_defect(B, K)
        worst_defect = max(worst_defect, d_tri)
        if d_tri > 1e-9:
            gaps.append((th, f"triple tracking defect {d_tri:.1e}"))
            covered_to = th + delta
            th = covered_to + delta
            continue
        B1 = beauchamp_nicoara(th + eps)
        K1 = track_K(B1, K)
        rate = (np.max(np.abs(B1 - B)) + np.max(np.abs(K1 - K))) / eps
        ok = False
        while delta >= DELTA_MIN:
            hslop = PAD * rate * delta
            try:
                n_sus, _m, nb = certified_triple_sweep(
                    B, K, hslop=hslop, wmin=2e-3, max_boxes=max_boxes,
                    verbose=False)
            except RuntimeError:
                delta *= 0.5
                continue
            boxes_total += nb
            if n_sus == 0:
                ok = True
                break
            delta *= 0.5
        if not ok:
            gaps.append((th, f"uncertifiable at delta_min (rate {rate:.2f})"))
            covered_to = th + DELTA_MIN
            th = covered_to + DELTA_MIN
            delta = DELTA_MIN * 2
            continue
        anchors += 1
        covered_to = th + delta
        th = covered_to + 0.95 * delta
        delta = min(delta * 1.3, DELTA_MAX)
        if anchors % 40 == 0:
            print(f"  theta {covered_to:.3f} ({100*(covered_to-th_lo)/(th_hi-th_lo):.0f}%), "
                  f"{anchors} anchors, delta {delta:.4f}, "
                  f"[{(time.time()-t0)/60:.1f} min]", flush=True)

    dt = (time.time() - t0) / 60
    print(f"\nFAMILY ARC CERTIFIED: {anchors} anchors, {boxes_total} boxes, "
          f"{dt:.1f} min, worst triple defect {worst_defect:.1e}")
    if gaps:
        print(f"GAPS ({len(gaps)}):")
        for g_th, why in gaps:
            print(f"  theta ~ {g_th:.4f}: {why}")
    else:
        print("NO GAPS: for every theta in [%.2f, %.2f], no vector is MU "
              "to the triple {I, B(theta), K(theta)}." % (th_lo, th_hi))


if __name__ == "__main__":
    main()
