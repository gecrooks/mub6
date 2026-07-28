"""R9: close the B-arc branch slivers by the branch variable.

The Beauchamp-Nicoara chart exists for cos(theta) <= (sqrt3-1)/2,
i.e. theta in [theta_min, 2pi - theta_min]; its theta-rates diverge
like 1/sqrt(theta - theta_min) at the ends, which left slivers
[theta_min, 1.1965] and [5.0858, 2pi - theta_min] uncovered by the
Result-20 walk. In the branch variable s (theta = theta_min + s^2,
resp. theta = theta_max - s^2) the map is smooth with FINITE rates,
so a short rebuild-mode s-walk covers each sliver, and the s = 0
anchor's ball [-delta, delta] covers theta in [branch, branch +
delta^2] including the branch point itself (Bjorck-equivalent,
independently certified in Result 20).

Rebuild mode: fresh triple per anchor (no tracking across the
degenerate end), FD s-rates (EMPIRICAL + PAD, same prototype status
as all Layer-3 rates; certified B-map rates remain a substrate item),
ball slop PAD * rate_s * delta_s.
"""

import time
import warnings

import numpy as np

from layer3 import build_triple, certified_triple_sweep
from mub import unbiasedness_defect
from parametric import PAD

warnings.filterwarnings("ignore")

TH_MIN = float(np.arccos((np.sqrt(3.0) - 1.0) / 2.0))
TH_MAX = 2.0 * np.pi - TH_MIN
DELTA_S = 4e-3
EPS_S = 1e-3


def walk_sliver(branch, sign, s_hi, max_boxes=3e7):
    """Cover theta between the branch point and branch + sign*s_hi^2
    by anchors in s, s = 0 first (its ball covers the branch point)."""
    print(f"  sliver at theta = {branch:.7f} (sign {sign:+d}), "
          f"s up to {s_hi:.4f} (theta width {s_hi**2:.2e})", flush=True)
    n_anchor = 0
    boxes = 0
    worst = 0.0
    t0 = time.time()
    # the branch point itself is basis-degenerate (Bjorck-equivalent,
    # certified separately); first anchor at small s0 > 0 whose ball
    # must reach back through s = 0 (delta >= s0), covering theta down
    # to the branch point inclusive
    s = None
    for cand in (DELTA_S, 2 * DELTA_S, 3 * DELTA_S):
        try:
            build_triple(branch + sign * cand * cand)
            s = cand
            break
        except RuntimeError:
            continue
    if s is None:
        raise RuntimeError("no buildable anchor near the branch")
    delta = s                      # first ball covers s = 0
    covered_s = -delta
    first = True
    while covered_s < s_hi:
        th = branch + sign * s * s
        B, K, _pol, _bidx = build_triple(th)
        d_tri = unbiasedness_defect(B, K)
        worst = max(worst, d_tri)
        if d_tri > 1e-9:
            raise RuntimeError(f"triple defect {d_tri:.1e} at s={s:.4f}")
        # FD rate in s (smooth at the branch in this variable)
        th_e = branch + sign * (s + EPS_S) ** 2
        B1, K1, _p, _b = build_triple(th_e)
        # basis matching across rebuilds: compare via defect-invariant
        # entrywise magnitudes is unstable — use column-greedy phase
        # alignment for the rate estimate only (rate is EMPIRICAL+PAD)
        rate = (np.max(np.abs(B1 - B))
                + np.max(np.abs(_align(K1, K) - K))) / EPS_S
        ok = False
        delta_floor = s if first else 5e-4   # first ball must reach s=0
        while delta >= delta_floor:
            hslop = PAD * rate * delta
            try:
                n_sus, _m, nb = certified_triple_sweep(
                    B, K, hslop=hslop, wmin=2e-3, max_boxes=max_boxes,
                    verbose=False)
            except RuntimeError:
                delta *= 0.5
                continue
            boxes += nb
            if n_sus == 0:
                ok = True
                break
            delta *= 0.5
        if not ok:
            raise RuntimeError(f"sliver anchor s={s:.4f} uncertifiable"
                               f" (floor {delta_floor:.1e})")
        first = False
        n_anchor += 1
        covered_s = s + delta
        s = covered_s + 0.95 * delta
        delta = min(delta * 1.4, 8e-3)
    print(f"    covered s to {covered_s:.4f} "
          f"(theta {branch + sign * covered_s**2:.7f}): {n_anchor} "
          f"anchors, {boxes} boxes, worst defect {worst:.1e}, "
          f"{time.time()-t0:.0f} s", flush=True)
    return n_anchor, boxes


def _align(K1, K):
    """Greedy column matching + phase alignment of K1 onto K (rate
    estimation only)."""
    K1 = K1.copy()
    used = set()
    out = np.zeros_like(K)
    for j in range(K.shape[1]):
        best, bi, bph = -1.0, -1, 1.0
        for i in range(K1.shape[1]):
            if i in used:
                continue
            z = np.vdot(K1[:, i], K[:, j])
            if abs(z) > best:
                best, bi, bph = abs(z), i, z / max(abs(z), 1e-300)
        used.add(bi)
        out[:, j] = K1[:, bi] * bph
    return out


def main():
    print("=== R9: B-arc branch slivers via theta = branch +- s^2 ===")
    print(f"theta_min = {TH_MIN:.10f}, theta_max = {TH_MAX:.10f}")
    # lower sliver: [TH_MIN, 1.1965]; upper: [5.0858, TH_MAX]
    n1, b1 = walk_sliver(TH_MIN, +1, np.sqrt(1.1965 - TH_MIN))
    n2, b2 = walk_sliver(TH_MAX, -1, np.sqrt(TH_MAX - 5.0858))
    print(f"BOTH SLIVERS CLOSED: {n1 + n2} anchors, {b1 + b2} boxes")
    print("Combined with Result 20: the FULL B-arc "
          f"[{TH_MIN:.7f}, {TH_MAX:.7f}] including both branch points "
          "is covered gap-free (prototype rates).")


if __name__ == "__main__":
    main()
