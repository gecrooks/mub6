"""Layer 3 over the Dita arc: for every x in the fundamental domain
[-1/8, 1/8] of the affine Dita family, certify that EVERY MU triple
{I, D(x), K} is strongly unextendible. Rebuild-mode anchors (the family
crosses the maximally-rich D(0) with 120 vectors / 10 bases, where any
tracking would lose branches); at anchors where no basis exists the
segment is triple-free and Layer 3 is vacuous (recorded)."""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from layer3 import certified_triple_sweep
from layer3_param import track_K
from mub import (bases_matrix, dita, find_bases, find_mu_vectors,
                 mu_vector_residuals, unbiasedness_defect)
from parametric import PAD

warnings.filterwarnings("ignore")


def triples_at(x, n_starts=9000, seed=0):
    D = dita(x)
    pool = find_mu_vectors([D], n_starts=n_starts, seed=seed)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([D],), method="lm",
                          xtol=3e-16, ftol=3e-16, gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-4 for u in pol):
            pol.append(w)
    pol = np.array(pol)
    bases = find_bases(pol, tol=1e-6)
    Ks = []
    for b in bases:
        K = bases_matrix(pol, b)
        U_, _, Vt_ = np.linalg.svd(K)
        Ks.append(U_ @ Vt_)
    return D, Ks, len(pol)


def main(x_lo=-0.125, x_hi=0.125, delta0=1.5e-3, delta_min=2e-4):
    t0 = time.time()
    print(f"=== Layer 3 over the Dita arc x in [{x_lo}, {x_hi}] ===")
    x = x_lo + delta0
    delta = delta0
    covered_to = x_lo
    anchors = 0
    gaps = []
    bare = []           # triple-free anchor segments (Layer 3 vacuous)
    eps = 5e-5
    max_bases = 0
    while covered_to < x_hi:
        D, Ks, nvec = triples_at(x, seed=anchors)
        max_bases = max(max_bases, len(Ks))
        if not Ks:
            bare.append((x - delta, x + delta))
            covered_to = x + delta
            x = covered_to + 0.95 * delta
            delta = min(delta * 1.3, 3e-3)
            continue
        ok_all = True
        for K in Ks:
            dtri = unbiasedness_defect(D, K)
            if dtri > 1e-9:
                ok_all = False
                gaps.append((x, f"triple defect {dtri:.1e}"))
                break
            D1 = dita(x + eps)
            K1 = track_K(D1, K)
            rate = (np.max(np.abs(D1 - D)) + np.max(np.abs(K1 - K))) / eps
            got = False
            d = delta
            while d >= delta_min:
                try:
                    n_sus, _m, _nb = certified_triple_sweep(
                        D, K, hslop=PAD * rate * d, wmin=2e-3,
                        max_boxes=1.2e7, verbose=False)
                except RuntimeError:
                    d *= 0.5
                    continue
                if n_sus == 0:
                    got = True
                    break
                d *= 0.5
            if not got:
                ok_all = False
                gaps.append((x, f"uncertifiable (rate {rate:.1f}, "
                                f"{len(Ks)} bases)"))
                break
            delta = d
        if ok_all:
            anchors += 1
            covered_to = x + delta
            x = covered_to + 0.95 * delta
            delta = min(delta * 1.3, 3e-3)
            if anchors % 15 == 0:
                print(f"  x {covered_to:+.4f} "
                      f"({100*(covered_to-x_lo)/(x_hi-x_lo):.0f}%), "
                      f"nvec={nvec} nbases={len(Ks)}, delta={delta:.2e} "
                      f"[{(time.time()-t0)/60:.1f} min]", flush=True)
        else:
            covered_to = x + delta_min
            x = covered_to + delta_min
            delta = max(delta, 4 * delta_min)
    print(f"\nDITA ARC DONE: {anchors} triple-bearing anchors certified, "
          f"max {max_bases} bases/anchor, {len(gaps)} gaps, "
          f"{(time.time()-t0)/60:.1f} min")
    if bare:
        lo = min(a for a, b in bare)
        hi = max(b for a, b in bare)
        print(f"triple-free anchor range: [{lo:+.4f}, {hi:+.4f}] "
              f"({len(bare)} anchors, Layer 3 vacuous there)")
    for g_x, why in gaps:
        print(f"  GAP x ~ {g_x:+.4f}: {why}")
    if not gaps:
        print("NO GAPS: every MU triple {I, D(x), K} on the arc is "
              "strongly unextendible.")


if __name__ == "__main__":
    main()
