"""Layer 3 in the generic X6^(2) interior: spot certificates at sampled
alpha, then a rebuild-mode chain line across the deltoid region D — the
2-D walk demonstrator. Full 2-D coverage of D is campaign-scale
(~4e5 anchors at these rates); tonight fixes the machinery and the
per-anchor cost."""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from layer3 import certified_triple_sweep
from layer3_param import track_K
from mub import (bases_matrix, find_bases, find_mu_vectors,
                 mu_vector_residuals, unbiasedness_defect)
from parametric import PAD
from szollosi import in_D, szollosi_map

warnings.filterwarnings("ignore")

CHOICE = ((0, 1), (0, 1))


def triples_at(alpha, n_starts=9000, seed=0):
    H, _ch, d = szollosi_map(alpha, choice=CHOICE)
    if d > 1e-9:
        raise RuntimeError(f"not Hadamard at alpha={alpha} (defect {d:.1e})")
    pool = find_mu_vectors([H], n_starts=n_starts, seed=seed)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([H],), method="lm",
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
    return H, Ks, len(pol)


def spot_certs():
    print("=== X-interior spot certificates ===")
    rng = np.random.default_rng(9)
    done = 0
    while done < 8:
        a = rng.uniform(-1.2, 1.2) + 1j * rng.uniform(-1.2, 1.2)
        if not in_D(a):
            continue
        try:
            H, Ks, nvec = triples_at(a, seed=done)
        except RuntimeError as e:
            print(f"  alpha {a:.3f}: {e}")
            continue
        stat = []
        for K in Ks:
            dtri = unbiasedness_defect(H, K)
            n_sus, _m, _nb = certified_triple_sweep(
                H, K, hslop=1e-9, wmin=2e-3, max_boxes=1.2e7, verbose=False)
            stat.append((dtri, n_sus))
        ok = all(s[1] == 0 and s[0] < 1e-9 for s in stat)
        print(f"  alpha {a.real:+.3f}{a.imag:+.3f}i: nvec={nvec} "
              f"bases={len(Ks)} {'ALL CERTIFIED' if ok and Ks else 'CHECK'}",
              flush=True)
        done += 1


def line_walk(c=0.15, t_lo=-0.7, t_hi=0.7, delta0=1.5e-3, delta_min=2e-4):
    print(f"=== chain line across D: alpha = t + {c}i, "
          f"t in [{t_lo}, {t_hi}] ===")
    t0c = time.time()
    t = t_lo + delta0
    delta = delta0
    covered = t_lo
    anchors = 0
    gaps = []
    eps = 5e-5
    while covered < t_hi:
        a = t + 1j * c
        if not in_D(a):
            gaps.append((t, "outside D"))
            covered = t + delta
            t = covered + delta
            continue
        try:
            H, Ks, nvec = triples_at(a, seed=anchors)
        except RuntimeError as e:
            gaps.append((t, str(e)))
            covered = t + delta_min
            t = covered + delta_min
            continue
        if not Ks:
            gaps.append((t, f"no basis ({nvec} vecs)"))
            covered = t + delta_min
            t = covered + delta_min
            continue
        ok_all = True
        for K in Ks:
            if unbiasedness_defect(H, K) > 1e-9:
                ok_all = False
                gaps.append((t, "triple defect"))
                break
            H1, _c1, _d1 = szollosi_map(a + eps, choice=CHOICE)
            K1 = track_K(H1, K)
            rate = (np.max(np.abs(H1 - H)) + np.max(np.abs(K1 - K))) / eps
            got = False
            d = delta
            while d >= delta_min:
                try:
                    n_sus, _m, _nb = certified_triple_sweep(
                        H, K, hslop=PAD * rate * d, wmin=2e-3,
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
                gaps.append((t, f"uncertifiable (rate {rate:.1f})"))
                break
            delta = d
        if ok_all:
            anchors += 1
            covered = t + delta
            t = covered + 0.95 * delta
            delta = min(delta * 1.3, 3e-3)
            if anchors % 20 == 0:
                print(f"  t {covered:+.3f} "
                      f"({100*(covered-t_lo)/(t_hi-t_lo):.0f}%), "
                      f"nvec={nvec} nbases={len(Ks)} delta={delta:.2e} "
                      f"[{(time.time()-t0c)/60:.1f} min]", flush=True)
        else:
            covered = t + delta_min
            t = covered + delta_min
            delta = max(delta, 4 * delta_min)
    print(f"\nX LINE DONE: {anchors} anchors, {len(gaps)} gaps, "
          f"{(time.time()-t0c)/60:.1f} min")
    for g_t, why in gaps[:12]:
        print(f"  GAP t ~ {g_t:+.4f}: {why}")


if __name__ == "__main__":
    spot_certs()
    print()
    line_walk()
