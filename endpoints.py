"""B-arc endpoint completion:
1. Certified triple sweep at Bjorck's C itself (upgrading Result 5's
   multistart observation to a certificate).
2. Rebuild-mode walks over the endpoint segments [1.1965, 1.26] and
   [5.02, 5.0855], where the family chart's sqrt-branch rates diverge;
   the residual sliver to theta_min = arccos((sqrt3-1)/2) ~ 1.19606
   (and its mirror 5.08712) is reported explicitly — the endpoints
   themselves are Bjorck-equivalent and covered by the C certificate.
"""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from layer3 import certified_triple_sweep
from layer3_patch import main as patch_main
from mub import (bases_matrix, bjorck_c, find_bases, find_mu_vectors,
                 mu_vector_residuals, unbiasedness_defect)

warnings.filterwarnings("ignore")


def certify_C():
    print("=== certified triple sweep at Bjorck's C ===")
    C = bjorck_c()
    pool = find_mu_vectors([C], n_starts=20000, seed=11)
    pol = []
    for v in pool:
        th = np.angle(v * np.sqrt(6))[1:]
        th = least_squares(mu_vector_residuals, th, args=([C],), method="lm",
                          xtol=3e-16, ftol=3e-16, gtol=3e-16).x
        w = np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6)
        if not any(np.max(np.abs(w - u)) < 1e-4 for u in pol):
            pol.append(w)
    pol = np.array(pol)
    bases = find_bases(pol, tol=1e-6)
    print(f"  {len(pol)} MU vectors, {len(bases)} basis(es)")
    for b in bases:
        K = bases_matrix(pol, b)
        U_, _, Vt_ = np.linalg.svd(K)
        K = U_ @ Vt_
        print(f"  triple defect: {unbiasedness_defect(C, K):.2e}")
        t0 = time.time()
        n_sus, _m, nb = certified_triple_sweep(C, K, hslop=1e-9,
                                               wmin=2e-3, verbose=False)
        print(f"  sweep: {nb} boxes, suspects {n_sus} "
              f"[{time.time()-t0:.1f} s]")
        if n_sus == 0:
            print("  ==> CERTIFIED: no vector is MU to the Bjorck triple "
                  "{I, C, K} (1e-9 balls).")


if __name__ == "__main__":
    certify_C()
    print()
    patch_main(th_lo=1.1965, th_hi=1.26, delta0=5e-4, delta_min=2e-5)
    print()
    patch_main(th_lo=5.02, th_hi=5.0855, delta0=5e-4, delta_min=2e-5)
