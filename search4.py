"""Global multistart search for sets of four MU bases {I, B1, B2, B3}.

Controls: d=5 and d=7 (four MU bases exist -> cost ~ 0 reachable).
Target:   d=6 (conjectured impossible -> cost floor > 0).

Usage: search4.py D N_STARTS [--fix-fourier]
"""

import sys

import numpy as np

from mub import avg_sq_distance, fourier, search_mub_set, unbiasedness_defect


def polish_unitary(M):
    """Project a near-unitary matrix onto U(d) (polar decomposition)."""
    U, _, Vh = np.linalg.svd(M)
    return U @ Vh


def main():
    d = int(sys.argv[1])
    n_starts = int(sys.argv[2])
    fix_fourier = "--fix-fourier" in sys.argv

    fixed = [fourier(d)] if fix_fourier else []
    m = 2 if fix_fourier else 3

    label = f"d={d}, {'{I,F,B2,B3}' if fix_fourier else '{I,B1,B2,B3}'}"
    print(f"Searching four MU bases: {label}, {n_starts} starts")

    best, mats, costs = search_mub_set(d, m, fixed=fixed, n_starts=n_starts,
                                       seed=42, verbose=True)
    print(f"\nbest cost (0.5*sum residuals^2): {best:.6e}")
    print(f"median cost over starts:          {np.median(costs):.6e}")
    print(f"starts within 2x of best:         {np.sum(costs < 2*best)}/{n_starts}")

    all_b = [np.eye(d)] + fixed + [polish_unitary(M) for M in mats]
    dbar = avg_sq_distance(all_b)
    worst = max(unbiasedness_defect(all_b[a], all_b[b])
                for a in range(4) for b in range(a + 1, 4))
    print(f"D-bar_4 of best solution (unitarized): {dbar:.7f}")
    print(f"worst pairwise unbiasedness defect:    {worst:.6f}")
    print("(exact four MU bases would give D-bar_4 = 1, defect = 0;")
    print(" known numerical max in d=6 is D-bar_4 = 0.9982917)")

    tag = f"d{d}" + ("_fixF" if fix_fourier else "")
    np.savez(f"best4_{tag}.npz", costs=costs,
             **{f"B{i+1}": b for i, b in enumerate(all_b[1:])})


if __name__ == "__main__":
    main()
