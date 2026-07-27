"""Reproduce Grassl's obstruction: the pair {I, F6} admits exactly 48 MU
vectors, which form 16 orthonormal bases, no two of which are mutually
unbiased => no MU quadruple contains {I, F6}."""

import numpy as np

from mub import (bases_matrix, bjorck_c, dephase, find_bases, find_mu_vectors,
                 fourier, mu_base_pairs, unbiasedness_defect)


def main():
    F = fourier(6)
    print("Finding vectors MU to both I and F6 (biunimodular sequences)...")
    vecs = find_mu_vectors([F], n_starts=4000, seed=1)
    print(f"  distinct MU vectors found: {len(vecs)}   (Grassl: 48)")

    # verify quality
    worst = 0.0
    for v in vecs:
        worst = max(worst, np.max(np.abs(np.abs(v) ** 2 - 1 / 6)),
                    np.max(np.abs(np.abs(F.conj().T @ v) ** 2 - 1 / 6)))
    print(f"  worst deviation from unbiasedness: {worst:.2e}")

    bases = find_bases(vecs)
    print(f"  orthonormal bases among them: {len(bases)}   (Grassl: 16)")

    pairs = mu_base_pairs(vecs, bases)
    print(f"  mutually unbiased base-pairs: {len(pairs)}   (Grassl: 0)")

    # how close do the best pairs get?
    best = np.inf
    for a in range(len(bases)):
        A = bases_matrix(vecs, bases[a])
        for b in range(a + 1, len(bases)):
            B = bases_matrix(vecs, bases[b])
            best = min(best, unbiasedness_defect(A, B))
    print(f"  smallest unbiasedness defect between any two bases: {best:.4f}")

    if len(pairs) == 0 and len(bases) > 0:
        print("\n=> CONFIRMED numerically: no four MU bases contain {I, F6}.")
    return vecs, bases, pairs


if __name__ == "__main__":
    main()
