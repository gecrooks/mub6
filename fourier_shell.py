"""Fourier-shell prototype (gap R1): quadruple-exclusion AT the F-locus.

Near the Fourier family, MU triples exist (16 bases through {I, F6}), so
the campaign certificate must exclude QUADRUPLES: among the certified-
complete MU vectors of {I, H}, every pair of orthonormal bases fails
mutual unbiasedness by a certified margin. At F6 the margin is the exact
"1/6 wall" (Grassl); this prototype

  1. re-runs the certified pointwise pipeline at F6,
  2. enumerates ALL 6-cliques of the certified possible-edge graph
     (candidate bases -- a certified superset of actual bases),
  3. certifies for every clique pair a lower bound on the max deviation
     of cross overlaps from 1/6 (the quadruple obstruction), and
  4. measures the drift rate of that margin along the Fourier family,
     giving the shell thickness the campaign can use.
"""

import time
import warnings

import numpy as np
from scipy.optimize import least_squares

from certify import certified_graph, sweep, verify_all, _uvec
from mub import fourier, fourier_family, find_mu_vectors, mu_vector_residuals

warnings.filterwarnings("ignore")


def all_cliques(n, edges, size=6):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    out = []

    def extend(cl, cand):
        if len(cl) == size:
            out.append(tuple(cl))
            return
        if len(cl) + len(cand) < size:
            return
        for v in sorted(cand):
            extend(cl + [v], cand & adj[v] & set(range(v + 1, n)))

    extend([], set(range(n)))
    return out


def clique_pair_margin(enclosures, c1, c2):
    """Certified lower bound on max_ij | |<u_i, v_j>|^2 - 1/6 | for the two
    candidate bases, over the enclosure boxes."""
    cs = np.array([enclosures[i][0] for i in c1])
    rs = np.array([enclosures[i][1] for i in c1])
    ds = np.array([enclosures[j][0] for j in c2])
    qs = np.array([enclosures[j][1] for j in c2])
    U1 = _uvec(cs)
    U2 = _uvec(ds)
    O = np.abs(U1.conj() @ U2.T)                      # moduli
    dev = np.abs(O ** 2 - 1.0 / 6.0)
    # overlap drift: |d<u,v>| <= (sum r_u + sum r_v)/6; deviation drift
    # |d(o^2)| <= 2 o drift + drift^2, o <= 1
    drift = (rs.sum(axis=1)[:, None] + qs.sum(axis=1)[None, :]) / 6.0
    dev_lo = dev - 2.0 * drift - drift ** 2
    return float(dev_lo.max())


def main():
    print("=== Fourier-shell prototype at F6 ===")
    F = fourier(6)
    t0 = time.time()
    sus_C, sus_W, stats = sweep(F, hslop=1e-11, wmin=1e-4, verbose=True)
    enclosures, ok = verify_all(F, sus_C, sus_W, hslop=1e-11)
    if not ok:
        raise SystemExit("verification failed")
    edges, margin = certified_graph(enclosures)
    print(f"  {len(enclosures)} roots, {len(edges)} edges "
          f"[{time.time()-t0:.0f} s]")

    cliques = all_cliques(len(enclosures), edges)
    print(f"  6-cliques of the certified graph: {len(cliques)} "
          f"(Grassl: 16 bases)")

    worst = np.inf
    for a in range(len(cliques)):
        for b in range(a + 1, len(cliques)):
            m = clique_pair_margin(enclosures, cliques[a], cliques[b])
            worst = min(worst, m)
    print(f"  min over clique pairs of certified max-deviation from 1/6: "
          f"{worst:.6f}   (the 1/6 wall: 0.166667)")
    if worst > 0:
        print("  ==> QUADRUPLE-EXCLUSION CERTIFIED at F6: no two MU bases "
              "through {I, F6} are mutually unbiased.")

    # shell thickness: drift rate of the wall along the Fourier family
    print("drift of the wall along F(a,b):")
    base_vecs = None
    for eps in (0.0, 1e-3, 3e-3, 1e-2):
        H = fourier_family(eps, eps / 2)
        pool = find_mu_vectors([H], n_starts=3000, seed=3)
        pol = []
        for v in pool:
            th = np.angle(v * np.sqrt(6))[1:]
            th = least_squares(mu_vector_residuals, th, args=([H],),
                               method="lm", xtol=3e-16, ftol=3e-16,
                               gtol=3e-16).x
            pol.append(np.exp(1j * np.concatenate(([0.0], th))) / np.sqrt(6))
        pol = np.array(pol)
        from mub import find_bases, bases_matrix, unbiasedness_defect
        bases = find_bases(pol, tol=1e-6)
        best = np.inf
        for a in range(len(bases)):
            A = bases_matrix(pol, bases[a])
            for b in range(a + 1, len(bases)):
                Bb = bases_matrix(pol, bases[b])
                best = min(best, unbiasedness_defect(A, Bb))
        print(f"  eps={eps:6g}: vecs={len(pol):3d} bases={len(bases):3d} "
              f"min pair-defect={best:.6f}")


if __name__ == "__main__":
    main()
