"""Computational lab for the MUB problem in dimension 6.

Conventions: bases are 6x6 complex matrices whose COLUMNS are the basis
vectors. A complex Hadamard matrix here is unitary with all |H_jk| = 1/sqrt(d)
(i.e. already normalized), so {I, H} is a pair of mutually unbiased bases.

Explicit order-6 Hadamard families follow Bengtsson-Bruzda-Ericsson-Larsson-
Tadej-Zyczkowski, quant-ph/0610161, and McNulty-Weigert, arXiv:2410.23997.
"""

import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(20260724)

# ----------------------------------------------------------------------------
# Hadamard matrix constructors (all return unitary matrices, entries 1/sqrt(6))
# ----------------------------------------------------------------------------

Q = np.exp(2j * np.pi / 6)          # primitive 6th root
W3 = np.exp(2j * np.pi / 3)         # primitive 3rd root
BJORCK_D = (1 - np.sqrt(3)) / 2 + 1j * np.sqrt(np.sqrt(3) / 2)  # d^2-(1-r3)d+1=0


def fourier(d=6):
    j, k = np.meshgrid(np.arange(d), np.arange(d), indexing="ij")
    return np.exp(2j * np.pi * j * k / d) / np.sqrt(d)


def fourier_family(x1, x2):
    """Two-parameter affine Fourier family F(x1,x2); F(0,0) = F6."""
    z1 = np.exp(2j * np.pi * x1)
    z2 = np.exp(2j * np.pi * x2)
    q = Q
    rows = [
        [1, 1, 1, 1, 1, 1],
        [1, q * z1, q**2 * z2, q**3, q**4 * z1, q**5 * z2],
        [1, q**2, q**4, 1, q**2, q**4],
        [1, q**3 * z1, z2, q**3, z1, q**3 * z2],
        [1, q**4, q**2, 1, q**4, q**2],
        [1, q**5 * z1, q**4 * z2, q**3, q**2 * z1, q * z2],
    ]
    return np.array(rows, dtype=complex) / np.sqrt(6)


def fourier_family_T(x1, x2):
    return fourier_family(x1, x2).T.copy()


def bjorck_c():
    """Bjorck's circulant Hadamard C6: row j is the first row shifted by j."""
    d = BJORCK_D
    r = np.array([1, 1j * d, -d, -1j, -d.conjugate(), 1j * d.conjugate()])
    C = np.empty((6, 6), dtype=complex)
    for j in range(6):
        C[j] = np.roll(r, j)
    return C / np.sqrt(6)


def dita(x):
    """One-parameter affine Dita family D(x); D(0) is Dita's matrix."""
    z = np.exp(2j * np.pi * x)
    zb = z.conjugate()
    i = 1j
    rows = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, i, -i, -i, i],
        [1, i, -1, i * z, -i * z, -i],
        [1, -i, i * zb, -1, i, -i * zb],
        [1, -i, -i * zb, i, -1, i * zb],
        [1, i, -i, -i * z, i * z, -1],
    ]
    return np.array(rows, dtype=complex) / np.sqrt(6)


def tao_s6():
    w = W3
    rows = [
        [1, 1, 1, 1, 1, 1],
        [1, 1, w, w, w**2, w**2],
        [1, w, 1, w**2, w**2, w],
        [1, w, w**2, 1, w, w**2],
        [1, w**2, w**2, w, 1, w],
        [1, w**2, w, w**2, w, 1],
    ]
    return np.array(rows, dtype=complex) / np.sqrt(6)


def beauchamp_nicoara(theta):
    """Hermitian (self-adjoint) family B(theta); needs cos(theta)<=(sqrt3-1)/2.

    Endpoints of the allowed range give (permutations of) Bjorck's C;
    theta = pi gives Dita's matrix.
    """
    y = np.exp(1j * theta)
    z = (1 + 2 * y - y**2) / (y * (-1 + 2 * y + y**2))
    disc = 1 + 2 * y + 2 * y**3 + y**4
    x = (1 + 2 * y + y**2 - np.sqrt(2) * np.sqrt(disc)) / (1 + 2 * y - y**2)
    t = x * y * z
    xb, yb, zb, tb = (v.conjugate() for v in (x, y, z, t))
    rows = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, -xb, -y, y, xb],
        [1, -x, 1, y, zb, -tb],
        [1, -yb, yb, -1, -tb, tb],
        [1, yb, z, -t, 1, -xb],
        [1, x, -t, t, -x, -1],
    ]
    return np.array(rows, dtype=complex) / np.sqrt(6)


# ----------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------

def is_unitary(U, tol=1e-9):
    return np.max(np.abs(U.conj().T @ U - np.eye(U.shape[0]))) < tol


def is_hadamard(H, tol=1e-9):
    d = H.shape[0]
    return is_unitary(H, tol) and np.max(np.abs(np.abs(H) - 1 / np.sqrt(d))) < tol


def unbiasedness_defect(A, B):
    """max_ij | |<a_i|b_j>|^2 - 1/d |  (0 iff bases A,B are MU)."""
    d = A.shape[0]
    G = np.abs(A.conj().T @ B) ** 2
    return np.max(np.abs(G - 1 / d))


def dephase(H):
    """Equivalence-normal form: make first row and column real positive."""
    H = H.copy()
    H = H / (H[:, 0:1] / np.abs(H[:, 0:1]))     # rows: kill phase of col 0
    H = H / (H[0:1, :] / np.abs(H[0:1, :]))     # cols: kill phase of row 0
    return H


def defect(H, tol=1e-6):
    """Defect of a Hadamard matrix (Tadej-Zyczkowski): dimension of the
    solution space of the linearized unitarity constraints for entrywise
    phase perturbations, minus the trivial 2d-1 dephasing directions.
    H_jk -> H_jk exp(i eps_jk): unitarity to first order requires, for j<k,
    sum_m H_jm conj(H_km) i(eps_jm - eps_km) = 0.  Count solution dim - (2d-1).
    """
    d = H.shape[0]
    rows = []
    for j in range(d):
        for k in range(j + 1, d):
            # complex equation -> 2 real rows in the d*d unknowns eps
            coef = np.zeros((d, d), dtype=complex)
            for m in range(d):
                c = 1j * H[j, m] * np.conjugate(H[k, m])
                coef[j, m] += c
                coef[k, m] -= c
            rows.append(coef.real.ravel())
            rows.append(coef.imag.ravel())
    Mst = np.array(rows)
    s = np.linalg.svd(Mst, compute_uv=False)
    null_dim = d * d - int(np.sum(s > tol * s[0]))
    return null_dim - (2 * d - 1)


# ----------------------------------------------------------------------------
# MU vectors: unit vectors unbiased to I and to a list of Hadamard bases
# ----------------------------------------------------------------------------

def _vec_from_phases(phases):
    v = np.exp(1j * np.concatenate(([0.0], phases))) / np.sqrt(6)
    return v


def mu_vector_residuals(phases, mats):
    v = _vec_from_phases(phases)
    res = []
    for H in mats:
        res.append(np.abs(H.conj().T @ v) ** 2 - 1 / 6)
    return np.concatenate(res)


def find_mu_vectors(mats, n_starts=2000, tol=1e-20, dedup_tol=1e-6, seed=None):
    """Find unit vectors MU to the identity basis and to every basis in mats.

    Vector ansatz: entries exp(i theta_j)/sqrt(6), theta_0 = 0 (global phase).
    Returns array of distinct solutions (each a length-6 complex vector).
    """
    local_rng = np.random.default_rng(seed) if seed is not None else rng
    sols = []
    for _ in range(n_starts):
        th0 = local_rng.uniform(0, 2 * np.pi, size=5)
        r = least_squares(mu_vector_residuals, th0, args=(mats,), method="lm",
                          xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if r.cost < tol:
            v = _vec_from_phases(r.x)
            for w in sols:
                if np.max(np.abs(v - w)) < dedup_tol:
                    break
            else:
                sols.append(v)
    return np.array(sols)


def orthogonal_pairs(vecs, tol=1e-8):
    n = len(vecs)
    G = np.abs(vecs.conj() @ vecs.T)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if G[a, b] < tol:
                adj[a].add(b)
                adj[b].add(a)
    return adj


def find_bases(vecs, tol=1e-8):
    """All orthonormal bases (6-cliques of the orthogonality graph)."""
    adj = orthogonal_pairs(vecs, tol)
    n = len(vecs)
    bases = []

    def extend(clique, candidates):
        if len(clique) == 6:
            bases.append(tuple(clique))
            return
        for v in sorted(candidates):
            extend(clique + [v], candidates & adj[v] & set(range(v + 1, n)))

    extend([], set(range(n)))
    return bases


def bases_matrix(vecs, base):
    return np.array([vecs[i] for i in base]).T


def mu_base_pairs(vecs, bases):
    """Which pairs of bases (as index tuples) are mutually unbiased?"""
    out = []
    for a in range(len(bases)):
        A = bases_matrix(vecs, bases[a])
        for b in range(a + 1, len(bases)):
            B = bases_matrix(vecs, bases[b])
            if unbiasedness_defect(A, B) < 1e-8:
                out.append((a, b))
    return out


# ----------------------------------------------------------------------------
# Random complex Hadamard matrices (generic points of the Hadamard variety)
# ----------------------------------------------------------------------------

def _hadamard_residuals(phases, d):
    E = np.exp(1j * phases.reshape(d - 1, d - 1))
    H = np.ones((d, d), dtype=complex)
    H[1:, 1:] = E
    G = H.conj().T @ H / d - np.eye(d)
    iu = np.triu_indices(d, 1)
    return np.concatenate([G[iu].real, G[iu].imag])


def random_hadamard(d=6, seed=None, max_tries=50):
    """Least-squares solve for a dephased complex Hadamard from a random start.

    Generic solutions for d=6 should be generic points of the Hadamard variety
    (conjecturally the four-parameter family G6^(4), defect 4).
    """
    local_rng = np.random.default_rng(seed) if seed is not None else rng
    for _ in range(max_tries):
        p0 = local_rng.uniform(0, 2 * np.pi, size=(d - 1) ** 2)
        r = least_squares(_hadamard_residuals, p0, args=(d,), method="lm",
                          xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if r.cost < 1e-22:
            E = np.exp(1j * r.x.reshape(d - 1, d - 1))
            H = np.ones((d, d), dtype=complex)
            H[1:, 1:] = E
            return H / np.sqrt(d)
    raise RuntimeError("no Hadamard found")


# ----------------------------------------------------------------------------
# Global search for sets of MU bases: {I, B1, ..., Bm} in dimension d
# ----------------------------------------------------------------------------

def _phase_param_sizes(d, m, first_dephased):
    """Every basis MU to I is a Hadamard: entries exp(i th)/sqrt(d).
    Column phases (row 0 of th) are free gauge -> fixed to 0 for all bases.
    The shared row-diagonal gauge dephases column 0 of the FIRST free basis
    (only when no fixed Hadamard already used that gauge)."""
    sizes = []
    for i in range(m):
        if i == 0 and first_dephased:
            sizes.append((d - 1) ** 2)
        else:
            sizes.append(d * (d - 1))
    return sizes


def _unpack_phases(params, d, m, first_dephased):
    mats, off = [], 0
    for i, n in enumerate(_phase_param_sizes(d, m, first_dephased)):
        chunk = params[off:off + n]
        off += n
        th = np.zeros((d, d))
        if i == 0 and first_dephased:
            th[1:, 1:] = chunk.reshape(d - 1, d - 1)
        else:
            th[1:, :] = chunk.reshape(d - 1, d)
        mats.append(np.exp(1j * th) / np.sqrt(d))
    return mats


def _mub_residuals(params, d, m, fixed, first_dephased):
    mats = _unpack_phases(params, d, m, first_dephased)
    res = []
    iu = np.triu_indices(d, 1)
    for M in mats:                       # unitarity (diagonal automatic)
        G = (M.conj().T @ M)[iu]
        res.append(G.real)
        res.append(G.imag)
    hs = fixed + mats                    # unbiasedness among Hadamard bases
    for a in range(len(hs)):
        for b in range(a + 1, len(hs)):
            P = np.abs(hs[a].conj().T @ hs[b]) ** 2 - 1 / d
            res.append(P.ravel())
    return np.concatenate(res)


def search_mub_set(d, m, fixed=(), n_starts=30, seed=None, verbose=False,
                   max_nfev=40000):
    """Try to find m free Hadamard bases forming an MU set {I, *fixed, B1..Bm}.

    Unbiasedness to I is exact by construction; the returned cost is
    0.5 * sum of squared residuals (unitarity + pairwise unbiasedness),
    which vanishes iff an exact MU set exists at the found point.
    """
    local_rng = np.random.default_rng(seed) if seed is not None else rng
    fixed = list(fixed)
    first_dephased = len(fixed) == 0
    nparams = sum(_phase_param_sizes(d, m, first_dephased))
    best, best_mats, costs = np.inf, None, []
    for s in range(n_starts):
        p0 = local_rng.uniform(0, 2 * np.pi, size=nparams)
        r = least_squares(_mub_residuals, p0, args=(d, m, fixed, first_dephased),
                          method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                          max_nfev=max_nfev)
        costs.append(r.cost)
        if r.cost < best:
            best = r.cost
            best_mats = _unpack_phases(r.x, d, m, first_dephased)
        if verbose:
            print(f"  start {s}: cost {r.cost:.6e}   (best {best:.6e})")
    return best, best_mats, np.array(costs)


def avg_sq_distance(bases):
    """Mean squared distance D-bar over all pairs of bases,
    D_{b,b'} = 1 - (1/(d-1)) sum_{ij} (|<b_i|c_j>|^2 - 1/d)^2.
    Equals 1 iff all pairs are MU; 0 for identical bases.
    (Bengtsson et al. / Raynal-Lu-Englert measure; known numerical max for
    four bases in d=6 is 0.9982917.)"""
    d = bases[0].shape[0]
    tot, cnt = 0.0, 0
    for a in range(len(bases)):
        for b in range(a + 1, len(bases)):
            P = np.abs(bases[a].conj().T @ bases[b]) ** 2
            tot += 1 - np.sum((P - 1 / d) ** 2) / (d - 1)
            cnt += 1
    return tot / cnt
