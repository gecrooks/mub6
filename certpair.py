"""Certified signed pair rates — the certified-grade version of the
collar tile's v3 primitive (NOTES_LP_BRIDGE 4.34).

The demo tile taxes pair overlaps with SAMPLED rates (FD S from
root_data2). This module produces enclosures instead:

  1. root enclosures via krawczyk_verify (existing, certified);
  2. J at the enclosure center with a Hessian-row pad over the ball
     (HESS_ROW = 11/18, rigorous);
  3. dg/dbeta from dual_karlsson (certified dH/dbeta partials),
     with a theta-ball pad (prototype grade: Lipschitz sample x
     PAD);
  4. S enclosure by midpoint-inverse residual bound:
       S = S0 + E,  |E|_inf <= kappa (|dJ| |S0| + w_dgdb)
       kappa = |J0^-1|_inf / (1 - |J0^-1|_inf |dJ|_inf);
  5. pair-rate enclosure d<u_i,u_j>/dbeta_l from S enclosures and
     the root balls.

Output: certified interval vs the FD signed rate, at a generic
breaking edge and at a wall breaking edge — the fatness ratio is
the feasibility number for the certified collar pass.
"""

import numpy as np

from certify import HESS_ROW, _g_and_J, krawczyk_verify
from collar_tile import _pool_phases, _uvecs
from dual import dual_karlsson
from interval import IV
from karlsson import karlsson_map
from mub import find_bases
from parametric import root_data2

PAD_PROTO = 2.0     # pad on the (sampled) theta-Lipschitz of dgdb


def _dgdb_enclosure(beta, c, r):
    """(5,3) center + width of dg/dbeta over the root ball c +- r.
    Center and beta-partials certified via dual_karlsson at the
    point beta; ball dependence padded (prototype grade)."""
    Hd = dual_karlsson(IV(beta[0]), IV(beta[1]), IV(beta[2]))
    u = np.exp(1j * np.concatenate(([0.0], c))) / np.sqrt(6.0)
    H0 = karlsson_map(*beta)
    s = u @ H0.conj()
    dgdb = np.zeros((5, 3))
    wid = np.zeros((5, 3))
    for l in range(3):
        dH = np.array([[complex(Hd[a][b].d[l].re.mid,
                                Hd[a][b].d[l].im.mid)
                        for b in range(6)] for a in range(6)])
        dHw = np.array([[max(Hd[a][b].d[l].re.width,
                             Hd[a][b].d[l].im.width)
                         for b in range(6)] for a in range(6)])
        sb = u @ dH.conj()
        dgdb[:, l] = (2.0 * np.real(np.conj(s) * sb))[1:6]
        # interval width of the dual partials + ball pad (sampled
        # Lipschitz of dgdb in theta ~ 2|dH| entrywise, PADded)
        lip = 2.0 * np.abs(dH).max() * r.sum()
        wid[:, l] = 2.0 * np.abs(s).max() * dHw.max() \
            + PAD_PROTO * lip
    return dgdb, wid


def certified_S(beta, th, r0=1e-4):
    """Enclosure of the root sensitivity S = -J^-1 dg/dbeta.
    Returns (S0, S_err_inf, c, r) or raises."""
    H = karlsson_map(*beta)
    ok, c, r, _R = krawczyk_verify(H, th.copy(),
                                   np.full(5, r0))[:4]
    if not ok:
        raise RuntimeError("krawczyk failed")
    _g, J0 = _g_and_J(H, c)
    Jinv = np.linalg.inv(J0)
    nJinv = np.abs(Jinv).sum(axis=1).max()
    dJ = HESS_ROW * r.sum()               # entry/row pad over ball
    if nJinv * dJ * 5 >= 1.0:
        raise RuntimeError(f"J too ill-conditioned: "
                           f"|Jinv| {nJinv:.1f} dJ {dJ:.2e}")
    dgdb, wid = _dgdb_enclosure(beta, c, r)
    S0 = -np.linalg.solve(J0, dgdb)
    kappa = nJinv / (1.0 - nJinv * dJ * 5)
    nS0 = np.abs(S0).sum(axis=0).max()
    err = kappa * (dJ * 5 * nS0 + wid.sum(axis=0).max())
    return S0, err, c, r


def pair_rate(beta, th_i, th_j, r0=1e-4):
    """Certified enclosure of |d<u_i,u_j>/dbeta_l| for l=0,1,2,
    plus the FD reference. Returns (rate0, rate_err, fd)."""
    Si, ei, ci, ri = certified_S(beta, th_i, r0)
    Sj, ej, cj, rj = certified_S(beta, th_j, r0)
    ui = np.exp(1j * np.concatenate(([0.0], ci))) / np.sqrt(6.0)
    uj = np.exp(1j * np.concatenate(([0.0], cj))) / np.sqrt(6.0)
    rate0 = np.empty(3)
    for l in range(3):
        d = 1j * (Sj[:, l] - Si[:, l])
        rate0[l] = abs((np.conj(ui[1:]) * uj[1:] * d).sum()) / 6.0
    # error: S errors + ball phase pads (|u| exact, phase in ball)
    smag = max(np.abs(Si).max(), np.abs(Sj).max())
    err = (ei + ej) / 6.0 * 5 + (ri.sum() + rj.sum()) * smag / 6.0
    # FD reference via root_data2
    fd = None
    try:
        Sfi = None
        Sfj = None
        for de in (1e-5, 2e-6):
            try:
                Sfi = root_data2(beta, th_i, delta=de)[0]
                Sfj = root_data2(beta, th_j, delta=de)[0]
                break
            except Exception:
                continue
        if Sfi is not None:
            fd = np.empty(3)
            for l in range(3):
                d = 1j * (Sfj[:, l] - Sfi[:, l])
                fd[l] = abs((np.conj(ui[1:]) * uj[1:]
                             * d).sum()) / 6.0
    except Exception:
        pass
    return rate0, err, fd


def main():
    for tag, beta in (("generic", (0.005, 1.0, 2.0)),
                      ("wall", (0.005, 1.1, float(np.pi / 3)))):
        ph = _pool_phases(beta, n_starts=4000)
        U = _uvecs(ph)
        bases = find_bases(U, tol=5e-2)
        bs = bases[0]
        B = np.array([U[k] for k in bs])
        G = np.abs(B.conj() @ B.T)
        iu = np.triu_indices(6, 1)
        k = int(np.argmax(G[iu]))
        a, b = bs[iu[0][k]], bs[iu[1][k]]
        try:
            r0, err, fd = pair_rate(beta, ph[a], ph[b])
            fds = (np.array2string(fd, precision=4)
                   if fd is not None else "n/a")
            print(f"{tag}: breaking edge o={G[iu][k]:.5f}\n"
                  f"  certified rate {np.array2string(r0, precision=4)}"
                  f" +- {err:.2e}\n  FD reference   {fds}", flush=True)
        except RuntimeError as e:
            print(f"{tag}: {e}", flush=True)


if __name__ == "__main__":
    main()
