"""Collar-tile chain (amortization measurement for the collar
campaign).

The collar tile's dominant cost is the fresh MU-root enumeration
(~5 s of ~6 s at 6000 starts). Chaining along beta2 warm-starts
each tile's pool from the previous tile's roots (re-polish, ~0.2 s)
plus a small top-up enumeration; a lost root or changed count
fails LOUDLY into a full re-enumeration (the certified pass anchors
enumeration trust in the fat-sweep coverage, as in the bulk).

Usage: chain(theta_lo, theta_hi, b2_0, b3, hf, n_steps, ...)
steps beta2 by 2*hf per tile (adjacent boxes tile the beta2 line).
"""

import time

import numpy as np

from collar_tile import (_polish, _pool_phases, _repolish_pool,
                         collar_tile)
from karlsson import karlsson_map
from mub import find_mu_vectors


def _warm_pool(beta, ph_prev, topup=500, seed=11):
    """Re-polish the previous roots at the new beta and top up with
    a short enumeration; returns (phases, n_new, full_fallback)."""
    H = karlsson_map(*beta)
    ph = []
    for t in _repolish_pool(beta, ph_prev):
        if not any(np.max(np.abs(np.exp(1j * t) - np.exp(1j * u)))
                   < 1e-5 for u in ph):
            ph.append(t)
    n_kept = len(ph)
    for v in find_mu_vectors([H], n_starts=topup, seed=seed):
        t = _polish(H, np.angle(v * np.sqrt(6))[1:])
        if not any(np.max(np.abs(np.exp(1j * t) - np.exp(1j * u)))
                   < 1e-5 for u in ph):
            ph.append(t)
    if n_kept < len(ph_prev) or len(ph) != len(ph_prev):
        return _pool_phases(beta), len(ph) - n_kept, True
    return np.array(ph), len(ph) - n_kept, False


def chain(theta_lo, theta_hi, b2_0, b3, hf, n_steps,
          adjacency="signed", hf3=None, curv=None):
    t00 = time.time()
    n_ok = n_fall = 0
    ph_prev = None
    for k in range(n_steps):
        b2 = b2_0 + 2.0 * hf * k
        t0 = time.time()
        if ph_prev is None:
            ph = _pool_phases((theta_lo, b2, b3))
            fell = False
        else:
            ph, _new, fell = _warm_pool((theta_lo, b2, b3), ph_prev)
        n_fall += fell
        kw = {} if curv is None else {"curv": curv}
        ok = collar_tile(theta_lo, theta_hi, b2, b3, hf,
                         adjacency=adjacency, hf3=hf3,
                         pool=ph, **kw)
        n_ok += bool(ok)
        ph_prev = ph
    dt = time.time() - t00
    print(f"COLLAR_CHAIN b3={b3:.4f} slab [{theta_lo:g},{theta_hi:g}]"
          f" {n_steps} steps: {n_ok} certified, {n_fall} pool "
          f"fallbacks, {dt/n_steps:.2f} s/tile amortized "
          f"[{dt:.0f}s]", flush=True)
    return n_ok


def main():
    # generic line and wall line, 20 tiles each
    chain(0.005, 0.015, 1.0, 2.0, hf=5e-3, n_steps=20, hf3=2.5e-4)
    chain(0.005, 0.015, 1.0, float(np.pi / 3), hf=5e-3, n_steps=20,
          hf3=2.5e-4)


if __name__ == "__main__":
    main()
