"""Array-module-generic certified triple sweep (NumPy or CuPy).

Rigor note: the sweep's certificate rests on the static rounding-error
lemma (Result 15), which requires only IEEE-754-compliant FP64 add/mul —
true on NVIDIA GPUs — so the CuPy path certifies to the same standard
as the NumPy path. The |s|-local ball slop matches layer3.py.
"""

import numpy as _np

from certify import SLOP


def certified_triple_sweep_xp(B, K, hslop=1e-9, wmin=2e-3, chunk=1_000_000,
                              max_boxes=4e8, xp=None, verbose=False):
    if xp is None:
        xp = _np
    stacks = [xp.asarray(B.conj()), xp.asarray(K.conj())]
    SQ6 = 2.0 * float(_np.sqrt(6.0))
    inv_sqrt6 = 1.0 / float(_np.sqrt(6.0))
    stack_C = [xp.full((1, 5), _np.pi)]
    stack_W = [xp.full((1, 5), _np.pi)]
    n_sus = 0
    total = 0
    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        total += int(C.shape[0])
        if total > max_boxes:
            raise RuntimeError(f"exceeded {max_boxes:g} boxes")
        u = xp.empty((C.shape[0], 6), dtype=xp.complex128)
        u[:, 0] = inv_sqrt6
        u[:, 1:] = xp.exp(1j * C) * inv_sqrt6
        sw = W.sum(axis=1)
        best = xp.full(C.shape[0], -_np.inf)
        for Hc in stacks:
            s = u @ Hc
            g = xp.abs(s) ** 2 - 1.0 / 6.0
            smod_w = xp.minimum(xp.abs(s) + sw[:, None] / 6.0, 1.0)
            L = 2.0 * smod_w / 6.0
            marg = xp.abs(g) - L * sw[:, None] - (SLOP + SQ6 * smod_w * hslop)
            best = xp.maximum(best, marg.max(axis=1))
        keep = best <= 0
        Ck, Wk = C[keep], W[keep]
        small = Wk.max(axis=1) <= wmin
        n_sus += int(small.sum())
        Cb, Wb = Ck[~small], Wk[~small]
        if Cb.shape[0]:
            j = xp.argmax(Wb, axis=1)
            rows = xp.arange(Cb.shape[0])
            Wn = Wb.copy()
            Wn[rows, j] = Wb[rows, j] / 2.0
            Cl, Cr = Cb.copy(), Cb.copy()
            Cl[rows, j] -= Wn[rows, j]
            Cr[rows, j] += Wn[rows, j]
            stack_C.append(xp.concatenate([Cl, Cr]))
            stack_W.append(xp.concatenate([Wn, Wn]))
    return n_sus, 0.0, total
