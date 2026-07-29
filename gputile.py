"""Array-module-generic zoned tile sweep (NumPy or CuPy) — the GPU port
of `parametric.zoned_sweep`, box-exact: same splitting order, same taxes,
same guard/oracle collection, so xp=numpy reproduces the CPU sweep's
(stuck, D0, total) exactly and the CuPy path certifies to the same
standard (IEEE-754 FP64, static rounding lemma of Result 15).

Oracles and caches cross the device boundary: oracle polylines are
uploaded per call (small), cache/stuck outputs are returned as host
float32 arrays (the chain machinery consumes them on CPU).
"""

import numpy as _np

from certify import SLOP
from trig_kernel import kexp_i


def _to_host(xp, a):
    return a.get() if hasattr(a, "get") else _np.asarray(a)


def zoned_sweep_xp(H0, roots, coef0, coef1, guards, far_tax,
                   zone_R=0.85, wmin=1e-4, chunk=1_000_000, max_boxes=6e7,
                   init_C=None, init_W=None, s_drift=0.0, oracles=None,
                   beta_rate=None, cache=None, beta_unit=0.0,
                   stuck_out=None, xp=None, fo=None):
    if xp is None:
        xp = _np
    fo_d = None
    if fo is not None:
        fo_d = dict(dH0c=[xp.asarray(a) for a in fo["dH0c"]],
                    WD=[xp.asarray(w) for w in fo["WD"]],
                    hv=[float(x) for x in fo["hv"]],
                    s_drift=float(fo["s_drift"]))
    R = xp.asarray(_np.array(roots))
    n = int(R.shape[0])
    Hc = xp.asarray(H0.conj())
    inv6 = 1.0 / float(_np.sqrt(6.0))
    coef0 = xp.asarray(_np.asarray(coef0, float))
    coef1 = xp.asarray(_np.asarray(coef1, float))
    guards_d = xp.asarray(_np.asarray(guards, float))
    beta_rate_d = None if beta_rate is None \
        else xp.asarray(_np.asarray(beta_rate, float))
    bu_vec = _np.ndim(beta_unit) > 0
    beta_unit_d = xp.asarray(_np.asarray(beta_unit, float)) if bu_vec \
        else float(beta_unit)
    orc = {}
    if oracles is not None:
        for i, o in oracles.items():
            if o is None:
                continue
            olist = o if isinstance(o, list) else [o]
            orc[i] = [dict(th0=xp.asarray(q["th0"]),
                           w=xp.asarray(q["w"]), Wc=xp.asarray(q["Wc"]),
                           tgrid=xp.asarray(q["tgrid"]),
                           Yc=xp.asarray(q["Yc"]),
                           rho_y=float(q["rho_y"]), T=float(q["T"]))
                      for q in olist]
    if init_C is None:
        stack_C = [xp.full((1, 5), _np.pi)]
        stack_W = [xp.full((1, 5), _np.pi)]
    else:
        stack_C = [xp.asarray(_np.asarray(init_C, float).reshape(-1, 5))]
        stack_W = [xp.asarray(_np.asarray(init_W, float).reshape(-1, 5))]
    D0 = xp.zeros((n, 5))
    stuck = 0
    total = 0
    TWO_PI = 2.0 * _np.pi

    def torus(d):
        return (d + _np.pi) % TWO_PI - _np.pi

    while stack_C:
        C = stack_C.pop()
        W = stack_W.pop()
        if len(C) > chunk:
            stack_C.append(C[chunk:])
            stack_W.append(W[chunk:])
            C, W = C[:chunk], W[:chunk]
        m = int(C.shape[0])
        total += m
        if total > max_boxes:
            raise RuntimeError(f"zoned sweep exceeded {max_boxes:g} boxes")
        u = xp.empty((m, 6), dtype=xp.complex128)
        u[:, 0] = inv6
        u[:, 1:] = kexp_i(C, xp) * inv6      # certified kernel (E_TRIG)
        s = u @ Hc
        g = xp.abs(s) ** 2 - 1.0 / 6.0
        sw = W.sum(axis=1)
        delta = xp.abs(torus(C[:, None, :] - R[None, :, :]))
        dist = delta.max(axis=2)
        reach = dist + sw[:, None]
        taxes = coef0[None, :] + coef1[None, :] * reach
        taxes = xp.where(dist <= zone_R, taxes, xp.inf)
        tax_box = xp.minimum(far_tax, taxes.min(axis=1))
        smod_w = xp.abs(s) + sw[:, None] / 6.0 + s_drift
        L = 2.0 * xp.minimum(smod_w, 1.0) / 6.0
        margin = xp.abs(g) - L * sw[:, None]
        if beta_rate_d is not None:
            beta_tax = beta_rate_d * xp.minimum(smod_w, 1.0)
            if fo_d is not None:
                smod1 = xp.minimum(smod_w, 1.0)
                t1 = xp.zeros_like(beta_tax)
                for j in range(3):
                    sb = u @ fo_d["dH0c"][j]
                    d0g = 2.0 * xp.real(xp.conj(s) * sb)
                    t1 += fo_d["hv"][j] * (xp.abs(d0g)
                                           + 2.0 * fo_d["s_drift"]
                                           * xp.abs(sb)
                                           + 2.0 * smod1
                                           * fo_d["WD"][j])
                beta_tax = xp.minimum(beta_tax, t1 + SLOP)
            tax_arr = xp.minimum(beta_tax, tax_box[:, None])
        else:
            beta_tax = None
            tax_arr = tax_box[:, None]
        excl = (margin > SLOP + tax_arr).any(axis=1)
        if cache is not None and bool(excl.any()):
            if beta_tax is not None:
                exb = margin - SLOP - beta_tax
                kb = xp.argmax(exb, axis=1)
                rows_i = xp.arange(m)
                eb = exb[rows_i, kb]
                cacheable = excl & (eb > 0)
            else:
                cacheable = xp.zeros(m, dtype=bool)
            if bool(cacheable.any()):
                sm = xp.minimum(smod_w[rows_i, kb], 1.0)
                bu = (beta_unit_d[kb[cacheable]] if bu_vec
                      else beta_unit_d)
                cache["C"].append(_to_host(
                    xp, C[cacheable]).astype(_np.float32))
                cache["W"].append(_to_host(
                    xp, W[cacheable]).astype(_np.float32))
                cache["E"].append(_to_host(
                    xp, eb[cacheable]).astype(_np.float32))
                cache["R"].append(_to_host(
                    xp, bu * sm[cacheable]).astype(_np.float32))
            rest = excl & ~cacheable
            if bool(rest.any()):
                cache["SC"].append(_to_host(xp, C[rest]).astype(_np.float32))
                cache["SW"].append(_to_host(xp, W[rest]).astype(_np.float32))
        keep = ~excl
        C, W, delta = C[keep], W[keep], delta[keep]
        if len(C):
            inside_hull = (delta + W[:, None, :]
                           <= guards_d[None, :, :]).all(axis=2)
            collected = xp.zeros(len(C), dtype=bool)
            for i in range(n):
                sel = inside_hull[:, i]
                if not bool(sel.any()):
                    continue
                if i in orc:
                    fine = xp.zeros(int(sel.sum()), dtype=bool)
                    for o in orc[i]:
                        d = torus(C[sel] - o["th0"])
                        tb = d @ o["w"]
                        yb = d @ o["Wc"]
                        tW = W[sel] @ xp.abs(o["w"])
                        yW = W[sel] @ xp.abs(o["Wc"])
                        yc = xp.stack(
                            [xp.interp(tb, o["tgrid"], o["Yc"][:, j])
                             for j in range(4)], axis=1)
                        fine = fine | (
                            (xp.abs(tb) + tW <= o["T"]) &
                            ((xp.abs(yb - yc) + yW
                              <= o["rho_y"]).all(axis=1)))
                    collected[xp.where(sel)[0][fine]] = True
                else:
                    collected = collected | sel
                    D0[i] = xp.maximum(
                        D0[i], (delta[sel, i, :] + W[sel]).max(axis=0))
            if bool(collected.any()):
                C, W = C[~collected], W[~collected]
        if len(C):
            small = W.max(axis=1) <= wmin
            n_small = int(small.sum())
            stuck += n_small
            if stuck_out is not None and n_small:
                stuck_out.append((_to_host(xp, C[small]).copy(),
                                  _to_host(xp, W[small]).copy()))
            Cb, Wb = C[~small], W[~small]
            if len(Cb):
                j = xp.argmax(Wb, axis=1)
                rows = xp.arange(len(Cb))
                Wn = Wb.copy()
                Wn[rows, j] = Wb[rows, j] / 2.0
                Cl, Cr = Cb.copy(), Cb.copy()
                Cl[rows, j] -= Wn[rows, j]
                Cr[rows, j] += Wn[rows, j]
                stack_C.append(xp.concatenate([Cl, Cr]))
                stack_W.append(xp.concatenate([Wn, Wn]))
    return stuck, _to_host(xp, D0), total


def main():
    """Validation: xp=numpy must reproduce parametric.zoned_sweep exactly
    on the reference certified tile's sweep inputs."""
    import time

    import numpy as np

    from parametric import certify_tile
    import parametric

    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    h = 3e-4

    # capture the sweep inputs by monkey-patching, then replay via xp
    captured = {}
    orig = parametric.zoned_sweep

    def spy(H0, roots, coef0, coef1, guards, far_tax, **kw):
        r = orig(H0, roots, coef0, coef1, guards, far_tax, **kw)
        if "stage" not in captured:          # first (stage A) call
            captured.update(stage=(H0, [np.array(t) for t in roots],
                                   np.array(coef0), np.array(coef1),
                                   np.array(guards), far_tax, dict(kw)),
                            result=r)
        return r

    parametric.zoned_sweep = spy
    try:
        res = certify_tile(beta, h, verbose=False, use_certified=True)
    finally:
        parametric.zoned_sweep = orig
    assert res["ok"], res
    H0, roots, c0, c1, gu, ft, kw = captured["stage"]
    kw.pop("stuck_out", None)
    stuck0, D00, tot0 = captured["result"]
    t0 = time.time()
    stuck1, D01, tot1 = zoned_sweep_xp(H0, roots, c0, c1, gu, ft,
                                       stuck_out=[], **kw)
    dt = time.time() - t0
    print(f"cpu reference: stuck={stuck0} total={tot0}")
    print(f"xp=numpy:      stuck={stuck1} total={tot1}  [{dt:.1f} s]")
    print(f"D0 max |diff| = {np.max(np.abs(D00 - D01)):.2e}")
    ok = stuck0 == stuck1 and tot0 == tot1 and np.array_equal(D00, D01)
    print("EXACT MATCH" if ok else "MISMATCH")


if __name__ == "__main__":
    main()
