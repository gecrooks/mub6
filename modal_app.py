"""Modal app: GPU-accelerated Layer-3 X-cover.

- bench(): sweep timings CPU-vs-GPU at several deltas on one anchor.
- line(c): one full chain-line walk (walker on container CPU, sweeps on
  GPU via monkeypatched certified_triple_sweep).
Run:  ./.venv/bin/python -m modal run modal_app.py            (bench)
      ./.venv/bin/python -m modal run modal_app.py --lines "0.4,0.45"
"""

import modal

app = modal.App("mub6-xcover")

image = (modal.Image.from_registry(
             "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.12")
         .pip_install("numpy", "scipy", "cupy-cuda12x")
         .add_local_python_source(
             "gpusweep", "layer3", "layer3_param", "layer3_x", "mub",
             "szollosi", "parametric", "certify", "fold", "karlsson",
             "interval", "dual", "rates", "xcover", "families",
             "degenerate", "gputile", "tm", "tmres", "cache", "campaign"))


def _gpu_sweep_fn():
    import cupy as cp

    from gpusweep import certified_triple_sweep_xp

    def sweep(B, K, hslop=1e-9, wmin=2e-3, chunk=1_000_000,
              max_boxes=4e8, verbose=False):
        return certified_triple_sweep_xp(B, K, hslop=hslop, wmin=wmin,
                                         chunk=chunk, max_boxes=max_boxes,
                                         xp=cp)
    return sweep


@app.function(gpu="A100", image=image, timeout=1800)
def bench():
    import time
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np

    import layer3
    from layer3_param import track_K
    from layer3_x import triples_at, CHOICE
    from parametric import PAD
    from szollosi import szollosi_map

    gpu_sweep = _gpu_sweep_fn()
    a = 0.5 + 0.131j
    H, Ks, _nv = triples_at(a, seed=1)
    K = Ks[0]
    eps = 5e-5
    H1, _c, _d = szollosi_map(a + eps, choice=CHOICE)
    K1 = track_K(H1, K)
    rate = (np.max(np.abs(H1 - H)) + np.max(np.abs(K1 - K))) / eps
    out = []
    for delta in (6e-3, 1.2e-2, 2e-2):
        hs = PAD * rate * delta
        t0 = time.time()
        n1, _m, b1 = layer3.certified_triple_sweep(
            H, K, hslop=hs, wmin=2e-3, max_boxes=4e7, verbose=False)
        t_cpu = time.time() - t0
        t0 = time.time()
        n2, _m2, b2 = gpu_sweep(H, K, hslop=hs)
        t_gpu = time.time() - t0
        t0 = time.time()
        n3, _m3, _b3 = gpu_sweep(H, K, hslop=hs)
        t_gpu2 = time.time() - t0
        out.append((delta, b1, t_cpu, n1, b2, t_gpu, t_gpu2, n3))
    return out


@app.function(gpu="A100", image=image, timeout=3600)
def line(c: float):
    import warnings
    warnings.filterwarnings("ignore")
    import xcover

    xcover.certified_triple_sweep = _gpu_sweep_fn()
    return xcover.walk_line(c)


@app.function(gpu="A100", image=image, timeout=3600)
def patch_line(c: float):
    import warnings
    warnings.filterwarnings("ignore")
    import xcover

    xcover.certified_triple_sweep = _gpu_sweep_fn()
    return xcover.walk_line(c, delta_min=2e-4, n_starts=20000,
                            rebuild_every=5)


@app.function(gpu="A100", image=image, timeout=1800, memory=16384)
def tile_bench():
    """Zoned TILE sweep (certify_tile stage A) CPU vs A100: capture the
    reference tile's sweep inputs via the spy pattern, then replay with
    xp=numpy (in-container reference) and xp=cupy (cold + warm)."""
    import time
    import warnings
    warnings.filterwarnings("ignore")
    import cupy as cp
    import numpy as np

    import parametric
    from gputile import zoned_sweep_xp
    from parametric import certify_tile

    beta = (5.978503016422594, 4.007534549834652, 1.6327649325136653)
    captured = {}
    orig = parametric.zoned_sweep

    def spy(H0, roots, coef0, coef1, guards, far_tax, **kw):
        r = orig(H0, roots, coef0, coef1, guards, far_tax, **kw)
        if "stage" not in captured:
            captured.update(stage=(H0, [np.array(t) for t in roots],
                                   np.array(coef0), np.array(coef1),
                                   np.array(guards), far_tax, dict(kw)),
                            result=r)
        return r

    parametric.zoned_sweep = spy
    try:
        t0 = time.time()
        res = certify_tile(beta, 3e-4, verbose=False, use_certified=True)
        t_tile = time.time() - t0
    finally:
        parametric.zoned_sweep = orig
    assert res["ok"], res
    print(f"TILE_OK {t_tile:.0f}s", flush=True)
    H0, roots, c0, c1, gu, ft, kw = captured["stage"]
    kw.pop("stuck_out", None)
    s0, _D, b0 = captured["result"]
    t0 = time.time()
    s1, _D1, b1 = zoned_sweep_xp(H0, roots, c0, c1, gu, ft, **kw)
    t_np = time.time() - t0
    print(f"XP_NP_OK {t_np:.1f}s", flush=True)
    t0 = time.time()
    s2, _D2, b2 = zoned_sweep_xp(H0, roots, c0, c1, gu, ft, xp=cp, **kw)
    t_gpu = time.time() - t0
    t0 = time.time()
    s3, _D3, b3 = zoned_sweep_xp(H0, roots, c0, c1, gu, ft, xp=cp, **kw)
    t_gpu2 = time.time() - t0
    out = dict(tile_s=t_tile, cpu=(s0, b0), xp_np=(s1, b1, t_np),
               gpu_cold=(s2, b2, t_gpu), gpu_warm=(s3, b3, t_gpu2))
    print("TILE_BENCH_RESULT", out, flush=True)   # lands in app logs
    return out


@app.function(gpu="A100", image=image, timeout=3600, memory=16384)
def campaign_line(spec: dict):
    """One campaign theta-line with the zoned sweep on GPU. Designed for
    .spawn(): ledger records and the summary go to app logs."""
    import io
    import json
    import warnings
    warnings.filterwarnings("ignore")
    import cupy as cp

    import cache as cache_mod
    import parametric
    from gputile import zoned_sweep_xp

    def sweep(*a, **kw):
        kw.setdefault("chunk", 1_000_000)
        return zoned_sweep_xp(*a, **kw, xp=cp)

    parametric.zoned_sweep = sweep      # certify_tile path
    cache_mod.zoned_sweep = sweep       # anchored_tile/chain_step path
    from campaign import run_line
    buf = io.StringIO()
    covered, n = run_line(spec["phi"], spec["lam"], spec["th_lo"],
                          spec["th_hi"], spec["h"], buf,
                          start_at=spec.get("start"))
    for rec in buf.getvalue().splitlines():
        print("LEDGER " + rec, flush=True)
    print("CAMPAIGN_LINE_RESULT " + json.dumps(
        dict(covered=covered, n=n, **spec)), flush=True)
    return dict(covered=covered, n=n, ledger=buf.getvalue())


@app.local_entrypoint()
def main(lines: str = "", patch: str = "", tile: bool = False,
         tile_spawn: bool = False, campaign: str = ""):
    if campaign:
        phi, lam, th_lo, th_hi, h = [float(x) for x in campaign.split(",")]
        call = campaign_line.spawn(dict(phi=phi, lam=lam, th_lo=th_lo,
                                        th_hi=th_hi, h=h))
        print(f"SPAWNED {call.object_id} — collect via app logs "
              f"(LEDGER / CAMPAIGN_LINE_RESULT)", flush=True)
        return
    if tile_spawn:
        call = tile_bench.spawn()
        print(f"SPAWNED {call.object_id} — result in app logs "
              f"(TILE_BENCH_RESULT)", flush=True)
        return
    if tile:
        r = tile_bench.remote()
        print(f"in-container tile: {r['tile_s']:.0f} s; "
              f"cpu ref {r['cpu']}")
        print(f"xp=numpy  {r['xp_np'][:2]}  {r['xp_np'][2]:.1f} s")
        print(f"gpu cold  {r['gpu_cold'][:2]}  {r['gpu_cold'][2]:.1f} s")
        print(f"gpu warm  {r['gpu_warm'][:2]}  {r['gpu_warm'][2]:.1f} s")
        return
    if patch:
        cs = [float(x) for x in patch.split(",")]
        for r in patch_line.map(cs):
            print(r, flush=True)
    elif lines:
        cs = [float(x) for x in lines.split(",")]
        for r in line.map(cs):
            print(r, flush=True)
    else:
        rows = bench.remote()
        print("delta | CPU boxes/time | GPU cold | GPU warm | suspects")
        for (d, b1, tc, n1, b2, tg, tg2, n3) in rows:
            print(f"{d:g} | {b1} boxes {tc:.2f}s | {tg:.2f}s | {tg2:.2f}s "
                  f"| cpu {n1} gpu {n3}")
