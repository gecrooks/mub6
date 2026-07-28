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
             "degenerate"))


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


@app.local_entrypoint()
def main(lines: str = "", patch: str = ""):
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
