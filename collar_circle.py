"""Full-circle collar chains: walk beta2 around the whole torus at
fixed beta3 — the collar campaign's scale test (warm pools, signed
tiles, wall-aware widths). Usage:
    python collar_circle.py <b3> <hf3> [curv3]
"""

import sys

import numpy as np

from collar_chain import chain


def main():
    b3 = float(sys.argv[1])
    hf3 = float(sys.argv[2])
    curv3 = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
    hf = 5e-3
    n = int(np.floor(2 * np.pi / (2 * hf)))
    chain(0.005, 0.01, 0.05, b3, hf=hf, n_steps=n, hf3=hf3,
          curv=(10.0, 10.0, curv3))


if __name__ == "__main__":
    main()
