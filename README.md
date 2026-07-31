# Certified exclusion for the dimension-6 MUB problem

A computer-assisted proof machine attacking **Zauner's MUB
conjecture**: that ℂ⁶ admits at most 3 mutually unbiased bases
(every other small dimension being a prime power, d = 6 is the first
open case, and it has been open since 1999).

**Target theorem (not yet proven — see status):** no four mutually
unbiased bases of ℂ⁶ contain a complex Hadamard matrix from the
H₂-reducible (Karlsson) family K₆⁽³⁾ — the three-parameter family
that existing LP/SDP/Gröbner methods provably cannot reach.

Everything here was built in a ~5-day sprint by Claude (Anthropic),
directed by Gavin Crooks, as an experiment in agentic mathematics.
The session record is honest to a fault: `RESULTS.md` documents 34
results *including* the dead ends, walls, and an $8 GPU survey whose
main finding was that the naive full campaign is a ~$10⁶ object.

## What is actually proven (and at what grade)

| Claim | Status |
|---|---|
| No MU quadruple contains the Björck/Beauchamp–Nicoară arc B(θ), full arc incl. branch points | **Certified, gap-free** (arc rates FD+PAD — "prototype grade") |
| Same for the Diţă arc and the full two-parameter Szöllősi region | **Certified, gap-free** (same grade) |
| Pointwise anchors: S₆, Björck's C, K₃ points, F₆ quadruple wall at exactly 1/6 | **Certified** |
| Tile certificates (3-param Karlsson boxes): machinery | **Assumption-free**: no PAD, no sampled constants, no libm, no enumeration trust (`PROOF_SKELETON.md` §4) |
| The target theorem over the full fundamental domain | **NOT proven** — two-regime program (PROOF_SKELETON.md §5): bulk coarse tiles validated 6/6 and chain-certified, ~$25–50k CPU at current amortization with mechanical multipliers pending; collar = thickened Layer-3 (ingredients built, certificates to come) |

## Read first

- `PROOF_SKELETON.md` — the full lemma chain, referee-first: the
  compactness/root-count argument (enumeration is never trusted;
  missed roots fail loudly), and the assumption ledger.
- `RESULTS.md` — the 34-result session record, in order, with
  failures included.
- `PROOF_ROADMAP.md` — gap list and open threads (mode-B branch-aware
  solves; second-order sweep taxes; the campaign economics).
- `NOTES_LP_BRIDGE.md` — the theory bet: the observed exact-1/6 walls
  as an LP-duality face; a ℤ₂-graded Delsarte LP from H₂-reducibility
  that could replace the compute campaign with one infeasibility
  certificate.

## Run a certificate

```
python -m venv .venv && .venv/bin/pip install numpy scipy mpmath
.venv/bin/python -c "
from parametric import certify_tile
r = certify_tile((5.978503016422594, 4.007534549834652,
                  1.6327649325136653), 3e-4, use_certified=True)
print(r['ok'])"          # ~30 s: one assumption-free tile theorem
```

`gputile.py` is the CuPy port (box-exact vs CPU, 41× on A100);
`campaign.py`/`dispatch.py` the (validated, unlaunched) campaign
driver; `trig_kernel.py` the self-contained certified sin/cos that
removed the last libm dependence.

## Provenance

Built with Claude (Claude Code) over 2026-07-24 → 07-29. Key prior
art: Karlsson's K₆⁽³⁾ parametrization (arXiv:1003.4177), the
McNulty–Weigert review (arXiv:2410.23997), Matolcsi–Weiner /
Jaming et al. on LP methods, Grassl, Brierley–Weigert, and
Bengtsson et al.'s Hadamard catalogue. Contact: Gavin Crooks
(gavincrooks@gmail.com).
