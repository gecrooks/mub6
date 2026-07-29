# CLAUDE.md — mub6

Certified-exclusion program for the dimension-6 MUB problem (Zauner's
MUB conjecture). This is the PRIMARY repo as of 2026-07-29 —
development previously lived in `~/Work/ops/MUB6` (now frozen; its
git history was subtree-split into this repo's lineage).

## Read first, in order
1. `README.md` — claims-vs-status table (what is proven, at what grade)
2. `PROOF_SKELETON.md` — the lemma chain and assumption ledger
3. `RESULTS.md` — the full session record (Results 1-34+)
4. `PROOF_ROADMAP.md` — open threads; `NOTES_LP_BRIDGE.md` — the
   two-regime/LP theory program (active research direction)

## Layout & conventions
- Bare repo: `.bare/` + `.git` pointer; worktree per branch
  (`main/` = the main worktree). Work happens in worktrees.
- **Never run git/gh commands that change remote state** (no push,
  no gh mutations). Prepare locally; Gavin pushes.
- Cloud (Modal, account `gavincrooks`, 10-GPU cap, A100s): ask
  before anything beyond a few dollars; cents-scale benchmarks fine.
  Dispatch via the DEPLOYED app + `.spawn()` (live `.remote()` calls
  die when the local client drops); containers need memory=16384.

## Run
```
.venv/bin/python -c "from parametric import certify_tile; ..."
```
venv has numpy/scipy/mpmath/modal. GPU sweeps: `gputile.py` (CuPy,
box-exact vs CPU). Long jobs: run backgrounded, logs to files.
Recurring gotcha: shell cwd resets between commands — always
`cd ~/Work/mub6/main` first.

## State pointers
- Tile certificates are assumption-free (no PAD/libm/sampling) —
  ledger in PROOF_SKELETON.md section 4.
- Best-territory h = 7e-4 (Result 34); 1e-3 walls on subdivision
  economics. Full-domain flat campaign ~$200-400k — superseded as a
  plan by the two-regime program (NOTES_LP_BRIDGE.md sections
  4.5-4.7): thin rigidity collar at the theta=0 branch face
  (= Fourier locus), clique-starved bulk, fat starvation tiles TBD.
