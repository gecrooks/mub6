# Parent coverage reuse benchmark

Run on 2026-08-02 at the existing rigorous reference point
`(5.9785030164, 4.0075345498, 1.6327649325)`.

## Measurements

- An isotropic parent at `h=6e-4` failed before producing a coverage
  artifact (36.9 seconds). It therefore has zero reusable children.
- The parent at `h=5e-4` certified rigorously with 48 root zones in
  28.9 seconds.
- It contains eight overlapping isotropic children at `h=2.6e-4`.
- Representative child indices 0, 4, and 7 were run both with the parent
  artifact and with independent coverage.
- All three independent children certified rigorously, taking 24.3, 24.7,
  and 26.3 seconds (25.1-second mean).
- All three reused children failed rigorously and fail-closed: each reached
  the pair engine in about 0.45 seconds, but all 48 child tubes failed.

There is consequently **no valid measured speedup or break-even count yet**.
The benchmark program suppresses those projections whenever either path does
not certify rigorously.

## Finding

The artifact proves that the parent's root structures exhaust the parent
parameter box, but it currently exports only the structures' root coordinates
at the parent center. Those coordinates are not roots at an offset child
center and are inadequate seeds for the child tube construction.

The next artifact version must carry a certified root-restriction map. For
each regular parent structure this should provide a child seed/enclosure from
the parent continuation data (for example its certified `S`, `Q`, residual,
and remainder bounds), and folds need the corresponding restricted oracle.
Numerical polishing may propose a child seed, but cannot establish the
parent-to-child association by itself. The child Krawczyk/tube proof must then
verify the restricted enclosure before pair bounds run.

The independent baseline suggests that coverage is indeed the dominant cost,
so the potential remains substantial once continuation is exported. The
failed `6e-4` parent and successful `5e-4` parent put the observed isotropic
parent ceiling for this point in `(5e-4, 6e-4)`.
