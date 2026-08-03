# Proof campaign ledger schema

Rigorous resume records use schema `mub6-ledger-v3-binary64-box`.

Every parameter box records:

- `beta_bits`: three IEEE-754 binary64 values, each encoded as exactly 16
  lowercase hexadecimal digits, matching `BALL_COVERAGE_SCHEMA.md`;
- `hv_bits`: the three half-widths in the same encoding;
- `box_bounds_bits`: the lower and upper binary64 results for all three
  componentwise `beta - hv` and `beta + hv` operations;
- `theta_interval_bits`: a checked alias of `box_bounds_bits[0]`;
- decimal JSON `beta` and `hv` mirrors for human inspection only.

The mirrors must decode to the identical bit patterns. NaNs, infinities,
noncanonical hex strings, mismatched mirrors, and recomputed endpoint
mismatches are invalid.

Rigorous frontier assembly:

1. rejects decimal-only legacy geometry;
2. groups transverse lines by the exact `beta_bits` values, with no decimal
   rounding;
3. deduplicates theta intervals by their endpoint bit patterns;
4. joins closed intervals only on exact overlap, with zero gap tolerance.

Historical decimal-only ledgers remain usable for experimental or sampled
analysis, but cannot advance a `RIGOROUS` resume frontier. This is intentional:
there is no canonical bit-level proposition for a future Lean rechecker to
recover from a rounded decimal convention.

The v3 all-axis bounds allow the bundle verifier to prove that every theta
interval also covers its manifest-assigned transverse cell; theta-line
completeness alone is not treated as three-dimensional coverage.

The canonical encoder and validator are in `ledger_bits.py`. New campaign,
parent-child, and rigorous-census writers all use it.
