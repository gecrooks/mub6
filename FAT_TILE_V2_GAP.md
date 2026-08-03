# Fat-tile to continuation-v2 gap

The landed `fat_tile.py` cannot currently produce a continuation-v2 artifact.
This is a fail-closed interface mismatch, not merely a serialization task.

## What the current result exports

- overall `SAMPLED_BOUND` result and evidence;
- root count, chromatic result, blanket drift and pair budget;
- aggregate sweep count and stage timings.

## Required data not exported

1. A complete `ParentCoverageArtifact`, including the parameter box, global
   sweep witness, stable root indices, and each zone's tube/fold kind.
2. For regular zones, certified `S[5,3]`, symmetric `Q[5,3,3]`, curve
   residual, Jacobian-remainder, and tube-radius bounds.
3. For folds, the certified tangent/transverse oracle frame and all eight
   octants when a split is used.

The current ball strategy may ultimately avoid `S/Q` continuation entirely.
If so, it should receive a separate versioned ball-coverage artifact whose
proof payload is the complete box-wise hull-to-ball association and blanket
drift bound. It should not populate continuation-v2 fields with sampled or
dummy values.

`fat_tile_continuation_adapter.py` is the integration boundary. It emits a
machine-readable readiness report with the missing proof classes. When a
complete parent artifact and typed continuation payloads are supplied, it
validates root counts, zone kinds, artifact grade, and digest before returning
v2. The fat tile's separate pair/coloring outcome does not override the
artifact's independently derived grade.
