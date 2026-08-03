# Ledger-to-certificate artifact binding

A rigorous ledger record authorizes resume only when all three layers agree:

1. its geometry passes `mub6-ledger-v3-binary64-box` validation;
2. `artifact_schema` and `artifact_id` resolve to an immutable object in the
   content-addressed artifact store;
3. that object is a complete full-tile certificate of the required grade,
   matches the ledger box bit-for-bit, and passes its Python kernel model.

The first registered full-tile schema is `mub6-ball-coverage-v1`. Parent
coverage v1 and continuation v2 are deliberately not accepted: they prove
dependencies of a child certificate, not the child pair/coloring conclusion.

Missing files, unknown schemas, digest changes, open sweep frontiers, weak
grades, box mismatches, and kernel failures all prevent frontier advancement.
Unbound historical records remain auditable but are not rigorous resume input.

`campaign_manifest.py` defines `mub6-campaign-manifest-v2-cells`. Its digest
commits to binary64 domain bounds, a Cartesian partition of the two transverse
axes and its representative line for each cell, the named symmetry quotient
and factor, the required grade, and the ordered names and SHA-256 digests of
every ledger shard. Shard order is significant. Artifact JSON files are
referenced transitively by the schema/digest pairs in those shards rather than
embedded in the manifest.

`campaign_bundle_verify.py` is the fail-closed Python reference rechecker. It
checks the manifest and shard hashes before parsing records, replays each
certificate artifact, requires every accepted tile to contain its complete
transverse cell, and proves exact theta coverage on every cell. Thus a set of
complete theta lines is not mistaken for coverage of the full 3D domain.
