# Campaign bundle reference verifier

`campaign_bundle_verify.py` is the executable reference contract for a
complete rigorous campaign. It is deliberately small and uses binary64 bit
patterns throughout so that its accepted input can later be mirrored by a
Lean rechecker.

A bundle consists of:

- a `mub6-campaign-manifest-v2-cells` manifest;
- the ordered JSONL ledger shards committed by that manifest; and
- a content-addressed store of full-tile certificate artifacts.

The manifest's transverse cells must be the exact Cartesian product of
adjacent phi and lambda partitions covering the declared domain. Each cell
names a representative line. For every cell, the verifier accepts only
rigorous ledger records whose artifact passes kernel replay, whose box matches
the ledger bit-for-bit, and whose phi/lambda box contains the entire cell. The
accepted theta intervals must then cover the declared theta domain without a
binary64 gap.

Run it with:

```sh
python campaign_bundle_verify.py MANIFEST.json \
  --base-directory CAMPAIGN_DIR \
  --artifact-store CAMPAIGN_DIR/certificate_artifacts
```

It emits a deterministic JSON report and exits nonzero for a malformed
manifest or record, changed or missing shard, missing/tampered/open/weak
artifact, box mismatch, unlisted line, transverse undercoverage, or theta gap.
The verifier does not infer coverage between sampled lines.
