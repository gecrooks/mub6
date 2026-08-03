# Certificate ledger storage & distribution

*Drafted 2026-08-02. Builds on `campaign_manifest.py`
(`mub6-campaign-manifest-v1`) and ARTIFACT_BINDING.md — this file
adds the hosting/distribution policy only; the commitment layer is
already specified there. For joint review at the schema freeze.*

## Size basis (measured, 4.82-era serializers)

| object | raw | gzipped |
|---|---|---|
| ball-coverage-v1 tile (48 roots) | ~31 KB | ~4 KB |
| continuation-v2-style tile | ~58 KB | ~26 KB |
| campaign ledger (2.5e5–4e6 tiles, realistic mix) | — | **~1–30 GB** (worst case ~100 GB) |

Too big for GitHub (100 MB/file cap, repos degrade past a few GB;
LFS caps at 2 GiB/file and bills egress). The ledger is also
*replayable* — regenerable from the campaign spec — so only the
frozen published copy needs archival guarantees.

## Principle: the repo hosts the commitment, not the data

`CampaignManifest.manifest_id` (SHA-256 over the canonical
manifest, which commits to domain bounds, symmetry quotient,
required grade, and the ordered shard name/digest list) is the
single hash the paper cites. Anyone can fetch any subset of shards
from any host and verify against the manifest with
`verify_shards()`; no host is trusted. A Merkle tree over shards
was considered and rejected: the flat manifest is only a few MB,
so partial verification never needs tree paths.

## Tier 0 — GitHub repo (small, permanent, versioned)

- code + Lean re-checker,
- the campaign manifest JSON (few MB),
- this policy, the schema docs, and the frozen schema versions.

Nothing else. No certificate payloads, no LFS.

## Tier 1 — working store during the campaign

- **Where:** Modal volume for the GPU workers' write path (already
  colocated); mirrored to **Cloudflare R2** for sharing
  (S3-compatible, ~$0.015/GB-month → full ledger ≈ $0.50/month,
  **zero egress fees** so collaborators pull free).
- **Ledger shards:** gzipped JSONL, ~1,000 records/shard
  (~4–25 MB), keyed by territory/stratum so a verifier streams one
  region at a time. Shard naming = the `LedgerShard.name` safe
  relative paths from the manifest; shard order is significant
  (ARTIFACT_BINDING.md).
- **Artifact store:** content-addressed as
  `artifacts/<schema>/<id[:2]>/<id>.json.gz` — immutable objects,
  referenced transitively from ledger records by schema/digest
  pairs, per the binding rules. Write-once; digest changes are
  tamper evidence, never updates.
- A worker holds only its in-flight batch in memory
  (tens–hundreds of MB); nothing ever needs the whole ledger
  resident.

## Tier 2 — publication (paper submission)

- **Zenodo** record: frozen ledger shards + manifest + the exact
  re-checker release. Free, CERN-operated, versioned, DOI-minted;
  50 GB/record (more on request; worst case splits across records
  of one Zenodo "community" without weakening the commitment —
  the manifest spans records).
- The paper cites: DOI + `manifest_id`. Verification workflow for
  a third party: clone repo (tier 0) → fetch shards (tier 2 or any
  mirror) → `verify_shards()` → run the re-checker.
- Optional cheap redundancy: mirror the frozen record to the
  Internet Archive / Academic Torrents. Redundant hosts cost
  nothing in trust because of the manifest.

## Open questions for the freeze

1. Shard keying: territory/stratum boundaries should match the
   frontier-ledger grouping (coverage-verifier's call).
2. Artifact-store GC: pre-freeze superseded artifacts may be
   pruned; anything referenced by a manifest-committed shard is
   permanent.
3. R2 bucket naming/ownership and whether the campaign writes to
   R2 directly or syncs from the Modal volume in batches
   (recommend: batch sync, hourly, sync failure never blocks
   workers).
