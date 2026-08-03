"""End-to-end Python reference verifier for rigorous campaign bundles."""

import argparse
from dataclasses import dataclass
import json
from pathlib import Path

from campaign_artifacts import TileArtifactStore, verify_tile_binding
from campaign_coverage import analyze_ledger_records, grade_accepted
from campaign_manifest import CampaignManifest
from ledger_bits import decode_box_bounds, decode_box_record, float_bits


@dataclass(frozen=True)
class BundleVerificationReport:
    ok: bool
    manifest_id: str | None
    shard_count: int
    record_count: int
    candidate_count: int
    verified_artifacts: int
    complete_cells: int
    expected_cells: int
    failures: tuple[str, ...]

    def as_dict(self):
        return {
            "ok": self.ok,
            "manifest_id": self.manifest_id,
            "shard_count": self.shard_count,
            "record_count": self.record_count,
            "candidate_count": self.candidate_count,
            "verified_artifacts": self.verified_artifacts,
            "complete_cells": self.complete_cells,
            "expected_cells": self.expected_cells,
            "failures": list(self.failures),
        }


def load_manifest(path):
    with open(path) as stream:
        return CampaignManifest.from_dict(json.load(stream))


def verify_bundle(manifest, base_directory, artifact_store):
    """Verify manifest, shards, artifacts, and exact 3D box coverage."""
    base = Path(base_directory)
    store = (artifact_store if isinstance(artifact_store, TileArtifactStore)
             else TileArtifactStore(artifact_store))
    failures = list(manifest.verify_shards(base))
    if not manifest.transverse_cells:
        failures.append("manifest has no transverse coverage cells")
    if failures:
        return BundleVerificationReport(
            False, manifest.manifest_id, len(manifest.ledger_shards), 0, 0,
            0, 0, len(manifest.transverse_cells), tuple(failures)
        )

    records = []
    for shard in manifest.ledger_shards:
        path = base / shard.name
        with path.open() as stream:
            for line_number, line in enumerate(stream, 1):
                try:
                    value = json.loads(line)
                except json.JSONDecodeError as error:
                    failures.append(f"{shard.name}:{line_number}: {error}")
                    continue
                if not isinstance(value, dict):
                    failures.append(f"{shard.name}:{line_number}: "
                                    "record is not an object")
                    continue
                records.append((shard.name, line_number, value))

    cells = {(float_bits(cell.line[0]), float_bits(cell.line[1])): cell
             for cell in manifest.transverse_cells}
    eligible = {key: [] for key in cells}
    candidates = verified = 0
    for shard, line_number, record in records:
        if not grade_accepted(record, manifest.required_grade):
            continue
        candidates += 1
        prefix = f"{shard}:{line_number}"
        binding = verify_tile_binding(record, store, manifest.required_grade)
        if not binding.accepted:
            failures.append(f"{prefix}: {binding.reason}")
            continue
        try:
            beta, _widths, _interval, _token = decode_box_record(record)
            bounds = decode_box_bounds(record)
        except (KeyError, TypeError, ValueError) as error:
            failures.append(f"{prefix}: {error}")
            continue
        key = (float_bits(beta[1]), float_bits(beta[2]))
        cell = cells.get(key)
        if cell is None:
            failures.append(f"{prefix}: line is absent from manifest cells")
            continue
        if bounds[1][0] > cell.phi_bounds[0] \
                or bounds[1][1] < cell.phi_bounds[1] \
                or bounds[2][0] > cell.lambda_bounds[0] \
                or bounds[2][1] < cell.lambda_bounds[1]:
            failures.append(f"{prefix}: tile does not cover transverse cell")
            continue
        verified += 1
        eligible[key].append(record)

    complete_cells = 0
    theta_lo, theta_hi = manifest.domain[0]
    for key, cell in cells.items():
        cell_records = eligible[key]
        if not cell_records:
            failures.append(f"cell {cell.line}: no verified tile records")
            continue
        reports = analyze_ledger_records(
            cell_records, theta_lo, theta_hi, manifest.required_grade
        )
        report = reports.get(cell.line)
        if report is None or not report.complete:
            gap_count = "missing" if report is None else len(report.gaps)
            failures.append(f"cell {cell.line}: theta coverage {gap_count}")
            continue
        complete_cells += 1

    return BundleVerificationReport(
        not failures, manifest.manifest_id, len(manifest.ledger_shards),
        len(records), candidates, verified, complete_cells, len(cells),
        tuple(failures),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    parser.add_argument("--base-directory", default=".")
    parser.add_argument("--artifact-store", default="certificate_artifacts")
    args = parser.parse_args()
    try:
        manifest = load_manifest(args.manifest)
        report = verify_bundle(
            manifest, args.base_directory, args.artifact_store
        )
    except (OSError, KeyError, TypeError, ValueError,
            json.JSONDecodeError) as error:
        report = BundleVerificationReport(
            False, None, 0, 0, 0, 0, 0, 0, (str(error),)
        )
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    raise SystemExit(0 if report.ok else 1)


if __name__ == "__main__":
    main()
