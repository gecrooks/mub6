import json
import math
from pathlib import Path
import tempfile
import unittest

from ball_coverage_artifact import (BallCoverageArtifact, DriftBoundClaim,
                                    SweepReplaySpec)
from campaign_artifacts import TileArtifactStore, bind_tile_artifact
from campaign_bundle_verify import verify_bundle
from campaign_manifest import TransverseCell, build_manifest
from certificate_result import CertificateGrade
from ledger_bits import bits_float, box_record, float_bits


def tile_artifact(frontier=True):
    beta = (0.5, 0.5, 0.5)
    hv = (0.5, 0.5, 0.5)
    drift = DriftBoundClaim("test certified rate", beta, hv, 0.0, 0.0)
    sweep = SweepReplaySpec(
        "test deterministic sweep", 0.1, 0.1, "test interval taxes",
        1, 1, frontier_complete=frontier,
        arithmetic_grade=CertificateGrade.RIGOROUS,
    )
    return BallCoverageArtifact(
        beta, hv, 0.01, ((0.0,) * 5,), drift, sweep, 0.01, (), (0,), 1
    )


class BundleFixture:
    def __init__(self, directory, *, frontier=True):
        self.base = Path(directory)
        self.store = TileArtifactStore(self.base / "artifacts")
        self.artifact = tile_artifact(frontier)
        self.store.save(self.artifact)
        record = box_record({"ok": True, "grade": "RIGOROUS"},
                            self.artifact.parameter_center,
                            self.artifact.parameter_half_widths)
        self.record = bind_tile_artifact(record, self.artifact)
        self.shard = self.base / "ledger.jsonl"
        self.write_record(self.record)
        self.rebuild_manifest()

    def write_record(self, record):
        self.record = record
        self.shard.write_text(json.dumps(record, sort_keys=True) + "\n")

    def rebuild_manifest(self):
        self.manifest = build_manifest(
            ((0.0, 1.0),) * 3, "identity test domain", 1,
            CertificateGrade.RIGOROUS, (self.shard,),
            base_directory=self.base,
            transverse_cells=(TransverseCell(
                (0.5, 0.5), (0.0, 1.0), (0.0, 1.0)
            ),),
        )


class CampaignBundleVerifyTests(unittest.TestCase):
    def test_complete_synthetic_bundle_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory)
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
        self.assertTrue(report.ok)
        self.assertEqual(report.complete_cells, 1)
        self.assertEqual(report.verified_artifacts, 1)

    def test_shard_tampering_fails_before_records_are_trusted(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory)
            fixture.shard.write_text(fixture.shard.read_text() + " \n")
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
        self.assertFalse(report.ok)
        self.assertIn("digest mismatch", report.failures[0])

    def test_missing_artifact_and_box_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory)
            fixture.store.path_for(fixture.artifact.artifact_id).unlink()
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
            self.assertFalse(report.ok)
            self.assertTrue(any("No such file" in x for x in report.failures))

        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory)
            wrong = box_record({"ok": True, "grade": "RIGOROUS",
                                "artifact_schema": fixture.record["artifact_schema"],
                                "artifact_id": fixture.record["artifact_id"]},
                               (0.4, 0.5, 0.5), (0.4, 0.5, 0.5))
            fixture.write_record(wrong)
            fixture.rebuild_manifest()
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
            self.assertFalse(report.ok)
            self.assertTrue(any("beta" in x for x in report.failures))

    def test_open_artifact_and_one_bit_endpoint_change_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory, frontier=False)
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
            self.assertFalse(report.ok)
            self.assertTrue(any("incomplete" in x or "below" in x
                                for x in report.failures))

        with tempfile.TemporaryDirectory() as directory:
            fixture = BundleFixture(directory)
            changed = json.loads(json.dumps(fixture.record))
            upper = bits_float(changed["box_bounds_bits"][0][1])
            changed["box_bounds_bits"][0][1] = float_bits(
                math.nextafter(upper, math.inf)
            )
            fixture.write_record(changed)
            fixture.rebuild_manifest()
            report = verify_bundle(fixture.manifest, fixture.base,
                                   fixture.store)
            self.assertFalse(report.ok)
            self.assertTrue(any("bounds bits" in x for x in report.failures))


if __name__ == "__main__":
    unittest.main()
