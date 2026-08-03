import json
from pathlib import Path
import tempfile
import unittest

from ball_coverage_artifact import (BallCoverageArtifact, DriftBoundClaim,
                                    SweepReplaySpec)
from campaign_artifacts import (TileArtifactStore, bind_tile_artifact,
                                verify_tile_binding)
from certificate_result import CertificateGrade
from ledger_bits import box_record


def artifact(beta=(0.5, 1.0, 2.0), hv=(0.1, 0.1, 0.1)):
    drift = DriftBoundClaim("test certified rate", beta, hv, 0.0, 0.0)
    sweep = SweepReplaySpec(
        "test deterministic sweep", 0.1, 0.1, "test interval taxes",
        1, 1, frontier_complete=True,
        arithmetic_grade=CertificateGrade.RIGOROUS,
    )
    return BallCoverageArtifact(
        beta, hv, 0.01, ((0.0,) * 5,), drift, sweep, 0.01, (), (0,), 1
    )


class CampaignArtifactTests(unittest.TestCase):
    def test_store_binding_and_kernel_verification(self):
        tile = artifact()
        record = box_record({"ok": True, "grade": "RIGOROUS"},
                            tile.parameter_center,
                            tile.parameter_half_widths)
        record = bind_tile_artifact(record, tile)
        with tempfile.TemporaryDirectory() as directory:
            store = TileArtifactStore(directory)
            path = store.save(tile)
            self.assertTrue(path.exists())
            report = verify_tile_binding(record, store)
        self.assertTrue(report.accepted)

    def test_missing_file_and_geometry_mismatch_fail_closed(self):
        tile = artifact()
        record = bind_tile_artifact(
            box_record({"ok": True, "grade": "RIGOROUS"},
                       tile.parameter_center, tile.parameter_half_widths), tile
        )
        with tempfile.TemporaryDirectory() as directory:
            report = verify_tile_binding(record, TileArtifactStore(directory))
        self.assertFalse(report.accepted)
        wrong = box_record({"ok": True, "grade": "RIGOROUS"},
                           (0.6, 1.0, 2.0), (0.1, 0.1, 0.1))
        with self.assertRaisesRegex(ValueError, "beta"):
            bind_tile_artifact(wrong, tile)

    def test_tampered_stored_artifact_is_rejected(self):
        tile = artifact()
        record = bind_tile_artifact(
            box_record({"ok": True, "grade": "RIGOROUS"},
                       tile.parameter_center, tile.parameter_half_widths), tile
        )
        with tempfile.TemporaryDirectory() as directory:
            store = TileArtifactStore(directory)
            path = store.save(tile)
            value = json.loads(path.read_text())
            value["chi_bound"] = 2
            path.write_text(json.dumps(value))
            report = verify_tile_binding(record, store)
        self.assertFalse(report.accepted)
        self.assertIn("digest", report.reason)
