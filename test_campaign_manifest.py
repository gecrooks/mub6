import json
from pathlib import Path
import tempfile
import unittest

from campaign_manifest import CampaignManifest, TransverseCell, build_manifest
from certificate_result import CertificateGrade


class CampaignManifestTests(unittest.TestCase):
    def test_manifest_round_trip_commits_to_ordered_shards(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            first, second = base / "a.jsonl", base / "b.jsonl"
            first.write_bytes(b"a\n")
            second.write_bytes(b"b\n")
            manifest = build_manifest(
                ((0.0, 1.0), (0.0, 2.0), (0.0, 3.0)),
                "Karlsson order-32 fundamental domain", 32,
                CertificateGrade.RIGOROUS, (first, second),
                base_directory=base,
            )
            restored = CampaignManifest.from_dict(
                json.loads(json.dumps(manifest.as_dict()))
            )
            self.assertEqual(restored, manifest)
            self.assertEqual(restored.verify_shards(base), ())
            reversed_manifest = build_manifest(
                manifest.domain, manifest.symmetry, 32,
                CertificateGrade.RIGOROUS, (second, first),
                base_directory=base,
            )
            self.assertNotEqual(manifest.manifest_id,
                                reversed_manifest.manifest_id)

    def test_shard_tampering_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            shard = base / "a.jsonl"
            shard.write_text("before\n")
            manifest = build_manifest(
                ((0.0, 1.0),) * 3, "test symmetry", 1,
                CertificateGrade.RIGOROUS, (shard,), base_directory=base,
            )
            shard.write_text("after\n")
            self.assertIn("digest mismatch", manifest.verify_shards(base)[0])

    def test_transverse_cells_must_be_exact_cartesian_partition(self):
        domain = ((0.0, 1.0), (0.0, 2.0), (0.0, 2.0))
        cells = tuple(
            TransverseCell((phi + 0.5, lam + 0.5),
                           (phi, phi + 1.0), (lam, lam + 1.0))
            for phi in (0.0, 1.0) for lam in (0.0, 1.0)
        )
        manifest = CampaignManifest(domain, "test", 1,
                                    CertificateGrade.RIGOROUS, (), cells)
        self.assertEqual(len(manifest.transverse_cells), 4)
        with self.assertRaisesRegex(ValueError, "Cartesian grid"):
            CampaignManifest(domain, "test", 1,
                             CertificateGrade.RIGOROUS, (), cells[:-1])
        with self.assertRaisesRegex(ValueError, "gap or overlap"):
            CampaignManifest(
                domain, "test", 1, CertificateGrade.RIGOROUS, (),
                (TransverseCell((0.25, 1.0), (0.0, 0.5), (0.0, 2.0)),
                 TransverseCell((1.25, 1.0), (1.0, 2.0), (0.0, 2.0))),
            )


if __name__ == "__main__":
    unittest.main()
