import io
import json
from pathlib import Path
import tempfile
import unittest

from certificate_result import CertificateGrade, CertificateResult, Evidence
from coverage_contract import RootCoverageZone, SweepCoverageWitness
from parent_coverage_campaign import (
    CoverageArtifactStore,
    certify_and_store_parent,
    completed_child_keys,
    partition_children,
    run_children,
)


def artifact(half_widths=(1.0, 1.0, 1.0)):
    zone = RootCoverageZone(
        kind="tube", center=(0.0,) * 5, guard_radii=(0.2,) * 5,
        collected_reach=(0.1,) * 5, enclosure_radii=(0.01,) * 5,
        handoff_complete=True, grade=CertificateGrade.RIGOROUS,
    )
    return SweepCoverageWitness(
        zones=(zone,), global_sweep_complete=True,
        arithmetic_grade=CertificateGrade.RIGOROUS,
        boxes_processed=12, parameter_center=(1.0, 2.0, 3.0),
        parameter_half_widths=half_widths,
    ).artifact()


class ParentCoverageCampaignTests(unittest.TestCase):
    def test_store_round_trip_is_content_addressed(self):
        parent = artifact()
        with tempfile.TemporaryDirectory() as directory:
            store = CoverageArtifactStore(directory)
            path = store.save(parent)
            self.assertEqual(path.name, f"{parent.artifact_id}.json")
            self.assertEqual(store.save(parent), path)
            self.assertEqual(store.load(parent.artifact_id), parent)

    def test_partition_is_contained_and_covers_expected_grid(self):
        parent = artifact((1.0, 0.5, 0.5))
        children = partition_children(parent, (0.5, 0.25, 0.5))
        self.assertEqual(len(children), 4)
        self.assertTrue(all(parent.restrict(*child).complete
                            for child in children))

    def test_parent_creation_persists_coverage_even_if_pairs_fail(self):
        parent = artifact()

        def certify(beta, widths, **kwargs):
            return {"ok": False, "coverage_artifact": parent}

        with tempfile.TemporaryDirectory() as directory:
            made = certify_and_store_parent(
                (1.0, 2.0, 3.0), (1.0, 1.0, 1.0),
                CoverageArtifactStore(directory), certifier=certify,
                verbose=False,
            )
            self.assertFalse(made["parent_result_ok"])
            self.assertEqual(made["coverage_grade"],
                             CertificateGrade.RIGOROUS)
            self.assertTrue(Path(made["path"]).exists())

    def test_resume_requires_same_artifact_grade_and_success(self):
        parent = artifact()
        base = {"mode": "parent-coverage-child",
                "coverage_artifact_id": parent.artifact_id,
                "beta": [1.0, 2.0, 3.0], "hv": [1.0, 1.0, 1.0],
                "ok": True, "grade": "RIGOROUS"}
        records = [base, dict(base, coverage_artifact_id="0" * 64),
                   dict(base, grade="SAMPLED_BOUND"), dict(base, ok=False)]
        self.assertEqual(len(completed_child_keys(
            records, parent.artifact_id)), 1)

    def test_runner_binds_records_and_skips_accepted_children(self):
        parent = artifact((0.5, 0.5, 0.5))
        calls = []

        def certify(*args, **kwargs):
            calls.append((args, kwargs))
            return CertificateResult.from_evidence(
                True, [Evidence("pair", CertificateGrade.RIGOROUS)]
            )

        ledger = io.StringIO()
        first = run_children(parent, 0.5, ledger, certify_child=certify,
                             verbose=False)
        records = [json.loads(line) for line in ledger.getvalue().splitlines()]
        second_ledger = io.StringIO()
        second = run_children(
            parent, 0.5, second_ledger, prior_records=records,
            certify_child=certify, verbose=False,
        )
        self.assertEqual(first["written"], 1)
        self.assertEqual(first["accepted"], 1)
        self.assertEqual(records[0]["coverage_artifact_id"],
                         parent.artifact_id)
        self.assertIs(calls[0][1]["coverage_artifact"], parent)
        self.assertEqual(second["skipped"], 1)
        self.assertEqual(second["written"], 0)


if __name__ == "__main__":
    unittest.main()
