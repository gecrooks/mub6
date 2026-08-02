"""Measure parent-coverage reuse against independent child certification."""

import argparse
import json
import math
import time

from certificate_result import CertificateGrade
from parent_coverage_campaign import partition_children


DEFAULT_BETA = (5.978503016422594, 4.007534549834652,
                1.6327649325136653)


def representative_indices(count, sample_count):
    if count < 1 or sample_count < 1:
        return ()
    if sample_count >= count:
        return tuple(range(count))
    return tuple(sorted(set(round(i * (count - 1) / (sample_count - 1))
                            for i in range(sample_count)))) \
        if sample_count > 1 else (count // 2,)


def summarize(parent_seconds, child_count, reused_seconds,
              independent_seconds):
    reused = tuple(float(value) for value in reused_seconds)
    independent = tuple(float(value) for value in independent_seconds)
    if not reused or len(reused) != len(independent):
        raise ValueError("paired nonempty child timings are required")
    reused_mean = sum(reused) / len(reused)
    independent_mean = sum(independent) / len(independent)
    saving = independent_mean - reused_mean
    projected_reuse = parent_seconds + child_count * reused_mean
    projected_independent = child_count * independent_mean
    return {
        "sample_count": len(reused),
        "reused_child_mean_seconds": reused_mean,
        "independent_child_mean_seconds": independent_mean,
        "marginal_child_speedup": (independent_mean / reused_mean
                                   if reused_mean else None),
        "coverage_saving_per_child_seconds": saving,
        "break_even_children": (max(1, math.ceil(parent_seconds / saving))
                                if saving > 0 else None),
        "projected_all_children_reuse_seconds": projected_reuse,
        "projected_all_children_independent_seconds": projected_independent,
        "projected_all_children_speedup": (
            projected_independent / projected_reuse
            if projected_reuse else None
        ),
    }


def benchmark(beta=DEFAULT_BETA, parent_h=5e-4, child_h=2.6e-4,
              sample_count=3, verbose=False):
    from parametric import certify_tile
    from rigor_tile import fully_rigorous_signed_tile

    started = time.time()
    parent_run = certify_tile(
        beta, (parent_h,) * 3, verbose=verbose, use_certified=True
    )
    parent_seconds = time.time() - started
    artifact = parent_run.get("coverage_artifact")
    if artifact is None or not artifact.complete:
        return {
            "ok": False,
            "stage": "parent",
            "reason": parent_run.get("reason",
                                     "parent coverage artifact unavailable"),
            "parent_seconds": parent_seconds,
        }
    if artifact.witness.evidence().grade < CertificateGrade.RIGOROUS:
        return {"ok": False, "stage": "parent",
                "reason": "parent coverage is not rigorous",
                "parent_seconds": parent_seconds}

    children = partition_children(artifact, (child_h,) * 3)
    indices = representative_indices(len(children), sample_count)
    samples = []
    for index in indices:
        center, _ = children[index]
        args = (center[0] - child_h, center[0] + child_h,
                center[1], center[2])
        started = time.time()
        reused = fully_rigorous_signed_tile(
            *args, hf=child_h, hf3=child_h, verbose=verbose,
            coverage_artifact=artifact,
        )
        reused_seconds = time.time() - started
        started = time.time()
        independent = fully_rigorous_signed_tile(
            *args, hf=child_h, hf3=child_h, verbose=verbose,
        )
        independent_seconds = time.time() - started
        samples.append({
            "child_index": index,
            "beta": list(center),
            "reused_ok": reused.ok,
            "reused_grade": reused.grade.name,
            "reused_seconds": reused_seconds,
            "independent_ok": independent.ok,
            "independent_grade": independent.grade.name,
            "independent_seconds": independent_seconds,
            "reused_roots": reused.metadata.get("n_roots"),
            "independent_roots": independent.metadata.get("n_roots"),
            "reused_tube_failures": reused.metadata.get("n_tube_fail"),
            "independent_tube_failures": independent.metadata.get(
                "n_tube_fail"
            ),
            "independent_coverage_seconds": independent.metadata.get(
                "coverage_seconds"
            ),
        })
    all_rigorous = all(
        sample["reused_ok"] and sample["independent_ok"]
        and sample["reused_grade"] == "RIGOROUS"
        and sample["independent_grade"] == "RIGOROUS"
        for sample in samples
    )
    result = {
        "ok": all_rigorous,
        "beta": list(beta),
        "parent_h": parent_h,
        "child_h": child_h,
        "artifact_id": artifact.artifact_id,
        "parent_seconds": parent_seconds,
        "child_count": len(children),
        "sample_indices": list(indices),
        "samples": samples,
        "comparison_valid": all_rigorous,
    }
    if all_rigorous:
        result.update(summarize(
            parent_seconds, len(children),
            [sample["reused_seconds"] for sample in samples],
            [sample["independent_seconds"] for sample in samples],
        ))
    else:
        result["reason"] = (
            "at least one paired child did not certify rigorously; "
            "speedup and break-even projection suppressed"
        )
        result["diagnostic_reused_mean_seconds"] = sum(
            sample["reused_seconds"] for sample in samples
        ) / len(samples)
        result["diagnostic_independent_mean_seconds"] = sum(
            sample["independent_seconds"] for sample in samples
        ) / len(samples)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta", default=",".join(map(str, DEFAULT_BETA)))
    parser.add_argument("--parent-h", type=float, default=5e-4)
    parser.add_argument("--child-h", type=float, default=2.6e-4)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    beta = tuple(float(value) for value in args.beta.split(","))
    if len(beta) != 3:
        parser.error("--beta must contain theta,phi,lambda")
    result = benchmark(beta, args.parent_h, args.child_h,
                       args.samples, args.verbose)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
