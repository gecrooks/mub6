"""Summarize rigorous-census JSONL output."""

import argparse
import json
import statistics


def summarize(records):
    points = {}
    failures = {}
    rigorous = []
    for record in records:
        beta = tuple(float(x) for x in record["beta"])
        hv = tuple(float(x) for x in record["hv"])
        point = points.setdefault(beta, {
            "beta": list(beta), "attempts": 0, "rigorous_passes": 0,
            "max_rigorous_h": None,
        })
        point["attempts"] += 1
        if record.get("accepted_rigorous"):
            point["rigorous_passes"] += 1
            h = min(hv)
            point["max_rigorous_h"] = max(
                h, point["max_rigorous_h"] or h
            )
            rigorous.append(record)
        else:
            result = record.get("result", {})
            reason = result.get("reason") or "unknown failure"
            failures[reason] = failures.get(reason, 0) + 1

    seconds = [float(record["seconds"]) for record in rigorous]
    boxes = [
        int(record.get("result", {}).get("metadata", {})
            .get("coverage_boxes", 0))
        for record in rigorous
    ]
    return {
        "records": len(records),
        "points": len(points),
        "rigorous_passes": len(rigorous),
        "pass_fraction": len(rigorous) / len(records) if records else 0.0,
        "seconds_median": statistics.median(seconds) if seconds else None,
        "seconds_max": max(seconds) if seconds else None,
        "coverage_boxes_median": statistics.median(boxes) if boxes else None,
        "coverage_boxes_max": max(boxes) if boxes else None,
        "failures": failures,
        "point_results": [points[key] for key in sorted(points)],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    with open(args.ledger) as stream:
        records = [json.loads(line) for line in stream if line.strip()]
    report = summarize(records)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(f"{report['rigorous_passes']}/{report['records']} rigorous "
          f"across {report['points']} points; median "
          f"{report['seconds_median']}s, median "
          f"{report['coverage_boxes_median']} coverage boxes")
    for point in report["point_results"]:
        print(f"  {tuple(point['beta'])}: {point['rigorous_passes']}/"
              f"{point['attempts']}, max h={point['max_rigorous_h']}")
    for reason, count in sorted(report["failures"].items()):
        print(f"  FAIL x{count}: {reason}")


if __name__ == "__main__":
    main()
