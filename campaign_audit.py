"""Audit a campaign JSONL ledger for graded, contiguous line coverage."""

import argparse
import json

from campaign_coverage import analyze_ledger_records
from certificate_result import CertificateGrade


def read_ledger(path):
    records = []
    malformed = 0
    with open(path) as stream:
        for line in stream:
            try:
                value = json.loads(line)
            except json.JSONDecodeError:
                malformed += 1
                continue
            if isinstance(value, dict):
                records.append(value)
            else:
                malformed += 1
    return records, malformed


def audit(path, theta_lo, theta_hi, required_grade):
    records, malformed = read_ledger(path)
    reports = analyze_ledger_records(
        records, theta_lo, theta_hi, required_grade
    )
    lines = [report.as_dict() for _, report in sorted(reports.items())]
    return {
        "ledger": str(path),
        "required_grade": CertificateGrade(required_grade).name,
        "theta_domain": [theta_lo, theta_hi],
        "records": len(records),
        "malformed_records": malformed,
        "line_count": len(lines),
        "complete_lines": sum(line["complete"] for line in lines),
        "lines_with_islands": sum(bool(line["islands"]) for line in lines),
        "lines": lines,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger")
    parser.add_argument("--theta-lo", type=float, default=0.0)
    parser.add_argument("--theta-hi", type=float, default=1.5707963267948966)
    parser.add_argument("--required-grade",
                        choices=[grade.name for grade in CertificateGrade],
                        default=CertificateGrade.RIGOROUS.name)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = audit(
        args.ledger, args.theta_lo, args.theta_hi,
        CertificateGrade[args.required_grade],
    )
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    print(f"{result['complete_lines']}/{result['line_count']} lines complete "
          f"at {result['required_grade']}; "
          f"{result['lines_with_islands']} with disconnected islands; "
          f"{result['malformed_records']} malformed records")
    for line in result["lines"]:
        status = "COMPLETE" if line["complete"] else \
            f"frontier={line['frontier']:.9g}, gaps={len(line['gaps'])}"
        print(f"  ({line['phi']:.9g}, {line['lambda']:.9g}): {status}; "
              f"islands={len(line['islands'])}, "
              f"rejected={line['rejected_records']}, "
              f"duplicates={line['duplicate_records']}")


if __name__ == "__main__":
    main()
