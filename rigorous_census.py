"""Small JSONL census of fully rigorous signed tiles."""

import argparse
import json
import time

from certificate_result import CertificateGrade
from rigor_tile import fully_rigorous_signed_tile


DEFAULT_POINTS = (
    (5.978503016422594, 4.007534549834652, 1.6327649325136653),
    (0.50025, 1.0, 2.0),
)


def parse_point(value):
    parts = tuple(float(part) for part in value.split(","))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("point must be theta,phi,lambda")
    return parts


def run_point(beta, half_widths):
    theta, phi, lam = beta
    h_theta, h_phi, h_lam = half_widths
    started = time.time()
    result = fully_rigorous_signed_tile(
        theta - h_theta, theta + h_theta, phi, lam,
        hf=h_phi, hf3=h_lam, verbose=False,
    )
    return {
        "beta": list(beta),
        "hv": list(half_widths),
        "ok": result.ok,
        "grade": result.grade.name,
        "accepted_rigorous": result.accepted_at(CertificateGrade.RIGOROUS),
        "seconds": time.time() - started,
        "result": result.as_dict(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--point", action="append", type=parse_point)
    parser.add_argument("--h", type=float, action="append",
                        help="half-width rung; repeat to survey a ladder")
    parser.add_argument("--ledger", default="rigorous_census.jsonl")
    args = parser.parse_args()
    points = tuple(args.point) if args.point else DEFAULT_POINTS
    half_widths = tuple(args.h) if args.h else (2.5e-4,)
    with open(args.ledger, "a") as ledger:
        for h in half_widths:
            for point in points:
                record = run_point(point, (h, h, h))
                ledger.write(json.dumps(record) + "\n")
                ledger.flush()
                print(f"{tuple(round(x, 6) for x in point)} h={h:g}: "
                      f"{'OK' if record['ok'] else 'FAIL'}"
                      f"[{record['grade']}] {record['seconds']:.1f}s",
                      flush=True)


if __name__ == "__main__":
    main()
