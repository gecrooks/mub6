"""Adaptive campaign driver (CPU reference): cover the fundamental
domain of the Karlsson parameter torus with certified tiles, chaining
where the territory is regular and dropping to standalone anisotropic
tiles across valley strata.

Fundamental domain (Result 26, order-32 quotient of [0,2pi)^3):
period pi in every axis reduces to [0,pi)^3; the point group
{id, (pi-t,pi-p,pi-l), (t,pi-p,pi-l), (pi-t,p,l)} then gives
    F = [0, pi/2] x [0, pi/2] x [0, pi].
Closed tiles may overhang F's boundary; over-coverage into mirror
copies is sound (the union still covers F), so no half-open
bookkeeping is needed for correctness — only for cost accounting.

Scheduler: theta-chain-lines on a (phi, lam) grid with transverse
spacing 1.8*h_base (tiles half-width >= h_base overlap between lines).
Along each line: margin-cached chain steps; on ANY chain failure
(valley emerged, tube failure, stuck patch) the driver switches to
standalone certified tiles with an ADAPTIVE h-vector (Result 27:
fatten the axis least coupled to the flagged roots, retry ladder
down to isotropic), then re-anchors a chain past the stratum.
Coverage: only-fatten policy (hv >= h_base in every axis) keeps the
h_base line grid sufficient in the transverse axes; theta progress is
tracked per tile as [theta - hv_t, theta + hv_t].

Ledger: append-only JSONL, one record per certified tile / failure;
resume skips lines whose certified theta-span already covers the
domain segment. PROTOTYPE rigor model as everywhere.
"""

import argparse
import json
import os
import time
import warnings

import numpy as np

from cache import anchored_tile, chain_step
from campaign_coverage import analyze_ledger_records, grade_accepted
from certificate_result import CertificateGrade
from karlsson import karlsson_map
from mub import find_mu_vectors
from parametric import _g_and_J, certify_tile, dg_dbeta, polish_root

warnings.filterwarnings("ignore")

PI = np.pi
DOMAIN = ((0.0, PI / 2), (0.0, PI / 2), (0.0, PI))
FATTEN_LADDER = (3.0, 2.0, 1.5, 1.0)


def coupling_profile(beta, sig_cut=0.06):
    """Per-axis max |u . dg/dbeta| over near-singular roots at beta
    (the anisotropy diagnostics of Result 27). Returns (profile, n_flag);
    profile is None when no root is near-singular."""
    H0 = karlsson_map(*beta)
    vecs = find_mu_vectors([H0], n_starts=4000, seed=99)
    roots = [polish_root(H0, np.angle(v * np.sqrt(6))[1:]) for v in vecs]
    prof = np.zeros(3)
    n_flag = 0
    for th in roots:
        _, J = _g_and_J(H0, th)
        U, sv, _ = np.linalg.svd(J)
        if sv[-1] < sig_cut:
            n_flag += 1
            s = np.abs(U[:, -1] @ dg_dbeta(beta, th))
            prof = np.maximum(prof, s)
    return (prof if n_flag else None), n_flag


def adaptive_tile(beta, h_base, ledger, verbose=True,
                  required_grade=CertificateGrade.RIGOROUS):
    """Standalone certified tile with the fatten-weakest-axis ladder.
    Returns hv on success, None on failure (after the full ladder)."""
    prof, n_flag = coupling_profile(beta)
    axis = int(np.argmin(prof)) if prof is not None else 0
    for fac in FATTEN_LADDER:
        hv = np.full(3, h_base)
        hv[axis] = fac * h_base
        r = certify_tile(beta, hv, verbose=False, use_certified=True)
        result = r.get("result")
        grade = (CertificateGrade.EXPERIMENTAL if result is None
                 else result.grade)
        rec = dict(mode="tile", beta=list(map(float, beta)),
                   hv=[float(x) for x in hv], ok=bool(r["ok"]),
                   grade=grade.name,
                   seconds=float(r["seconds"]), n_flag=n_flag,
                   reason=(None if r["ok"] else
                           str(r.get("reason", "certificate failed"))),
                   evidence=(None if result is None else
                             result.as_dict()["dependencies"]))
        ledger.write(json.dumps(rec) + "\n")
        ledger.flush()
        if grade_accepted(rec, required_grade):
            if verbose:
                print(f"    tile th={beta[0]:.6f} hv=({hv[0]:g},{hv[1]:g},"
                      f"{hv[2]:g}) OK [{r['seconds']:.0f} s]", flush=True)
            return hv
        if r["ok"] and verbose:
            print(f"    tile th={beta[0]:.6f} rejected by grade policy: "
                  f"{grade.name} < {CertificateGrade(required_grade).name}",
                  flush=True)
        elif verbose:
            print(f"    tile th={beta[0]:.6f} fac={fac:g} FAIL "
                  f"({r.get('reason', 'certificate failed')})", flush=True)
    return None


def run_line(phi, lam, th_lo, th_hi, h_base, ledger, start_at=None,
             verbose=True, required_grade=CertificateGrade.RIGOROUS):
    """Certify one theta-line. Returns (covered_to, n_tiles) — covered_to
    is the certified theta frontier (>= th_hi on full success)."""
    th = th_lo + h_base if start_at is None else start_at
    n_tiles = 0
    state = None
    while th - h_base < th_hi:
        if state is not None:
            res = chain_step(state, (th, phi, lam), verbose=verbose)
            if res.get("ok"):
                ledger.write(json.dumps(dict(
                    mode="chain", beta=[float(th), float(phi), float(lam)],
                    hv=[h_base] * 3, ok=True,
                    grade=CertificateGrade.EXPERIMENTAL.name,
                    seconds=float(res["seconds"]))) + "\n")
                ledger.flush()
                n_tiles += 1
                th += 1.6 * h_base
                continue
            if verbose:
                print(f"    chain broke at th={th:.6f} "
                      f"({res.get('reason')}) -> adaptive", flush=True)
            state = None
        elif required_grade <= CertificateGrade.EXPERIMENTAL:
            # try to open a chain at th; valley territory raises
            try:
                state = anchored_tile((th, phi, lam), h_base,
                                      verbose=verbose, use_certified=True)
                ledger.write(json.dumps(dict(
                    mode="anchor", beta=[float(th), float(phi), float(lam)],
                    hv=[h_base] * 3, ok=True,
                    grade=CertificateGrade.EXPERIMENTAL.name,
                    seconds=float(state["anchor_seconds"]))) + "\n")
                ledger.flush()
                n_tiles += 1
                th += 1.6 * h_base
                continue
            except RuntimeError as e:
                if verbose:
                    print(f"    anchor th={th:.6f}: {e} -> adaptive",
                          flush=True)
        # standalone adaptive tile across the stratum
        hv = adaptive_tile((th, phi, lam), h_base, ledger, verbose=verbose,
                           required_grade=required_grade)
        if hv is None:
            return th - h_base, n_tiles          # frontier stalled
        n_tiles += 1
        th += 1.6 * hv[0]
    return th_hi, n_tiles


def load_frontiers(path, required_grade=CertificateGrade.RIGOROUS,
                   theta_lo=0.0, theta_hi=PI / 2):
    """Resume support: per-(phi,lam) certified theta frontier from the
    ledger (max certified theta+hv_t of contiguously-OK records; the
    conservative re-check is one overlapping tile at resume)."""
    if not os.path.exists(path):
        return {}
    records = []
    with open(path) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                records.append({"invalid_json": line.rstrip("\n")})
    reports = analyze_ledger_records(
        records, theta_lo, theta_hi, required_grade
    )
    return {key: report.frontier for key, report in reports.items()
            if report.frontier > theta_lo}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=3e-4)
    ap.add_argument("--phi0", type=float, default=None)
    ap.add_argument("--lam0", type=float, default=None)
    ap.add_argument("--th-lo", type=float, default=None)
    ap.add_argument("--th-hi", type=float, default=None)
    ap.add_argument("--n-lines", type=int, default=1)
    ap.add_argument("--ledger", default="campaign_ledger.jsonl")
    ap.add_argument("--required-grade",
                    choices=[grade.name for grade in CertificateGrade],
                    default=CertificateGrade.RIGOROUS.name)
    args = ap.parse_args()

    h = args.h
    spacing = 1.8 * h
    phi0 = args.phi0 if args.phi0 is not None else DOMAIN[1][0] + h
    lam0 = args.lam0 if args.lam0 is not None else DOMAIN[2][0] + h
    th_lo = args.th_lo if args.th_lo is not None else DOMAIN[0][0]
    th_hi = args.th_hi if args.th_hi is not None else DOMAIN[0][1]

    required_grade = CertificateGrade[args.required_grade]
    frontiers = load_frontiers(args.ledger, required_grade, th_lo, th_hi)
    t0 = time.time()
    total = 0
    with open(args.ledger, "a") as ledger:
        for il in range(args.n_lines):
            phi = phi0 + (il % 1) * spacing        # v1: vary lam only
            lam = lam0 + il * spacing
            key = (round(phi, 9), round(lam, 9))
            start = frontiers.get(key)
            if start is not None and start >= th_hi:
                print(f"line (phi={phi:.6f}, lam={lam:.6f}): already done")
                continue
            print(f"line (phi={phi:.6f}, lam={lam:.6f}) "
                  f"th in [{th_lo:.6f}, {th_hi:.6f}]"
                  + (f" resume at {start:.6f}" if start else ""),
                  flush=True)
            covered, n = run_line(phi, lam, th_lo, th_hi, h, ledger,
                                  start_at=start,
                                  required_grade=required_grade)
            total += n
            status = "COMPLETE" if covered >= th_hi else \
                f"STALLED at {covered:.6f}"
            print(f"  line {status}: {n} tiles", flush=True)
    print(f"\n{total} tiles, {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
