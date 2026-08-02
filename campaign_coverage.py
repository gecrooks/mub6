"""Grade-aware interval accounting for campaign ledgers."""

from dataclasses import dataclass

from certificate_result import CertificateGrade


@dataclass(frozen=True)
class LineCoverage:
    key: tuple[float, float]
    frontier: float
    merged_intervals: tuple[tuple[float, float], ...]
    gaps: tuple[tuple[float, float], ...]
    islands: tuple[tuple[float, float], ...]
    accepted_records: int
    rejected_records: int
    invalid_records: int
    duplicate_records: int

    @property
    def complete(self):
        return not self.gaps

    def as_dict(self):
        return {
            "phi": self.key[0],
            "lambda": self.key[1],
            "frontier": self.frontier,
            "complete": self.complete,
            "merged_intervals": [list(x) for x in self.merged_intervals],
            "gaps": [list(x) for x in self.gaps],
            "islands": [list(x) for x in self.islands],
            "accepted_records": self.accepted_records,
            "rejected_records": self.rejected_records,
            "invalid_records": self.invalid_records,
            "duplicate_records": self.duplicate_records,
        }


def grade_accepted(record, required=CertificateGrade.RIGOROUS):
    name = record.get("grade", CertificateGrade.EXPERIMENTAL.name)
    try:
        grade = CertificateGrade[name]
    except (KeyError, TypeError):
        return False
    return bool(record.get("ok")) and grade >= CertificateGrade(required)


def _record_interval(record):
    beta = record["beta"]
    hv = record["hv"]
    if len(beta) != 3 or len(hv) != 3:
        raise ValueError("beta and hv must have three components")
    center = float(beta[0])
    radius = float(hv[0])
    phi = round(float(beta[1]), 9)
    lam = round(float(beta[2]), 9)
    if radius < 0:
        raise ValueError("negative tile half-width")
    values = (center, radius, phi, lam)
    if not all(value == value and abs(value) != float("inf")
               for value in values):
        raise ValueError("non-finite tile record")
    return (phi, lam), (center - radius, center + radius)


def _merge(intervals, atol):
    merged = []
    for lo, hi in sorted(intervals):
        if not merged or lo > merged[-1][1] + atol:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    return tuple((float(lo), float(hi)) for lo, hi in merged)


def analyze_ledger_records(records, theta_lo, theta_hi,
                           required_grade=CertificateGrade.RIGOROUS,
                           *, atol=1e-12):
    """Return per-line coverage connected to ``theta_lo``.

    Accepted tiles beyond a gap are retained as disconnected islands but
    cannot advance ``frontier``. Records may arrive in any order.
    """
    if theta_hi < theta_lo:
        raise ValueError("theta_hi must be at least theta_lo")
    grouped = {}
    rejected = {}
    invalid = {}
    seen = {}
    duplicates = {}
    for record in records:
        try:
            key, interval = _record_interval(record)
        except (KeyError, TypeError, ValueError, OverflowError):
            raw_beta = record.get("beta", [None, None, None]) \
                if isinstance(record, dict) else [None, None, None]
            try:
                key = (round(float(raw_beta[1]), 9),
                       round(float(raw_beta[2]), 9))
            except (TypeError, ValueError, IndexError):
                key = (float("nan"), float("nan"))
            invalid[key] = invalid.get(key, 0) + 1
            continue
        if not grade_accepted(record, required_grade):
            rejected[key] = rejected.get(key, 0) + 1
            continue
        token = (round(interval[0], 15), round(interval[1], 15))
        key_seen = seen.setdefault(key, set())
        if token in key_seen:
            duplicates[key] = duplicates.get(key, 0) + 1
        else:
            key_seen.add(token)
            grouped.setdefault(key, []).append(interval)

    keys = set(grouped) | set(rejected) | set(invalid)
    reports = {}
    for key in keys:
        merged = _merge(grouped.get(key, ()), atol)
        frontier = float(theta_lo)
        for lo, hi in merged:
            if lo <= frontier + atol and hi >= frontier - atol:
                frontier = max(frontier, hi)
            elif lo > frontier + atol:
                break
        frontier = min(frontier, float(theta_hi))

        clipped = []
        for lo, hi in merged:
            lo_c, hi_c = max(lo, theta_lo), min(hi, theta_hi)
            if hi_c >= lo_c - atol:
                clipped.append((lo_c, hi_c))
        gaps = []
        cursor = float(theta_lo)
        for lo, hi in clipped:
            if lo > cursor + atol:
                gaps.append((cursor, lo))
            cursor = max(cursor, hi)
        if cursor < theta_hi - atol:
            gaps.append((cursor, float(theta_hi)))
        islands = tuple((lo, hi) for lo, hi in merged
                        if lo > frontier + atol)
        reports[key] = LineCoverage(
            key=key,
            frontier=frontier,
            merged_intervals=merged,
            gaps=tuple(gaps),
            islands=islands,
            accepted_records=len(grouped.get(key, ())),
            rejected_records=rejected.get(key, 0),
            invalid_records=invalid.get(key, 0),
            duplicate_records=duplicates.get(key, 0),
        )
    return reports
