"""Grade-aware interval accounting for campaign ledgers."""

from dataclasses import dataclass

from certificate_result import CertificateGrade
from ledger_bits import decode_box_record, float_bits


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


def _record_interval(record, *, require_bits):
    beta, _widths, interval, token = decode_box_record(
        record, require_bits=require_bits
    )
    key = (beta[1], beta[2])
    key_bits = (float_bits(beta[1]), float_bits(beta[2]))
    return key, key_bits, interval, token


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
                           *, atol=None):
    """Return per-line coverage connected to ``theta_lo``.

    Accepted tiles beyond a gap are retained as disconnected islands but
    cannot advance ``frontier``. Records may arrive in any order.
    """
    if theta_hi < theta_lo:
        raise ValueError("theta_hi must be at least theta_lo")
    rigorous = CertificateGrade(required_grade) >= CertificateGrade.RIGOROUS
    if rigorous and atol not in (None, 0, 0.0):
        raise ValueError("rigorous resume forbids gap-bridging tolerance")
    atol = 0.0 if rigorous else (1e-12 if atol is None else float(atol))
    grouped = {}
    rejected = {}
    invalid = {}
    seen = {}
    duplicates = {}
    display_keys = {}
    for record in records:
        accepted = grade_accepted(record, required_grade)
        try:
            key, key_bits, interval, token = _record_interval(
                record,
                require_bits=(accepted and rigorous),
            )
        except (KeyError, TypeError, ValueError, OverflowError):
            raw_beta = record.get("beta", [None, None, None]) \
                if isinstance(record, dict) else [None, None, None]
            try:
                key = (float(raw_beta[1]), float(raw_beta[2]))
                key_bits = (float_bits(key[0]), float_bits(key[1]))
                display_keys[key_bits] = key
            except (TypeError, ValueError, IndexError, OverflowError):
                key = (float("nan"), float("nan"))
                key_bits = ("invalid", "invalid")
                display_keys[key_bits] = key
            invalid[key_bits] = invalid.get(key_bits, 0) + 1
            continue
        display_keys[key_bits] = key
        if not accepted:
            rejected[key_bits] = rejected.get(key_bits, 0) + 1
            continue
        key_seen = seen.setdefault(key_bits, set())
        if token in key_seen:
            duplicates[key_bits] = duplicates.get(key_bits, 0) + 1
        else:
            key_seen.add(token)
            grouped.setdefault(key_bits, []).append(interval)

    keys = set(grouped) | set(rejected) | set(invalid)
    reports = {}
    for key_bits in keys:
        key = display_keys.get(key_bits, key_bits)
        merged = _merge(grouped.get(key_bits, ()), atol)
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
            accepted_records=len(grouped.get(key_bits, ())),
            rejected_records=rejected.get(key_bits, 0),
            invalid_records=invalid.get(key_bits, 0),
            duplicate_records=duplicates.get(key_bits, 0),
        )
    return reports
