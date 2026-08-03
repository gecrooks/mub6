"""Canonical IEEE-754 binary64 encoding for proof-campaign ledgers."""

import math
import re
import struct


LEDGER_SCHEMA = "mub6-ledger-v2-binary64"
_BITS = re.compile(r"^[0-9a-f]{16}$")


def float_bits(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("ledger floats must be finite")
    bits = struct.unpack(">Q", struct.pack(">d", value))[0]
    return f"{bits:016x}"


def bits_float(value):
    if not isinstance(value, str) or _BITS.fullmatch(value) is None:
        raise ValueError("binary64 value must be 16 lowercase hex digits")
    result = struct.unpack(">d", struct.pack(">Q", int(value, 16)))[0]
    if not math.isfinite(result):
        raise ValueError("ledger floats must be finite")
    return result


def _vector_bits(values, name):
    values = tuple(float(value) for value in values)
    if len(values) != 3:
        raise ValueError(f"{name} must have three components")
    return values, [float_bits(value) for value in values]


def box_record(record, beta, half_widths):
    """Return a record with exact bits and matching decimal display mirrors."""
    beta, beta_bits = _vector_bits(beta, "beta")
    widths, width_bits = _vector_bits(half_widths, "hv")
    if any(width < 0 or float_bits(width) == "8000000000000000"
           for width in widths):
        raise ValueError("half-widths must be nonnegative")
    theta_interval = (beta[0] - widths[0], beta[0] + widths[0])
    return {
        **record,
        "ledger_schema": LEDGER_SCHEMA,
        "beta": list(beta),
        "beta_bits": beta_bits,
        "hv": list(widths),
        "hv_bits": width_bits,
        "theta_interval_bits": [float_bits(value)
                                for value in theta_interval],
    }


def decode_box_record(record, *, require_bits=True):
    """Decode and cross-check a ledger box.

    Returns ``(beta, half_widths, theta_interval, interval_bit_token)``.
    Legacy decimal-only records are accepted only when ``require_bits=False``.
    """
    if not isinstance(record, dict):
        raise ValueError("ledger record must be an object")
    if record.get("ledger_schema") != LEDGER_SCHEMA:
        if require_bits:
            raise ValueError("rigorous resume requires binary64 ledger schema")
        beta = tuple(float(value) for value in record["beta"])
        widths = tuple(float(value) for value in record["hv"])
        if len(beta) != 3 or len(widths) != 3:
            raise ValueError("beta and hv must have three components")
        if not all(math.isfinite(value) for value in beta + widths) \
                or any(width < 0 or
                       float_bits(width) == "8000000000000000"
                       for width in widths):
            raise ValueError("invalid legacy box")
        interval = (beta[0] - widths[0], beta[0] + widths[0])
        return beta, widths, interval, tuple(float_bits(x) for x in interval)

    beta_bits = record.get("beta_bits")
    width_bits = record.get("hv_bits")
    interval_bits = record.get("theta_interval_bits")
    if not all(isinstance(values, list) for values in
               (beta_bits, width_bits, interval_bits)) \
            or len(beta_bits) != 3 or len(width_bits) != 3 \
            or len(interval_bits) != 2:
        raise ValueError("incomplete binary64 box fields")
    beta = tuple(bits_float(value) for value in beta_bits)
    widths = tuple(bits_float(value) for value in width_bits)
    interval = tuple(bits_float(value) for value in interval_bits)
    if any(width < 0 or float_bits(width) == "8000000000000000"
           for width in widths):
        raise ValueError("half-widths must be nonnegative")
    if [float_bits(value) for value in record.get("beta", ())] != beta_bits:
        raise ValueError("beta decimal mirror does not match beta_bits")
    if [float_bits(value) for value in record.get("hv", ())] != width_bits:
        raise ValueError("hv decimal mirror does not match hv_bits")
    computed = [float_bits(beta[0] - widths[0]),
                float_bits(beta[0] + widths[0])]
    if computed != interval_bits:
        raise ValueError("theta interval bits do not match beta/hv bits")
    return beta, widths, interval, tuple(interval_bits)
