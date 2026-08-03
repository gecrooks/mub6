import unittest

from ledger_bits import (
    LEDGER_SCHEMA,
    bits_float,
    box_record,
    decode_box_record,
    float_bits,
)


class LedgerBitsTests(unittest.TestCase):
    def test_binary64_round_trip_preserves_signed_zero(self):
        for value in (0.0, -0.0, 0.1, -3.25, 1e-300):
            self.assertEqual(float_bits(bits_float(float_bits(value))),
                             float_bits(value))
        self.assertNotEqual(float_bits(0.0), float_bits(-0.0))

    def test_box_record_includes_exact_interval_bits(self):
        record = box_record({"ok": True}, (0.1, 1.0, 2.0),
                            (0.03, 0.04, 0.05))
        beta, widths, interval, token = decode_box_record(record)
        self.assertEqual(record["ledger_schema"], LEDGER_SCHEMA)
        self.assertEqual(beta, (0.1, 1.0, 2.0))
        self.assertEqual(widths, (0.03, 0.04, 0.05))
        self.assertEqual(token, tuple(record["theta_interval_bits"]))
        self.assertEqual(interval, (0.1 - 0.03, 0.1 + 0.03))

    def test_decimal_mirror_or_endpoint_tampering_is_rejected(self):
        record = box_record({}, (0.1, 1.0, 2.0), (0.03, 0.04, 0.05))
        record["beta"][0] = 0.2
        with self.assertRaisesRegex(ValueError, "mirror"):
            decode_box_record(record)
        record = box_record({}, (0.1, 1.0, 2.0), (0.03, 0.04, 0.05))
        record["theta_interval_bits"][0] = float_bits(0.0)
        with self.assertRaisesRegex(ValueError, "interval bits"):
            decode_box_record(record)

    def test_noncanonical_hex_and_legacy_rigorous_input_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "lowercase"):
            bits_float("0X0000000000000000")
        with self.assertRaisesRegex(ValueError, "requires binary64"):
            decode_box_record({"beta": [0.0] * 3, "hv": [0.1] * 3})
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            box_record({}, (0.0, 0.0, 0.0), (-0.0, 0.1, 0.1))


if __name__ == "__main__":
    unittest.main()
