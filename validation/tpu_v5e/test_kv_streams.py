#!/usr/bin/env python3
"""Regression tests for DRAM-only KV predictions and comparison input guards."""
import csv
from pathlib import Path
import tempfile
import unittest

from kv_streams import MIB, kv_bytes, parse_probe, stream_time
from compare_stack import load_kv_times


class KVStreamsTest(unittest.TestCase):
    def test_geometry_counts_both_k_and_v(self):
        self.assertEqual(kv_bytes(8, 8192), 256 * MIB)

    def test_probe_uses_clock_ticks_not_rounded_display_time(self):
        output = 'tCK=0.625ns\nreqs=32768 ticks=2238 time=1.4us (read_bytes=1048576 write_bytes=0)'
        self.assertAlmostEqual(parse_probe(output, MIB), 1.39875)
        with self.assertRaises(ValueError):
            parse_probe(output, 2 * MIB)
        with self.assertRaises(ValueError):
            parse_probe(output.replace('write_bytes=0', 'write_bytes=32'), MIB)

    def test_large_stream_retains_startup_once(self):
        calls = []
        def probe(size):
            calls.append(size)
            return 0.5 + size / 750e3
        us, method = stream_time(4096 * MIB, 32 * MIB, probe)
        self.assertAlmostEqual(us, 0.5 + 4096 * MIB / 750e3)
        self.assertEqual(method, 'extrapolated')
        self.assertEqual(calls, [16 * MIB, 32 * MIB])
        self.assertEqual(stream_time(MIB, 32 * MIB, probe)[1], 'exact')

    def test_reject_rounding_and_bad_slope(self):
        for size in [0, -MIB, MIB + 32]:
            with self.assertRaises(ValueError):
                stream_time(size, 32 * MIB, lambda n: 1)
        with self.assertRaises(ValueError):
            stream_time(64 * MIB, 32 * MIB, lambda n: 1)

    def test_comparison_rejects_bad_geometry_duplicate_and_nan(self):
        good = dict(B=1, L=256, kv_bytes=MIB, kv_us=1.4, method='exact')
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'kv.csv'
            def write(rows):
                with p.open('w', newline='') as f:
                    writer = csv.DictWriter(f, list(good))
                    writer.writeheader()
                    writer.writerows(rows)
            write([good])
            self.assertEqual(load_kv_times(p), {(1, 256): (1.4, 'exact')})
            for rows in [[good, good], [dict(good, kv_bytes=2*MIB)], [dict(good, kv_us='nan')],
                         [dict(good, kv_us=0)], [dict(good, method='measured')]]:
                write(rows)
                with self.assertRaises(ValueError):
                    load_kv_times(p)


if __name__ == '__main__':
    unittest.main()
