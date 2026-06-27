#!/usr/bin/env python3
import csv
import json
import math
import pathlib
import sys


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    return 1


def read_json(path):
    with pathlib.Path(path).open() as f:
        return json.load(f)


def read_csv(path):
    with pathlib.Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def main(argv):
    if len(argv) < 3:
        return fail("usage: check_persistence_outputs.py <out-dir> <expected-status>")
    out = pathlib.Path(argv[1])
    expected = argv[2]
    summary_path = out / "captop_persistence_summary.json"
    if not summary_path.exists():
        return fail(f"missing {summary_path}")
    data = read_json(summary_path)
    if data.get("status") != expected:
        return fail(f"status {data.get('status')} != {expected}")
    for name in ["persistence_pairs.csv", "barcode_summary.csv", "betti_curve.csv", "diagram_dim0.csv", "diagram_dim1.csv", "diagram_dim2.csv", "captop_persistence_summary.json", "captop_persistence_report.txt"]:
        if not (out / name).exists():
            return fail(f"missing {name}")
    pairs = read_csv(out / "persistence_pairs.csv")
    for row in pairs:
        if row.get("death_type") not in {"finite", "infinite_death", "essential_unpaired"}:
            return fail("invalid death_type")
        if int(row.get("coefficient_field", "0")) < 2:
            return fail("invalid coefficient field in pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
