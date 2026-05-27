#!/usr/bin/env python3

from __future__ import annotations

import math
import re
from pathlib import Path

RUST_DIR = Path("~/repos/semi-analytic-afterglow-rs").expanduser()
F90_DIR = Path("~/repos/semi-analytic-afterglow").expanduser()

RUST_OUT = RUST_DIR / "saa_out.dat"
F90_OUT = F90_DIR / "saa_out.dat"

RUST_SEDS = RUST_DIR / "saa_SEDs.dat"
F90_SEDS = F90_DIR / "saa_SEDs.dat"


FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

RAD_RE = re.compile(
    rf"Rad transfer complete for i_tm =\s*(\d+);\s*"
    rf"Gam_0 =\s*({FLOAT_RE});\s*"
    rf"t_obs =\s*({FLOAT_RE})"
)


def close(a: float, b: float, *, rtol: float = 1e-8, atol: float = 1e-12) -> bool:
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def parse_rad_transfer(path: Path) -> list[tuple[int, float, float]]:
    rows = []

    for line in path.read_text().splitlines():
        match = RAD_RE.search(line)
        if match:
            rows.append(
                (
                    int(match.group(1)),
                    float(match.group(2)),
                    float(match.group(3)),
                )
            )

    rows.sort(key=lambda row: row[0])
    return rows


def compare_rad_transfer() -> None:
    rust = parse_rad_transfer(RUST_OUT)
    f90 = parse_rad_transfer(F90_OUT)

    print("=== saa_out.dat Rad transfer comparison ===")
    print(f"Rust rows: {len(rust)}")
    print(f"F90 rows:  {len(f90)}")

    if len(rust) != len(f90):
        print("FAILED: row counts differ")
        return

    failures = []

    for r_row, f_row in zip(rust, f90):
        r_i, r_gam, r_tobs = r_row
        f_i, f_gam, f_tobs = f_row

        if r_i != f_i or not close(r_gam, f_gam) or not close(r_tobs, f_tobs):
            failures.append((r_row, f_row))

    if not failures:
        print("PASS: all Rad transfer values are close enough")
    else:
        print(f"FAILED: {len(failures)} mismatches")
        for rust_row, f90_row in failures[:20]:
            print(f"Rust: {rust_row}")
            print(f"F90:  {f90_row}")
            print()


def parse_seds(path: Path) -> tuple[list[float], list[list[float]]]:
    photon_energies: list[float] = []
    rows: list[list[float]] = []

    in_table = False

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("i_tm"):
            in_table = True
            continue

        parts = line.split()

        if not in_table:
            for part in parts:
                try:
                    photon_energies.append(float(part))
                except ValueError:
                    pass
            continue

        try:
            values = [float(part) for part in parts]
        except ValueError:
            continue

        if len(values) >= 6:
            rows.append(values)

    return photon_energies, rows


def compare_seds() -> None:
    rust_ph, rust_rows = parse_seds(RUST_SEDS)
    f90_ph, f90_rows = parse_seds(F90_SEDS)

    print()
    print("=== saa_SEDs.dat comparison ===")
    print(f"Rust photon energies: {len(rust_ph)}")
    print(f"F90 photon energies:  {len(f90_ph)}")
    print(f"Rust table rows:      {len(rust_rows)}")
    print(f"F90 table rows:       {len(f90_rows)}")

    if len(rust_ph) != len(f90_ph):
        print("FAILED: photon energy counts differ")
    else:
        max_ph_diff = max(abs(a - b) for a, b in zip(rust_ph, f90_ph))
        print(f"Max photon energy abs diff: {max_ph_diff:.6e}")

    if len(rust_rows) != len(f90_rows):
        print("FAILED: table row counts differ")
        return

    max_abs_diff = 0.0
    worst = None
    mismatch_count = 0

    for idx, (r_row, f_row) in enumerate(zip(rust_rows, f90_rows), start=1):
        width = min(len(r_row), len(f_row))

        for col in range(width):
            diff = abs(r_row[col] - f_row[col])

            if diff > max_abs_diff:
                max_abs_diff = diff
                worst = (idx, col + 1, r_row, f_row, diff)

            if not close(r_row[col], f_row[col], rtol=1e-6, atol=1e-10):
                mismatch_count += 1

    print(f"Max table abs diff: {max_abs_diff:.6e}")
    print(f"Mismatch count:     {mismatch_count}")

    if worst:
        row_idx, col_idx, r_row, f_row, diff = worst
        print()
        print("Worst mismatch:")
        print(f"  table row: {row_idx}")
        print(f"  column:    {col_idx}")
        print(f"  diff:      {diff:.6e}")
        print(f"  Rust:      {r_row}")
        print(f"  F90:       {f_row}")

    if mismatch_count == 0:
        print("PASS: all SED table values are close enough")
    else:
        print("FAILED: SED table values differ")


def main() -> None:
    for path in [RUST_OUT, F90_OUT, RUST_SEDS, F90_SEDS]:
        if not path.exists():
            raise FileNotFoundError(path)

    compare_rad_transfer()
    compare_seds()


if __name__ == "__main__":
    main()
