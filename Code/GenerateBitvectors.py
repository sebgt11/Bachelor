"""
generate_bitvectors.py

Reads dataset.csv and writes bitvectors.csv, converting each set of integers
to a bitvector string of length u.

Example: elements "4 5 6 9 10" with u=10  →  bitvector "0001110011"
         position:                              1234567890

Output: bitvectors.csv
  Columns:
    id        – matches dataset.csv
    u         – universe size
    bitvector – string of 0/1 of length u
"""

import csv

INPUT_FILE = "Bachelor/Code/dataset.csv"
OUTPUT_FILE = "Bachelor/Code/bitvectors.csv"


def set_to_bitvector(elements: list[int], u: int) -> str:
    bv = ['0'] * u
    for e in elements:
        bv[e - 1] = '1'
    return "".join(bv)


def main():
    with open(INPUT_FILE, newline="") as f:
        rows = list(csv.DictReader(f))

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "u", "bitvector"])
        writer.writeheader()
        for row in rows:
            u        = int(row["u"])
            elements = [int(x) for x in row["elements"].split()] if row["elements"].strip() else []
            writer.writerow({
                "id":        row["id"],
                "u":         u,
                "bitvector": set_to_bitvector(elements, u),
            })

    print(f"Written {len(rows):,} rows to '{OUTPUT_FILE}'")

    # Spot-check
    print("\nSpot-checks:")
    checks = [
        ([4, 5, 6, 9, 10], 10, "0001110011"),
        ([4, 5, 8, 10],    10, "0001100101"),
    ]
    for elems, u, expected in checks:
        got = set_to_bitvector(elems, u)
        print(f"  {elems} → '{got}'  {'✓' if got == expected else f'✗ expected {expected}'}")


if __name__ == "__main__":
    main()