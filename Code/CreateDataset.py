"""
generate_dataset.py

Generates a large CSV of random integer sets drawn from a universe [1..u].
Each row is one set; elements are stored as a sorted, space-separated list.

Output: dataset.csv
  Columns:
    id        – unique row id
    u         – universe size  (bitvector has positions 1..u)
    size      – number of elements in the set |S|
    elements  – sorted space-separated element list, e.g. "4 5 8 10"
"""

import csv
import random

SEED = 42
random.seed(SEED)

OUTPUT_FILE = "Bachelor/Code/dataset.csv"

# ── (universe_size, set_size_range, num_samples) ──────────────────────────────
CONFIGS = [
    (16,   (1,  8),    500),
    (32,   (1,  16),   500),
    (64,   (1,  32),   1_000),
    (128,  (2,  64),   1_000),
    (256,  (2,  128),  1_000),
    (512,  (1,  50),   1_000),
    (512,  (100, 256), 500),
    (1024, (1,  100),  2_000),
    (1024, (200, 512), 500),
    (4096, (1,  200),  2_000),
    (4096, (500, 2048),500),
]

EDGE_CASES = [
    (16,  []),
    (16,  [8]),
    (16,  list(range(1, 17))),
    (16,  list(range(1, 17, 2))),
    (32,  list(range(2, 33, 2))),
    (10,  [4, 5, 8, 10]),           # the user's example
    (20,  [1, 3, 7, 12, 15, 19]),
    (100, list(range(1, 101, 10))),
]


def random_set(u: int, lo: int, hi: int) -> list[int]:
    size = random.randint(lo, min(hi, u))
    return sorted(random.sample(range(1, u + 1), size))


def main():
    rows = []
    uid  = 0

    for u, (lo, hi), n_samples in CONFIGS:
        for _ in range(n_samples):
            elems = random_set(u, lo, hi)
            rows.append({"id": uid, "u": u,
                         "size": len(elems),
                         "elements": " ".join(map(str, elems))})
            uid += 1

    for u, elems in EDGE_CASES:
        rows.append({"id": uid, "u": u,
                     "size": len(elems),
                     "elements": " ".join(map(str, elems))})
        uid += 1

    fieldnames = ["id", "u", "size", "elements"]
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written {len(rows):,} rows to '{OUTPUT_FILE}'")

    example = next(r for r in rows if r["elements"] == "4 5 8 10")
    print(f"\nSpot-check – user's example:")
    print(f"  id={example['id']}  u={example['u']}  "
          f"elements=[{example['elements']}]")


if __name__ == "__main__":
    main()