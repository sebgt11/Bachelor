import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

CSV_PATH = "results.csv"
OUT_DIR = "figures"
SUMMARY_PATH = "summary_table.csv"

B_FIXED = 16            # block size used for everything except the block-size study
COL_BASE = "#7f7f7f"
COL_CEB = "#1f77b4"
COL_NEST = "#1f77b4"
COL_INDEP = "#d62728"

FLOAT_FIELDS = ["baseline_bpe", "ceb_bpe", "H_bpe", "L_bpe", "saving_pct",
                "baseline_total_bpe", "ceb_total_bpe", "total_saving_pct",
                "ceb_build_ms", "base_build_ms", "mean_depth",
                "ceb_retrieve_us", "ceb_rank_us",
                "base_retrieve_us", "base_rank_us"]
INT_FIELDS = ["b", "contraction", "max_depth", "universe"]


def load_rows(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for k in FLOAT_FIELDS:
                r[k] = float(r[k])
            for k in INT_FIELDS:
                r[k] = int(r[k])
            rows.append(r)
    return rows


def sel(rows, **cond):
    out = []
    for r in rows:
        if all(r[k] == v for k, v in cond.items()):
            out.append(r)
    return out


def agg(rows, field):
    vals = [r[field] for r in rows]
    return float(np.mean(vals)), float(np.std(vals))


def save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


# RQ1  space against the entropy bounds, favourable and pathological data
def plot_space_comparison(rows):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    datasets = ["nested", "independent"]
    titles = {"nested": "Nested (favourable)", "independent": "Independent (unfavourable)"}
    series = [("H_bpe", "Independent entropy H", "#bbbbbb"),
              ("baseline_bpe", "Baseline (measured)", COL_BASE),
              ("L_bpe", "Containment entropy L", "#9ecae1"),
              ("ceb_bpe", "Containment (measured)", COL_CEB)]
    x = np.arange(len(datasets))
    w = 0.2
    for i, (field, label, color) in enumerate(series):
        means, stds = [], []
        for d in datasets:
            m, s = agg(sel(rows, dataset=d, b=B_FIXED, contraction=1), field)
            means.append(m); stds.append(s)
        ax.bar(x + (i - 1.5) * w, means, w, yerr=stds, capsize=3,
               label=label, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels([titles[d] for d in datasets])
    ax.set_ylabel("Bits per element")
    ax.set_title("Space cost versus information-theoretic bounds (b = %d)" % B_FIXED)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    save(fig, "space_comparison.png")


# RQ1b  encoding size against total size including the query index, both methods
def plot_space_with_index(rows):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    datasets = ["nested", "independent"]
    titles = {"nested": "Nested (favourable)", "independent": "Independent (unfavourable)"}
    series = [("baseline_bpe", "Baseline, encoding", COL_BASE, None),
              ("baseline_total_bpe", "Baseline, with index", COL_BASE, "//"),
              ("ceb_bpe", "Containment, encoding", COL_CEB, None),
              ("ceb_total_bpe", "Containment, with index", COL_CEB, "//")]
    x = np.arange(len(datasets))
    w = 0.2
    for i, (field, label, color, hatch) in enumerate(series):
        means, stds = [], []
        for d in datasets:
            m, s = agg(sel(rows, dataset=d, b=B_FIXED, contraction=1), field)
            means.append(m); stds.append(s)
        ax.bar(x + (i - 1.5) * w, means, w, yerr=stds, capsize=3,
               label=label, color=color, hatch=hatch, edgecolor="white")
    ax.set_xticks(x)
    ax.set_xticklabels([titles[d] for d in datasets])
    ax.set_ylabel("Bits per element")
    ax.set_title("Encoding size versus total size with the query index (b = %d)" % B_FIXED)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    save(fig, "space_with_index.png")


# RQ2  contraction shortens the hierarchy and the queries
def plot_contraction_effect(rows):
    fig, (axd, axt) = plt.subplots(1, 2, figsize=(10, 4.5))
    labels = ["Without\ncontraction", "With\ncontraction"]
    x = np.arange(2)
    w = 0.35

    for col, field, color, lbl in [(0, "max_depth", "#9ecae1", "max depth"),
                                   (1, "mean_depth", "#1f77b4", "mean depth")]:
        means, stds = [], []
        for con in (0, 1):
            m, s = agg(sel(rows, dataset="nested", b=B_FIXED, contraction=con), field)
            means.append(m); stds.append(s)
        axd.bar(x + (col - 0.5) * w, means, w, yerr=stds, capsize=3,
                color=color, label=lbl)
    axd.set_xticks(x); axd.set_xticklabels(labels)
    axd.set_ylabel("Hierarchy depth (edges to root)")
    axd.set_title("Depth")
    axd.legend(); axd.grid(axis="y", alpha=0.3)

    for col, field, color, lbl in [(0, "ceb_retrieve_us", "#9ecae1", "retrieve"),
                                   (1, "ceb_rank_us", "#1f77b4", "rank")]:
        means, stds = [], []
        for con in (0, 1):
            m, s = agg(sel(rows, dataset="nested", b=B_FIXED, contraction=con), field)
            means.append(m); stds.append(s)
        axt.bar(x + (col - 0.5) * w, means, w, yerr=stds, capsize=3,
                color=color, label=lbl)
    axt.set_xticks(x); axt.set_xticklabels(labels)
    axt.set_ylabel("Time per query (microseconds)")
    axt.set_title("Query time")
    axt.legend(); axt.grid(axis="y", alpha=0.3)

    fig.suptitle("Effect of parent contraction (nested data, b = %d)" % B_FIXED)
    save(fig, "contraction_effect.png")


# RQ3  query time of the containment structure against the baseline
def plot_query_time_comparison(rows):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ops = [("retrieve", "ceb_retrieve_us", "base_retrieve_us"),
           ("rank", "ceb_rank_us", "base_rank_us")]
    x = np.arange(len(ops))
    w = 0.35
    base_m, base_s, ceb_m, ceb_s = [], [], [], []
    sub = sel(rows, dataset="nested", b=B_FIXED, contraction=1)
    for _, ceb_f, base_f in ops:
        bm, bs = agg(sub, base_f); cm, cs = agg(sub, ceb_f)
        base_m.append(bm); base_s.append(bs); ceb_m.append(cm); ceb_s.append(cs)
    ax.bar(x - w / 2, base_m, w, yerr=base_s, capsize=4, label="Baseline", color=COL_BASE)
    ax.bar(x + w / 2, ceb_m, w, yerr=ceb_s, capsize=4, label="Containment", color=COL_CEB)
    ax.set_xticks(x); ax.set_xticklabels([o[0] for o in ops])
    ax.set_ylabel("Time per query (microseconds)")
    ax.set_title("Query time, containment versus baseline (nested, b = %d)" % B_FIXED)
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    save(fig, "query_time_comparison.png")


# RQ3  build time of the two solutions
def plot_build_time(rows):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    datasets = ["nested", "independent"]
    titles = {"nested": "Nested (favourable)", "independent": "Independent (unfavourable)"}
    x = np.arange(len(datasets))
    w = 0.35
    base_m, base_s, ceb_m, ceb_s = [], [], [], []
    for d in datasets:
        sub = sel(rows, dataset=d, b=B_FIXED, contraction=1)
        bm, bs = agg(sub, "base_build_ms")
        cm, cs = agg(sub, "ceb_build_ms")
        base_m.append(bm); base_s.append(bs); ceb_m.append(cm); ceb_s.append(cs)
    ax.bar(x - w / 2, base_m, w, yerr=base_s, capsize=4, label="Baseline", color=COL_BASE)
    ax.bar(x + w / 2, ceb_m, w, yerr=ceb_s, capsize=4, label="Containment", color=COL_CEB)
    ax.set_xticks(x)
    ax.set_xticklabels([titles[d] for d in datasets])
    ax.set_ylabel("Build time (milliseconds)")
    ax.set_title("Construction time, containment versus baseline (b = %d)" % B_FIXED)
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    save(fig, "build_time.png")


# Block size study  one figure two panels
def plot_block_size(rows):
    bs = sorted({r["b"] for r in rows})
    fig, (axs, axt) = plt.subplots(1, 2, figsize=(11, 4.5))

    for field, color, lbl, mk in [("baseline_bpe", COL_BASE, "Baseline", "s"),
                                  ("ceb_bpe", COL_CEB, "Containment", "o")]:
        m = [agg(sel(rows, dataset="nested", b=b, contraction=1), field)[0] for b in bs]
        s = [agg(sel(rows, dataset="nested", b=b, contraction=1), field)[1] for b in bs]
        axs.errorbar(bs, m, yerr=s, marker=mk, capsize=4, color=color, label=lbl)
    axs.set_xscale("log", base=2); axs.set_xticks(bs)
    axs.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    axs.set_xlabel("Block size b"); axs.set_ylabel("Bits per element")
    axs.set_title("Space"); axs.legend(); axs.grid(alpha=0.3)

    for field, color, lbl, mk in [("ceb_retrieve_us", "#1f77b4", "retrieve", "o"),
                                  ("ceb_rank_us", "#ff7f0e", "rank", "^")]:
        m = [agg(sel(rows, dataset="nested", b=b, contraction=1), field)[0] for b in bs]
        s = [agg(sel(rows, dataset="nested", b=b, contraction=1), field)[1] for b in bs]
        axt.errorbar(bs, m, yerr=s, marker=mk, capsize=4, color=color, label=lbl)
    axt.set_xscale("log", base=2); axt.set_xticks(bs)
    axt.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    axt.set_xlabel("Block size b"); axt.set_ylabel("Time per query (microseconds)")
    axt.set_title("Query time"); axt.legend(); axt.grid(alpha=0.3)

    fig.suptitle("Effect of block size b (nested data, with contraction)")
    save(fig, "block_size.png")


def write_summary(rows):
    fields = ["baseline_bpe", "ceb_bpe", "H_bpe", "L_bpe", "saving_pct",
              "baseline_total_bpe", "ceb_total_bpe", "total_saving_pct",
              "max_depth", "mean_depth", "ceb_build_ms", "base_build_ms",
              "ceb_retrieve_us", "ceb_rank_us", "base_retrieve_us", "base_rank_us"]
    groups = [("nested", B_FIXED, 1, "nested, contraction on"),
              ("nested", B_FIXED, 0, "nested, contraction off"),
              ("independent", B_FIXED, 1, "independent, contraction on")]
    with open(SUMMARY_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "metric", "mean", "std"])
        for d, b, con, name in groups:
            sub = sel(rows, dataset=d, b=b, contraction=con)
            for fld in fields:
                m, s = agg(sub, fld)
                w.writerow([name, fld, "%.4f" % m, "%.4f" % s])
    print("wrote", SUMMARY_PATH)


def main():
    rows = load_rows(CSV_PATH)
    plot_space_comparison(rows)
    plot_space_with_index(rows)
    plot_contraction_effect(rows)
    plot_query_time_comparison(rows)
    plot_build_time(rows)
    plot_block_size(rows)
    write_summary(rows)


if __name__ == "__main__":
    main()
