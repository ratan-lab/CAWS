#!/usr/bin/env python3
"""Mean signal profile around region reference points, one panel per condition.

Reads a deepTools computeMatrix (reference-point) file whose sample labels are the
sampleIDs, and stats/samples.tsv for each sample's display label and condition.
Panels share the y axis so conditions can be compared; IgG controls are drawn last.

usage: plot_region_profile.py <matrix.gz> <samples.tsv> <region> <ylabel> <out.pdf>
"""

import gzip
import itertools
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]
LINESTYLES = ["-", "--", ":", "-."]


def mean_profiles(matrix_file, chunk=20000):
    """Per-sample mean over regions, streamed so large region sets fit in memory."""
    with gzip.open(matrix_file, "rt") as f:
        meta = json.loads(f.readline().strip()[1:])  # first line: '@' + JSON
        total, count = None, None
        while True:
            lines = list(itertools.islice(f, chunk))
            if not lines:
                break
            values = np.array([line.rstrip("\n").split("\t")[6:] for line in lines], dtype=float)
            s, n = np.nansum(values, axis=0), np.sum(~np.isnan(values), axis=0)
            total, count = (s, n) if total is None else (total + s, count + n)
    if total is None:
        return meta, None
    with np.errstate(invalid="ignore", divide="ignore"):
        return meta, total / count


def main(matrix_file, samples_file, region, ylabel, out_file):
    meta, profile = mean_profiles(matrix_file)
    sample_ids = meta["sample_labels"]
    bounds = meta["sample_boundaries"]

    with open(samples_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    info = {r["sampleID"]: r for r in rows}

    # experimental conditions in samplesheet order, IgG controls last
    conditions = list(dict.fromkeys(info[s]["condition"] for s in sample_ids))
    conditions = [c for c in conditions if c != "Control"] + (["Control"] if "Control" in conditions else [])

    ncols = min(4, len(conditions))
    nrows = -(-len(conditions) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.6 * nrows), sharey=True,
                             squeeze=False, layout="constrained")
    upstream, downstream = meta["upstream"][0], meta["downstream"][0]

    for ax, condition in zip(axes.flat, conditions):
        members = [i for i, s in enumerate(sample_ids) if info[s]["condition"] == condition]
        for k, i in enumerate(members):
            if profile is None:
                continue
            y = profile[bounds[i]:bounds[i + 1]]
            x = np.linspace(-upstream, downstream, len(y)) / 1000
            ax.plot(x, y, color=OKABE_ITO[k % len(OKABE_ITO)],
                    linestyle=LINESTYLES[(k // len(OKABE_ITO)) % len(LINESTYLES)],
                    linewidth=1.5, label=info[sample_ids[i]]["label"])
        ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
        ax.set_title("IgG controls" if condition == "Control" else condition, fontsize=10)
        ax.set_xlabel(f"Distance from {region} (kb)", fontsize=8)
        ax.tick_params(labelsize=8)
        # legend under the panel, so long labels never cover the curves
        if profile is not None:
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=7, frameon=False,
                      ncol=1 if len(members) <= 4 else 2)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=8)
    for ax in list(axes.flat)[len(conditions):]:
        ax.set_visible(False)
    fig.suptitle(region, fontweight="bold")
    if profile is None:
        fig.text(0.5, 0.5, "no regions with signal", ha="center")
    fig.savefig(out_file)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
