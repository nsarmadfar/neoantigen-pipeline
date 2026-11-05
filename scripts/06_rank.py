#!/usr/bin/env python
"""Filter and rank the candidates, and draw the summary figures.

A mutant peptide is only a plausible vaccine target if three things hold at
once: it is predicted to be presented by one of the patient's own HLA alleles,
the gene it comes from is actually transcribed, and the mutant version binds at
least as well as the normal one. Each of those cuts away a large part of the
list, so the filters are reported as a funnel rather than as a single number.

Usage:
  python scripts/06_rank.py <settings.yaml> <expression.tsv> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# Chart colours: slots 1-3 of a palette checked for colour-vision deficiency
# over all pairs, which is what a scatter plot needs.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK_SOFT = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK_SOFT,
    "xtick.color": INK_SOFT, "ytick.color": INK_SOFT,
    "axes.edgecolor": "#d9d8d2", "grid.color": "#ebeae4",
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "semibold",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})


def save(fig, out_dir, name):
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(Path(out_dir) / f"{name}.{suffix}", bbox_inches="tight")
    plt.close(fig)


def funnel_figure(steps, out_dir):
    """How many peptide/allele pairs survive each filter."""
    labels = [label for label, _ in steps][::-1]
    values = [value for _, value in steps][::-1]
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(labels) + 1.4))
    bars = ax.barh(labels, values, color=BLUE, height=0.58)
    ax.bar_label(bars, labels=[f"{v:,}" for v in values], padding=6,
                 color=INK, fontsize=9.5)
    ax.set_xlim(0, max(values) * 1.18 if max(values) else 1)
    ax.set_xlabel("mutant peptide / HLA allele pairs")
    ax.set_title("What each filter removes")
    ax.xaxis.grid(True, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    save(fig, out_dir, "filter_funnel")


def binding_figure(df, out_dir):
    """Mutant against normal binding affinity, for substitutions."""
    paired = df.dropna(subset=["wt_affinity", "affinity"])
    paired = paired[paired["wt_peptide"] != ""]
    if paired.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    limits = [
        min(paired["affinity"].min(), paired["wt_affinity"].min()) * 0.6,
        max(paired["affinity"].max(), paired["wt_affinity"].max()) * 1.6,
    ]
    ax.plot(limits, limits, color="#b9b8b1", linewidth=1, zorder=1)
    ax.axhline(500, color="#b9b8b1", linewidth=1, linestyle=(0, (4, 3)), zorder=1)

    groups = [("passes every filter", AQUA), ("presented, but filtered out", ORANGE),
              ("not presented", BLUE)]
    for label, colour in groups:
        subset = paired[paired["group"] == label]
        if subset.empty:
            continue
        ax.scatter(subset["wt_affinity"], subset["affinity"], s=26, color=colour,
                   alpha=0.75, linewidths=0.6, edgecolors=SURFACE, label=label,
                   zorder=3 if label.startswith("passes") else 2)

    best = paired[paired["group"] == groups[0][0]].nsmallest(5, "presentation_percentile")
    for _, row in best.iterrows():
        ax.annotate(f"{row['gene_name']} {row['peptide']}",
                    (row["wt_affinity"], row["affinity"]),
                    textcoords="offset points", xytext=(7, -3),
                    fontsize=8, color=INK)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(limits); ax.set_ylim(limits)
    ax.set_xlabel("normal peptide, predicted affinity (nM)")
    ax.set_ylabel("mutant peptide, predicted affinity (nM)")
    ax.set_title("Does the mutation improve binding?")
    ax.text(0.03, 0.97, "below the diagonal:\nmutant binds better",
            transform=ax.transAxes, va="top", fontsize=8.5, color=INK_SOFT)
    ax.grid(True, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    save(fig, out_dir, "mutant_vs_normal_affinity")


def allele_figure(candidates, alleles, out_dir):
    """How many candidates each HLA allele accounts for."""
    counts = candidates["allele"].value_counts().reindex(alleles).fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(6.4, 0.56 * len(counts) + 1.5))
    bars = ax.barh(counts.index[::-1], counts.values[::-1], color=BLUE, height=0.56)
    ax.bar_label(bars, labels=[f"{v:,}" for v in counts.values[::-1]],
                 padding=6, color=INK, fontsize=9.5)
    ax.set_xlim(0, max(counts.max() * 1.18, 1))
    ax.set_xlabel("candidate peptides")
    ax.set_title("Candidates by presenting allele")
    ax.xaxis.grid(True, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    save(fig, out_dir, "candidates_per_allele")


def main(settings_path, expression_path, results_dir):
    config = yaml.safe_load(Path(settings_path).read_text())
    cuts = config["filters"]
    results = Path(results_dir)

    peptides = pd.read_csv(results / "peptides.tsv", sep="\t").fillna({"wt_peptide": ""})
    mutant = pd.read_csv(results / "binding_mutant.tsv", sep="\t")
    normal = pd.read_csv(results / "binding_normal.tsv", sep="\t")
    expression = pd.read_csv(expression_path, sep="\t")

    if peptides.empty or mutant.empty:
        print("Nothing to rank.")
        return

    df = peptides.merge(mutant, on="peptide", how="inner")
    df = df.merge(
        normal[["peptide", "allele", "affinity", "presentation_percentile"]].rename(
            columns={"peptide": "wt_peptide",
                     "affinity": "wt_affinity",
                     "presentation_percentile": "wt_presentation_percentile"}),
        on=["wt_peptide", "allele"], how="left",
    )

    tpm = dict(zip(expression["gene_id"], expression["tpm"]))
    df["tpm"] = df["gene_id"].map(tpm)

    # How much better the mutant binds than the normal peptide. Above 1 means
    # better. Frameshifts and indels have no normal counterpart, so this is
    # left empty for them and they are not filtered on it.
    df["agretopicity"] = df["wt_affinity"] / df["affinity"]

    presented = df["presentation_percentile"] <= cuts["max_presentation_percentile"]
    expressed = df["tpm"].fillna(0) >= cuts["min_tpm"]
    improved = df["agretopicity"].isna() | (df["agretopicity"] >= cuts["min_agretopicity"])

    df["group"] = np.where(
        presented & expressed & improved, "passes every filter",
        np.where(presented, "presented, but filtered out", "not presented"))

    candidates = df[presented & expressed & improved].sort_values(
        ["presentation_percentile", "affinity"]).reset_index(drop=True)
    candidates.insert(0, "rank", candidates.index + 1)

    steps = [
        ("all mutant peptides scored", len(df)),
        (f"presented (percentile <= {cuts['max_presentation_percentile']:g})",
         int(presented.sum())),
        (f"and gene expressed (>= {cuts['min_tpm']:g} TPM)",
         int((presented & expressed).sum())),
        (f"and binds better than normal (>= {cuts['min_agretopicity']:g}x)",
         len(candidates)),
    ]

    columns = [
        "rank", "gene_name", "peptide", "wt_peptide", "allele",
        "presentation_percentile", "presentation_score", "affinity",
        "wt_affinity", "agretopicity", "processing_score", "tpm", "vaf",
        "consequence", "protein_change", "variant_id", "transcript",
        "position_in_protein", "length",
    ]
    candidates[columns].to_csv(results / "candidates.tsv", sep="\t", index=False,
                               float_format="%.4g")
    df.drop(columns=["wt_protein", "mt_protein"], errors="ignore").to_csv(
        results / "all_scored_peptides.tsv.gz", sep="\t", index=False,
        float_format="%.4g", compression="gzip")

    alleles = sorted(df["allele"].unique())
    funnel_figure(steps, results)
    binding_figure(df, results)
    allele_figure(candidates, alleles, results)

    with open(results / "summary.tsv", "w") as handle:
        handle.write("item\tvalue\n")
        for label, value in steps:
            handle.write(f"{label}\t{value}\n")
        handle.write(f"distinct candidate peptides\t{candidates['peptide'].nunique()}\n")
        handle.write(f"genes with a candidate\t{candidates['gene_name'].nunique()}\n")
        handle.write(f"variants with a candidate\t{candidates['variant_id'].nunique()}\n")

    for label, value in steps:
        print(f"  {label}: {value:,}")
    print(f"{candidates['peptide'].nunique()} distinct candidate peptides "
          f"in {candidates['gene_name'].nunique()} genes")
    if not candidates.empty:
        print("\nTop 10:")
        preview = candidates.head(10)[
            ["rank", "gene_name", "peptide", "allele",
             "presentation_percentile", "affinity", "tpm"]]
        print(preview.to_string(index=False))
    print(f"\nDone. Results are in {results}/")


if __name__ == "__main__":
    main(*sys.argv[1:4])
