#!/usr/bin/env python3
"""Assemble claim-led, Nature-style vector figures from audited result tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd


BLUE = "#1F5A94"
BLUE_LIGHT = "#91B7D9"
TEAL = "#2A8C82"
RED = "#B85C5C"
GREY = "#767676"
LIGHT_GREY = "#D9D9D9"
BLACK = "#222222"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6.6,
        "ytick.labelsize": 6.6,
        "legend.fontsize": 6.6,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def panel_label(ax: plt.Axes, label: str, x: float = -0.14, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def deterministic_jitter(n: int, width: float) -> np.ndarray:
    """Return display-only offsets without sampling or changing observations."""
    if n <= 1:
        return np.zeros(n)
    return np.linspace(-width, width, n)


def save_bundle(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def clean_subset_label(value: str) -> str:
    replacements = {
        "Myeloid_c4_DCs_pDC_IRF7": "IRF7+ pDC",
        "T_cells_c6_IFIT1": "IFIT1+ T cells",
        "Myeloid_c9_Macrophage_2_CXCL10": "CXCL10+ macrophages",
        "Myeloid_c5_Macrophage_3_SIGLEC1": "SIGLEC1+ macrophages",
        "Endothelial CXCL12": "CXCL12+ endothelial",
        "Endothelial ACKR1": "ACKR1+ endothelial",
        "B cells Naive": "Naive B cells",
        "Myeloid_c2_LAM2_APOE": "APOE+ macrophages",
    }
    return replacements.get(value, value.replace("_", " "))


def figure3(project: Path, out_dir: Path) -> dict[str, object]:
    source = project / "results" / "single_cell"
    summary = pd.read_csv(source / "GSE176078_TNBC_IFI44L_celltype_major_summary.tsv", sep="\t")
    donor = pd.read_csv(source / "GSE176078_TNBC_IFI44L_by_donor_celltype_major.tsv", sep="\t")
    # Use the same donor-compartment eligibility rule as the paired tests.
    # Otherwise low-cell groups are displayed despite being excluded from inference.
    donor = donor.loc[donor["n_cells"] >= 20].copy()
    tests = pd.read_csv(source / "GSE176078_TNBC_author_major_paired_tests.tsv", sep="\t")
    subset = pd.read_csv(source / "GSE176078_TNBC_IFI44L_celltype_subset_summary.tsv", sep="\t")

    order = ["Myeloid", "Endothelial", "T-cells", "B-cells", "PVL", "CAFs", "Cancer Epithelial"]
    summary = summary.loc[summary["celltype_major"].isin(order)].set_index("celltype_major").loc[order].reset_index()

    fig = plt.figure(figsize=(7.2, 5.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.92, 1.38], height_ratios=[1.0, 1.05])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # a: expression/detection landscape
    y = np.arange(len(summary))[::-1]
    sizes = 25 + 230 * summary["median_detection_fraction"].to_numpy()
    colors = [RED if x == "Cancer Epithelial" else BLUE for x in summary["celltype_major"]]
    ax_a.scatter(summary["median_mean_logCP10k"], y, s=sizes, c=colors, alpha=0.9,
                 edgecolor="white", linewidth=0.6)
    ax_a.set_yticks(y, summary["celltype_major"])
    ax_a.set_xlabel("Median IFI44L expression\n(mean ln(CP10K + 1))")
    ax_a.set_title("Major cell compartments", loc="left", fontweight="bold")
    ax_a.set_xlim(left=0)
    ax_a.text(0.98, 0.04, "Dot size: median detection", transform=ax_a.transAxes,
              ha="right", va="bottom", fontsize=6.6, color=GREY)
    panel_label(ax_a, "a", x=-0.27)

    # b: donor-by-compartment heatmap
    heat_order = ["Myeloid", "Endothelial", "T-cells", "CAFs", "PVL", "Cancer Epithelial"]
    heat = donor.loc[donor["celltype_major"].isin(heat_order)].pivot(
        index="celltype_major", columns="donor_id", values="IFI44L_mean_logCP10k"
    ).reindex(heat_order)
    donors = sorted(heat.columns)
    heat = heat[donors]
    cmap = mpl.colormaps["Blues"].copy()
    cmap.set_bad("#F2F2F2")
    im = ax_b.imshow(heat.to_numpy(), aspect="auto", cmap=cmap, vmin=0,
                     vmax=np.nanquantile(heat.to_numpy(), 0.95))
    ax_b.set_yticks(np.arange(len(heat_order)), heat_order)
    ax_b.set_xticks(np.arange(len(donors)), [x.replace("CID", "") for x in donors],
                    rotation=45, rotation_mode="anchor", ha="right")
    ax_b.tick_params(length=0)
    ax_b.set_xlabel("Patient")
    ax_b.set_title("Patient-level expression", loc="left", fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.035, pad=0.02)
    cbar.set_label("Mean ln(CP10K + 1)")
    cbar.outline.set_linewidth(0.5)
    panel_label(ax_b, "b", x=-0.18)

    # c: paired donor differences against malignant epithelium
    classes = ["Myeloid", "Endothelial", "T-cells", "CAFs", "PVL"]
    comparison_q = []
    for i, cell_class in enumerate(classes):
        wide_unfiltered = donor.loc[donor["celltype_major"].isin([cell_class, "Cancer Epithelial"])].pivot(
            index="donor_id", columns="celltype_major", values="IFI44L_mean_logCP10k"
        )
        n_before = len(wide_unfiltered)
        wide = wide_unfiltered.dropna()
        n_after = len(wide)
        assert n_after <= n_before
        diff = wide[cell_class] - wide["Cancer Epithelial"]
        ax_c.scatter(i + deterministic_jitter(len(diff), 0.09), diff, color=BLUE, s=18,
                     alpha=0.78, edgecolor="white", linewidth=0.4)
        ax_c.plot([i - 0.18, i + 0.18], [np.median(diff)] * 2, color=BLACK, lw=1.4)
        row = tests.loc[(tests["cell_class"] == cell_class) & (tests["endpoint"] == "IFI44L_mean_logCP10k")]
        if len(row):
            comparison_q.append((i, float(row.iloc[0]["fdr_bh"])))
    if comparison_q and np.allclose([q for _, q in comparison_q], comparison_q[0][1]):
        ax_c.text(0.03, 0.98, f"All five comparisons: q = {comparison_q[0][1]:.3f}",
                  transform=ax_c.transAxes, ha="left", va="top", fontsize=6.6)
    else:
        for i, q in comparison_q:
            ax_c.text(i, 1.45, f"q =\n{q:.3f}", ha="center", va="bottom", fontsize=6.6)
    ax_c.axhline(0, color=GREY, lw=0.8, ls="--")
    ax_c.set_xticks(range(len(classes)), ["Myeloid", "Endo-\nthelial", "T cells", "CAFs", "PVL"])
    ax_c.set_xlim(-0.65, len(classes) - 0.35)
    ax_c.set_ylabel("Paired difference versus\ncancer epithelium")
    ax_c.set_ylim(-0.25, 1.65)
    ax_c.set_title("Within-patient comparison", loc="left", fontweight="bold")
    panel_label(ax_c, "c", x=-0.27)

    # d: fine subsets with highest median expression
    top = subset.loc[subset["n_donors"] >= 3].head(8).iloc[::-1].copy()
    labels = [clean_subset_label(x) for x in top["celltype_subset"]]
    yd = np.arange(len(top))
    ax_d.hlines(yd, 0, top["median_mean_logCP10k"], color=BLUE_LIGHT, lw=2)
    ax_d.scatter(top["median_mean_logCP10k"], yd, s=20 + 110 * top["median_detection_fraction"],
                 color=TEAL, edgecolor="white", linewidth=0.5)
    ax_d.set_yticks(yd, labels)
    ax_d.set_xlabel("Median IFI44L expression\n(mean ln(CP10K + 1))")
    ax_d.set_xlim(0, 1.82)
    ax_d.set_title("IFI44L-high cellular states", loc="left", fontweight="bold")
    for yi, (_, row) in enumerate(top.iterrows()):
        ax_d.text(1.74, yi, f"n={int(row['n_donors'])}", ha="right", va="center", fontsize=6.6, color=GREY)
    panel_label(ax_d, "d", x=-0.18)

    base = out_dir / "Figure3_single_cell_compartment"
    save_bundle(fig, base)
    return {"figure": "Figure 3", "base": str(base), "n_patients": 10, "n_cells": 42512,
            "minimum_cells_per_donor_compartment": 20,
            "paired_donor_counts": {cell_class: int(tests.loc[(tests["cell_class"] == cell_class) &
                (tests["endpoint"] == "IFI44L_mean_logCP10k"), "n_paired_donors"].iloc[0])
                for cell_class in classes},
            "inference_unit": "patient", "panels": ["a", "b", "c", "d"]}


def figure4(project: Path, out_dir: Path) -> dict[str, object]:
    source = project / "results" / "spatial" / "GSE210616"
    sections = pd.read_csv(source / "section_qc_and_IFI44L_spatial_autocorrelation.tsv", sep="\t")
    corr = pd.read_csv(source / "patient_IFI44L_signature_correlations.tsv", sep="\t")
    partial = pd.read_csv(source / "patient_partial_spearman.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 4.65), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.08, 1.42], height_ratios=[1, 1])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    # a: two sections per patient, showing reproducible spatial aggregation
    patient_order = sections.groupby("patient_id")["IFI44L_moran_I"].median().sort_values().index.tolist()
    for i, patient in enumerate(patient_order):
        values = sections.loc[sections["patient_id"] == patient].sort_values("section")
        ax_a.plot(values["IFI44L_moran_I"], [i] * len(values), color=LIGHT_GREY, lw=1, zorder=1)
        for _, row in values.iterrows():
            marker = "o" if int(row["section"]) == 1 else "s"
            color = BLUE if row["moran_fdr_bh"] < 0.05 and row["IFI44L_moran_I"] > 0 else GREY
            ax_a.scatter(row["IFI44L_moran_I"], i, marker=marker, s=25, color=color,
                         edgecolor="white", linewidth=0.5, zorder=2)
    ax_a.axvline(0, color=GREY, lw=0.8, ls="--")
    ax_a.set_yticks(range(len(patient_order)), patient_order)
    ax_a.set_xlabel("IFI44L Moran's I")
    ax_a.set_ylabel("Patient")
    ax_a.set_title("Spatial autocorrelation across 28 sections", loc="left", fontweight="bold")
    ax_a.scatter([], [], marker="o", color=BLUE, label="Section 1")
    ax_a.scatter([], [], marker="s", color=BLUE, label="Section 2")
    ax_a.legend(loc="lower right")
    panel_label(ax_a, "a", x=-0.24)

    # b: patient-level correlations for seven prespecified signatures
    order = ["IFN_response_excluding_IFI44L", "Myeloid", "T_cell", "Cytotoxicity", "Checkpoint", "Endothelial", "Epithelial"]
    labels = ["IFN response", "Myeloid", "T cell", "Cytotoxicity", "Checkpoint", "Endothelial", "Epithelial"]
    for i, sig in enumerate(order):
        values_unfiltered = corr.loc[corr["signature"] == sig, "spearman_rho"]
        n_before = len(values_unfiltered)
        values = values_unfiltered.dropna().to_numpy()
        n_after = len(values)
        assert n_after <= n_before
        color = RED if sig == "Epithelial" else (TEAL if sig == "IFN_response_excluding_IFI44L" else BLUE)
        ax_b.scatter(i + deterministic_jitter(len(values), 0.10), values, color=color, s=15,
                     alpha=0.72, edgecolor="white", linewidth=0.35)
        ax_b.plot([i - 0.20, i + 0.20], [np.median(values)] * 2, color=BLACK, lw=1.4)
    ax_b.axhline(0, color=GREY, lw=0.8, ls="--")
    ax_b.set_xticks(range(len(order)), labels, rotation=35, rotation_mode="anchor", ha="right")
    ax_b.set_ylabel("Within-tumour Spearman ρ")
    ax_b.set_title("Patient-level spatial co-localization", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.15)

    # c: partial correlations preserve IFN association but not a large checkpoint effect
    specs = [
        ("IFN_adjusted_for_lineage_modules", "IFN adjusted for lineage", TEAL),
        ("Checkpoint_adjusted_for_IFN_and_lineage_modules", "Checkpoint adjusted\nfor IFN and lineage", GREY),
    ]
    for i, (analysis, label, color) in enumerate(specs):
        values_unfiltered = partial.loc[partial["analysis"] == analysis, "partial_spearman_rho"]
        n_before = len(values_unfiltered)
        values = values_unfiltered.dropna().to_numpy()
        n_after = len(values)
        assert n_after <= n_before
        ax_c.scatter(i + deterministic_jitter(len(values), 0.07), values, color=color, s=20,
                     alpha=0.78, edgecolor="white", linewidth=0.4)
        ax_c.plot([i - 0.18, i + 0.18], [np.median(values)] * 2, color=BLACK, lw=1.5)
    ax_c.axhline(0, color=GREY, lw=0.8, ls="--")
    ax_c.set_xticks([0, 1], [x[1] for x in specs])
    ax_c.set_ylabel("Partial Spearman ρ")
    ax_c.set_title("Composition-adjusted associations", loc="left", fontweight="bold")
    panel_label(ax_c, "c", x=-0.15)

    base = out_dir / "Figure4_spatial_IFI44L_ecology_quantitative"
    save_bundle(fig, base)
    return {"figure": "Figure 4 quantitative draft", "base": str(base), "n_patients": 14,
            "n_sections": 28, "n_spots_after_qc": int(sections["n_spots_after_qc"].sum()),
            "inference_unit": "patient (two sections averaged)", "panels": ["a", "b", "c"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    set_style()
    out_dir = args.project / "figures" / "manuscript_nature_style"
    audits = [figure3(args.project, out_dir), figure4(args.project, out_dir)]
    (out_dir / "figure_build_audit.json").write_text(json.dumps(audits, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
