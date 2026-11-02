#!/usr/bin/env python3
"""Build claim-led manuscript Figures 3, 4 and 7 from audited result tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


BLUE = "#2B6CA3"
BLUE_LIGHT = "#A9C8E2"
TEAL = "#23877D"
RED = "#B64E55"
ORANGE = "#D28A35"
GREY = "#747474"
LIGHT_GREY = "#D8D8D8"
PALE_GREY = "#F3F3F3"
BLACK = "#222222"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2,
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


def save_bundle(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure1(project: Path, out_dir: Path) -> dict[str, object]:
    source = project / "results" / "bulk"
    assoc = pd.read_csv(source / "continuous_IFI44L_gene_associations.tsv", sep="\t")
    gsea = pd.read_csv(source / "GSEA_GO_BP_all.tsv", sep="\t")
    immune = pd.read_csv(source / "IFI44L_immune_signature_correlations.tsv", sep="\t")
    checkpoint = pd.read_csv(source / "IFI44L_checkpoint_partial_correlations.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.03, 0.97])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # a: continuous IFI44L association landscape
    dat = assoc.loc[assoc["gene"] != "IFI44L"].copy()
    dat["mlog10fdr"] = -np.log10(dat["fdr_bh"].clip(lower=np.finfo(float).tiny))
    sig = dat["fdr_bh"] < 0.05
    highlight = sig & (dat["beta_per_IFI44L_SD"].abs() >= 0.585)
    ax_a.scatter(dat.loc[~sig, "beta_per_IFI44L_SD"], dat.loc[~sig, "mlog10fdr"],
                 s=4, color=LIGHT_GREY, alpha=0.45, linewidth=0)
    ax_a.scatter(dat.loc[sig & ~highlight, "beta_per_IFI44L_SD"], dat.loc[sig & ~highlight, "mlog10fdr"],
                 s=5, color=BLUE_LIGHT, alpha=0.55, linewidth=0)
    ax_a.scatter(dat.loc[highlight, "beta_per_IFI44L_SD"], dat.loc[highlight, "mlog10fdr"],
                 s=8, color=BLUE, alpha=0.8, linewidth=0)
    label_positions = {
        "CMPK2": (0.80, 73.0),
        "OAS2": (1.34, 69.0),
        "IFIT1": (1.39, 63.0),
    }
    for gene, position in label_positions.items():
        row = dat.loc[dat["gene"] == gene]
        if len(row):
            x, y = row.iloc[0][["beta_per_IFI44L_SD", "mlog10fdr"]]
            ax_a.annotate(gene, (x, y), xytext=position, textcoords="data", fontsize=5.2,
                          color=BLACK, ha="center", va="center",
                          arrowprops={"arrowstyle": "-", "lw": 0.35, "color": GREY,
                                      "shrinkA": 4, "shrinkB": 3})
    ax_a.axhline(-np.log10(0.05), color=GREY, lw=0.7, ls="--")
    ax_a.axvline(0, color=GREY, lw=0.6)
    ax_a.set_xlabel("Expression change per IFI44L s.d.")
    ax_a.set_ylabel("−log10 FDR")
    ax_a.set_title("Continuous transcriptome associations", loc="left", fontweight="bold")
    panel_label(ax_a, "a", x=-0.19)

    # b: non-redundant, prespecified biological-process view
    pathway_map = {
        "defense response to virus": "Defence response to virus",
        "regulation of innate immune response": "Innate immune response",
        "regulation of T cell activation": "T-cell activation",
        "positive regulation of cytokine production": "Cytokine production",
        "cytoplasmic translation": "Cytoplasmic translation",
        "oxidative phosphorylation": "Oxidative phosphorylation",
        "mitochondrial ATP synthesis coupled electron transport": "Mitochondrial electron transport",
    }
    rows = []
    for description, label in pathway_map.items():
        hit = gsea.loc[gsea["Description"] == description]
        if len(hit):
            row = hit.iloc[0].copy()
            row["display"] = label
            rows.append(row)
    pathways = pd.DataFrame(rows).sort_values("NES")
    y = np.arange(len(pathways))
    colours = [BLUE if value > 0 else ORANGE for value in pathways["NES"]]
    ax_b.hlines(y, 0, pathways["NES"], color=colours, lw=2)
    ax_b.scatter(pathways["NES"], y, s=22, color=colours, edgecolor="white", linewidth=0.5)
    ax_b.axvline(0, color=GREY, lw=0.7)
    ax_b.set_yticks(y, pathways["display"])
    ax_b.set_xlabel("Normalized enrichment score")
    ax_b.set_title("GO biological-process enrichment", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.25)

    # c: immune-state correlations
    order = ["IFN_ALPHA_RESPONSE", "ANTIGEN_PRESENTATION", "IFN_GAMMA_RESPONSE", "CHECKPOINT_AXIS",
             "T_CELL_INFLAMED", "CYTOLYTIC_ACTIVITY", "CD8_T_CELL", "EPITHELIAL_CONTENT"]
    labels = ["IFN-α response", "Antigen presentation", "IFN-γ response", "Checkpoint axis",
              "T-cell inflamed", "Cytolytic activity", "CD8 T cell", "Epithelial content"]
    plot = immune.set_index("signature").loc[order].iloc[::-1]
    y = np.arange(len(plot))
    colours = [RED if name == "EPITHELIAL_CONTENT" else TEAL for name in plot.index]
    ax_c.hlines(y, 0, plot["spearman_rho"], color=colours, lw=2)
    ax_c.scatter(plot["spearman_rho"], y, s=24, color=colours, edgecolor="white", linewidth=0.5)
    ax_c.axvline(0, color=GREY, lw=0.7)
    ax_c.set_yticks(y, labels[::-1])
    ax_c.set_xlim(-0.15, 1.0)
    ax_c.set_xlabel("Spearman ρ with IFI44L")
    ax_c.set_title("Prespecified immune signatures", loc="left", fontweight="bold")
    panel_label(ax_c, "c", x=-0.29)

    # d: raw and inflammation-adjusted checkpoint associations
    checkpoint = checkpoint.sort_values("raw_spearman_rho").reset_index(drop=True)
    y = np.arange(len(checkpoint))
    for i, row in checkpoint.iterrows():
        ax_d.plot([row["adjusted_spearman_rho"], row["raw_spearman_rho"]], [i, i], color=LIGHT_GREY, lw=1.2)
    ax_d.scatter(checkpoint["raw_spearman_rho"], y, s=22, color=BLUE, label="Raw")
    adj_colours = [TEAL if q < 0.05 else GREY for q in checkpoint["adjusted_fdr_bh"]]
    ax_d.scatter(checkpoint["adjusted_spearman_rho"], y, s=22, color=adj_colours, marker="s", label="Adjusted")
    ax_d.axvline(0, color=GREY, lw=0.7)
    ax_d.set_yticks(y, checkpoint["gene"])
    ax_d.set_xlim(-0.30, 0.78)
    ax_d.set_xlabel("Spearman ρ")
    ax_d.set_title("Checkpoint associations after state adjustment", loc="left", fontweight="bold")
    ax_d.legend(loc="lower right")
    panel_label(ax_d, "d", x=-0.19)

    base = out_dir / "Figure1_TCGA_IFI44L_inflammatory_state"
    save_bundle(fig, base)
    return {"figure": "Figure 1", "base": str(base), "cohort": "TCGA TNBC", "n": 187,
            "inference": "sample-level continuous analysis", "panels": ["a", "b", "c", "d"]}


def km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(time)
    time, event = time[order], event[order]
    event_times = np.unique(time[event == 1])
    survival = 1.0
    xs, ys = [0.0], [1.0]
    for t in event_times:
        at_risk = np.sum(time >= t)
        events = np.sum((time == t) & (event == 1))
        xs.extend([float(t), float(t)])
        ys.extend([survival, survival * (1 - events / at_risk)])
        survival = ys[-1]
    xs.append(float(np.max(time)))
    ys.append(survival)
    return np.asarray(xs), np.asarray(ys)


def figure2(project: Path, out_dir: Path) -> dict[str, object]:
    bulk = project / "results" / "bulk"
    source = project / "results" / "validation_GSE58812"
    tcga = pd.read_csv(bulk / "IFI44L_immune_signature_correlations.tsv", sep="\t")
    val = pd.read_csv(source / "GSE58812_IFI44L_immune_signature_correlations.tsv", sep="\t")
    checkpoints = pd.read_csv(source / "GSE58812_IFI44L_checkpoint_correlations.tsv", sep="\t")
    cox = pd.read_csv(source / "GSE58812_continuous_IFI44L_Cox.tsv", sep="\t")
    clinical = pd.read_csv(source / "GSE58812_clinical_IFI44L.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.08, 0.9, 1.12])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[0, 2])
    ax_e = fig.add_subplot(gs[1, 2])

    order = ["IFN_ALPHA_RESPONSE", "IFN_GAMMA_RESPONSE", "ANTIGEN_PRESENTATION", "CHECKPOINT_AXIS",
             "T_CELL_INFLAMED", "CYTOLYTIC_ACTIVITY", "CD8_T_CELL", "EPITHELIAL_CONTENT"]
    labels = ["IFN-α response", "IFN-γ response", "Antigen presentation", "Checkpoint axis",
              "T-cell inflamed", "Cytolytic activity", "CD8 T cell", "Epithelial content"]
    t = tcga.set_index("signature").loc[order, "spearman_rho"]
    v = val.set_index("signature").loc[order, "spearman_rho"]
    y = np.arange(len(order))[::-1]
    for yi, a, b in zip(y, t, v):
        ax_a.plot([a, b], [yi, yi], color=LIGHT_GREY, lw=1.2)
    ax_a.scatter(t, y, s=24, color=BLUE, label="TCGA (n=187)")
    ax_a.scatter(v, y, s=24, color=TEAL, marker="s", label="GSE58812 (n=107)")
    ax_a.axvline(0, color=GREY, lw=0.7)
    ax_a.set_yticks(y, labels)
    ax_a.set_xlim(-0.4, 1.0)
    ax_a.set_xlabel("Spearman ρ with IFI44L")
    ax_a.set_title("Immune-state replication", loc="left", fontweight="bold")
    ax_a.legend(loc="lower right")
    panel_label(ax_a, "a", x=-0.28)

    checkpoints = checkpoints.sort_values("spearman_rho").reset_index(drop=True)
    yb = np.arange(len(checkpoints))
    ax_b.hlines(yb, 0, checkpoints["spearman_rho"], color=BLUE_LIGHT, lw=2)
    ax_b.scatter(checkpoints["spearman_rho"], yb, s=22,
                 color=[BLUE if q < 0.05 else GREY for q in checkpoints["fdr_bh"]])
    ax_b.axvline(0, color=GREY, lw=0.7)
    ax_b.set_yticks(yb, checkpoints["gene"])
    ax_b.set_xlim(0, 0.7)
    ax_b.set_xlabel("Spearman ρ")
    ax_b.set_title("Checkpoint genes", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.30)

    cox_order = [
        ("Metastasis-free survival", "Unadjusted", "MFS, unadjusted"),
        ("Metastasis-free survival", "Age-adjusted", "MFS, age-adjusted"),
        ("Overall survival", "Unadjusted", "OS, unadjusted"),
        ("Overall survival", "Age-adjusted", "OS, age-adjusted"),
    ]
    rows = []
    for endpoint, model, label in cox_order:
        row = cox.loc[(cox["endpoint"] == endpoint) & (cox["model"] == model)].iloc[0].copy()
        row["display"] = label
        rows.append(row)
    forest = pd.DataFrame(rows).iloc[::-1].reset_index(drop=True)
    yc = np.arange(len(forest))
    ax_c.hlines(yc, forest["ci95_low"], forest["ci95_high"], color=BLACK, lw=1.1)
    ax_c.scatter(forest["hazard_ratio_per_IFI44L_SD"], yc, s=25, color=TEAL, zorder=2)
    ax_c.axvline(1, color=GREY, lw=0.7, ls="--")
    ax_c.set_yticks(yc, forest["display"])
    ax_c.set_xlim(0.35, 1.10)
    ax_c.set_xlabel("Hazard ratio per IFI44L s.d. (95% CI)")
    ax_c.set_title("Continuous survival models", loc="left", fontweight="bold")
    panel_label(ax_c, "c", x=-0.30)

    for ax, endpoint, time_col, event_col, label in [
        (ax_d, "MFS", "mfs_days", "metastasis_event", "Metastasis-free survival"),
        (ax_e, "OS", "os_days", "death_event", "Overall survival"),
    ]:
        for group, colour in [("High", TEAL), ("Low", GREY)]:
            sub = clinical.loc[clinical["IFI44L_group"] == group]
            x, surv = km_curve(sub[time_col].to_numpy(float) / 365.25, sub[event_col].to_numpy(int))
            ax.step(x, surv, where="post", color=colour, lw=1.4, label=f"{group} (n={len(sub)})")
        ax.set_xlim(left=0)
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("Survival probability")
        ax.set_title(label, loc="left", fontweight="bold")
        ax.legend(loc="lower left")
    panel_label(ax_d, "d", x=-0.20)
    panel_label(ax_e, "e", x=-0.20)

    base = out_dir / "Figure2_GSE58812_replication_survival"
    save_bundle(fig, base)
    return {"figure": "Figure 2", "base": str(base), "cohort": "GSE58812", "n": 107,
            "events": {"MFS": 31, "OS": 29}, "primary_survival_inference": "continuous Cox models",
            "median_split": "visual sensitivity analysis", "panels": ["a", "b", "c", "d", "e"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    set_style()
    out_dir = args.project / "figures" / "manuscript_nature_style"
    audits = [figure1(args.project, out_dir), figure2(args.project, out_dir)]
    audit_path = out_dir / "figure_build_audit_figure1_2_v7.json"
    audit_path.write_text(json.dumps(audits, indent=2), encoding="utf-8")
    print(audit_path)


if __name__ == "__main__":
    main()
