#!/usr/bin/env python3
"""Build v11 public-multiomics Figures 2 and 5.

The v11 figures deliberately exclude local patient/IHC material, cell-culture
perturbation assays and the low-powered treatment-response boundary analysis. TCGA OS/PFI results are displayed whether
statistically significant or not.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm


BLUE = "#2B6CA3"
BLUE_LIGHT = "#A9C8E2"
TEAL = "#23877D"
RED = "#B64E55"
ORANGE = "#D28A35"
GREY = "#747474"
LIGHT_GREY = "#D8D8D8"
BLACK = "#222222"
PALE_BLUE = "#F1F6FA"
PALE_TEAL = "#EFF8F6"
PALE_GREY = "#F5F5F5"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7.2,
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


def panel_label(ax: plt.Axes, label: str, x: float = -0.15, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom")


def save_bundle(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Kaplan-Meier with Greenwood log-log 95% confidence limits.

    Ties include all observations at a time in its risk set before event/censor
    removal. This is cross-checked against survival::survfit(conf.type='log-log').
    """
    order = np.argsort(time)
    time, event = time[order], event[order]
    event_times = np.unique(time[event == 1])
    survival = 1.0
    xs, ys, lows, highs = [0.0], [1.0], [1.0], [1.0]
    greenwood = 0.0
    z = float(norm.ppf(0.975))
    for value in event_times:
        at_risk = np.sum(time >= value)
        events = np.sum((time == value) & (event == 1))
        survival *= 1 - events / at_risk
        if events < at_risk:
            greenwood += events / (at_risk * (at_risk - events))
        if survival <= 0:
            low, high = 0.0, 0.0
        elif survival >= 1:
            low, high = 1.0, 1.0
        else:
            loglog = np.log(-np.log(survival))
            se_loglog = np.sqrt(greenwood) / abs(np.log(survival))
            low = float(np.exp(-np.exp(loglog + z * se_loglog)))
            high = float(np.exp(-np.exp(loglog - z * se_loglog)))
        xs.append(float(value)); ys.append(survival)
        lows.append(low); highs.append(high)
    xs.append(float(np.max(time)))
    ys.append(survival)
    lows.append(lows[-1]); highs.append(highs[-1])
    return np.asarray(xs), np.asarray(ys), np.asarray(lows), np.asarray(highs)


def logrank_p(frame: pd.DataFrame, time_col: str, event_col: str, group_col: str) -> float:
    groups = list(frame[group_col].dropna().unique())
    if len(groups) != 2:
        raise ValueError("Log-rank helper requires exactly two groups")
    times = np.sort(frame.loc[frame[event_col] == 1, time_col].unique())
    observed = 0.0
    expected = 0.0
    variance = 0.0
    group = groups[0]
    for value in times:
        at_risk_all = np.sum(frame[time_col] >= value)
        at_risk_group = np.sum((frame[time_col] >= value) & (frame[group_col] == group))
        events_all = np.sum((frame[time_col] == value) & (frame[event_col] == 1))
        events_group = np.sum(
            (frame[time_col] == value) & (frame[event_col] == 1) & (frame[group_col] == group)
        )
        if at_risk_all <= 1:
            continue
        observed += events_group
        expected += events_all * at_risk_group / at_risk_all
        variance += (
            at_risk_group
            * (at_risk_all - at_risk_group)
            * events_all
            * (at_risk_all - events_all)
            / (at_risk_all**2 * (at_risk_all - 1))
        )
    if variance <= 0:
        return float("nan")
    statistic = (observed - expected) ** 2 / variance
    return float(chi2.sf(statistic, 1))


def draw_km(
    ax: plt.Axes,
    risk_ax: plt.Axes,
    frame: pd.DataFrame,
    time_col: str,
    event_col: str,
    group_col: str,
    title: str,
) -> tuple[float, list[dict], list[dict], list[dict]]:
    max_years = float(frame[time_col].max()) / 365.25
    step_years = 2 if max_years <= 10 else 5
    ticks = np.arange(0, max_years + 1e-8, step_years)
    km_audit_rows, censor_audit_rows, risk_audit_rows = [], [], []
    for group_index, (group, colour) in enumerate([("High", TEAL), ("Low", GREY)]):
        subset = frame.loc[frame[group_col] == group]
        time = subset[time_col].to_numpy(float) / 365.25
        event = subset[event_col].to_numpy(int)
        x, survival, low, high = km_curve(
            time, event,
        )
        ax.fill_between(x, low, high, step="post", color=colour, alpha=0.12, linewidth=0)
        ax.step(x, survival, where="post", color=colour, lw=1.4, label=f"{group} (n={len(subset)})")
        censor_times = time[event == 0]
        censor_survival = survival[np.maximum(np.searchsorted(x, censor_times, side="right") - 1, 0)]
        ax.scatter(censor_times, censor_survival, marker="+", s=10, linewidths=0.55, color=colour)
        for value, probability, lower, upper in zip(x, survival, low, high):
            km_audit_rows.append({"endpoint": title, "group": group, "time_years": value,
                "survival": probability, "ci95_low": lower, "ci95_high": upper})
        for value, probability in zip(censor_times, censor_survival):
            censor_audit_rows.append({"endpoint": title, "group": group,
                "time_years": value, "survival_at_censor": probability})
        for value in ticks:
            n_risk = int(np.sum(time >= value))
            risk_ax.text(value, 1 - group_index, str(n_risk), color=colour, ha="center", va="center", fontsize=6.6)
            risk_audit_rows.append({"endpoint": title, "group": group, "time_years": value, "n_at_risk": n_risk})
        risk_ax.text(-0.055, 1 - group_index, group, transform=risk_ax.get_yaxis_transform(),
            color=colour, ha="right", va="center", fontsize=6.6)
    p_value = logrank_p(frame, time_col, event_col, group_col)
    ax.set_xlim(0, max_years * 1.02)
    ax.set_xticks(ticks)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("")
    ax.set_ylabel("Survival probability")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.text(0.98, 0.97, f"Log-rank P = {p_value:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.6)
    ax.legend(loc="lower left")
    risk_ax.set_xlim(ax.get_xlim())
    risk_ax.set_xticks(ticks)
    risk_ax.set_ylim(-0.6, 2.0)
    risk_ax.set_yticks([])
    risk_ax.tick_params(axis="x", length=0, labelbottom=False)
    for spine in risk_ax.spines.values(): spine.set_visible(False)
    risk_ax.text(0, 1.70, "Number at risk", fontsize=6.6, va="center")
    risk_ax.set_xlabel("Time (years)", labelpad=0.5)
    return p_value, km_audit_rows, censor_audit_rows, risk_audit_rows


def build_figure1(project: Path, output_dir: Path) -> dict[str, object]:
    source = project / "results" / "bulk"
    associations = pd.read_csv(source / "continuous_IFI44L_gene_associations.tsv", sep="\t")
    gsea = pd.read_csv(source / "GSEA_GO_BP_all.tsv", sep="\t")
    immune = pd.read_csv(source / "IFI44L_immune_signature_correlations.tsv", sep="\t")
    checkpoints = pd.read_csv(source / "IFI44L_checkpoint_partial_correlations.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.03, 0.97])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    data = associations.loc[associations["gene"] != "IFI44L"].copy()
    data["mlog10fdr"] = -np.log10(data["fdr_bh"].clip(lower=np.finfo(float).tiny))
    significant = data["fdr_bh"] < 0.05
    ax_a.scatter(data.loc[~significant, "beta_per_IFI44L_SD"], data.loc[~significant, "mlog10fdr"], s=4, color=LIGHT_GREY, alpha=0.45, linewidth=0)
    ax_a.scatter(data.loc[significant, "beta_per_IFI44L_SD"], data.loc[significant, "mlog10fdr"], s=5, color=BLUE, alpha=0.7, linewidth=0)
    for gene, offset in {"CMPK2": (-23, 18), "OAS2": (10, 10), "IFIT1": (19, -8)}.items():
        row = data.loc[data["gene"] == gene]
        if len(row):
            x, y = row.iloc[0][["beta_per_IFI44L_SD", "mlog10fdr"]]
            ax_a.annotate(gene, (x, y), xytext=offset, textcoords="offset points", fontsize=6.6, color=BLACK, ha="center", va="center", arrowprops={"arrowstyle": "-", "lw": 0.35, "color": GREY, "shrinkA": 4, "shrinkB": 3})
    ax_a.axhline(-np.log10(0.05), color=GREY, lw=0.7, ls="--")
    ax_a.axvline(0, color=GREY, lw=0.6)
    ax_a.set_xlabel("Expression coefficient per IFI44L s.d.\n(archived units)")
    ax_a.set_ylabel("-log10 FDR")
    ax_a.set_title("Continuous transcriptome associations", loc="left", fontweight="bold")
    panel_label(ax_a, "a", x=-0.19)

    pathway_map = {
        "defense response to virus": "Defence response to virus",
        "regulation of innate immune response": "Innate immune response",
        "regulation of T cell activation": "T-cell activation",
        "positive regulation of cytokine production": "Cytokine production",
        "cytoplasmic translation": "Cytoplasmic translation",
        "oxidative phosphorylation": "Oxidative phosphorylation",
        "mitochondrial ATP synthesis coupled electron transport": "Mitochondrial electron transport",
    }
    pathway_rows = []
    for description, display in pathway_map.items():
        hit = gsea.loc[gsea["Description"] == description]
        if len(hit):
            row = hit.iloc[0].copy()
            row["display"] = display
            pathway_rows.append(row)
    pathways = pd.DataFrame(pathway_rows).sort_values("NES")
    y = np.arange(len(pathways))
    colours = [BLUE if value > 0 else ORANGE for value in pathways["NES"]]
    ax_b.hlines(y, 0, pathways["NES"], color=colours, lw=2)
    ax_b.scatter(pathways["NES"], y, s=22, color=colours, edgecolor="white", linewidth=0.5)
    ax_b.axvline(0, color=GREY, lw=0.7)
    ax_b.set_yticks(y, pathways["display"])
    ax_b.set_xlabel("Normalized enrichment score")
    ax_b.set_title("GO biological-process enrichment", loc="left", fontweight="bold")
    panel_label(ax_b, "b", x=-0.25)

    order = ["IFN_ALPHA_RESPONSE", "ANTIGEN_PRESENTATION", "IFN_GAMMA_RESPONSE", "CHECKPOINT_AXIS", "T_CELL_INFLAMED", "CYTOLYTIC_ACTIVITY", "CD8_T_CELL", "EPITHELIAL_CONTENT"]
    labels = ["IFN-alpha response", "Antigen presentation", "IFN-gamma response", "Checkpoint axis", "T-cell inflamed", "Cytolytic activity", "CD8 T cell", "Epithelial content"]
    plot = immune.set_index("signature").loc[order].iloc[::-1]
    y = np.arange(len(plot))
    colours = [RED if name == "EPITHELIAL_CONTENT" else TEAL for name in plot.index]
    ax_c.hlines(y, 0, plot["spearman_rho"], color=colours, lw=2)
    ax_c.scatter(plot["spearman_rho"], y, s=24, color=colours, edgecolor="white", linewidth=0.5)
    ax_c.axvline(0, color=GREY, lw=0.7)
    ax_c.set_yticks(y, labels[::-1])
    ax_c.set_xlim(-0.15, 1.0)
    ax_c.set_xlabel("Spearman rho with IFI44L")
    ax_c.set_title("Prespecified immune signatures", loc="left", fontweight="bold")
    panel_label(ax_c, "c", x=-0.29)

    checkpoints = checkpoints.sort_values("raw_spearman_rho").reset_index(drop=True)
    y = np.arange(len(checkpoints))
    for index, row in checkpoints.iterrows():
        ax_d.plot([row["adjusted_spearman_rho"], row["raw_spearman_rho"]], [index, index], color=LIGHT_GREY, lw=1.2)
    ax_d.scatter(checkpoints["raw_spearman_rho"], y, s=22, color=BLUE, label="Raw")
    adjusted_colours = [TEAL if value < 0.05 else GREY for value in checkpoints["adjusted_fdr_bh"]]
    ax_d.scatter(checkpoints["adjusted_spearman_rho"], y, s=22, color=adjusted_colours, marker="s", label="Adjusted")
    ax_d.axvline(0, color=GREY, lw=0.7)
    ax_d.set_yticks(y, checkpoints["gene"])
    ax_d.set_xlim(-0.30, 0.78)
    ax_d.set_xlabel("Spearman rho")
    ax_d.set_title("Checkpoint associations after state adjustment", loc="left", fontweight="bold")
    ax_d.legend(loc="lower right")
    panel_label(ax_d, "d", x=-0.19)

    base = output_dir / "Figure1"
    save_bundle(fig, base)
    return {"figure": "Figure 1", "base": str(base), "cohort": "TCGA receptor-defined TNBC",
            "n": int(immune["n"].iloc[0]), "panels": ["a", "b", "c", "d"]}


def build_figure2(project: Path, output_dir: Path) -> dict[str, object]:
    bulk = project / "results" / "bulk"
    validation = project / "results" / "validation_GSE58812"
    tcga_cor = pd.read_csv(bulk / "IFI44L_immune_signature_correlations.tsv", sep="\t")
    gse_cor = pd.read_csv(validation / "GSE58812_IFI44L_immune_signature_correlations.tsv", sep="\t")
    gse_cox = pd.read_csv(validation / "GSE58812_continuous_IFI44L_Cox.tsv", sep="\t")
    gse_clinical = pd.read_csv(validation / "GSE58812_clinical_IFI44L.tsv", sep="\t")
    tcga_cox = pd.read_csv(bulk / "TCGA_TNBC_continuous_IFI44L_Cox.tsv", sep="\t")
    tcga_clinical = pd.read_csv(bulk / "TCGA_TNBC_survival_merged.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 8.0), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.22, 1.22])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    curve_axes, risk_axes = [], []
    for row, col in [(1, 0), (1, 1), (2, 0), (2, 1)]:
        nested = grid[row, col].subgridspec(2, 1, height_ratios=[4.0, 1.0], hspace=0.02)
        curve_axes.append(fig.add_subplot(nested[0, 0]))
        risk_axes.append(fig.add_subplot(nested[1, 0]))
    ax_c, ax_d, ax_e, ax_f = curve_axes

    order = [
        "IFN_ALPHA_RESPONSE", "IFN_GAMMA_RESPONSE", "ANTIGEN_PRESENTATION", "CHECKPOINT_AXIS",
        "T_CELL_INFLAMED", "CYTOLYTIC_ACTIVITY", "CD8_T_CELL", "EPITHELIAL_CONTENT",
    ]
    labels = [
        "IFN-alpha response", "IFN-gamma response", "Antigen presentation", "Checkpoint axis",
        "T-cell inflamed", "Cytolytic activity", "CD8 T cell", "Epithelial content",
    ]
    t = tcga_cor.set_index("signature").loc[order, "spearman_rho"]
    g = gse_cor.set_index("signature").loc[order, "spearman_rho"]
    y = np.arange(len(order))[::-1]
    for yi, first, second in zip(y, t, g):
        ax_a.plot([first, second], [yi, yi], color=LIGHT_GREY, lw=1.2)
    ax_a.scatter(t, y, s=24, color=BLUE, label=f"TCGA (n={int(tcga_cor['n'].iloc[0])})")
    ax_a.scatter(g, y, s=24, color=TEAL, marker="s", label=f"GSE58812 (n={int(gse_cor['n'].iloc[0])})")
    ax_a.axvline(0, color=GREY, lw=0.7)
    ax_a.set_yticks(y, labels)
    ax_a.set_xlim(-0.4, 1.0)
    ax_a.set_xlabel("Spearman rho with IFI44L")
    ax_a.set_title("Immune-state replication", loc="left", fontweight="bold")
    ax_a.legend(loc="lower right")
    panel_label(ax_a, "a", x=-0.29)

    forest_rows: list[dict[str, object]] = []
    gse_specs = [
        ("Metastasis-free survival", "GSE58812 MFS"),
        ("Overall survival", "GSE58812 OS"),
    ]
    for endpoint, label in gse_specs:
        for model in ["Unadjusted", "Age-adjusted"]:
            row = gse_cox.loc[(gse_cox["endpoint"] == endpoint) & (gse_cox["model"] == model)].iloc[0]
            forest_rows.append({
                "label": label,
                "model": model,
                "cohort": "GSE58812",
                "hr": row["hazard_ratio_per_IFI44L_SD"],
                "low": row["ci95_low"],
                "high": row["ci95_high"],
            })
    for endpoint in ["OS", "PFI"]:
        for model in ["Unadjusted", "Age-adjusted"]:
            row = tcga_cox.loc[(tcga_cox["endpoint"] == endpoint) & (tcga_cox["model"] == model)].iloc[0]
            forest_rows.append({
                "label": f"TCGA {endpoint}",
                "model": model,
                "cohort": "TCGA",
                "hr": row["hazard_ratio_per_IFI44L_SD"],
                "low": row["ci95_low"],
                "high": row["ci95_high"],
            })
    forest = pd.DataFrame(forest_rows)
    label_order = ["GSE58812 MFS", "GSE58812 OS", "TCGA OS", "TCGA PFI"]
    y_positions = np.arange(len(label_order))[::-1]
    for position, label in zip(y_positions, label_order):
        subset = forest.loc[forest["label"] == label].reset_index(drop=True)
        colour = TEAL if subset.loc[0, "cohort"] == "GSE58812" else BLUE
        for index, row in subset.iterrows():
            offset = 0.10 if index == 0 else -0.10
            marker = "o" if index == 0 else "s"
            ax_b.hlines(position + offset, row["low"], row["high"], color=colour, lw=1.0)
            ax_b.scatter(row["hr"], position + offset, s=23, color=colour, marker=marker, zorder=2)
    ax_b.axvline(1, color=GREY, lw=0.7, ls="--")
    ax_b.set_xscale("log")
    ax_b.set_xlim(0.35, 1.60)
    ax_b.set_xticks([0.4, 0.6, 0.8, 1.0, 1.4], ["0.4", "0.6", "0.8", "1.0", "1.4"])
    ax_b.set_yticks(y_positions, label_order)
    ax_b.set_ylim(-0.95, 3.35)
    ax_b.set_xlabel("Hazard ratio per IFI44L s.d. (95% CI)")
    ax_b.set_title("Continuous survival models", loc="left", fontweight="bold")
    ax_b.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLACK, markersize=4.5, label="Unadjusted"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=BLACK, markersize=4.5, label="Adjusted sensitivity"),
    ], loc="lower left", fontsize=6.6, borderpad=0.2)
    panel_label(ax_b, "b", x=-0.28)

    km_results = [
        draw_km(ax_c, risk_axes[0], gse_clinical, "mfs_days", "metastasis_event", "IFI44L_group", "GSE58812 MFS"),
        draw_km(ax_d, risk_axes[1], gse_clinical, "os_days", "death_event", "IFI44L_group", "GSE58812 OS"),
        draw_km(ax_e, risk_axes[2], tcga_clinical, "OS_days", "OS_event", "median_group", "TCGA OS"),
        draw_km(ax_f, risk_axes[3], tcga_clinical, "PFI_days", "PFI_event", "median_group", "TCGA PFI"),
    ]
    gse_mfs_p, gse_os_p, tcga_os_p, tcga_pfi_p = [result[0] for result in km_results]
    for index, name in [(1, "KM_confidence"), (2, "KM_censoring"), (3, "KM_risk_counts")]:
        pd.DataFrame([row for result in km_results for row in result[index]]).to_csv(
            output_dir / f"Figure2_{name}.tsv", sep="\t", index=False)
    panel_label(ax_c, "c", x=-0.20)
    panel_label(ax_d, "d", x=-0.20)
    panel_label(ax_e, "e", x=-0.20)
    panel_label(ax_f, "f", x=-0.20)

    base = output_dir / "Figure2_cross_cohort_replication_survival"
    save_bundle(fig, base)
    return {
        "figure": "Figure 2",
        "base": str(base),
        "panels": ["a", "b", "c", "d", "e", "f"],
        "GSE58812": {"n": len(gse_clinical), "MFS_events": int(gse_clinical["metastasis_event"].sum()),
            "OS_events": int(gse_clinical["death_event"].sum()), "MFS_logrank_p": gse_mfs_p, "OS_logrank_p": gse_os_p},
        "TCGA": {"n": len(tcga_clinical), "OS_events": int(tcga_clinical["OS_event"].sum()),
            "PFI_events": int(tcga_clinical["PFI_event"].sum()), "OS_logrank_p": tcga_os_p, "PFI_logrank_p": tcga_pfi_p},
        "primary_inference": "continuous Cox models; median-split curves are graphical sensitivity analyses",
        "KM_confidence": "Greenwood log-log 95% confidence intervals; censor marks and risk counts shown",
    }


def add_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    lines: list[str],
    edge: str,
    face: str,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.0, edgecolor=edge, facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height - 0.055, title, ha="center", va="top", fontsize=7.2, fontweight="bold", color=edge)
    ax.text(x + 0.024, y + height - 0.135, "\n".join(lines), ha="left", va="top", fontsize=6.6, linespacing=1.35, color=BLACK)


def build_figure5(output_dir: Path) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(7.2, 4.1), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.97, "a", fontsize=8.5, fontweight="bold", va="top")
    ax.text(0.08, 0.95, "Compartment-aware synthesis of public TNBC transcriptomic evidence", fontsize=8.2, fontweight="bold", color=BLUE)

    add_box(ax, (0.04, 0.57), 0.25, 0.27, "Bulk state", [
        "IFI44L tracks IFN and", "antigen-presentation programmes", "across receptor-defined", "TNBC tumours"
    ], BLUE, PALE_BLUE)
    add_box(ax, (0.375, 0.57), 0.25, 0.27, "Cellular source", [
        "highest transcript abundance", "in myeloid, endothelial", "and T-cell compartments;", "lower in malignant epithelium"
    ], TEAL, PALE_TEAL)
    add_box(ax, (0.71, 0.57), 0.25, 0.27, "Spatial niches", [
        "focal IFI44L aggregation", "within multicellular IFN-", "and immune-enriched", "Visium neighbourhoods"
    ], TEAL, PALE_TEAL)
    add_box(ax, (0.18, 0.15), 0.27, 0.26, "Outcome context", [
        "cohort-specific survival", "estimates are interpreted with", "event counts and available", "clinical covariates"
    ], BLUE, PALE_BLUE)
    add_box(ax, (0.55, 0.15), 0.27, 0.26, "Evidence boundary", [
        "tissue-state marker, not", "a causal mechanism, independent", "prognostic biomarker or", "immunotherapy predictor"
    ], GREY, PALE_GREY)

    for start, end, colour in [
        ((0.29, 0.705), (0.375, 0.705), BLUE),
        ((0.625, 0.705), (0.71, 0.705), TEAL),
        ((0.50, 0.57), (0.33, 0.41), BLUE),
        ((0.835, 0.57), (0.69, 0.41), TEAL),
        ((0.45, 0.29), (0.55, 0.29), GREY),
    ]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, linewidth=1.0, color=colour))
    ax.text(
        0.50, 0.06,
        "Interpretation: IFI44L marks a compartment-dependent IFN-inflamed tissue state; outcome associations require further validation.",
        ha="center", va="bottom", fontsize=6.6, fontweight="bold",
    )

    base = output_dir / "Figure5_public_multiomics_conceptual_synthesis"
    save_bundle(fig, base)
    return {
        "figure": "Figure 5",
        "base": str(base),
        "type": "descriptive conceptual synthesis",
        "causal_claim": False,
        "validated_biomarker_claim": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    project = args.project.resolve()
    set_style()
    output_dir = project / "figures" / "v11_public_multiomics"
    output_dir.mkdir(parents=True, exist_ok=True)

    source_dir = project / "figures" / "manuscript_nature_style"
    copied = []
    for number, source_name in [
        (3, "Figure3_single_cell_compartment"),
        (4, "Figure4_spatial_IFI44L_ecology"),
    ]:
        for suffix in [".pdf", ".png", ".svg", ".tiff"]:
            source = source_dir / f"{source_name}{suffix}"
            destination = output_dir / f"Figure{number}{suffix}"
            if not source.exists():
                raise FileNotFoundError(source)
            shutil.copy2(source, destination)
        copied.append({"figure": number, "source": source_name})

    reports = [build_figure1(project, output_dir), build_figure2(project, output_dir), build_figure5(output_dir)]
    audit = {
        "version": "v11_public_multiomics",
        "excluded": ["local patient/IHC", "cell-culture assays", "exploratory treatment-response boundary analysis"],
        "copied_unchanged": copied,
        "rebuilt": reports,
    }
    (output_dir / "figure_build_audit_v11.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
