#!/usr/bin/env python3
"""Donor-aware IFI44L single-cell analysis in the 10 TNBC tumors from GSE176078."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path

import anndata as ad
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-codex-ifi44l")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon
import seaborn as sns


DEFAULT_PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = DEFAULT_PROJECT / "source_metadata" / "GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad"

GENE_SETS = {
    "interferon_response_excluding_IFI44L": [
        "IFI44", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2",
        "OAS3", "OASL", "RSAD2", "IFI6", "IFIH1", "DDX58", "DDX60", "HERC5",
        "HERC6", "CMPK2", "XAF1", "SAMD9L", "STAT1", "IRF7",
    ],
    "cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB"],
    "T_cell_inflamed": ["CD3D", "CD3E", "CD8A", "CXCL9", "CXCL10", "CXCL11", "IFNG", "GZMA", "PRF1"],
    "checkpoint_axis": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9", "TOX"],
}

MARKERS = [
    "IFI44L", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CD274", "PDCD1", "PDCD1LG2",
    "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9", "EPCAM", "KRT5",
    "KRT8", "KRT14", "KRT18", "KRT19", "EGFR", "CD68", "LST1", "C1QA", "C1QB",
    "C1QC", "MRC1", "MS4A1", "CD79A", "COL1A1", "COL1A2", "PECAM1", "VWF",
]

EPITHELIAL_TYPES = {
    "mammary gland epithelial cell",
    "luminal epithelial cell of mammary gland",
    "basal-myoepithelial cell of mammary gland",
}


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    good = np.isfinite(p)
    if good.any():
        vals = p[good]
        order = np.argsort(vals)
        ranked = vals[order]
        q = ranked * len(vals) / np.arange(1, len(vals) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        restored = np.empty_like(q)
        restored[order] = q
        out[np.flatnonzero(good)] = restored
    return pd.Series(out, index=pvalues.index)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def compartment(cell_type: str) -> str:
    if cell_type in EPITHELIAL_TYPES:
        return "Epithelial"
    if any(x in cell_type.lower() for x in ["t cell", "natural killer", "b cell", "plasmablast"]):
        return "Lymphoid"
    if any(x in cell_type.lower() for x in ["macrophage", "monocyte", "myeloid"]):
        return "Myeloid"
    return "Stromal/vascular"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--min-cells", type=int, default=20)
    args = parser.parse_args()

    out = args.project / "results" / "single_cell"
    fig = args.project / "figures" / "single_cell"
    logs = args.project / "logs"
    out.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    tnbc = pd.read_csv(out / "Wu2021_TNBC_donors.tsv", sep="\t")["donor_id"].tolist()
    source = ad.read_h5ad(args.h5ad, backed="r")
    required_obs = {"donor_id", "cell_type"}
    missing_obs = required_obs.difference(source.obs.columns)
    if missing_obs:
        raise RuntimeError(f"Missing required obs columns: {sorted(missing_obs)}")

    if source.raw is None:
        raise RuntimeError("Source H5AD has no .raw matrix; integer-count analysis cannot proceed")
    raw_var = source.raw.var
    feature_name = raw_var["feature_name"].astype(str) if "feature_name" in raw_var else pd.Series(raw_var.index, index=raw_var.index)
    all_targets = sorted(set(MARKERS).union(*GENE_SETS.values()))
    present_mask = feature_name.isin(all_targets).to_numpy()
    present = feature_name.loc[present_mask].tolist()
    missing = sorted(set(all_targets).difference(present))
    if "IFI44L" not in present:
        raise RuntimeError("IFI44L is absent from the source H5AD")

    cell_mask = source.obs["donor_id"].astype(str).isin(tnbc).to_numpy()
    # In this Discover asset .X is normalized; integer UMI counts are stored in .raw.X.
    # Compute library size across every raw gene before taking the targeted gene slice.
    full_matrix = source.raw[cell_mask, :].X
    full_library_sum = np.asarray(full_matrix.sum(axis=1)).ravel()
    library_sum_source = "row sum across all genes in source .raw.X"
    # anndata.Raw boolean var slicing can retain all columns in backed mode; slice the
    # materialized sparse matrix explicitly to guarantee the intended dimensions.
    full_matrix = full_matrix.tocsr() if sp.issparse(full_matrix) else sp.csr_matrix(full_matrix)
    counts = full_matrix[:, present_mask].tocsr()
    selected_obs = source.obs.loc[cell_mask].copy()
    selected_var = raw_var.loc[present_mask].copy()
    data = ad.AnnData(X=counts.copy(), obs=selected_obs, var=selected_var)
    data.var_names = feature_name.loc[present_mask].to_numpy()
    data.var_names_make_unique()
    data.obs["donor_id"] = data.obs["donor_id"].astype(str)
    data.obs["cell_type"] = data.obs["cell_type"].astype(str)
    data.obs["compartment"] = data.obs["cell_type"].map(compartment)

    if counts.data.size and (counts.data < 0).any():
        raise RuntimeError("Negative values found in matrix expected to contain raw counts")
    if counts.data.size and not np.allclose(counts.data, np.round(counts.data)):
        raise RuntimeError("Non-integer values found in source .raw.X")
    raw_sum = full_library_sum
    if (~np.isfinite(raw_sum)).any() or (raw_sum <= 0).any():
        raise RuntimeError("Invalid raw_sum values prevent whole-library normalization")
    logcp10k = counts.multiply((1e4 / raw_sum)[:, None]).tocsr()
    logcp10k.data = np.log1p(logcp10k.data)
    data.layers["counts"] = counts.copy()
    data.X = logcp10k

    expr = pd.DataFrame.sparse.from_spmatrix(logcp10k, index=data.obs_names, columns=data.var_names)
    count_df = pd.DataFrame.sparse.from_spmatrix(counts, index=data.obs_names, columns=data.var_names)
    for score_name, genes in GENE_SETS.items():
        use = [g for g in genes if g in expr.columns]
        data.obs[score_name] = np.asarray(expr[use].mean(axis=1)).ravel() if use else np.nan

    # Small reusable TNBC-only object: normalized X plus raw targeted counts.
    data.write_h5ad(out / "GSE176078_TNBC_targeted_IFI44L_immune.h5ad", compression="gzip")

    meta_cols = ["donor_id", "cell_type", "compartment"] + list(GENE_SETS)
    frame = data.obs[meta_cols].copy()
    frame["IFI44L_logCP10k"] = expr["IFI44L"].sparse.to_dense().to_numpy()
    frame["IFI44L_detected"] = (count_df["IFI44L"].sparse.to_dense().to_numpy() > 0).astype(int)

    donor_type = (
        frame.groupby(["donor_id", "cell_type", "compartment"], observed=True)
        .agg(
            n_cells=("IFI44L_detected", "size"),
            IFI44L_detection_fraction=("IFI44L_detected", "mean"),
            IFI44L_mean_logCP10k=("IFI44L_logCP10k", "mean"),
            **{f"{name}_mean": (name, "mean") for name in GENE_SETS},
        )
        .reset_index()
    )
    donor_type.to_csv(out / "GSE176078_TNBC_IFI44L_by_donor_celltype.tsv", sep="\t", index=False)

    eligible = donor_type.loc[donor_type["n_cells"] >= args.min_cells].copy()
    celltype_summary = (
        eligible.groupby(["cell_type", "compartment"], observed=True)
        .agg(
            n_donors=("donor_id", "nunique"),
            total_cells=("n_cells", "sum"),
            median_detection_fraction=("IFI44L_detection_fraction", "median"),
            iqr_detection_fraction=("IFI44L_detection_fraction", lambda x: x.quantile(.75) - x.quantile(.25)),
            median_mean_logCP10k=("IFI44L_mean_logCP10k", "median"),
            iqr_mean_logCP10k=("IFI44L_mean_logCP10k", lambda x: x.quantile(.75) - x.quantile(.25)),
        )
        .reset_index()
        .sort_values("median_mean_logCP10k", ascending=False)
    )
    celltype_summary.to_csv(out / "GSE176078_TNBC_IFI44L_celltype_summary.tsv", sep="\t", index=False)

    # Preserve the author's biologically richer annotation hierarchy.
    author_levels = [x for x in ["celltype_major", "celltype_minor", "celltype_subset"] if x in data.obs]
    author_outputs = {}
    for level in author_levels:
        author_frame = data.obs[["donor_id", level]].copy()
        author_frame["IFI44L_logCP10k"] = frame["IFI44L_logCP10k"].to_numpy()
        author_frame["IFI44L_detected"] = frame["IFI44L_detected"].to_numpy()
        grouped = (
            author_frame.groupby(["donor_id", level], observed=True)
            .agg(
                n_cells=("IFI44L_detected", "size"),
                IFI44L_detection_fraction=("IFI44L_detected", "mean"),
                IFI44L_mean_logCP10k=("IFI44L_logCP10k", "mean"),
            )
            .reset_index()
        )
        grouped.to_csv(out / f"GSE176078_TNBC_IFI44L_by_donor_{level}.tsv", sep="\t", index=False)
        grouped_eligible = grouped.loc[grouped["n_cells"] >= args.min_cells].copy()
        summary = (
            grouped_eligible.groupby(level, observed=True)
            .agg(
                n_donors=("donor_id", "nunique"),
                total_cells=("n_cells", "sum"),
                median_detection_fraction=("IFI44L_detection_fraction", "median"),
                median_mean_logCP10k=("IFI44L_mean_logCP10k", "median"),
            )
            .reset_index()
            .sort_values("median_mean_logCP10k", ascending=False)
        )
        summary.to_csv(out / f"GSE176078_TNBC_IFI44L_{level}_summary.tsv", sep="\t", index=False)
        author_outputs[level] = (grouped_eligible, summary)

    # Paired tests at the author's major-cell-class level, with malignant epithelium
    # as the biologically relevant reference.
    author_major_tests = []
    if "celltype_major" in author_outputs:
        major_eligible, _ = author_outputs["celltype_major"]
        major_reference = "Cancer Epithelial"
        for cell_class in sorted(major_eligible["celltype_major"].astype(str).unique()):
            if cell_class == major_reference:
                continue
            for endpoint in ["IFI44L_mean_logCP10k", "IFI44L_detection_fraction"]:
                wide = major_eligible.loc[major_eligible["celltype_major"].isin([major_reference, cell_class])].pivot(
                    index="donor_id", columns="celltype_major", values=endpoint
                ).dropna()
                if len(wide) < 5:
                    continue
                diff = wide[cell_class] - wide[major_reference]
                pvalue = float(wilcoxon(diff, alternative="two-sided", method="auto").pvalue) if np.any(diff != 0) else 1.0
                author_major_tests.append({
                    "cell_class": cell_class,
                    "reference": major_reference,
                    "endpoint": endpoint,
                    "n_paired_donors": len(wide),
                    "median_paired_difference": float(np.median(diff)),
                    "wilcoxon_p": pvalue,
                })
    author_major_tests = pd.DataFrame(author_major_tests)
    if len(author_major_tests):
        author_major_tests["fdr_bh"] = author_major_tests.groupby("endpoint")["wilcoxon_p"].transform(bh_adjust)
    author_major_tests.to_csv(out / "GSE176078_TNBC_author_major_paired_tests.tsv", sep="\t", index=False)

    # Paired donor-level sensitivity tests against the dominant epithelial category.
    reference = "mammary gland epithelial cell"
    tests = []
    for cell_type in sorted(eligible["cell_type"].unique()):
        if cell_type == reference:
            continue
        for endpoint in ["IFI44L_mean_logCP10k", "IFI44L_detection_fraction"]:
            wide = eligible.loc[eligible["cell_type"].isin([reference, cell_type])].pivot(
                index="donor_id", columns="cell_type", values=endpoint
            ).dropna()
            if len(wide) < 5:
                continue
            diff = wide[cell_type] - wide[reference]
            try:
                pvalue = float(wilcoxon(diff, alternative="two-sided", method="auto").pvalue) if np.any(diff != 0) else 1.0
            except ValueError:
                pvalue = np.nan
            tests.append({
                "cell_type": cell_type,
                "reference": reference,
                "endpoint": endpoint,
                "n_paired_donors": len(wide),
                "median_paired_difference": float(np.median(diff)),
                "wilcoxon_p": pvalue,
            })
    tests = pd.DataFrame(tests)
    if len(tests):
        tests["fdr_bh"] = tests.groupby("endpoint")["wilcoxon_p"].transform(bh_adjust)
    tests.to_csv(out / "GSE176078_TNBC_celltype_paired_tests.tsv", sep="\t", index=False)

    # Aggregate malignant epithelial cells within each donor; donor, not cell, is the inference unit.
    if "celltype_major" in data.obs:
        epi_mask = data.obs["celltype_major"].astype(str).eq("Cancer Epithelial").to_numpy()
        epithelial_definition = "author celltype_major == Cancer Epithelial"
    else:
        epi_mask = data.obs["cell_type"].isin(EPITHELIAL_TYPES).to_numpy()
        epithelial_definition = "standardized epithelial cell types"
    epi = data.obs.loc[epi_mask, ["donor_id"] + list(GENE_SETS)].copy()
    epi_expr = expr.loc[epi_mask]
    for gene in ["IFI44L", "CD274", "PDCD1LG2", "LAG3", "CXCL10", "IFNG"]:
        if gene in epi_expr:
            epi[f"{gene}_mean_logCP10k"] = epi_expr[gene].sparse.to_dense().to_numpy()
    epi_donor = epi.groupby("donor_id", observed=True).mean(numeric_only=True).reset_index()
    epi_counts = data.obs.loc[epi_mask].groupby("donor_id", observed=True).size().rename("n_epithelial_cells")
    epi_donor = epi_donor.merge(epi_counts, on="donor_id", how="left")
    epi_donor.to_csv(out / "GSE176078_TNBC_cancer_epithelial_donor_summaries.tsv", sep="\t", index=False)

    correlations = []
    xcol = "IFI44L_mean_logCP10k"
    for ycol in [c for c in epi_donor.columns if c not in {"donor_id", "n_epithelial_cells", xcol}]:
        xy = epi_donor[[xcol, ycol]].dropna()
        if len(xy) >= 6 and xy[ycol].nunique() > 1:
            rho, pvalue = spearmanr(xy[xcol], xy[ycol])
            correlations.append({"outcome": ycol, "n_donors": len(xy), "spearman_rho": rho, "p_value": pvalue})
    correlations = pd.DataFrame(correlations)
    if len(correlations):
        correlations["fdr_bh"] = bh_adjust(correlations["p_value"])
    correlations.to_csv(out / "GSE176078_TNBC_cancer_epithelial_IFI44L_correlations.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    plot_summary = celltype_summary.loc[celltype_summary["n_donors"] >= 3].copy().sort_values("median_mean_logCP10k")
    plot_summary["plot_label"] = plot_summary.apply(lambda r: f"{r['cell_type']} (n={int(r['n_donors'])})", axis=1)
    fig1, ax = plt.subplots(figsize=(8.2, max(5.5, .34 * len(plot_summary))))
    sizes = 40 + 900 * plot_summary["median_detection_fraction"].clip(0, 1)
    colors = plot_summary["compartment"].map({"Epithelial": "#D73027", "Lymphoid": "#4575B4", "Myeloid": "#7B3294", "Stromal/vascular": "#1A9850"})
    ax.scatter(plot_summary["median_mean_logCP10k"], plot_summary["plot_label"], s=sizes, c=colors, alpha=.82, edgecolor="black", linewidth=.35)
    ax.set_xlabel("Median donor-level IFI44L mean log(CP10K + 1)")
    ax.set_ylabel("")
    ax.set_title("IFI44L localization across TNBC cell types (donor-level summaries)")
    color_map = {"Epithelial": "#D73027", "Lymphoid": "#4575B4", "Myeloid": "#7B3294", "Stromal/vascular": "#1A9850"}
    color_handles = [Line2D([0], [0], marker="o", color="w", label=k, markerfacecolor=v, markeredgecolor="black", markersize=7) for k, v in color_map.items()]
    legend1 = ax.legend(handles=color_handles, title="Compartment", loc="lower right", frameon=True)
    ax.add_artist(legend1)
    size_handles = [plt.scatter([], [], s=40 + 900 * v, color="grey", alpha=.7, edgecolor="black", linewidth=.3, label=f"{v:.0%}") for v in (.1, .3, .5)]
    ax.legend(handles=size_handles, title="Median detection", loc="center right", frameon=True)
    fig1.tight_layout()
    fig1.savefig(fig / "GSE176078_TNBC_IFI44L_celltype_dotplot.pdf")
    fig1.savefig(fig / "GSE176078_TNBC_IFI44L_celltype_dotplot.png", dpi=400)
    plt.close(fig1)

    heat = eligible.pivot(index="cell_type", columns="donor_id", values="IFI44L_mean_logCP10k")
    heat = heat.loc[heat.mean(axis=1).sort_values(ascending=False).index]
    fig2, ax = plt.subplots(figsize=(10.8, max(5.8, .34 * len(heat))))
    sns.heatmap(heat, cmap="mako", linewidths=.25, linecolor="white", ax=ax, cbar_kws={"label": "Mean log(CP10K + 1)"})
    ax.set_xlabel("TNBC donor")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    fig2.tight_layout()
    fig2.savefig(fig / "GSE176078_TNBC_IFI44L_donor_celltype_heatmap.pdf")
    fig2.savefig(fig / "GSE176078_TNBC_IFI44L_donor_celltype_heatmap.png", dpi=400)
    plt.close(fig2)

    if "celltype_subset" in author_outputs:
        _, subset_summary = author_outputs["celltype_subset"]
        subset_plot = subset_summary.loc[subset_summary["n_donors"] >= 3].head(30).sort_values("median_mean_logCP10k")
        subset_plot["plot_label"] = subset_plot.apply(lambda r: f"{r['celltype_subset']} (n={int(r['n_donors'])})", axis=1)
        fig_sub, ax = plt.subplots(figsize=(8.8, max(6.2, .30 * len(subset_plot))))
        ax.scatter(
            subset_plot["median_mean_logCP10k"], subset_plot["plot_label"],
            s=35 + 800 * subset_plot["median_detection_fraction"].clip(0, 1),
            c="#6A3D9A", alpha=.82, edgecolor="black", linewidth=.3,
        )
        ax.set_xlabel("Median donor-level IFI44L mean log(CP10K + 1)")
        ax.set_ylabel("")
        ax.set_title("Top IFI44L-expressing author-defined subsets in TNBC")
        subset_size_handles = [plt.scatter([], [], s=35 + 800 * v, color="#6A3D9A", alpha=.7, edgecolor="black", linewidth=.3, label=f"{v:.0%}") for v in (.2, .5, .8)]
        ax.legend(handles=subset_size_handles, title="Median detection", loc="lower right", frameon=True)
        fig_sub.tight_layout()
        fig_sub.savefig(fig / "GSE176078_TNBC_IFI44L_author_subset_top30.pdf")
        fig_sub.savefig(fig / "GSE176078_TNBC_IFI44L_author_subset_top30.png", dpi=400)
        plt.close(fig_sub)

    panel_outcomes = [
        "interferon_response_excluding_IFI44L",
        "CD274_mean_logCP10k",
        "PDCD1LG2_mean_logCP10k",
        "CXCL10_mean_logCP10k",
    ]
    panel_outcomes = [x for x in panel_outcomes if x in epi_donor]
    fig3, axes = plt.subplots(1, len(panel_outcomes), figsize=(4.1 * len(panel_outcomes), 3.8), squeeze=False)
    for ax, ycol in zip(axes.ravel(), panel_outcomes):
        sns.regplot(data=epi_donor, x=xcol, y=ycol, ci=None, scatter_kws={"s": 45}, line_kws={"color": "#B2182B"}, ax=ax)
        rho, pvalue = spearmanr(epi_donor[xcol], epi_donor[ycol])
        qmatch = correlations.loc[correlations["outcome"].eq(ycol), "fdr_bh"]
        qvalue = float(qmatch.iloc[0]) if len(qmatch) else np.nan
        ax.set_title(f"Spearman ρ={rho:.2f}, P={pvalue:.3g}, q={qvalue:.3g}")
        ax.set_xlabel("Cancer epithelial IFI44L")
        ax.set_ylabel(ycol.replace("_mean_logCP10k", "").replace("_", " "))
    fig3.suptitle("Exploratory cancer-epithelial associations (8 TNBC donors)", y=.99, fontsize=11)
    fig3.tight_layout(rect=(0, 0, 1, .92))
    fig3.savefig(fig / "GSE176078_TNBC_cancer_epithelial_IFI44L_associations.pdf")
    fig3.savefig(fig / "GSE176078_TNBC_cancer_epithelial_IFI44L_associations.png", dpi=400)
    plt.close(fig3)

    audit = {
        "source_h5ad": str(args.h5ad),
        "source_sha256": sha256(args.h5ad),
        "n_tnbc_donors": int(data.obs["donor_id"].nunique()),
        "n_tnbc_cells": int(data.n_obs),
        "n_cell_types": int(data.obs["cell_type"].nunique()),
        "target_genes_present": present,
        "target_genes_missing": missing,
        "normalization": "log1p(raw count / cell raw_sum * 10000)",
        "library_sum_source": library_sum_source,
        "cancer_epithelial_definition": epithelial_definition,
        "minimum_cells_per_donor_celltype_for_inference": args.min_cells,
        "inference_unit": "donor",
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
    }
    (out / "GSE176078_TNBC_single_cell_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (logs / "single_cell_python_versions.txt").write_text(
        f"python={platform.python_version()}\nscanpy={sc.__version__}\nanndata={ad.__version__}\n",
        encoding="utf-8",
    )
    print(json.dumps({k: audit[k] for k in ["n_tnbc_donors", "n_tnbc_cells", "n_cell_types", "target_genes_missing"]}, indent=2))


if __name__ == "__main__":
    main()
