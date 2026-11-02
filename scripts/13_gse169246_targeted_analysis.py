#!/usr/bin/env python3
"""Patient-level targeted analysis of the GSE169246 CD45+ scRNA-seq matrix."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon


SIGNATURES = {
    "IFN_response_excluding_IFI44L": ["IFI44", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "IFI6", "IFIH1", "DDX58", "STAT1", "IRF7"],
    "B_cell": ["CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD37", "CD74", "HLA-DRA"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "CD8A", "CD8B", "IL7R"],
    "NK_cytotoxic": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMA", "GZMB"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "S100A8", "S100A9"],
    "CXCL13_Tfh_like": ["CXCL13", "PDCD1", "ICOS", "TOX", "CXCR5", "CD4"],
    "Checkpoint": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
}


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    good = np.isfinite(p)
    if good.any():
        vals = p[good]
        order = np.argsort(vals)
        q = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        restored = np.empty_like(q)
        restored[order] = np.clip(q, 0, 1)
        out[np.flatnonzero(good)] = restored
    return pd.Series(out, index=values.index)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    return float((sum(a > b for a in x for b in y) - sum(a < b for a in x for b in y)) / (len(x) * len(y)))


def prepare_targets(features_path: Path, output_path: Path) -> pd.DataFrame:
    features = pd.read_csv(features_path, header=None, names=["gene"], compression="gzip", dtype=str)
    wanted = ["IFI44L"] + sorted(set().union(*SIGNATURES.values()))
    records = []
    for gene in wanted:
        hits = np.flatnonzero(features["gene"].eq(gene).to_numpy())
        if len(hits):
            records.append({"matrix_row_1based": int(hits[0] + 1), "gene": gene})
    target_rows = pd.DataFrame(records)
    if "IFI44L" not in set(target_rows["gene"]):
        raise RuntimeError("IFI44L not found in feature file")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    target_rows.to_csv(output_path, sep="\t", index=False, header=False)
    return target_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--prepare-targets-only", action="store_true")
    args = parser.parse_args()

    result_dir = args.project / "results" / "single_cell" / "GSE169246"
    figure_dir = args.project / "figures" / "single_cell" / "GSE169246"
    result_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    targets_path = result_dir / "target_rows.tsv"
    targets = prepare_targets(args.source_dir / "GSE169246_TNBC_RNA.feature.tsv.gz", targets_path)
    print(f"targets_present={len(targets)}")
    if args.prepare_targets_only:
        return

    metadata = pd.read_csv(result_dir / "GSE169246_cell_metadata.tsv.gz", sep="\t", low_memory=False)
    n_cells = len(metadata)
    library_path = result_dir / "library_sizes_u64.bin"
    counts_path = result_dir / "target_counts_u32.bin"
    expected_library = n_cells * np.dtype("uint64").itemsize
    expected_counts = len(targets) * n_cells * np.dtype("uint32").itemsize
    if library_path.stat().st_size != expected_library or counts_path.stat().st_size != expected_counts:
        raise RuntimeError("Targeted binary sizes do not match metadata/target dimensions")
    library = np.memmap(library_path, dtype="uint64", mode="r", shape=(n_cells,))
    counts = np.memmap(counts_path, dtype="uint32", mode="r", shape=(len(targets), n_cells))
    if np.any(library == 0):
        raise RuntimeError("Zero library sizes found")
    log_expr = np.log1p(counts.astype(np.float64) * (1e4 / library)[None, :])
    gene_to_row = {gene: i for i, gene in enumerate(targets["gene"])}

    cell_values = pd.DataFrame({"sample_id": metadata["sample_id"], "IFI44L": log_expr[gene_to_row["IFI44L"]]})
    cell_values["IFI44L_detected"] = (counts[gene_to_row["IFI44L"]] > 0).astype(np.int8)
    for signature, genes in SIGNATURES.items():
        rows = [gene_to_row[g] for g in genes if g in gene_to_row]
        cell_values[signature] = log_expr[rows].mean(axis=0)

    sample_values = (
        cell_values.groupby("sample_id", as_index=False)
        .agg(
            n_cells=("IFI44L", "size"),
            IFI44L_mean_logCP10k=("IFI44L", "mean"),
            IFI44L_detection_fraction=("IFI44L_detected", "mean"),
            **{f"{name}_mean": (name, "mean") for name in SIGNATURES},
        )
        .merge(metadata.drop_duplicates("sample_id").drop(columns=["cell_index_1based", "cell_barcode"]), on="sample_id", how="left", validate="one_to_one")
    )
    sample_values.to_csv(result_dir / "GSE169246_sample_targeted_expression.tsv", sep="\t", index=False)

    endpoints = ["IFI44L_mean_logCP10k", "IFI44L_detection_fraction", "IFN_response_excluding_IFI44L_mean", "B_cell_mean", "CXCL13_Tfh_like_mean"]
    baseline_tests = []
    anti = sample_values.loc[
        sample_values["treatment_arm"].str.contains("Anti-PD-L1", na=False)
        & sample_values["timepoint"].eq("Pre")
    ].copy()
    for source, group in anti.groupby("sample_source"):
        for endpoint in endpoints:
            x = group.loc[group["response_binary"].eq("Responder"), endpoint].dropna().to_numpy()
            y = group.loc[group["response_binary"].eq("Non-responder"), endpoint].dropna().to_numpy()
            if len(x) >= 2 and len(y) >= 2:
                u, p = mannwhitneyu(x, y, alternative="two-sided")
                baseline_tests.append({
                    "sample_source": source, "endpoint": endpoint,
                    "n_responders": len(x), "n_nonresponders": len(y),
                    "median_responders": np.median(x), "median_nonresponders": np.median(y),
                    "cliffs_delta_responder_minus_nonresponder": cliffs_delta(x, y),
                    "mannwhitney_u": u, "p_value": p,
                })
    baseline_tests = pd.DataFrame(baseline_tests)
    baseline_tests["fdr_bh"] = bh_adjust(baseline_tests["p_value"])
    baseline_tests.to_csv(result_dir / "GSE169246_antiPDL1_baseline_response_tests.tsv", sep="\t", index=False)

    paired_rows = []
    key = ["patient_id", "sample_source", "treatment_arm", "response_binary"]
    for endpoint in endpoints:
        wide = sample_values.loc[sample_values["timepoint"].isin(["Pre", "Post"])].pivot_table(index=key, columns="timepoint", values=endpoint, aggfunc="first").dropna()
        wide["change_Post_minus_Pre"] = wide["Post"] - wide["Pre"]
        for keys, group in wide.reset_index().groupby(["treatment_arm", "sample_source", "response_binary"]):
            values = group["change_Post_minus_Pre"].to_numpy()
            if len(values) >= 2:
                statistic, p = wilcoxon(values, alternative="two-sided", method="auto") if not np.allclose(values, 0) else (0, 1)
            else:
                statistic, p = np.nan, np.nan
            paired_rows.append({
                "endpoint": endpoint, "treatment_arm": keys[0], "sample_source": keys[1], "response_binary": keys[2],
                "n_paired_patients": len(values), "median_change_Post_minus_Pre": np.median(values),
                "wilcoxon_statistic": statistic, "p_value": p,
            })
    paired = pd.DataFrame(paired_rows)
    paired["fdr_bh"] = bh_adjust(paired["p_value"])
    paired.to_csv(result_dir / "GSE169246_paired_pre_post_changes.tsv", sep="\t", index=False)

    plot = anti.loc[anti["sample_source"].eq("tumor")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    for ax, endpoint, ylabel in [
        (axes[0], "IFI44L_mean_logCP10k", "IFI44L mean log(CP10K+1)"),
        (axes[1], "IFN_response_excluding_IFI44L_mean", "IFN module mean"),
    ]:
        groups = ["Responder", "Non-responder"]
        rng = np.random.default_rng(20260902)
        for i, label in enumerate(groups):
            vals = plot.loc[plot["response_binary"].eq(label), endpoint].to_numpy()
            ax.scatter(i + rng.uniform(-.08, .08, len(vals)), vals, s=42, color="#315A7D")
            ax.plot([i - .18, i + .18], [np.median(vals)] * 2, color="#B23A48", lw=2)
        ax.set_xticks([0, 1], groups)
        ax.set_ylabel(ylabel)
        ax.set_title("Baseline tumor CD45+ cells")
    fig.suptitle("GSE169246 anti-PD-L1 + chemotherapy cohort (patient-level exploratory)")
    fig.tight_layout()
    fig.savefig(figure_dir / "GSE169246_baseline_tumor_IFI44L_IFN_by_response.pdf")
    fig.savefig(figure_dir / "GSE169246_baseline_tumor_IFI44L_IFN_by_response.png", dpi=400)
    plt.close(fig)

    audit = {
        "dataset": "GSE169246",
        "matrix_shape": [int(len(targets)), int(n_cells)],
        "normalization": "log1p(raw count / cell total UMI * 10000)",
        "inference_unit": "patient/sample; no cell-level inferential tests",
        "scope": "CD45+ immune cells; not a whole-tumor cell atlas",
        "response_definition": "Responder=CR/PR; Non-responder=SD/PD from original supplementary clinical table",
        "python": platform.python_version(),
        "targets": targets["gene"].tolist(),
    }
    (result_dir / "analysis_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print("Baseline response tests")
    print(baseline_tests.to_string(index=False))
    print("\nPaired changes")
    print(paired.to_string(index=False))


if __name__ == "__main__":
    main()
