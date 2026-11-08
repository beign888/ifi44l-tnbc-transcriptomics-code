#!/usr/bin/env python3
"""Patient-aware spatial validation of IFI44L in 28 TNBC Visium sections.

Primary endpoints are prespecified as within-section correlations of IFI44L
with immune/vascular, interferon, and epithelial module scores. Two sections
from the same tumor are averaged before patient-level inference.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon


SIGNATURES = {
    "IFN_response_excluding_IFI44L": [
        "IFI44", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2",
        "OAS3", "OASL", "RSAD2", "IFI6", "IFIH1", "DDX58", "STAT1", "IRF7",
    ],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "C1QA", "C1QB", "C1QC"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "CD8A", "IL7R"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG", "RAMP2"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "EGFR"],
    "Cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB"],
    "Checkpoint": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
}


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    good = np.isfinite(p)
    if good.any():
        vals = p[good]
        order = np.argsort(vals)
        ranked = vals[order]
        q = ranked * len(vals) / np.arange(1, len(vals) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        restored = np.empty_like(q)
        restored[order] = np.clip(q, 0, 1)
        out[np.flatnonzero(good)] = restored
    return pd.Series(out, index=values.index)


def read_positions(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        pos = pd.read_csv(
            handle,
            header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres"],
        )
    return pos.set_index("barcode")


def find_one(sample_dir: Path, pattern: str) -> Path:
    # Removable macOS volumes may contain AppleDouble metadata sidecars.
    hits = [path for path in sample_dir.glob(pattern)
            if path.is_file() and not path.name.startswith("._") and path.name != ".DS_Store"]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one {pattern} under {sample_dir}, found {len(hits)}")
    return hits[0]


def moran_i(values: np.ndarray, rows: np.ndarray, cols: np.ndarray, permutations: int, seed: int) -> tuple[float, float, int]:
    coord_to_i = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    offsets = [(0, 2), (1, 1), (1, -1)]
    edge_i, edge_j = [], []
    for i, (r, c) in enumerate(zip(rows, cols)):
        for dr, dc in offsets:
            j = coord_to_i.get((int(r + dr), int(c + dc)))
            if j is not None:
                edge_i.append(i)
                edge_j.append(j)
    edge_i = np.asarray(edge_i, dtype=int)
    edge_j = np.asarray(edge_j, dtype=int)
    if len(edge_i) == 0 or np.nanstd(values) == 0:
        return np.nan, np.nan, len(edge_i)
    z = values - np.mean(values)
    denominator = np.sum(z * z)
    observed = len(values) / (2 * len(edge_i)) * (2 * np.sum(z[edge_i] * z[edge_j])) / denominator
    rng = np.random.default_rng(seed)
    permuted = np.empty(permutations)
    for k in range(permutations):
        zp = rng.permutation(z)
        permuted[k] = len(values) / (2 * len(edge_i)) * (2 * np.sum(zp[edge_i] * zp[edge_j])) / denominator
    p = (1 + np.sum(np.abs(permuted) >= abs(observed))) / (permutations + 1)
    return float(observed), float(p), len(edge_i)


def load_section(sample_dir: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    matrix_path = find_one(sample_dir, "*_filtered_feature_bc_matrix.h5")
    positions_path = find_one(sample_dir, "*_tissue_positions_list.csv.gz")
    adata = sc.read_10x_h5(matrix_path, gex_only=True)
    adata.var_names_make_unique()
    positions = read_positions(positions_path)
    shared = adata.obs_names.intersection(positions.index)
    adata = adata[shared].copy()
    positions = positions.loc[shared]

    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    total = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    mt = np.asarray(counts[:, mt_mask].sum(axis=1)).ravel()
    pct_mt = np.divide(mt, total, out=np.zeros_like(total, dtype=float), where=total > 0) * 100
    keep = (total >= 500) & (n_genes >= 200) & (pct_mt < 20) & positions["in_tissue"].eq(1).to_numpy()
    adata = adata[keep].copy()
    positions = positions.iloc[np.flatnonzero(keep)].copy()
    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    total = np.asarray(counts.sum(axis=1)).ravel()

    needed = sorted({"IFI44L"}.union(*SIGNATURES.values()))
    present = [g for g in needed if g in adata.var_names]
    selected = counts[:, adata.var_names.isin(present)].astype(np.float64)
    selected = selected.multiply((1e4 / total)[:, None]).tocsr()
    selected.data = np.log1p(selected.data)
    selected_genes = adata.var_names[adata.var_names.isin(present)].tolist()
    expression = pd.DataFrame(selected.toarray(), index=adata.obs_names, columns=selected_genes)

    out = positions.copy()
    out["IFI44L"] = expression["IFI44L"].to_numpy() if "IFI44L" in expression else 0.0
    for name, genes in SIGNATURES.items():
        use = [g for g in genes if g in expression]
        out[name] = expression[use].mean(axis=1).to_numpy() if use else np.nan
    qc = {
        "n_spots_before_qc": int(len(keep)),
        "n_spots_after_qc": int(keep.sum()),
        "median_umi_after_qc": float(np.median(total)),
        "IFI44L_detection_fraction": float((out["IFI44L"] > 0).mean()),
        "genes_present": present,
        "genes_missing": sorted(set(needed).difference(present)),
    }
    return out, qc


def signed_rank_against_zero(values: pd.Series) -> tuple[float, float]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if len(clean) == 0:
        return np.nan, np.nan
    if np.allclose(clean, 0):
        return 0.0, 1.0
    test = wilcoxon(clean, alternative="two-sided", method="auto")
    return float(test.statistic), float(test.pvalue)


def partial_spearman(frame: pd.DataFrame, x: str, y: str, controls: list[str]) -> float:
    """Partial Spearman correlation by residualizing ranked variables."""
    use = frame[[x, y] + controls].dropna().rank(method="average")
    if len(use) < 20 or use[x].nunique() < 2 or use[y].nunique() < 2:
        return np.nan
    design = np.column_stack([np.ones(len(use)), use[controls].to_numpy(float)])
    x_resid = use[x].to_numpy(float) - design @ np.linalg.lstsq(design, use[x].to_numpy(float), rcond=None)[0]
    y_resid = use[y].to_numpy(float) - design @ np.linalg.lstsq(design, use[y].to_numpy(float), rcond=None)[0]
    if np.std(x_resid) == 0 or np.std(y_resid) == 0:
        return np.nan
    return float(np.corrcoef(x_resid, y_resid)[0, 1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=499)
    args = parser.parse_args()

    source = args.data_dir
    meta = pd.read_csv(source / "GSE210616_primary28_samples.tsv", sep="\t")
    result_dir = args.project / "results" / "spatial" / "GSE210616"
    figure_dir = args.project / "figures" / "spatial" / "GSE210616"
    log_dir = args.project / "logs"
    result_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    section_frames: dict[str, pd.DataFrame] = {}
    section_rows, correlation_rows, partial_rows = [], [], []
    for index, row in meta.iterrows():
        gsm = row["gsm"]
        frame, qc = load_section(source / "files" / gsm)
        section_frames[gsm] = frame
        moran, moran_p, n_edges = moran_i(
            frame["IFI44L"].to_numpy(), frame["array_row"].to_numpy(), frame["array_col"].to_numpy(),
            args.permutations, seed=20260902 + index,
        )
        section_rows.append({
            "gsm": gsm, "patient_id": row["patient_id"], "section": int(row["section"]),
            **{k: v for k, v in qc.items() if k not in {"genes_present", "genes_missing"}},
            "IFI44L_mean_logCP10k": float(frame["IFI44L"].mean()),
            "IFI44L_moran_I": moran, "moran_permutation_p": moran_p, "n_spatial_edges": n_edges,
        })
        for signature in SIGNATURES:
            valid = frame[["IFI44L", signature]].dropna()
            rho, p = spearmanr(valid["IFI44L"], valid[signature]) if valid[signature].nunique() > 1 else (np.nan, np.nan)
            correlation_rows.append({
                "gsm": gsm, "patient_id": row["patient_id"], "section": int(row["section"]),
                "signature": signature, "n_spots": len(valid), "spearman_rho": rho, "spot_level_p_descriptive": p,
            })
        partial_specs = [
            (
                "IFN_adjusted_for_lineage_modules",
                "IFN_response_excluding_IFI44L",
                ["Myeloid", "T_cell", "Endothelial", "Epithelial"],
            ),
            (
                "Checkpoint_adjusted_for_IFN_and_lineage_modules",
                "Checkpoint",
                ["IFN_response_excluding_IFI44L", "Myeloid", "T_cell", "Endothelial", "Epithelial"],
            ),
        ]
        for analysis, outcome, controls in partial_specs:
            partial_rows.append({
                "gsm": gsm,
                "patient_id": row["patient_id"],
                "section": int(row["section"]),
                "analysis": analysis,
                "outcome": outcome,
                "controls": ";".join(controls),
                "partial_spearman_rho": partial_spearman(frame, "IFI44L", outcome, controls),
            })

    sections = pd.DataFrame(section_rows)
    sections["moran_fdr_bh"] = bh_adjust(sections["moran_permutation_p"])
    correlations = pd.DataFrame(correlation_rows)
    patient_corr = correlations.groupby(["patient_id", "signature"], as_index=False)["spearman_rho"].mean()
    tests = []
    for signature, group in patient_corr.groupby("signature"):
        statistic, p = signed_rank_against_zero(group["spearman_rho"])
        tests.append({
            "signature": signature, "n_patients": group["patient_id"].nunique(),
            "median_patient_rho": group["spearman_rho"].median(),
            "iqr_patient_rho": group["spearman_rho"].quantile(.75) - group["spearman_rho"].quantile(.25),
            "wilcoxon_statistic": statistic, "wilcoxon_p": p,
        })
    tests = pd.DataFrame(tests)
    tests["fdr_bh"] = bh_adjust(tests["wilcoxon_p"])

    partials = pd.DataFrame(partial_rows)
    patient_partials = partials.groupby(["patient_id", "analysis", "outcome", "controls"], as_index=False)["partial_spearman_rho"].mean()
    partial_tests = []
    for analysis, group in patient_partials.groupby("analysis"):
        statistic, p = signed_rank_against_zero(group["partial_spearman_rho"])
        partial_tests.append({
            "analysis": analysis,
            "n_patients": group["patient_id"].nunique(),
            "median_patient_partial_rho": group["partial_spearman_rho"].median(),
            "iqr_patient_partial_rho": group["partial_spearman_rho"].quantile(.75) - group["partial_spearman_rho"].quantile(.25),
            "wilcoxon_statistic": statistic,
            "wilcoxon_p": p,
        })
    partial_tests = pd.DataFrame(partial_tests)
    partial_tests["fdr_bh"] = bh_adjust(partial_tests["wilcoxon_p"])

    sections.to_csv(result_dir / "section_qc_and_IFI44L_spatial_autocorrelation.tsv", sep="\t", index=False)
    correlations.to_csv(result_dir / "section_IFI44L_signature_correlations.tsv", sep="\t", index=False)
    patient_corr.to_csv(result_dir / "patient_IFI44L_signature_correlations.tsv", sep="\t", index=False)
    tests.to_csv(result_dir / "patient_level_signature_tests.tsv", sep="\t", index=False)
    partials.to_csv(result_dir / "section_partial_spearman.tsv", sep="\t", index=False)
    patient_partials.to_csv(result_dir / "patient_partial_spearman.tsv", sep="\t", index=False)
    partial_tests.to_csv(result_dir / "patient_level_partial_spearman_tests.tsv", sep="\t", index=False)

    order = ["IFN_response_excluding_IFI44L", "Myeloid", "T_cell", "Endothelial", "Cytotoxicity", "Checkpoint", "Epithelial"]
    plot = patient_corr.copy()
    plot["signature"] = pd.Categorical(plot["signature"], categories=order, ordered=True)
    plot = plot.sort_values("signature")
    fig, ax = plt.subplots(figsize=(8.8, 4.9))
    rng = np.random.default_rng(20260902)
    for i, signature in enumerate(order):
        vals = plot.loc[plot["signature"].eq(signature), "spearman_rho"].to_numpy()
        ax.scatter(i + rng.uniform(-.12, .12, len(vals)), vals, s=27, alpha=.66, color="#315A7D")
        ax.plot([i - .24, i + .24], [np.median(vals)] * 2, color="#B23A48", lw=2.4)
    ax.axhline(0, color="black", lw=.8, ls="--")
    ax.set_xticks(range(len(order)), [x.replace("_", "\n") for x in order], rotation=0)
    ax.set_ylabel("Within-tumor Spearman rho with IFI44L\n(mean of two sections)")
    ax.set_title("Patient-level spatial colocalization of IFI44L in primary TNBC (n=14)")
    fig.tight_layout()
    fig.savefig(figure_dir / "GSE210616_IFI44L_patient_signature_correlations.pdf")
    fig.savefig(figure_dir / "GSE210616_IFI44L_patient_signature_correlations.png", dpi=400)
    plt.close(fig)

    # Display all patients without outcome-driven selection: section 1 for each tumor.
    section1 = meta.loc[meta["section"].eq(1)].sort_values("patient_id")
    fig, axes = plt.subplots(4, 4, figsize=(12.5, 11.2))
    axes = axes.ravel()
    vmax = np.quantile(np.concatenate([section_frames[g]["IFI44L"].to_numpy() for g in section1["gsm"]]), .99)
    for ax, (_, row) in zip(axes, section1.iterrows()):
        frame = section_frames[row["gsm"]]
        sca = ax.scatter(frame["pxl_col_fullres"], frame["pxl_row_fullres"], c=frame["IFI44L"], s=4, cmap="magma", vmin=0, vmax=vmax)
        ax.invert_yaxis()
        ax.set_title(str(row["patient_id"]), fontsize=9)
        ax.axis("off")
    for ax in axes[len(section1):]:
        ax.axis("off")
    fig.suptitle("IFI44L spatial expression across 14 primary TNBC tumors (section 1)")
    cbar = fig.colorbar(sca, ax=axes.tolist(), fraction=.018, pad=.01)
    cbar.set_label("IFI44L log(CP10K + 1)")
    fig.savefig(figure_dir / "GSE210616_IFI44L_all14_section1_spatial.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / "GSE210616_IFI44L_all14_section1_spatial.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Prespecified representative: section closest to cohort-median detection fraction.
    median_detection = sections["IFI44L_detection_fraction"].median()
    representative = sections.iloc[(sections["IFI44L_detection_fraction"] - median_detection).abs().argsort().iloc[0]]
    gsm = representative["gsm"]
    frame = section_frames[gsm]
    sample_dir = source / "files" / gsm
    image_path = find_one(sample_dir, "*_tissue_hires_image.png.gz")
    scale_path = find_one(sample_dir, "*_scalefactors_json.json.gz")
    with gzip.open(image_path, "rb") as handle:
        image = Image.open(handle).convert("RGB")
    with gzip.open(scale_path, "rt") as handle:
        scale = json.load(handle)["tissue_hires_scalef"]
    x = frame["pxl_col_fullres"] * scale
    y = frame["pxl_row_fullres"] * scale
    panels = ["IFI44L", "Myeloid", "Endothelial", "Epithelial"]
    fig, axes = plt.subplots(1, 4, figsize=(15.2, 4.1))
    for ax, feature in zip(axes, panels):
        ax.imshow(image)
        values = frame[feature]
        scatter = ax.scatter(x, y, c=values, s=5.5, cmap="magma", alpha=.78, linewidths=0)
        fig.colorbar(scatter, ax=ax, fraction=.046, pad=.02)
        ax.set_title(feature.replace("_", " "))
        ax.axis("off")
    fig.suptitle(f"Prespecified median-detection section: {representative['patient_id']} / {gsm}")
    fig.tight_layout()
    fig.savefig(figure_dir / "GSE210616_representative_HnE_IFI44L_modules.pdf")
    fig.savefig(figure_dir / "GSE210616_representative_HnE_IFI44L_modules.png", dpi=350)
    plt.close(fig)

    audit = {
        "dataset": "GSE210616",
        "primary_design": "28 Visium sections from 14 primary TNBC tumors",
        "n_sections": int(len(sections)),
        "n_patients": int(sections["patient_id"].nunique()),
        "qc": "in_tissue=1, total UMI>=500, detected genes>=200, mitochondrial percent<20",
        "inference_unit": "patient; two sections averaged before signed-rank tests",
        "spot_level_p_values": "descriptive only; not used as independent patient evidence",
        "moran_permutations": args.permutations,
        "representative_selection": "section closest to cohort median IFI44L detection fraction",
        "python": platform.python_version(),
        "scanpy": sc.__version__,
    }
    (result_dir / "analysis_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (log_dir / "spatial_GSE210616_python_versions.txt").write_text(
        f"python={platform.python_version()}\nscanpy={sc.__version__}\n", encoding="utf-8"
    )
    print(sections[["gsm", "patient_id", "n_spots_after_qc", "IFI44L_detection_fraction", "IFI44L_moran_I", "moran_fdr_bh"]].to_string(index=False))
    print("\nPatient-level signature tests")
    print(tests.sort_values("fdr_bh").to_string(index=False))
    print("\nPatient-level partial Spearman tests")
    print(partial_tests.sort_values("fdr_bh").to_string(index=False))


if __name__ == "__main__":
    main()
