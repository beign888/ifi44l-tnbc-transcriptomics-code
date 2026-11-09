#!/usr/bin/env python3
"""Exploratory donor-level IFI44L/IFN-response associations in fixed author lineages."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
from importlib.metadata import version
import itertools
import math

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import spearmanr, rankdata

MAJORS = ["B-cells", "CAFs", "Cancer Epithelial", "Endothelial", "Myeloid",
          "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"]
IFN_GENES = ["IFI44", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2",
             "OAS3", "OASL", "RSAD2", "IFI6", "IFIH1", "DDX58", "DDX60", "HERC5",
             "HERC6", "CMPK2", "XAF1", "SAMD9L", "STAT1", "IRF7"]
SCORE = "interferon_response_excluding_IFI44L"
LOCKED_PLAN_SHA256 = "f23a3904dad46db7690136ebf99d80ce91cc41d3b7400aa0f89cae3b83458cc9"
EXACT_AMENDMENT_SHA256 = "f5ac9f4db3dd75a3a1a26d8d33a421b694cd0ebad9f1a082a29fa6832551d664"
SOURCE_AMENDMENT_SHA256 = "8455a75c8dd79061eccc69be7111d5f0a08e0f42f33a6782fa5cdf4cac28a54d"
ORIGINAL_PLAN_SHA256 = "78e3af138cedbaa9b9e88ce3a3b3026491ad717f1e6f7a64a09c69112e37c1ab"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dense_column(matrix, column: int) -> np.ndarray:
    values = matrix[:, column]
    return np.asarray(values.toarray() if sparse.issparse(values) else values).ravel()


def bh_fixed_nine(p_values: np.ndarray) -> np.ndarray:
    """Missing tests occupy the nine-test family as P=1, then display as missing."""
    if len(p_values) != 9:
        raise ValueError("The multiplicity family must contain exactly nine categories")
    valid = np.isfinite(p_values)
    filled = np.where(valid, p_values, 1.0)
    order = np.argsort(filled, kind="stable")
    adjusted = np.minimum.accumulate((filled[order] * 9 / np.arange(1, 10))[::-1])[::-1]
    out = np.empty(9)
    out[order] = np.minimum(adjusted, 1.0)
    out[~valid] = np.nan
    return out


def format_p(value: float) -> str:
    return "NA" if not np.isfinite(value) else (f"{value:.2g}" if value < .001 else f"{value:.3f}")


def exact_permutation(x: np.ndarray, y: np.ndarray) -> dict:
    """Enumerate donor-label permutations, retaining multiplicity in tied ranks."""
    n = len(x)
    if n > 10:
        raise RuntimeError("Exact enumeration was planned only for n <= 10")
    rank_x = np.rint(2 * rankdata(x, method="average") - (n + 1)).astype(np.int64)
    rank_y = np.rint(2 * rankdata(y, method="average") - (n + 1)).astype(np.int64)
    observed_dot = int(rank_x @ rank_y)
    extreme, total = 0, 0
    permutations = itertools.permutations(rank_y.tolist())
    while True:
        block = list(itertools.islice(permutations, 50000))
        if not block:
            break
        values = np.asarray(block, dtype=np.int64)
        dots = values @ rank_x
        extreme += int(np.count_nonzero(np.abs(dots) >= abs(observed_dot)))
        total += len(block)
    if total != math.factorial(n) or extreme < 1:
        raise RuntimeError("Invalid exhaustive permutation count")
    return {"p_value_exact_donor_permutation": extreme / total, "n_permutations_exact": total,
            "n_extreme_two_sided": extreme, "observed_centered_rank_dot": observed_dot,
            "n_unique_IFN_values": len(np.unique(x)), "n_unique_IFI44L_values": len(np.unique(y)),
            "IFN_doubled_centered_average_ranks": json.dumps(rank_x.tolist()),
            "IFI44L_doubled_centered_average_ranks": json.dumps(rank_y.tolist())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--plan-dir", type=Path, default=None, help="Directory containing all locked plans and amendments; default: PROJECT/audit_v16")
    args = parser.parse_args()
    project = args.project.resolve()
    plan_dir = args.plan_dir.resolve() if args.plan_dir is not None else project / "audit_v16"
    out = project / "results/single_cell_context_v16"
    fig_dir = project / "figures/v16_supplementary"
    out.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    plan = plan_dir / "SINGLE_CELL_CONTEXT_PLAN.md"
    if sha256(plan) != LOCKED_PLAN_SHA256:
        raise RuntimeError("The pre-fit analysis plan differs from its locked SHA256")
    original_plan = plan_dir / "SINGLE_CELL_CONTEXT_PLAN_v1_before_availability_check.md"
    source_amendment = plan_dir / "SINGLE_CELL_CONTEXT_SOURCE_AVAILABILITY_AMENDMENT.md"
    if sha256(original_plan) != ORIGINAL_PLAN_SHA256 or sha256(source_amendment) != SOURCE_AMENDMENT_SHA256:
        raise RuntimeError("The original plan or source-availability amendment differs from its locked SHA256")
    exact_amendment = plan_dir / "SINGLE_CELL_CONTEXT_EXACT_PERMUTATION_AMENDMENT.md"
    if sha256(exact_amendment) != EXACT_AMENDMENT_SHA256:
        raise RuntimeError("The diagnostic-driven exact-permutation amendment differs from its locked SHA256")
    input_path = project / "results/single_cell/GSE176078_TNBC_targeted_IFI44L_immune.h5ad"
    baseline_path = project / "results/single_cell/GSE176078_TNBC_IFI44L_by_donor_celltype_major.tsv"
    data = ad.read_h5ad(input_path)
    if data.shape != (42512, 63):
        raise RuntimeError(f"Unexpected targeted H5AD dimensions {data.shape}")
    required_obs = {"donor_id", "celltype_major", SCORE}
    if not required_obs.issubset(data.obs) or "counts" not in data.layers:
        raise RuntimeError("Missing required author annotations, existing module or counts layer")
    unknown = set(data.obs["celltype_major"].astype(str)) - set(MAJORS)
    if unknown or data.obs[list(required_obs)].isna().any().any():
        raise RuntimeError(f"Missing annotations or unexpected major categories: {unknown}")
    missing = set(["IFI44L"] + IFN_GENES) - set(data.var_names)
    if "IFI44L" not in data.var_names or missing != {"DDX58"}:
        raise RuntimeError(f"Gene availability differs from the pre-fit amendment: {missing}")
    x = data.X.tocsr() if sparse.issparse(data.X) else np.asarray(data.X)
    x_values = x.data if sparse.issparse(x) else x
    if not np.isfinite(x_values).all() or (x_values < 0).any():
        raise RuntimeError("Invalid normalized expression values")
    genes = data.var_names.tolist()
    ifi = dense_column(x, genes.index("IFI44L")).astype(float)
    detected = dense_column(data.layers["counts"], genes.index("IFI44L")) > 0
    module_x = x[:, [genes.index(g) for g in IFN_GENES if g in genes]].astype(np.float64)
    reconstructed = np.asarray(module_x.mean(axis=1)).ravel()
    existing_score = data.obs[SCORE].to_numpy(dtype=float)
    score_max_error = float(np.max(np.abs(reconstructed - existing_score)))
    if not np.isfinite(existing_score).all() or score_max_error > 1e-6:
        raise RuntimeError(f"Existing 21-of-22-gene module disagreement: {score_max_error}")
    frame = data.obs[["donor_id", "celltype_major"]].astype(str).copy()
    frame["IFI44L_mean_logCP10k"] = ifi
    frame["IFI44L_detection_fraction"] = detected.astype(float)
    frame["IFN_response_mean_logCP10k"] = existing_score
    observed = frame.groupby(["donor_id", "celltype_major"], observed=True).agg(
        n_cells=("IFI44L_mean_logCP10k", "size"),
        IFI44L_mean_logCP10k=("IFI44L_mean_logCP10k", "mean"),
        IFI44L_detection_fraction=("IFI44L_detection_fraction", "mean"),
        IFN_response_mean_logCP10k=("IFN_response_mean_logCP10k", "mean"),
    ).reset_index()
    keys = ["donor_id", "celltype_major"]
    baseline = pd.read_csv(baseline_path, sep="\t")
    if baseline.duplicated(keys).any() or observed.duplicated(keys).any():
        raise RuntimeError("Duplicated donor-major baseline or new aggregate")
    checked = observed.merge(baseline, on=keys, how="outer", suffixes=("_new", "_baseline"), indicator=True)
    if not checked["_merge"].eq("both").all():
        raise RuntimeError("Baseline and new donor-major aggregate keys differ")
    baseline_errors = {}
    for field in ["n_cells", "IFI44L_mean_logCP10k", "IFI44L_detection_fraction"]:
        error = float((checked[f"{field}_new"] - checked[f"{field}_baseline"]).abs().max())
        baseline_errors[field] = error
        if error > (0 if field == "n_cells" else 1e-6):
            raise RuntimeError(f"Baseline disagreement in {field}: {error}")
    donors = sorted(frame["donor_id"].unique())
    if len(donors) != 10 or int(observed["n_cells"].sum()) != 42512:
        raise RuntimeError("Incorrect donor or retained-cell count")
    grid = pd.MultiIndex.from_product([donors, MAJORS], names=keys).to_frame(index=False)
    all_rows = grid.merge(observed, on=keys, how="left", validate="one_to_one")
    all_rows["n_cells"] = all_rows["n_cells"].fillna(0).astype(int)
    all_rows["eligible"] = (all_rows["n_cells"] >= 20) & np.isfinite(all_rows["IFI44L_mean_logCP10k"]) & np.isfinite(all_rows["IFN_response_mean_logCP10k"])
    all_rows["exclusion_reason"] = np.select(
        [all_rows["eligible"], all_rows["n_cells"].eq(0), all_rows["n_cells"].lt(20)],
        ["", "absent donor-major aggregate", "fewer than 20 cells"], default="nonfinite expression")
    eligible = all_rows.loc[all_rows["eligible"]].copy()
    tests = []
    for major in MAJORS:
        d = eligible.loc[eligible["celltype_major"].eq(major)]
        n = len(d)
        status, rho, p_value = "tested", np.nan, np.nan
        if n < 5:
            status = "insufficient_donors"
        elif d["IFI44L_mean_logCP10k"].nunique() < 2 or d["IFN_response_mean_logCP10k"].nunique() < 2:
            status = "constant_input"
        else:
            rho, p_value = spearmanr(d["IFN_response_mean_logCP10k"], d["IFI44L_mean_logCP10k"], alternative="two-sided")
        tests.append({"celltype_major": major, "n_eligible_donors": n, "n_cells_in_eligible_aggregates": int(d["n_cells"].sum()),
                      "spearman_rho": float(rho), "p_value_asymptotic_two_sided": float(p_value), "status": status,
                      "multiplicity_family_size": 9})
    tests = pd.DataFrame(tests)
    tests["fdr_bh_full_nine"] = bh_fixed_nine(tests["p_value_asymptotic_two_sided"].to_numpy())
    permutation_records = []
    for major in MAJORS:
        d = eligible.loc[eligible["celltype_major"].eq(major)].sort_values("donor_id")
        status = tests.loc[tests["celltype_major"].eq(major), "status"].iloc[0]
        record = {"celltype_major": major, "donor_order": json.dumps(d["donor_id"].tolist()),
                  "permutation_status": "not_tested", "p_value_exact_donor_permutation": np.nan,
                  "n_permutations_exact": 0, "n_extreme_two_sided": 0,
                  "observed_centered_rank_dot": np.nan, "n_unique_IFN_values": d["IFN_response_mean_logCP10k"].nunique(),
                  "n_unique_IFI44L_values": d["IFI44L_mean_logCP10k"].nunique(),
                  "IFN_doubled_centered_average_ranks": "", "IFI44L_doubled_centered_average_ranks": ""}
        if status == "tested":
            record.update(exact_permutation(d["IFN_response_mean_logCP10k"].to_numpy(), d["IFI44L_mean_logCP10k"].to_numpy()))
            record["permutation_status"] = "exhaustive_donor_label_sensitivity"
        permutation_records.append(record)
    permutation_table = pd.DataFrame(permutation_records)
    permutation_table["fdr_bh_exact_full_nine"] = bh_fixed_nine(permutation_table["p_value_exact_donor_permutation"].to_numpy())
    tests = tests.merge(permutation_table.drop(columns=["donor_order", "IFN_doubled_centered_average_ranks", "IFI44L_doubled_centered_average_ranks"]),
                        on="celltype_major", how="left", validate="one_to_one", sort=False)
    def save_table(df: pd.DataFrame, name: str) -> None:
        df.to_csv(out / name, sep="\t", index=False, na_rep="NA", float_format="%.17g")
    save_table(observed, "observed_donor_major_aggregates.tsv")
    save_table(all_rows, "all_donor_major_aggregates.tsv")
    save_table(eligible, "eligible_donor_major_aggregates.tsv")
    save_table(all_rows.loc[~all_rows["eligible"]], "excluded_donor_major_aggregates.tsv")
    save_table(eligible, "Supplementary_Figure_S3_sample_points.tsv")
    save_table(tests, "within_lineage_IFI44L_IFN_response_spearman.tsv")
    save_table(permutation_table, "exact_donor_label_permutation_sensitivity.tsv")
    save_table(pd.DataFrame({"gene_order": np.arange(1, 23), "gene": IFN_GENES,
                            "present_in_normalized_X": [g in genes for g in IFN_GENES],
                            "module": "broad IFN-response; IFI44L excluded; nominal22, available21; unweighted available-gene mean"}), "IFN_response_22_module_genes.tsv")
    save_table(checked.drop(columns="_merge"), "baseline_donor_major_comparison.tsv")
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5, "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3, "axes.linewidth": .75, "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "ifi44l-v16-S3",
        "legend.frameon": False, "savefig.facecolor": "white"})
    palette = ["#2B6CA3", "#23877D", "#B64E55", "#D28A35", "#756BB1", "#A6761D", "#1B9E77", "#E377C2", "#6B6B6B", "#17BECF"]
    donor_colors = dict(zip(donors, palette))
    fig, axes = plt.subplots(3, 3, figsize=(7.3, 7.6))
    fig.subplots_adjust(left=.10, right=.98, top=.875, bottom=.19, hspace=.78, wspace=.48)
    for idx, (major, ax) in enumerate(zip(MAJORS, axes.ravel())):
        d = eligible.loc[eligible["celltype_major"].eq(major)]
        test = tests.loc[tests["celltype_major"].eq(major)].iloc[0]
        for donor in donors:
            point = d.loc[d["donor_id"].eq(donor)]
            if len(point):
                ax.scatter(point["IFN_response_mean_logCP10k"], point["IFI44L_mean_logCP10k"], s=23, color=donor_colors[donor], alpha=.9, edgecolors="white", linewidths=.35)
        ax.set_title(major, loc="left", pad=27)
        stat = f"n = {int(test['n_eligible_donors'])} donors"
        if test["status"] == "tested":
            stat += f"; rho = {test['spearman_rho']:.2f}\nExact permutation FDR = {format_p(test['fdr_bh_exact_full_nine'])}"
        else:
            stat += "\nNot tested (<5 donors)" if test["status"] == "insufficient_donors" else "\nNot tested (constant input)"
        ax.text(0, 1.015, stat, transform=ax.transAxes, va="bottom", fontsize=6.4, linespacing=1.3)
        ax.text(-.24, 1.24, chr(97 + idx), transform=ax.transAxes, fontweight="bold", fontsize=9)
        ax.set_xlabel("IFN-response score (21/22 genes)")
        ax.set_ylabel("IFI44L mean ln(CP10K + 1)")
        ax.tick_params(width=.65, length=2.5)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=donor_colors[d], markeredgewidth=0, markersize=4, label=d) for d in donors]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.54, .07), ncol=5, fontsize=6.1, columnspacing=1.25, handletextpad=.25)
    fig.suptitle("IFI44L and IFN-response expression within author-defined cell lineages", y=.992, fontsize=9)
    fig.text(.50, .028, "One point per eligible donor (≥20 cells per lineage); exact donor-label permutation sensitivity; BH family: nine lineages.\nIFN response: nominal 22 genes, 21 available; IFI44L excluded. Distinct from the bulk IFN-alpha score.", ha="center", va="center", fontsize=6.4, linespacing=1.5)
    stem = fig_dir / "Supplementary_Figure_S3_within_lineage_IFN_context"
    fig.savefig(stem.with_suffix(".pdf"), metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(stem.with_suffix(".svg"), metadata={"Date": None})
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    plt.close(fig)
    caption = ("Supplementary Figure S3. Donor-level IFI44L and IFN-response expression within author-defined major cell lineages. "
        "Each point denotes one TNBC donor with at least 20 retained cells in that lineage; donors have consistent colours across panels. "
        "Axes show the mean whole-library-normalized ln(CP10K + 1) IFI44L expression and the mean pre-existing IFN-response score (nominal 22-gene definition; 21 available genes, with DDX58 absent), "
        "which excludes IFI44L. This broad single-cell module differs from the 19-gene bulk IFN-alpha module. "
        "Two-sided asymptotic Spearman tests were initially performed only with at least five eligible donors and variable inputs. "
        "After diagnosing an extreme asymptotic P value for a perfect ordering among five donors, a uniform exact donor-label permutation sensitivity "
        "was added for every testable lineage, enumerating all n! permutations, including tied-label multiplicities. Panels display its two-sided "
        "BH-adjusted P values; both original asymptotic and exact results are preserved in source tables. This sensitivity was diagnostic-driven, not preregistered. "
        "Each BH adjustment retains all nine prespecified lineage comparisons in its denominator; untestable comparisons remain NA. "
        "Panels with fewer than five donors are descriptive. No cells were treated as independent biological replicates. "
        "These are exploratory within-cohort associations and do not identify an absolute contribution to bulk RNA, a signalling source or causal effects.\n")
    (out / "Supplementary_Figure_S3_caption.txt").write_text(caption, encoding="utf-8")
    artifact_paths = sorted(out.glob("*.tsv")) + [out / "Supplementary_Figure_S3_caption.txt"] + [stem.with_suffix(ext) for ext in [".pdf", ".svg", ".png"]]
    audit = {"status": "PASS", "analysis": "exploratory donor-level within-lineage IFI44L/IFN-response association",
        "plan_sha256": sha256(plan), "script_sha256": sha256(Path(__file__)),
        "original_pre_availability_plan_sha256": sha256(original_plan),
        "inputs": {str(p.relative_to(project)): sha256(p) for p in [input_path, baseline_path]},
        "input_dimensions": list(data.shape), "n_unique_donors": len(donors), "prespecified_major_classes": MAJORS,
        "minimum_cells_per_donor_major": 20, "minimum_donors_per_test": 5, "inference_unit": "donor within a fixed lineage",
        "multiplicity": "BH, fixed nine-test family; unavailable tests use internal P=1 and displayed NA",
        "broad_IFN_response_module_genes": IFN_GENES, "module_missing_genes": sorted(missing),
        "module_defined_gene_count": 22, "module_available_gene_count": 21,
        "source_availability_amendment_sha256": sha256(source_amendment),
        "exact_permutation_amendment_sha256": sha256(exact_amendment),
        "exact_permutation_method": "diagnostic-driven sensitivity after asymptotic inspection; exhaustive n! donor labels; abs centered doubled-average-rank dot >= observed absolute dot; P=extreme/total; no Monte Carlo +1; tied labels retain multiplicity",
        "exact_permutation_block_size": 50000, "exact_permutation_seed": None,
        "normalization": "existing X: ln(CP10K + 1), full source raw library before targeted slicing; no re-normalization",
        "module_reconstruction_max_abs_difference": score_max_error, "baseline_max_abs_differences": baseline_errors,
        "observed_aggregates": len(observed), "complete_grid_aggregates": len(all_rows), "eligible_aggregates": len(eligible),
        "excluded_or_absent_aggregates": int((~all_rows["eligible"]).sum()),
        "tested_lineages": int(tests["status"].eq("tested").sum()),
        "results": json.loads(tests.to_json(orient="records")),
        "versions": {"python": platform.python_version(), "anndata": version("anndata"), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": mpl.__version__},
        "outputs_sha256": {str(p.relative_to(project)): sha256(p) for p in artifact_paths},
        "limits": ["same existing ten donors, not an independent validation", "asymptotic P values may be imprecise for small samples",
                   "association only; no absolute bulk source contributions or causal interpretation"]}
    (out / "SINGLE_CELL_CONTEXT_AUDIT.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(tests.to_string(index=False))
    print(f"PASS: {len(eligible)}/90 eligible donor-major aggregates; module error={score_max_error:.3g}; baseline errors={baseline_errors}")
    print(f"Audit: {out / 'SINGLE_CELL_CONTEXT_AUDIT.json'}")


if __name__ == "__main__":
    main()
