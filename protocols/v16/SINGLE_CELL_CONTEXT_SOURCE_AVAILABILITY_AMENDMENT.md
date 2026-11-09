# Pre-fit amendment: IFN-response module source availability

Date: 2026-11-09. Parent authorization: report nominal definition and available genes exactly; retain existing score, thresholds and nine-comparison family.

The initial script stopped before aggregation/correlation fitting because `DDX58` was absent from `results/single_cell/GSE176078_TNBC_targeted_IFI44L_immune.h5ad`. The original `04_single_cell_tnbc_analysis.py` defines 22 IFN-response genes but calculates its per-cell score from the available subset. The H5AD retains 21 of these genes. No value for DDX58 will be imputed or assigned zero.

Revision 2 reconstructs the arithmetic mean using those 21 available genes only, validates agreement with the unchanged existing obs score within 1e-6, reports all 22 genes and their availability, and labels figures and captions “nominal 22-gene definition; 21 available”. All nine author categories, cell/donor thresholds, asymptotic Spearman inference and fixed nine-test BH remain unchanged. This amendment was made before any correlation was fitted or examined.

The original plan and receipt are preserved as `SINGLE_CELL_CONTEXT_PLAN_v1_before_availability_check.md` and `SINGLE_CELL_CONTEXT_PLAN_LOCK_v1.json`. The stopped run is preserved as `SINGLE_CELL_CONTEXT_pre_fit_failure_v1.log`. New primary plan SHA256 is `f23a3904dad46db7690136ebf99d80ce91cc41d3b7400aa0f89cae3b83458cc9`.
