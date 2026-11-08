#!/usr/bin/env python3
"""Audit inputs, software and expected outputs for the non-experimental workflow."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess


PYTHON_PACKAGES = {
    "anndata": "anndata",
    "cellxgene-census": "cellxgene-census",
    "h5py": "h5py",
    "matplotlib": "matplotlib",
    "numpy": "numpy",
    "openpyxl": "openpyxl",
    "pandas": "pandas",
    "Pillow": "pillow",
    "PyMuPDF": "pymupdf",
    "scanpy": "scanpy",
    "scipy": "scipy",
    "seaborn": "seaborn",
    "statsmodels": "statsmodels",
    "zarr": "zarr",
}

R_PACKAGES = [
    "data.table", "limma", "GSVA", "ggplot2", "clusterProfiler", "org.Hs.eg.db",
    "GEOquery", "Biobase", "AnnotationDbi", "hgu133plus2.db", "survival", "survminer", "readxl",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, required: bool, deep_hash: bool) -> dict[str, object]:
    exists = path.is_file()
    record: dict[str, object] = {
        "path": str(path), "required": required, "exists": exists,
        "size_bytes": path.stat().st_size if exists else None,
    }
    if exists and deep_hash:
        record["sha256"] = sha256(path)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--deep-hash", action="store_true")
    parser.add_argument(
        "--require-outputs",
        action="store_true",
        help="Fail unless every declared final output is present; use after a complete pipeline run.",
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    raw_root = (args.raw_root or project.parent).resolve()

    inputs = [
        (raw_root / "数据/差异富集分析/3 差异分析/symbol.BRCA.txt", True),
        (raw_root / "数据/差异富集分析/3 差异分析/三阴.txt", True),
        (raw_root / "数据/差异富集分析/1 乳腺癌分型/3_乳腺癌分型.xlsx", True),
        (project / "source_metadata/GSE58812/GSE58812_series_matrix.txt.gz", True),
        (project / "source_metadata/Wu2021_NatGenet_Supplementary_Tables.xlsx", True),
        (project / "source_metadata/GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad", True),
        (project / "source_metadata/spatial/GSE210616/GSE210616_primary28_samples.tsv", True),
        (project / "source_metadata/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx", True),
    ]
    input_records = [file_record(path, required, args.deep_hash) for path, required in inputs]

    spatial_files = project / "source_metadata/spatial/GSE210616/files"
    spatial_counts = {
        "sample_directories": len([p for p in spatial_files.glob("GSM*") if p.is_dir()]),
        "matrix_h5": len(list(spatial_files.glob("GSM*/*_filtered_feature_bc_matrix.h5"))),
        "positions": len(list(spatial_files.glob("GSM*/*_tissue_positions_list.csv.gz"))),
        "images": len(list(spatial_files.glob("GSM*/*_tissue_hires_image.png.gz"))),
        "scalefactors": len(list(spatial_files.glob("GSM*/*_scalefactors_json.json.gz"))),
    }

    python_versions = {"python": platform.python_version()}
    for label, distribution in PYTHON_PACKAGES.items():
        try:
            python_versions[label] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            python_versions[label] = None

    rscript = shutil.which("Rscript")
    r_versions: dict[str, str | None] = {"Rscript": rscript}
    if rscript:
        expr = (
            "cat('R=',R.version$major,'.',R.version$minor,'\\n',sep='');"
            "p<-c(" + ",".join(repr(x) for x in R_PACKAGES) + ");"
            "for(x in p){v<-tryCatch(as.character(packageVersion(x)),error=function(e) NA);"
            "cat(x,'=',v,'\\n',sep='')}"
        )
        completed = subprocess.run([rscript, "-e", expr], text=True, capture_output=True)
        if completed.returncode == 0:
            for line in completed.stdout.splitlines():
                if "=" in line:
                    key, value = line.split("=", 1)
                    r_versions[key] = None if value == "NA" else value
        else:
            r_versions["error"] = completed.stderr.strip()

    expected_outputs = [
        project / "results/bulk/TCGA_TNBC_definition_audit.tsv",
        project / "results/bulk/TCGA_TNBC_sample_definition_manifest.tsv",
        project / "results/bulk/continuous_IFI44L_gene_associations.tsv",
        project / "results/bulk/TCGA_TNBC_continuous_IFI44L_Cox.tsv",
        project / "results/bulk/TCGA_TNBC_median_split_logrank_sensitivity.tsv",
        project / "results/bulk/TCGA_TNBC_survival_audit.tsv",
        project / "results/validation_GSE58812/GSE58812_continuous_IFI44L_Cox.tsv",
        project / "results/single_cell/GSE176078_TNBC_IFI44L_celltype_major_summary.tsv",
        project / "results/spatial/GSE210616/patient_level_signature_tests.tsv",
        project / "figures/v11_public_multiomics/Figure1.pdf",
        project / "figures/v11_public_multiomics/Figure2_cross_cohort_replication_survival.pdf",
        project / "figures/v11_public_multiomics/Figure3.pdf",
        project / "figures/v11_public_multiomics/Figure4.pdf",
        project / "figures/v11_public_multiomics/Figure5_public_multiomics_conceptual_synthesis.pdf",
    ]
    output_records = [file_record(path, True, args.deep_hash) for path in expected_outputs]

    missing_inputs = [r["path"] for r in input_records if r["required"] and not r["exists"]]
    missing_python = [k for k, v in python_versions.items() if k != "python" and v is None]
    missing_r = [k for k, v in r_versions.items() if k not in {"Rscript", "R", "error"} and v is None]
    spatial_ready = all(value == 28 for value in spatial_counts.values())
    report = {
        "project": str(project),
        "raw_root": str(raw_root),
        "inputs": input_records,
        "spatial_expected_28_each": spatial_counts,
        "python_versions": python_versions,
        "r_versions": r_versions,
        "expected_outputs": output_records,
        "status": {
            "inputs_ready": not missing_inputs and spatial_ready,
            "python_ready": not missing_python,
            "r_ready": bool(rscript) and not missing_r,
            "outputs_present": all(r["exists"] for r in output_records),
            "missing_inputs": missing_inputs,
            "missing_python_packages": missing_python,
            "missing_r_packages": missing_r,
        },
    }
    output = args.json_out or project / "reproducibility/preflight_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report["status"], ensure_ascii=False, indent=2))
    ready = (
        report["status"]["inputs_ready"]
        and report["status"]["python_ready"]
        and report["status"]["r_ready"]
    )
    if args.require_outputs:
        ready = ready and report["status"]["outputs_present"]
    if not ready:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
