#!/usr/bin/env python3
"""Run the complete non-experimental IFI44L analysis and figure pipeline."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import os
import shutil
import subprocess
import sys


def run_step(name: str, command: list[str], cwd: Path, log_dir: Path, env: dict[str, str]) -> None:
    log_path = log_dir / f"reproduce_{name}.log"
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {name}")
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        raise RuntimeError(f"Step failed: {name}. Inspect {log_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--stage", choices=["preflight", "figures", "public-analyses", "all"], default="preflight")
    parser.add_argument("--python", type=Path)
    parser.add_argument("--permutations", type=int, default=499)
    parser.add_argument(
        "--run-processing-sensitivity", action="store_true",
        help="Also reproduce the documented HiSeqV2 processing check from Additional file 5.",
    )
    parser.add_argument(
        "--run-census-audit",
        action="store_true",
        help="Also query the frozen CELLxGENE Census version used for the metadata audit.",
    )
    args = parser.parse_args()
    if args.run_processing_sensitivity and args.stage == "figures":
        parser.error("Processing sensitivity requires preflight, public-analyses or all stage")
    project = args.project.resolve()
    raw_root = (args.raw_root or project.parent).resolve()
    # Preserve a virtual-environment interpreter symlink. Path.resolve() would
    # replace it with the system interpreter and silently lose the venv site-packages.
    python = (args.python or Path(sys.executable)).absolute()
    scripts = project / "scripts"
    source = project / "source_metadata"
    results = project / "results"
    log_dir = project / "logs" / "reproduction"
    log_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(project / ".matplotlib-cache")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    sensitivity_work = source / "TCGA_documented_processing"
    if args.run_processing_sensitivity:
        provenance_path = sensitivity_work / "SELFCONTAINED_SUBSET_PROVENANCE.json"
        if not provenance_path.is_file():
            raise FileNotFoundError("Stage the TCGA_processing_sensitivity folder from Additional file 5 at " + str(sensitivity_work))
        provenance = json.loads(provenance_path.read_text())
        subset_path = sensitivity_work / provenance["subset_file"]
        if hashlib.sha256(subset_path.read_bytes()).hexdigest() != provenance["subset_sha256"]:
            raise RuntimeError("Documented-processing input differs from the supplied SHA-256")
        for name, expected in provenance["baseline_files"].items():
            if hashlib.sha256((sensitivity_work / name).read_bytes()).hexdigest() != expected["sha256"]:
                raise RuntimeError("Processing-sensitivity baseline input differs: " + name)

    original_matrix = raw_root / "数据/差异富集分析/3 差异分析/symbol.BRCA.txt"
    frozen_matrix = source / "TCGA_frozen/TCGA_BRCA_TNBC_frozen.tsv.gz"
    if not original_matrix.exists() and not frozen_matrix.exists():
        run_step("00_restore_frozen_tcga", [str(python), str(scripts / "00_restore_frozen_tcga_input.py"),
                 "--project", str(project)], raw_root, log_dir, env)

    preflight = [str(python), str(scripts / "17_reproducibility_preflight.py"),
                 "--project", str(project), "--raw-root", str(raw_root)]
    run_step("00_preflight", preflight, raw_root, log_dir, env)
    if args.stage == "preflight":
        return

    if args.stage in {"public-analyses", "all"}:
        rscript = shutil.which("Rscript")
        if not rscript:
            raise RuntimeError("Rscript is unavailable")
        analysis_steps = [
            ("01_tcga_definition_audit", [str(python), str(scripts / "21_audit_tcga_tnbc_definition.py"),
             "--project", str(project), "--raw-root", str(raw_root)], raw_root),
            ("02_tcga_bulk", [rscript, str(scripts / "01_bulk_tnbc_reanalysis.R"),
             "--project", str(project), "--raw-root", str(raw_root)], raw_root),
            ("02b_tcga_survival", [rscript, str(scripts / "06_tcga_tnbc_survival.R"),
             str(project)], raw_root),
            ("03_gse58812", [rscript, str(scripts / "05_gse58812_validation.R"),
             "--project", str(project)], raw_root),
            ("04_wu_supplement", [str(python), str(scripts / "03_parse_wu2021_supplement.py"),
             "--project", str(project)], raw_root),
            ("04_gse176078_single_cell", [str(python), str(scripts / "04_single_cell_tnbc_analysis.py"),
             "--project", str(project), "--h5ad", str(source / "GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad")], raw_root),
            ("05_gse210616_spatial", [str(python), str(scripts / "09_gse210616_spatial_analysis.py"),
             "--project", str(project), "--data-dir", str(source / "spatial/GSE210616"),
             "--permutations", str(args.permutations)], raw_root),
        ]
        if args.run_census_audit:
            analysis_steps.insert(
                3,
                ("03b_gse176078_census_audit", [str(python), str(scripts / "02_census_gse176078_audit.py"),
                 "--project", str(project)], raw_root),
            )
        for name, command, cwd in analysis_steps:
            run_step(name, command, cwd, log_dir, env)
        if args.run_processing_sensitivity:
            run_step("06_documented_processing_sensitivity", [
                rscript, str(scripts / "22_documented_tcga_processing_sensitivity.R"),
                "--work-dir", str(sensitivity_work),
                "--output-dir", str(results / "processing_sensitivity"),
            ], raw_root, log_dir, env)

    if args.stage in {"figures", "all", "public-analyses"}:
        for name, script in [
            ("11_figures_single_cell_spatial", "14_make_nature_style_figures.py"),
            ("12_figure_spatial_final", "16_make_spatial_main_figure.py"),
            ("13_figures_v11_public_multiomics", "15_make_bulk_validation_figures_v11.py"),
        ]:
            run_step(name, [str(python), str(scripts / script), "--project", str(project)], raw_root, log_dir, env)

    run_step("99_final_preflight", preflight + ["--require-outputs"], raw_root, log_dir, env)
    print(f"Reproduction complete. Logs: {log_dir}")


if __name__ == "__main__":
    main()
