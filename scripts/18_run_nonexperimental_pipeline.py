#!/usr/bin/env python3
"""Run the complete non-experimental IFI44L analysis and figure pipeline."""

from __future__ import annotations

import argparse
from datetime import datetime
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
        "--run-census-audit",
        action="store_true",
        help="Also query the frozen CELLxGENE Census version used for the metadata audit.",
    )
    args = parser.parse_args()
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
