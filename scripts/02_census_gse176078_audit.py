#!/usr/bin/env python3
"""Audit the versioned CELLxGENE representation of GSE176078."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import cellxgene_census
import pandas as pd


CENSUS_VERSION = "2025-11-08"
DATASET_ID = "9fddb063-056d-4202-8b8a-4b0ee531d3ce"
OBS_COLUMNS = [
    "soma_joinid",
    "dataset_id",
    "donor_id",
    "cell_type",
    "cell_type_ontology_term_id",
    "disease",
    "disease_ontology_term_id",
    "tissue",
    "tissue_general",
    "assay",
    "sex",
    "development_stage",
    "suspension_type",
    "is_primary_data",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    project = args.project.resolve()
    out = project / "results" / "single_cell"
    log = project / "logs"
    out.mkdir(parents=True, exist_ok=True)
    log.mkdir(parents=True, exist_ok=True)

    value_filter = f"dataset_id == '{DATASET_ID}' and is_primary_data == True"
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs = cellxgene_census.get_obs(
            census,
            "homo_sapiens",
            value_filter=value_filter,
            column_names=OBS_COLUMNS,
        )

    obs.to_csv(out / "GSE176078_census_cell_metadata.tsv.gz", sep="\t", index=False)
    summary_frames = []
    for column in ["donor_id", "cell_type", "disease", "tissue", "assay", "suspension_type"]:
        counts = (
            obs[column]
            .astype("string")
            .fillna("<missing>")
            .value_counts(dropna=False)
            .rename_axis("value")
            .reset_index(name="n_cells")
        )
        counts.insert(0, "field", column)
        summary_frames.append(counts)
    pd.concat(summary_frames, ignore_index=True).to_csv(
        out / "GSE176078_census_metadata_summary.tsv", sep="\t", index=False
    )

    donor_celltype = (
        obs.groupby(["donor_id", "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    donor_celltype.to_csv(out / "GSE176078_cells_by_donor_and_type.tsv", sep="\t", index=False)

    audit = {
        "census_version": CENSUS_VERSION,
        "dataset_id": DATASET_ID,
        "n_cells": int(len(obs)),
        "n_donors": int(obs["donor_id"].nunique(dropna=True)),
        "n_cell_types": int(obs["cell_type"].nunique(dropna=True)),
        "diseases": sorted(obs["disease"].dropna().astype(str).unique().tolist()),
        "tissues": sorted(obs["tissue"].dropna().astype(str).unique().tolist()),
        "python": platform.python_version(),
        "cellxgene_census": cellxgene_census.__version__,
    }
    (out / "GSE176078_census_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
