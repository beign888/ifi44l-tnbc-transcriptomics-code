#!/usr/bin/env python3
"""Extract auditable donor/subtype and cell-annotation tables from Wu et al. 2021."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    source = (
        args.source
        or project / "source_metadata" / "Wu2021_NatGenet_Supplementary_Tables.xlsx"
    ).resolve()
    out = project / "results" / "single_cell"
    out.mkdir(parents=True, exist_ok=True)

    raw = pd.read_excel(source, sheet_name="Supplementary Table 6", header=None)
    blocks = []
    for start in (0, 7, 14):
        block = raw.iloc[4:, start : start + 6].copy()
        block.columns = [
            "donor_id",
            "author_cell_type",
            "cell_number",
            "cell_proportion",
            "clinical_subtype",
            "classification_tier",
        ]
        block = block.dropna(subset=["donor_id", "author_cell_type"])
        block["donor_id"] = block["donor_id"].astype(str).str.strip()
        block["author_cell_type"] = block["author_cell_type"].astype(str).str.strip()
        block["clinical_subtype"] = block["clinical_subtype"].astype(str).str.strip()
        block["classification_tier"] = block["classification_tier"].astype(str).str.strip()
        block["cell_number"] = pd.to_numeric(block["cell_number"], errors="coerce").astype("Int64")
        block["cell_proportion"] = pd.to_numeric(block["cell_proportion"], errors="coerce")
        blocks.append(block)

    annotations = pd.concat(blocks, ignore_index=True).drop_duplicates()
    annotations.to_csv(out / "Wu2021_author_cell_annotations.tsv", sep="\t", index=False)
    conflicts = (
        annotations.groupby("donor_id", observed=True)["clinical_subtype"]
        .nunique()
        .loc[lambda values: values > 1]
    )
    if len(conflicts):
        raise RuntimeError(f"Conflicting clinical subtypes: {conflicts.to_dict()}")

    donors = (
        annotations[["donor_id", "clinical_subtype"]]
        .drop_duplicates()
        .sort_values(["clinical_subtype", "donor_id"])
    )
    donors.to_csv(out / "Wu2021_donor_clinical_subtypes.tsv", sep="\t", index=False)
    tnbc = donors.loc[donors["clinical_subtype"].eq("TNBC")].copy()
    tnbc.to_csv(out / "Wu2021_TNBC_donors.tsv", sep="\t", index=False)

    audit = {
        "source_file": str(source),
        "source_sheet": "Supplementary Table 6",
        "n_donors": int(donors["donor_id"].nunique()),
        "n_tnbc_donors": int(tnbc["donor_id"].nunique()),
        "tnbc_donors": tnbc["donor_id"].tolist(),
        "annotation_tiers": sorted(annotations["classification_tier"].unique().tolist()),
    }
    (out / "Wu2021_supplement_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
