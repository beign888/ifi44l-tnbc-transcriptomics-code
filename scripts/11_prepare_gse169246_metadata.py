#!/usr/bin/env python3
"""Build auditable patient/sample/cell metadata for GSE169246."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import openpyxl
import pandas as pd


def clinical_table(path: Path) -> pd.DataFrame:
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True)["Clinical information"]
    rows = []
    treatment = ""
    for row in ws.iter_rows(min_row=4, values_only=True):
        if not row[1]:
            continue
        if row[0]:
            treatment = str(row[0]).strip()
        efficacy = str(row[17]).strip() if row[17] is not None else ""
        rows.append(
            {
                "patient_id": str(row[1]).strip(),
                "treatment_arm": treatment,
                "biopsied_lesion": str(row[2]).strip() if row[2] is not None else "",
                "age": row[5],
                "PD_L1_reported": row[7],
                "TIL_percent_reported": row[8],
                "clinical_efficacy": efficacy,
                "response_binary": "Responder" if efficacy in {"CR", "PR"} else ("Non-responder" if efficacy in {"SD", "PD"} else "Unknown"),
                "time_to_progression_months_reported": row[18],
            }
        )
    return pd.DataFrame(rows)


def barcode_table(path: Path) -> pd.DataFrame:
    records = []
    with gzip.open(path, "rt") as handle:
        for cell_index, line in enumerate(handle):
            barcode = line.strip()
            match = re.search(r"\.(Pre|Post|Prog)_(P\d+)_([bt])$", barcode)
            if not match:
                raise RuntimeError(f"Unrecognized barcode suffix: {barcode}")
            timepoint, patient, source = match.groups()
            records.append(
                {
                    "cell_index_1based": cell_index + 1,
                    "cell_barcode": barcode,
                    "sample_id": f"{timepoint}_{patient}_{source}",
                    "timepoint": timepoint,
                    "patient_id": patient,
                    "sample_source": {"b": "blood", "t": "tumor"}[source],
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--barcodes", type=Path, required=True)
    parser.add_argument("--clinical", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    clinical = clinical_table(args.clinical)
    cells = barcode_table(args.barcodes).merge(clinical, on="patient_id", how="left", validate="many_to_one")
    samples = (
        cells.groupby(
            ["sample_id", "patient_id", "timepoint", "sample_source", "treatment_arm", "clinical_efficacy", "response_binary"],
            dropna=False,
        )
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["treatment_arm", "patient_id", "sample_source", "timepoint"])
    )
    clinical.to_csv(args.output_dir / "GSE169246_patient_clinical.tsv", sep="\t", index=False)
    samples.to_csv(args.output_dir / "GSE169246_sample_inventory.tsv", sep="\t", index=False)
    cells.to_csv(args.output_dir / "GSE169246_cell_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"patients={clinical.patient_id.nunique()} samples={samples.sample_id.nunique()} cells={len(cells)}")
    print(samples.groupby(["treatment_arm", "sample_source", "timepoint", "response_binary"], dropna=False)["n_cells"].agg(["count", "sum"]).to_string())


if __name__ == "__main__":
    main()
