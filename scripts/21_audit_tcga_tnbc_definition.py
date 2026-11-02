#!/usr/bin/env python3
"""Reconstruct and audit the archived receptor-defined TCGA TNBC sample list."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from openpyxl import load_workbook


def sheet_ids(workbook, sheet: str) -> set[str]:
    values = set()
    for row in workbook[sheet].iter_rows(min_row=2, values_only=True):
        if row[0] is not None and str(row[0]).strip():
            values.add(str(row[0]).strip())
    return values


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--raw-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    project = args.project.resolve()
    raw_root = args.raw_root.resolve()

    workbook_path = raw_root / "数据/差异富集分析/1 乳腺癌分型/3_乳腺癌分型.xlsx"
    list_path = raw_root / "数据/差异富集分析/3 差异分析/三阴.txt"
    matrix_path = raw_root / "数据/差异富集分析/3 差异分析/symbol.BRCA.txt"
    output_dir = project / "results/bulk"

    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    er_negative = sheet_ids(workbook, "ER-")
    pr_negative = sheet_ids(workbook, "pr-")
    her2_negative = sheet_ids(workbook, "HER2-")
    basal_like = sheet_ids(workbook, "Basal-like")
    reconstructed = er_negative & pr_negative & her2_negative
    listed = {line.strip() for line in list_path.read_text(encoding="utf-8").splitlines() if line.strip()}

    with matrix_path.open("r", encoding="utf-8") as handle:
        matrix_samples = set(handle.readline().rstrip("\n").split("\t")[1:])

    if reconstructed != listed:
        raise RuntimeError("Archived TNBC list does not equal ER-negative ∩ PR-negative ∩ HER2-negative")

    manifest = []
    for sample_id in sorted(reconstructed):
        manifest.append({
            "sample_id": sample_id,
            "ER_negative": True,
            "PR_negative": True,
            "HER2_negative": True,
            "PAM50_basal_like_in_archived_workbook": sample_id in basal_like,
            "present_in_expression_matrix": sample_id in matrix_samples,
            "included_in_analysis": sample_id in matrix_samples,
        })
    write_tsv(
        output_dir / "TCGA_TNBC_sample_definition_manifest.tsv",
        manifest,
        ["sample_id", "ER_negative", "PR_negative", "HER2_negative",
         "PAM50_basal_like_in_archived_workbook", "present_in_expression_matrix", "included_in_analysis"],
    )

    audit = [
        {"metric": "definition", "value": "intersection_of_archived_ER_negative_PR_negative_HER2_negative_lists"},
        {"metric": "equivocal_or_unknown_included", "value": "false"},
        {"metric": "ER_negative_n", "value": len(er_negative)},
        {"metric": "PR_negative_n", "value": len(pr_negative)},
        {"metric": "HER2_negative_n", "value": len(her2_negative)},
        {"metric": "reconstructed_TNBC_n", "value": len(reconstructed)},
        {"metric": "archived_TNBC_list_n", "value": len(listed)},
        {"metric": "exact_list_match", "value": str(reconstructed == listed).lower()},
        {"metric": "PAM50_basal_like_overlap_n", "value": len(reconstructed & basal_like)},
        {"metric": "PAM50_used_for_inclusion", "value": "false"},
        {"metric": "expression_matrix_samples_n", "value": len(matrix_samples)},
        {"metric": "TNBC_matched_to_expression_matrix_n", "value": len(reconstructed & matrix_samples)},
        {"metric": "TNBC_missing_from_expression_matrix_n", "value": len(reconstructed - matrix_samples)},
    ]
    write_tsv(output_dir / "TCGA_TNBC_definition_audit.tsv", audit, ["metric", "value"])
    print(f"Reconstructed {len(reconstructed)} receptor-defined TNBC samples; "
          f"{len(reconstructed & matrix_samples)} matched the expression matrix.")


if __name__ == "__main__":
    main()
