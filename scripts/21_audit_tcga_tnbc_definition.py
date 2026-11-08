#!/usr/bin/env python3
"""Reconstruct and audit the archived receptor-defined TCGA TNBC sample list."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from openpyxl import load_workbook

CLINICAL_SHA256 = "973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc"
CLINICAL_COMMIT = "4c975725c54df7e31e0610fad22a742ad01cdfa3"


def reconcile_clinical(project: Path, raw_root: Path, strict_ihc: bool) -> None:
    """Freeze receptor eligibility from clinical fields before reading outcomes."""
    clinical = project / "source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt"
    if not clinical.is_file():
        raise FileNotFoundError(f"Required pinned clinical file: {clinical}; see PUBLIC_DATA_ACCESS.md")
    if hashlib.sha256(clinical.read_bytes()).hexdigest() != CLINICAL_SHA256:
        raise RuntimeError("Clinical source SHA-256 mismatch; do not silently use a different version")
    with clinical.open(encoding="utf-8") as stream:
        patients = {r["PATIENT_ID"]: r for r in csv.DictReader(
            (line for line in stream if not line.startswith("#")), delimiter="\t")}
    matrix = raw_root / "数据/差异富集分析/3 差异分析/symbol.BRCA.txt"
    if matrix.exists():
        with matrix.open(encoding="utf-8") as stream:
            samples = stream.readline().rstrip("\r\n").split("\t")[1:]
    else:
        frozen_header = project / "source_metadata/TCGA_frozen/TCGA_BRCA_expression_header_samples.tsv"
        with frozen_header.open(encoding="utf-8") as stream:
            samples = [r["sample_id"] for r in csv.DictReader(stream, delimiter="\t")]
    legacy_path = raw_root / "数据/差异富集分析/3 差异分析/三阴.txt"
    legacy = set(legacy_path.read_text(encoding="utf-8").split()) if legacy_path.exists() else set()
    rows = []
    for sample in samples:
        patient = sample[:12]
        record = patients.get(patient, {})
        er = record.get("ER_STATUS_BY_IHC", "unavailable")
        pr = record.get("PR_STATUS_BY_IHC", "unavailable")
        ihc = record.get("IHC_HER2", "unavailable")
        fish = record.get("HER2_FISH_STATUS", "unavailable")
        score = record.get("HER2_IHC_SCORE", "unavailable")
        hormone_negative = er == "Negative" and pr == "Negative"
        positive = ihc == "Positive" or fish == "Positive" or score == "3"
        negative_evidence = ihc == "Negative" or fish == "Negative"
        strict = hormone_negative and not positive and ihc == "Negative"
        included = hormone_negative and not positive and negative_evidence
        if strict_ihc:
            included = strict
        reason = ("included" if included else "ER_or_PR_not_explicitly_negative" if not hormone_negative
                  else "HER2_positive_or_IHC_3plus" if positive else "HER2_not_resolved_negative")
        rows.append({
            "sample_id": sample, "patient_id": patient,
            "ER_status_by_IHC": er, "PR_status_by_IHC": pr,
            "HER2_IHC_status": ihc, "HER2_FISH_status": fish, "HER2_IHC_score": score,
            "ER_negative": er == "Negative", "PR_negative": pr == "Negative",
            "HER2_negative": not positive and negative_evidence,
            "strict_IHC_negative_eligible": strict,
            "included_in_analysis": included, "exclusion_reason": reason,
            "present_in_expression_matrix": True,
            "in_legacy_archived_list": sample in legacy,
        })
    output = project / "results/bulk"
    write_tsv(output / "TCGA_TNBC_sample_definition_manifest.tsv", rows, list(rows[0]))
    selected = [r for r in rows if r["included_in_analysis"]]
    if len({r["patient_id"] for r in selected}) != len(selected):
        raise RuntimeError("Multiple expression samples per eligible patient require explicit resolution")
    audit = [
        {"metric": "definition", "value": "clinical_ER_PR_negative_and_resolved_HER2_negative_no_positive_or_3plus"},
        {"metric": "clinical_source_commit", "value": CLINICAL_COMMIT},
        {"metric": "clinical_source_sha256", "value": CLINICAL_SHA256},
        {"metric": "expression_matrix_samples_n", "value": len(samples)},
        {"metric": "clinical_patients_matched_n", "value": sum(r["patient_id"] in patients for r in rows)},
        {"metric": "ER_PR_both_negative_n", "value": sum(r["ER_negative"] and r["PR_negative"] for r in rows)},
        {"metric": "TNBC_matched_to_expression_matrix_n", "value": len(selected)},
        {"metric": "strict_IHC_negative_n", "value": sum(r["strict_IHC_negative_eligible"] for r in rows)},
        {"metric": "strict_IHC_sensitivity_mode", "value": strict_ihc},
        {"metric": "PAM50_used_for_inclusion", "value": False},
        {"metric": "legacy_matched_n", "value": sum(r["in_legacy_archived_list"] for r in rows)},
        {"metric": "legacy_HER2_positive_or_3plus_n", "value": sum(
            r["in_legacy_archived_list"] and (r["HER2_IHC_status"] == "Positive" or
            r["HER2_FISH_status"] == "Positive" or r["HER2_IHC_score"] == "3") for r in rows)},
    ]
    write_tsv(output / "TCGA_TNBC_definition_audit.tsv", audit, ["metric", "value"])
    output.mkdir(parents=True, exist_ok=True)
    (output / "TCGA_TNBC_selected_samples.txt").write_text(
        "\n".join(r["sample_id"] for r in selected) + "\n", encoding="utf-8")
    print(json.dumps({"eligible_n": len(selected), "strict_IHC_n": sum(
        r["strict_IHC_negative_eligible"] for r in rows), "rule_frozen_without_outcome_access": True}))


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
    parser.add_argument("--legacy-archived", action="store_true", help="Audit historical lists only; not the current cohort")
    parser.add_argument("--strict-ihc", action="store_true", help="Require explicit HER2 IHC negativity as sensitivity")
    args = parser.parse_args()
    project = args.project.resolve()
    raw_root = args.raw_root.resolve()
    if not args.legacy_archived:
        reconcile_clinical(project, raw_root, args.strict_ihc)
        return

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
        {"metric": "equivocal_or_unknown_included", "value": "not_reconstructable_from_archived_lists"},
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
