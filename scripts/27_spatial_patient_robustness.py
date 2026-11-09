#!/usr/bin/env python3
"""Patient-resampling sensitivity of two preserved spatial partial correlations.

This script never reads spots or refits correlations. It resamples the existing
14 patient summaries jointly across the two comparisons, preserving their
within-patient pairing. Percentile intervals are nominal exploratory intervals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


COMPARISONS = (
    "IFN_adjusted_for_lineage_modules",
    "Checkpoint_adjusted_for_IFN_and_lineage_modules",
)


def quantile_linear(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--bootstrap-replicates", type=int, default=20000)
    args = parser.parse_args()
    source = args.input or args.project / "results/spatial/GSE210616/patient_partial_spearman.tsv"
    output = args.output_dir or args.project / "results/spatial/GSE210616/robustness"
    if args.bootstrap_replicates < 1000:
        raise ValueError("At least 1000 bootstrap replicates are required")
    raw = source.read_bytes()
    with source.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = list(reader.fieldnames or [])
        original = list(reader)
    required = {"patient_id", "analysis", "outcome", "controls", "partial_spearman_rho"}
    if not required.issubset(columns):
        raise ValueError(f"Missing required columns: {sorted(required-set(columns))}")
    if set(row["analysis"] for row in original) != set(COMPARISONS):
        raise ValueError("Input must contain exactly the two specified comparisons")
    patients = sorted({row["patient_id"] for row in original})
    if len(patients) != 14 or len(original) != 28:
        raise ValueError("Expected 14 patients with exactly two rows each")
    records = {}
    for row in original:
        key = (row["patient_id"], row["analysis"])
        if key in records:
            raise ValueError(f"Duplicate patient/comparison: {key}")
        value = float(row["partial_spearman_rho"])
        if not math.isfinite(value) or not -1 <= value <= 1:
            raise ValueError(f"Invalid partial correlation for {key}")
        records[key] = value
    if any((patient, comparison) not in records for patient in patients for comparison in COMPARISONS):
        raise ValueError("A paired patient/comparison value is missing")
    vectors = {name: [records[(patient, name)] for patient in patients] for name in COMPARISONS}
    full_medians = {name: statistics.median(values) for name, values in vectors.items()}
    output.mkdir(parents=True, exist_ok=True)
    # Preserve the exact source bytes rather than reserializing their decimals.
    (output / "actual_patient_partial_correlations.tsv").write_bytes(raw)
    loo = []
    for omitted_index, omitted_patient in enumerate(patients):
        for comparison in COMPARISONS:
            values = [value for index, value in enumerate(vectors[comparison]) if index != omitted_index]
            estimate = statistics.median(values)
            loo.append({
                "analysis": comparison,
                "omitted_patient_id": omitted_patient,
                "n_patients_retained": len(values),
                "full_patient_median_rho": full_medians[comparison],
                "leave_one_patient_median_rho": estimate,
                "difference_from_full_median": estimate-full_medians[comparison],
            })
    write_tsv(output / "actual_patient_LOO.tsv", list(loo[0]), loo)
    generator = random.Random(args.seed)
    bootstrap = []
    bootstrap_values = {name: [] for name in COMPARISONS}
    for replicate in range(1, args.bootstrap_replicates+1):
        indices = [generator.randrange(len(patients)) for _ in patients]
        row = {"bootstrap_replicate": replicate}
        for comparison in COMPARISONS:
            estimate = statistics.median([vectors[comparison][index] for index in indices])
            bootstrap_values[comparison].append(estimate)
            row[comparison] = estimate
        row["resampled_patient_ids"] = ";".join(patients[index] for index in indices)
        bootstrap.append(row)
    write_tsv(output / "actual_patient_bootstrap_medians.tsv", list(bootstrap[0]), bootstrap)
    summary = []
    for comparison in COMPARISONS:
        leave_out = [row["leave_one_patient_median_rho"] for row in loo if row["analysis"] == comparison]
        lower = quantile_linear(bootstrap_values[comparison], .025)
        upper = quantile_linear(bootstrap_values[comparison], .975)
        summary.append({
            "analysis": comparison,
            "n_independent_patients": len(patients),
            "median_patient_partial_rho": full_medians[comparison],
            "bootstrap_percentile_CI95_low": lower,
            "bootstrap_percentile_CI95_high": upper,
            "bootstrap_replicates": args.bootstrap_replicates,
            "seed": args.seed,
            "leave_one_patient_median_min": min(leave_out),
            "leave_one_patient_median_max": max(leave_out),
            "patients_with_positive_partial_rho": sum(value > 0 for value in vectors[comparison]),
            "leave_one_patient_positive_medians_n": sum(value > 0 for value in leave_out),
            "interval_scope": "nominal exploratory percentile interval; not simultaneous or external validation",
        })
    write_tsv(output / "patient_median_bootstrap_CI_summary.tsv", list(summary[0]), summary)
    audit = {
        "status": "PASS_PATIENT_LEVEL_RESAMPLING",
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_rows": len(original),
        "independent_patients": len(patients),
        "comparisons": list(COMPARISONS),
        "patients": patients,
        "resampling_unit": "patient; the same 14 sampled patient indices used jointly for both comparisons per replicate",
        "bootstrap_replicates": args.bootstrap_replicates,
        "random_generator": "Python random.Random / Mersenne Twister",
        "seed": args.seed,
        "CI_method": "2.5th and 97.5th bootstrap percentiles with linear interpolation",
        "leave_one_patient_out_runs": len(patients),
        "new_spot_level_models_fitted": False,
        "new_p_values_or_multiplicity_claims": False,
        "input_mutated": False,
        "summary": summary,
        "limits": [
            "Conditional uncertainty of medians of already-computed patient correlations; no spot-level measurement/model uncertainty propagated.",
            "Fourteen patients make bootstrap distributions discrete and intervals imprecise; percentile intervals are exploratory, not simultaneous coverage claims.",
            "Leave-one-patient-out summarizes influence; it is not cross-validation, independent replication or a new biological experiment.",
            "No absolute cell attribution, single-cell co-expression, causal regulation, treatment prediction or independent prognosis is established.",
        ],
    }
    (output / "SPATIAL_PATIENT_ROBUSTNESS_AUDIT.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps({"status": audit["status"], "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
