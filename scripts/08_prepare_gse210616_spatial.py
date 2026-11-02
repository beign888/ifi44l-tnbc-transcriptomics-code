#!/usr/bin/env python3
"""Parse GSE210616 SOFT metadata and optionally download the 28 primary Visium sections.

The series contains later-added records beyond the 28-section primary design.
This script selects the first 28 samples explicitly listed in GEO, representing
14 patients with two sections each, and excludes very large .cloupe files.
"""

from __future__ import annotations

import argparse
import gzip
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd


KEEP_SUFFIXES = (
    "_filtered_feature_bc_matrix.h5",
    "_scalefactors_json.json.gz",
    "_tissue_hires_image.png.gz",
    "_tissue_positions_list.csv.gz",
)


def parse_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"gsm": line.split("=", 1)[1].strip(), "urls": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
                match = re.search(r"patient\s+(\d+),\s*section\s+(\d+)", str(current["title"]), re.I)
                if match:
                    current["patient_id"] = f"P{int(match.group(1)):02d}"
                    current["section"] = int(match.group(2))
            elif current is not None and line.startswith("!Sample_supplementary_file_"):
                url = line.split("=", 1)[1].strip()
                if url.endswith(KEEP_SUFFIXES):
                    current["urls"].append(url.replace("ftp://", "https://"))
        if current:
            records.append(current)

    rows = []
    # GEO overall design specifies 28 sections from 14 primary tumors.
    for record in records[:28]:
        for url in record["urls"]:
            rows.append(
                {
                    "gsm": record["gsm"],
                    "patient_id": record.get("patient_id", ""),
                    "section": record.get("section", ""),
                    "title": record.get("title", ""),
                    "file_name": url.rsplit("/", 1)[-1],
                    "url": url,
                }
            )
    frame = pd.DataFrame(rows)
    expected = 28 * len(KEEP_SUFFIXES)
    if len(frame) != expected:
        raise RuntimeError(f"Expected {expected} files for 28 sections, found {len(frame)}")
    if frame["patient_id"].nunique() != 14:
        raise RuntimeError("The selected GEO records do not resolve to 14 patients")
    return frame


def download(frame: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for i, row in frame.iterrows():
        sample_dir = output_dir / str(row["gsm"])
        sample_dir.mkdir(parents=True, exist_ok=True)
        target = sample_dir / str(row["file_name"])
        if target.exists() and target.stat().st_size > 0:
            print(f"[{i + 1}/{len(frame)}] exists {target.name}")
            continue
        tmp = target.with_suffix(target.suffix + ".part")
        print(f"[{i + 1}/{len(frame)}] download {row['gsm']} {target.name}")
        for attempt in range(1, 4):
            try:
                req = urllib.request.Request(str(row["url"]), headers={"User-Agent": "IFI44L-reanalysis/1.0"})
                with urllib.request.urlopen(req, timeout=120) as response, tmp.open("wb") as out:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        out.write(block)
                tmp.replace(target)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()

    frame = parse_soft(args.soft)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "GSE210616_primary28_manifest.tsv", sep="\t", index=False)
    sample_meta = frame[["gsm", "patient_id", "section", "title"]].drop_duplicates()
    sample_meta.to_csv(args.output_dir / "GSE210616_primary28_samples.tsv", sep="\t", index=False)
    print(sample_meta.to_string(index=False))
    if args.download:
        download(frame, args.output_dir / "files")


if __name__ == "__main__":
    main()
