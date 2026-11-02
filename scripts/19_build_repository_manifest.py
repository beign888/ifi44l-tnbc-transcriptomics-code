#!/usr/bin/env python3
"""Build the repository SHA-256 manifest from the current tracked-file snapshot."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXCLUDED_PARTS = {".git", "__pycache__", ".venv"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    repository = args.repository.resolve()
    output = repository / "FILE_MANIFEST_SHA256.tsv"
    files = [
        path
        for path in repository.rglob("*")
        if path.is_file()
        and path != output
        and not any(part in EXCLUDED_PARTS for part in path.relative_to(repository).parts)
        and not path.name.startswith("._")
        and path.name != ".DS_Store"
    ]
    rows = ["relative_path\tsize_bytes\tsha256"]
    for path in sorted(files):
        relative = path.relative_to(repository).as_posix()
        rows.append(f"{relative}\t{path.stat().st_size}\t{sha256(path)}")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"manifest_files={len(files)} output={output}")


if __name__ == "__main__":
    main()
