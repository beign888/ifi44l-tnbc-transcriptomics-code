#!/usr/bin/env python3
"""Restore the exact TCGA subset from the two journal supporting-data archives."""
from pathlib import Path
import argparse, hashlib, json


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    a = p.parse_args()
    directory = a.project / 'source_metadata/TCGA_frozen'
    receipt = json.loads((directory / 'TCGA_frozen_input_provenance.json').read_text())
    output = directory / 'TCGA_BRCA_TNBC_frozen.tsv.gz'
    if not output.exists():
        parts = [directory / 'TCGA_BRCA_TNBC_frozen.tsv.gz.part01',
                 directory / 'TCGA_BRCA_TNBC_frozen.tsv.gz.part02']
        if not all(p.is_file() for p in parts):
            raise FileNotFoundError('Unzip both Additional files 3 and 4 into source_metadata first')
        with output.open('xb') as stream:
            for part in parts:
                with part.open('rb') as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b''):
                        stream.write(chunk)
    sha = hashlib.sha256(output.read_bytes()).hexdigest()
    if sha != receipt['frozen_subset_sha256']:
        raise RuntimeError('Restored TCGA input SHA-256 differs from the frozen public-data subset')
    print(json.dumps({'restored_input': output.name, 'sha256': sha, 'status': 'PASS'}))


if __name__ == '__main__':
    main()
