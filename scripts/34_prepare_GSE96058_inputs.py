#!/usr/bin/env python3
"""Acquire official SCAN-B inputs and freeze an exact all-gene clinical-TNBC subset.

Expression is streamed as decimal strings; this helper performs no IFI44L or
outcome association analysis. It does not alter the source log2(FPKM + 0.1).
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
from urllib.request import Request, urlopen

BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/'
EXPRESSION = 'GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz'
MATRICES = [f'GSE96058-{platform}_series_matrix.txt.gz' for platform in ['GPL11154', 'GPL18573']]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def download(url, path):
    if path.exists():
        return {'file': path.name, 'url': url, 'action': 'existing complete file retained',
                'bytes': path.stat().st_size, 'sha256': sha(path)}
    partial = path.with_name(path.name + '.part')
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with urlopen(Request(url, headers={'User-Agent': 'IFI44L-public-data-reanalysis/1.0'}), timeout=90) as response:
        if response.status != 200:
            raise RuntimeError(f'Expected complete HTTP200 response: {url}: {response.status}')
        total = response.headers.get('Content-Length')
        written = 0
        with partial.open('wb') as target:
            for block in iter(lambda: response.read(2**20), b''):
                target.write(block)
                written += len(block)
                if written % (64 * 2**20) < 2**20:
                    print(json.dumps({'download': path.name, 'bytes': written, 'content_length': total}), flush=True)
        if total is not None and written != int(total):
            raise RuntimeError(f'Incomplete download {written}/{total}: {url}')
        headers = {k: response.headers.get(k) for k in ['Content-Length', 'Last-Modified', 'ETag', 'Content-Type']}
    partial.replace(path)
    return {'file': path.name, 'url': url, 'action': 'actual HTTP200 full acquisition',
            'started_utc': started, 'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'bytes': written, 'sha256': sha(path), 'headers': headers}


def metadata(path):
    rows = []
    with gzip.open(path, 'rt', newline='') as handle:
        for line in handle:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                rows.append(next(csv.reader([line], delimiter='\t')))
    titles = next(row[1:] for row in rows if row[0] == '!Sample_title')
    records = [{'platform_id': re.search(r'GPL\d+', path.name).group()} for _ in titles]
    processing = set()
    for row in rows:
        key, values = row[0], row[1:]
        if key in ['!Sample_title', '!Sample_geo_accession', '!Sample_source_name_ch1']:
            for rec, val in zip(records, values):
                rec[key.replace('!Sample_', '')] = val
        elif key == '!Sample_characteristics_ch1':
            for rec, val in zip(records, values):
                name, content = val.split(': ', 1)
                rec[name] = content
        elif key == '!Sample_relation':
            for rec, val in zip(records, values):
                if val.startswith('BioSample: '):
                    rec['BioSample'] = val.split(': ', 1)[1]
        elif key == '!Sample_data_processing':
            processing.update(values)
    return records, sorted(processing)


def write_tsv(path, records, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({k: ('TRUE' if v else 'FALSE') if isinstance(v, bool) else v
                         for k, v in record.items()} for record in records)


def restore_parts(work):
    """Reassemble the exact source gzip offline; do not recompress or delete parts."""
    provenance = json.loads((work / 'GSE96058_frozen_input_provenance.json').read_text())
    target = work / provenance['subset_file']
    expected_sha = provenance['subset_sha256']
    if expected_sha != '2e49a251c2f9129658a2924e78372ebaf4fb936466e1553d7a7e1e6410e5b5ba':
        raise RuntimeError('The frozen v16 subset SHA differs; do not restore another dataset')
    expected_bytes = next(r['bytes'] for r in provenance['files'] if r['file'] == target.name)
    if target.exists():
        if target.stat().st_size != expected_bytes or sha(target) != expected_sha:
            raise RuntimeError('Existing subset fails identity checks; it was not overwritten')
        print(json.dumps({'status': 'PASS_EXISTING_FROZEN_SUBSET', 'bytes': expected_bytes,
                          'sha256': expected_sha, 'network_used': False}), flush=True)
        return
    parts = [work / (target.name + suffix) for suffix in ['.part01', '.part02']]
    for part in parts:
        if not part.exists():
            raise FileNotFoundError(f'Missing offline lossless volume: {part}')
    partial = target.with_name(target.name + '.tmp')
    with partial.open('wb') as handle:
        for part in parts:
            with part.open('rb') as source:
                shutil.copyfileobj(source, handle, 2**20)
    if partial.stat().st_size != expected_bytes or sha(partial) != expected_sha:
        raise RuntimeError('Reassembled subset failed byte-count/SHA256; temporary file retained, no final replacement')
    partial.replace(target)
    record = {'status': 'PASS_OFFLINE_EXACT_GZIP_RESTORATION', 'bytes': expected_bytes,
              'sha256': expected_sha, 'network_used': False, 'parts_retained': True,
              'parts': [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in parts]}
    (work / 'GSE96058_OFFLINE_RESTORATION_RECEIPT.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', default='.')
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--restore-parts', action='store_true', help='Offline exact-gzip concatenation only; never download or recompress')
    parser.add_argument('--analysis-plan')
    parser.add_argument('--expected-plan-sha256', default='f174e513f1c3091a19d0ab3233512048c64e0b02f564a46101fd324ab550a7ae')
    parser.add_argument('--analysis-plan-amendment')
    parser.add_argument('--expected-amendment-sha256', default='06b06ed2d8879e2b1aecc2df094eca9baabcd317c1686d55c0801d2029f8db22')
    args = parser.parse_args()
    project = Path(args.project).resolve()
    work = project / 'source_metadata/GSE96058'
    work.mkdir(parents=True, exist_ok=True)
    if args.restore_parts:
        if args.download:
            raise ValueError('--restore-parts is offline and cannot be combined with --download')
        restore_parts(work)
        return
    plan = Path(args.analysis_plan) if args.analysis_plan else project / 'protocols/v16/GSE96058_ANALYSIS_PLAN.md'
    if not plan.exists() and not args.analysis_plan:
        plan = project / 'audit_v16/GSE96058_ANALYSIS_PLAN.md'
    if sha(plan) != args.expected_plan_sha256:
        raise RuntimeError('Frozen analysis plan differs before source preparation')
    fixed_plan = work / 'GSE96058_ANALYSIS_PLAN.md'
    if fixed_plan.exists() and sha(fixed_plan) != sha(plan):
        raise RuntimeError('Do not overwrite a different fixed plan')
    shutil.copyfile(plan, fixed_plan)
    amendment = Path(args.analysis_plan_amendment) if args.analysis_plan_amendment else project / 'protocols/v16/GSE96058_PLATFORM_AMENDMENT.md'
    if not amendment.exists() and not args.analysis_plan_amendment:
        amendment = project / 'audit_v16/GSE96058_PLATFORM_AMENDMENT.md'
    if sha(amendment) != args.expected_amendment_sha256:
        raise RuntimeError('Frozen platform amendment differs before source preparation')
    shutil.copyfile(amendment, work / 'GSE96058_PLATFORM_AMENDMENT.md')
    downloads = []
    for name in MATRICES + [EXPRESSION]:
        url = BASE + ('matrix/' if name in MATRICES else 'suppl/') + name
        path = work / name
        if args.download:
            downloads.append(download(url, path))
        elif not path.exists():
            raise FileNotFoundError(f'Run with --download to obtain {url}')
    acquisition = work / 'GSE96058_download_receipt.json'
    if downloads:
        previous = json.loads(acquisition.read_text()) if acquisition.exists() else []
        acquisition.write_text(json.dumps(previous + downloads, ensure_ascii=False, indent=2) + '\n')
    records, processing = [], set()
    for name in MATRICES:
        rec, proc = metadata(work / name)
        records += rec
        processing.update(proc)
    assert len(records) == 3409 and len({r['geo_accession'] for r in records}) == 3409
    for r in records:
        r['sample_id'] = r['geo_accession']
        r['patient_source_id'] = re.sub('repl$', '', r['title'], flags=re.I)
        r['is_technical_replicate'] = bool(re.search('repl$', r['title'], re.I))
        r['eligible_receptor_TNBC'] = (not r['is_technical_replicate']) and all(r[k] == '0' for k in ['er status', 'pgr status', 'her2 status'])
        reasons = []
        if r['is_technical_replicate']:
            reasons.append('explicit source technical replicate')
        for key in ['er status', 'pgr status', 'her2 status']:
            if r[key] != '0':
                reasons.append(f'{key}={r[key]} (positive or missing clinical call)')
        r['selection_reason'] = '; '.join(reasons) if reasons else 'included primary clinical three-negative case'
    primary = [r for r in records if not r['is_technical_replicate']]
    selected = [r for r in records if r['eligible_receptor_TNBC']]
    assert len(primary) == 3273 and len(selected) == 143
    assert len({r['patient_source_id'] for r in primary}) == 3273
    by_title = {r['title']: r for r in primary}
    pairs = []
    clinical_fields = ['age at diagnosis', 'tumor size', 'lymph node status', 'er status', 'pgr status', 'her2 status', 'nhg', 'overall survival days', 'overall survival event']
    for r in records:
        if r['is_technical_replicate']:
            original = by_title[r['patient_source_id']]
            assert all(r[k] == original[k] for k in clinical_fields)
            pairs.append({'primary_title': original['title'], 'replicate_title': r['title'],
                          'primary_sample_id': original['sample_id'], 'replicate_sample_id': r['sample_id'],
                          'same_recorded_clinical_fields': True,
                          'same_deposited_BioSample_accession': r.get('BioSample') == original.get('BioSample'),
                          'policy': 'retain original primary; exclude explicit technical replicate without expression/outcome selection'})
    assert len(pairs) == 136
    fields = list(dict.fromkeys(k for r in records for k in r))
    write_tsv(work / 'GSE96058_all3409_metadata.tsv', records, fields)
    write_tsv(work / 'GSE96058_TNBC143_metadata.tsv', selected, fields)
    write_tsv(work / 'GSE96058_technical_replicate_pairs.tsv', pairs, list(pairs[0]))
    (work / 'GSE96058_source_data_processing.txt').write_text('\n\n'.join(sorted(processing)) + '\n')
    module_source = project / 'results/bulk/immune_signature_definitions.tsv'
    if not module_source.exists():
        module_source = project / 'source_metadata/GSE21653/baseline_module_definitions.tsv'
    shutil.copyfile(module_source, work / 'baseline_module_definitions.tsv')

    # Stream the full official CSV; never load all3409 columns into memory.
    selected = sorted(selected, key=lambda r: int(r['title'][1:]))
    subset = work / 'GSE96058_TNBC143_frozen_all_genes.tsv.gz'
    expr = work / EXPRESSION
    gene_ids = set()
    rows = 0
    with gzip.open(expr, 'rt', newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader)
        assert len(header) == 3410 and len(set(header[1:])) == 3409
        assert set(header[1:]) == {r['title'] for r in records}
        (work / 'GSE96058_complete_source_header.tsv').write_text('\t'.join(header) + '\n')
        indexes = [header.index(r['title']) for r in selected]
        with subset.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding='utf-8', newline='') as text:
                    writer = csv.writer(text, delimiter='\t', lineterminator='\n')
                    writer.writerow(['gene'] + [r['sample_id'] for r in selected])
                    for row in reader:
                        assert len(row) == len(header), (rows, len(row), len(header))
                        gene = row[0]
                        assert gene and gene not in gene_ids
                        gene_ids.add(gene)
                        values = [row[i] for i in indexes]
                        assert all(math.isfinite(float(x)) for x in values)
                        writer.writerow([gene] + values)
                        rows += 1
    assert rows == 30865 and 'IFI44L' in gene_ids
    linkage = [{'expression_title': r['title'], 'sample_id': r['sample_id'],
                'patient_source_id': r['patient_source_id'], 'platform_id': r['platform_id'],
                'source_column_index_0based': i} for r, i in zip(selected, indexes)]
    write_tsv(work / 'GSE96058_expression_linkage.tsv', linkage, list(linkage[0]))
    files = [p for p in work.iterdir() if p.is_file() and p.name != 'GSE96058_frozen_input_provenance.json']
    file_records = [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in sorted(files)]
    provenance = {'status': 'ACTUAL_FULL_SOURCE_AND_EXACT_ALL_GENE_143_SUBSET_PREPARED',
                  'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'full_source_file': EXPRESSION, 'full_source_sha256': sha(expr),
                  'subset_file': subset.name, 'subset_sha256': sha(subset),
                  'metadata_sha256': sha(work / 'GSE96058_all3409_metadata.tsv'),
                  'linkage_sha256': sha(work / 'GSE96058_expression_linkage.tsv'),
                  'module_definitions_sha256': sha(work / 'baseline_module_definitions.tsv'),
                  'analysis_plan_sha256': sha(plan), 'gene_rows': rows,
                  'platform_amendment_sha256': sha(amendment),
                  'full_source_expression_samples': 3409, 'primary_patients': 3273,
                  'excluded_technical_replicates': 136, 'eligible_clinical_TNBC': 143,
                  'source_expression_scale': 'log2(FPKM+0.1), source transcript FPKM summed by gene symbol; no further transform',
                  'source_values': 'exact decimal strings preserved for every gene and selected sample; gzip CRC checked by reading to EOF',
                  'phenotype': 'clinical er status=pgr status=her2 status=0; no MGC/SGC/PAM50 replacement; source diagnostic binary calls, no retrospectively imposed modern cutoff',
                  'BioSample_caveat': 'Official design describes same biological case/RNA source; deposited BioSample accessions differ, so paired source titles and matching clinical fields establish replicate exclusion policy.',
                  'no_analysis': 'No IFI44L association, outcome model, cutpoint or phenotype selection based on expression performed by preparation helper.',
                  'downloads': downloads, 'files': file_records}
    (work / 'GSE96058_frozen_input_provenance.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: provenance[k] for k in ['status', 'full_source_sha256', 'subset_sha256', 'gene_rows', 'eligible_clinical_TNBC']}), flush=True)


if __name__ == '__main__':
    main()
