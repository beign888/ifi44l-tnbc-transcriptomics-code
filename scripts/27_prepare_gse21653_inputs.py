#!/usr/bin/env python3
"""Acquire pinned public GSE21653 inputs and preserve an exact all-probe TNBC subset.

No expression/outcome associations are calculated. The original RMA values are
preserved; missing clinical values and zero-time events are not repaired.
"""
from pathlib import Path
from urllib.request import Request, urlopen
from datetime import datetime, timezone
import argparse, csv, gzip, hashlib, io, json, shutil

SOURCES = {
    'GSE21653_series_matrix.txt.gz': (
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21653/matrix/GSE21653_series_matrix.txt.gz',
        'f13e3bf972bd92a6fd53a1fbe42ad7c84caf4651df3a5e836f70dd3b22b5e76f', 200_000_000),
    'GSE21653_GEO_GSM_brief_SOFT.txt': (
        'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21653&targ=gsm&view=brief&form=text',
        'c3c156c50f968b08f8eabca8f077f1d68428eda1958cc2c01e4099cadef47d13', 5_000_000)
}
MODULE_SHA = '8bc3ff195ccaef81d53e6f4543d44f584e16ade8a3b0e04fb8231c3068562bb3'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda:stream.read(8*1024*1024),b''):h.update(part)
    return h.hexdigest()

def acquire(path,url,expected,limit):
    temporary=path.with_name(path.name+'.download')
    request=Request(url,headers={'User-Agent':'IFI44L-public-data-reproduction/1.0'})
    total=0
    with urlopen(request,timeout=60) as response,temporary.open('wb') as target:
        advertised=response.headers.get('Content-Length')
        if advertised and int(advertised)>limit:raise ValueError('File exceeds acquisition bound')
        while True:
            chunk=response.read(1024*1024)
            if not chunk:break
            total+=len(chunk)
            if total>limit:raise ValueError('File exceeds acquisition bound')
            target.write(chunk)
    if advertised and total!=int(advertised):raise ValueError('Incomplete HTTP body; retain .download for diagnosis')
    if sha(temporary)!=expected:raise ValueError('Source differs from frozen SHA; do not silently use changed inputs')
    temporary.replace(path)

def write_tsv(path,rows,fields):
    with path.open('w',newline='') as target:
        writer=csv.DictWriter(target,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--source-dir',type=Path)
    parser.add_argument('--module-definitions',type=Path)
    parser.add_argument('--download',action='store_true')
    args=parser.parse_args()
    work=args.project/'source_metadata/GSE21653';work.mkdir(parents=True,exist_ok=True)
    if args.source_dir:
        for filename in SOURCES:
            source=args.source_dir/filename
            if source.resolve()!=(work/filename).resolve():shutil.copyfile(source,work/filename)
    if args.download:
        for filename,(url,expected,limit) in SOURCES.items():
            path=work/filename
            if not path.exists() or sha(path)!=expected:acquire(path,url,expected,limit)
    for filename,(_,expected,_) in SOURCES.items():
        if sha(work/filename)!=expected:raise ValueError('Pinned source hash mismatch: '+filename)
    module_path=work/'baseline_module_definitions.tsv'
    if args.module_definitions and args.module_definitions.resolve()!=module_path.resolve():
        shutil.copyfile(args.module_definitions,module_path)
    if sha(module_path)!=MODULE_SHA:raise ValueError('Module definitions differ from frozen primary analysis')

    text=(work/'GSE21653_GEO_GSM_brief_SOFT.txt').read_text()
    rows=[]
    for block in text.split('^SAMPLE = ')[1:]:
        row={'sample_id':block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_title = '):row['source_case_id']=line.split(' = ',1)[1]
            if line.startswith('!Sample_characteristics_ch1 = '):
                key,value=line.split(' = ',1)[1].split(': ',1);row[key]=value
            if line.startswith('!Sample_data_processing = '):row['data_processing']=line.split(' = ',1)[1]
        row['eligible_receptor_TNBC']=all(row.get(k)=='0' for k in ['er ihc','pr ihc','erbb2'])
        rows.append(row)
    assert len(rows)==266 and len({r['sample_id'] for r in rows})==266
    selected=[r for r in rows if r['eligible_receptor_TNBC']]
    assert len(selected)==87 and len({r['source_case_id'] for r in selected})==87
    assert {r['data_processing'] for r in rows}=={'RMA using a quantile normalization'}
    fields=['sample_id','source_case_id']+sorted(set().union(*(r.keys() for r in rows))-{'sample_id','source_case_id'})
    write_tsv(work/'GSE21653_all266_metadata.tsv',rows,fields)
    write_tsv(work/'GSE21653_TNBC87_metadata.tsv',selected,fields)
    (work/'GSE21653_TNBC87_selected_samples.txt').write_text(''.join(r['sample_id']+'\n' for r in selected))

    subset=work/'GSE21653_TNBC87_frozen_all_probes.tsv.gz'
    selected_ids=[r['sample_id'] for r in selected]
    with gzip.open(work/'GSE21653_series_matrix.txt.gz','rt',newline='') as source:
        for line in source:
            if line.rstrip()=='!series_matrix_table_begin':break
        reader=csv.reader(source,delimiter='\t');header=next(reader)
        assert header[0]=='ID_REF' and len(header)==267 and len(set(header[1:]))==266
        assert set(header[1:])=={r['sample_id'] for r in rows}
        indexes=[header.index(s) for s in selected_ids]
        write_tsv(work/'GSE21653_all266_expression_header.tsv',[{'sample_id':s} for s in header[1:]],['sample_id'])
        probes=[]
        with subset.open('wb') as target,gzip.GzipFile(fileobj=target,filename='',mode='wb',mtime=0) as compressed:
            with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as textout:
                writer=csv.writer(textout,delimiter='\t',lineterminator='\n');writer.writerow(['ID_REF']+selected_ids)
                for values in reader:
                    if values and values[0]=='!series_matrix_table_end':break
                    assert len(values)==267
                    probes.append(values[0]);writer.writerow([values[0]]+[values[j] for j in indexes])
        # Read to gzip EOF to validate its trailer/CRC rather than accepting a prefix.
        for _ in source:pass
    assert len(probes)==54613 and len(set(probes))==54613 and '204439_at' in probes
    receipt={
        'dataset':'GSE21653','source_url':SOURCES['GSE21653_series_matrix.txt.gz'][0],
        'full_source_bytes':(work/'GSE21653_series_matrix.txt.gz').stat().st_size,
        'full_source_sha256':sha(work/'GSE21653_series_matrix.txt.gz'),
        'subset_file':subset.name,'subset_bytes':subset.stat().st_size,'subset_sha256':sha(subset),
        'probe_rows':len(probes),'source_sample_columns':266,'selected_sample_columns':87,
        'value_policy':'Exact source decimal strings; all probe rows; selected columns in frozen clinical-list order; no extra transformation',
        'processing_as_recorded':'RMA using a quantile normalization',
        'additional_log_transform':False,'IFI44L_fixed_probe':'204439_at',
        'metadata_sha256':sha(work/'GSE21653_GEO_GSM_brief_SOFT.txt'),
        'all266_metadata_sha256':sha(work/'GSE21653_all266_metadata.tsv'),
        'TNBC87_metadata_sha256':sha(work/'GSE21653_TNBC87_metadata.tsv'),
        'selected_sample_list_sha256':sha(work/'GSE21653_TNBC87_selected_samples.txt'),
        'module_definitions_sha256':sha(module_path),
        'sample_identity_limit':'Unique source case labels are not cross-repository patient identity verification',
        'prepared_utc':datetime.now(timezone.utc).isoformat()
    }
    (work/'GSE21653_frozen_input_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
