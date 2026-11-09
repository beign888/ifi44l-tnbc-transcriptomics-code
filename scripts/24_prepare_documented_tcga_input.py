#!/usr/bin/env python3
"""Stage the full documented header and an exact all-gene primary-TNBC subset.

The public reproduction bundle need not redistribute all 1,218 columns: a fixed
full header plus pinned clinical metadata reconstructs eligibility, and the exact
selected expression subset supplies the downstream values without unknown inputs.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,importlib.util,json,shutil
from pathlib import Path

SOURCE_SHA256='263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746'
CLINICAL_SHA256='973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def copy_if_distinct(source,destination):
    if source.resolve()!=destination.resolve():shutil.copyfile(source,destination)


def main():
    p=argparse.ArgumentParser();p.add_argument('--project',type=Path,required=True)
    p.add_argument('--expression-matrix',type=Path,required=True);p.add_argument('--source-receipt',type=Path,required=True)
    p.add_argument('--metadata-file',type=Path,required=True);p.add_argument('--clinical-file',type=Path,required=True)
    a=p.parse_args();project=a.project.resolve();project.mkdir(parents=True,exist_ok=True)
    actual=sha(a.expression_matrix);source_receipt=json.loads(a.source_receipt.read_text());metadata=json.loads(a.metadata_file.read_text())
    if actual!=SOURCE_SHA256 or source_receipt['sha256']!=actual:raise RuntimeError('Full documented matrix hash does not match its pinned acquisition receipt')
    if sha(a.clinical_file)!=CLINICAL_SHA256:raise RuntimeError('Pinned clinical source hash mismatch')
    if metadata.get('version')!='2017-10-13':raise RuntimeError('Unexpected documented dataset version')
    directory=project/'source_metadata/TCGA_documented_primary';directory.mkdir(parents=True,exist_ok=True)
    clinical=project/'source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt';clinical.parent.mkdir(parents=True,exist_ok=True)
    copy_if_distinct(a.clinical_file,clinical)
    copy_if_distinct(a.source_receipt,directory/'REACQUISITION_RECEIPT.json');copy_if_distinct(a.metadata_file,directory/'HiSeqV2.json')
    with gzip.open(a.expression_matrix,'rt') as f:header=f.readline().rstrip('\r\n').split('\t')
    sample_ids=header[1:]
    if len(sample_ids)!=1218 or len(set(sample_ids))!=len(sample_ids):raise RuntimeError('Unexpected or duplicated full documented sample header')
    helper=Path(__file__).with_name('21_audit_tcga_tnbc_definition.py')
    spec=importlib.util.spec_from_file_location('documented_tnbc_definition',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    header_file=directory/'TCGA_HiSeqV2_full_sample_header.tsv'
    with header_file.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source_column_index','sample_id','patient_id','sample_type_code'])
        for i,sample in enumerate(sample_ids,1):w.writerow([i,sample,*module.parse_sample(sample)])
    provenance=dict(schema_version=1,dataset='TCGA.BRCA.sampleMap/HiSeqV2',source_version='2017-10-13',
        source_matrix_file=a.expression_matrix.name,source_matrix_sha256=actual,source_matrix_size_bytes=a.expression_matrix.stat().st_size,
        source_url='https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz',
        units='log2(RSEM normalized count + 1)',additional_transformation='none; use supplied documented log2 matrix',
        metadata_file='HiSeqV2.json',metadata_sha256=sha(directory/'HiSeqV2.json'),
        full_header_file=header_file.name,full_header_sha256=sha(header_file),full_header_samples=len(sample_ids),
        source_gene_rows=None,clinical_source_commit=module.CLINICAL_COMMIT,clinical_source_sha256=CLINICAL_SHA256,
        subset_file='TCGA_BRCA_documented_primary_subset.tsv.gz',subset_sha256=None,
        subset_value_policy='original decimal strings, original gene order and full-header sample order; no expression transformation',
        one_sample_per_patient_rule='lexicographically first complete type01 source barcode within patient; no expression/outcome tie-break')
    receipt=directory/'TCGA_documented_input_provenance.json';receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    rows=module.reconcile_clinical(project,False,a.expression_matrix,receipt)
    selected=[r['sample_id'] for r in rows if r['primary_TNBC_eligible']];positions=[header.index(s) for s in selected]
    subset=directory/provenance['subset_file'];genes=set();gene_rows=0
    with subset.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as out:
            out.write(('\t'.join([header[0]]+selected)+'\n').encode())
            with gzip.open(a.expression_matrix,'rt') as f:
                next(f)
                for line in f:
                    cells=line.rstrip('\r\n').split('\t')
                    if len(cells)!=len(header):raise RuntimeError('Malformed full-matrix row')
                    genes.add(cells[0]);gene_rows+=1
                    out.write(('\t'.join([cells[0]]+[cells[i] for i in positions])+'\n').encode())
    if gene_rows!=20530 or len(genes)!=gene_rows:raise RuntimeError('Unexpected or duplicated input genes')
    provenance.update(source_gene_rows=gene_rows,subset_gene_rows=gene_rows,subset_sample_columns=len(selected),
        subset_sha256=sha(subset),subset_size_bytes=subset.stat().st_size,primary_cohort_n=len(selected),
        strict_IHC_cohort_n=sum(r['strict_IHC_negative_eligible'] for r in rows))
    receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(status='PASS_DOCUMENTED_FULL_HEADER_AND_EXACT_SUBSET',project=str(project),primary_n=len(selected),strict_n=provenance['strict_IHC_cohort_n'],subset_sha256=provenance['subset_sha256'],subset_bytes=provenance['subset_size_bytes'])))


if __name__=='__main__':main()
