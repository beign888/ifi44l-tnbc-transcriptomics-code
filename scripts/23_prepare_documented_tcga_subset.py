#!/usr/bin/env python3
"""Preserve exact selected source strings and baseline values for standalone use."""
from pathlib import Path
from urllib.request import urlopen
import argparse,csv,gzip,hashlib,io,json,shutil

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda:stream.read(8*1024*1024),b''):h.update(part)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--project',type=Path,required=True);p.add_argument('--work-dir',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args()
    bulk=a.project/'results/bulk';w=a.work_dir
    source=w/'TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz'
    assert sha(source)=='263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746'
    copies={'TCGA_TNBC_selected_samples.txt':'baseline_TCGA151_selected_samples.txt','immune_signature_definitions.tsv':'baseline_module_definitions.tsv','tnbc_sample_metadata.tsv':'baseline_TCGA_IFI44L.tsv','immune_signature_ssGSEA_scores_wide.tsv':'baseline_TCGA_ssGSEA_scores_wide.tsv'}
    for old,new in copies.items():shutil.copyfile(bulk/old,w/new)
    selected=(w/'baseline_TCGA151_selected_samples.txt').read_text().splitlines()
    assert len(selected)==151 and len(set(selected))==151 and len({s[:12] for s in selected})==151
    output=w/'documented_HiSeqV2_matched151_expression.tsv.gz'
    with gzip.open(source,'rt',newline='') as stream,output.open('wb') as target:
        reader=csv.reader(stream,delimiter='\t');header=next(reader)
        # All 151 match exactly. Do not manufacture aliases or average samples.
        indexes=[header.index(sample) for sample in selected]
        with gzip.GzipFile(fileobj=target,filename='',mode='wb',mtime=0) as compressed:
            with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as text:
                writer=csv.writer(text,delimiter='\t',lineterminator='\n');writer.writerow([header[0]]+selected)
                genes=[]
                for row in reader:
                    genes.append(row[0]);writer.writerow([row[0]]+[row[index] for index in indexes])
    receipt=dict(source_file=source.name,source_bytes=source.stat().st_size,source_sha256=sha(source),source_url='https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz',subset_file=output.name,subset_bytes=output.stat().st_size,subset_sha256=sha(output),gene_rows=len(genes),sample_columns=len(selected),duplicate_gene_rows=len(genes)-len(set(genes)),sample_order='frozen clinical151 list order',value_policy='source decimal strings preserved exactly; no transformation or filtering at subset preparation',baseline_files={n:dict(bytes=(w/n).stat().st_size,sha256=sha(w/n)) for n in copies.values()},unit='log2(norm_count+1), documented RSEM normalized counts',source_version='2017-10-13',independent_patient_validation=False,original_archived_expression_provenance_recovered=False)
    assert output.stat().st_size<18_000_000
    (w/'SELFCONTAINED_SUBSET_PROVENANCE.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n');print(json.dumps(receipt,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
