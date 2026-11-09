#!/usr/bin/env python3
"""Define TNBC from the complete, documented TCGA expression-column universe.

Only receptor/sample-type metadata are evaluated. The archived matrix and
historical eligibility lists are never read by this version.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,re
from collections import defaultdict
from pathlib import Path

CLINICAL_SHA256='973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc'
CLINICAL_COMMIT='4c975725c54df7e31e0610fad22a742ad01cdfa3'
SOURCE_SHA256='263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def write_tsv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def parse_sample(sample):
    match=re.fullmatch(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-(\d{2})([A-Z])?(?:-.*)?',sample)
    return (match[1],match[2]) if match else ('','unparsed')


def full_samples(project,expression_matrix=None,source_receipt=None):
    directory=project/'source_metadata/TCGA_documented_primary'
    receipt_path=source_receipt or directory/'TCGA_documented_input_provenance.json'
    provenance=json.loads(receipt_path.read_text())
    if provenance['source_matrix_sha256']!=SOURCE_SHA256:
        raise RuntimeError('Documented source receipt has an unexpected matrix SHA-256')
    header=directory/provenance['full_header_file']
    if sha(header)!=provenance['full_header_sha256']:
        raise RuntimeError('Frozen full-expression header SHA-256 mismatch')
    with header.open() as f:samples=[r['sample_id'] for r in csv.DictReader(f,delimiter='\t')]
    if len(samples)!=provenance['full_header_samples'] or len(set(samples))!=len(samples):
        raise RuntimeError('Full expression header is incomplete or duplicated')
    if expression_matrix is None:
        subset_path=directory/provenance['subset_file']
        full_path=directory/provenance['source_matrix_file']
        expression_matrix=subset_path if subset_path.is_file() else full_path
    if expression_matrix.is_file():
        observed=sha(expression_matrix)
        if observed==SOURCE_SHA256:
            opener=gzip.open if expression_matrix.suffix=='.gz' else open
            with opener(expression_matrix,'rt') as f:actual=f.readline().rstrip('\r\n').split('\t')[1:]
            if actual!=samples:raise RuntimeError('Supplied full-matrix columns differ from frozen full header')
        elif observed==provenance.get('subset_sha256'):
            # The subset may contain primary and strict samples, but never defines
            # the receptor-screening universe. The full header remains authoritative.
            pass
        else:raise RuntimeError('Supplied expression matrix is neither the pinned full matrix nor its exact supplied subset')
    else:
        raise FileNotFoundError('Stage the documented full matrix or exact subset; archived matrices are not a fallback')
    return samples,provenance


def reconcile_clinical(project,strict_ihc=False,expression_matrix=None,source_receipt=None,clinical_file=None):
    clinical=clinical_file or project/'source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt'
    if sha(clinical)!=CLINICAL_SHA256:raise RuntimeError('Pinned clinical source SHA-256 mismatch')
    with clinical.open() as f:records=list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))
    ids=[r['PATIENT_ID'] for r in records]
    if len(ids)!=len(set(ids)):raise RuntimeError('Duplicated clinical patient identifiers')
    patients={r['PATIENT_ID']:r for r in records}
    samples,provenance=full_samples(project,expression_matrix,source_receipt)
    by_patient=defaultdict(list)
    for sample in samples:
        patient,kind=parse_sample(sample)
        if kind=='01':by_patient[patient].append(sample)
    chosen={patient:sorted(barcodes)[0] for patient,barcodes in by_patient.items()}
    rows=[]
    for index,sample in enumerate(samples):
        patient,kind=parse_sample(sample);r=patients.get(patient,{})
        er=r.get('ER_STATUS_BY_IHC','unavailable');pr=r.get('PR_STATUS_BY_IHC','unavailable')
        ihc=r.get('IHC_HER2','unavailable');fish=r.get('HER2_FISH_STATUS','unavailable');score=r.get('HER2_IHC_SCORE','unavailable')
        hp=er=='Negative' and pr=='Negative'
        positive=ihc=='Positive' or fish=='Positive' or score.strip() in {'3','3+','3.0'}
        negative=ihc=='Negative' or fish=='Negative'
        selected_primary=kind=='01' and chosen.get(patient)==sample
        main=selected_primary and hp and not positive and negative
        strict=selected_primary and hp and not positive and ihc=='Negative'
        included=strict if strict_ihc else main
        reason=('included' if included else 'not_primary_tumour_sample_type_01' if kind!='01'
                else 'additional_primary_sample_same_patient' if not selected_primary
                else 'clinical_patient_not_available' if not r
                else 'ER_or_PR_not_explicitly_negative' if not hp
                else 'HER2_positive_or_IHC_3plus' if positive
                else 'HER2_not_explicitly_IHC_negative_for_strict' if strict_ihc and ihc!='Negative'
                else 'HER2_not_resolved_negative')
        rows.append(dict(source_column_index=index+1,sample_id=sample,patient_id=patient,sample_type_code=kind,
            is_primary_tumour=kind=='01',clinical_patient_matched=bool(r),primary_samples_for_patient=len(by_patient.get(patient,[])),
            deterministic_selected_primary_sample=chosen.get(patient,''),selected_primary_sample=selected_primary,
            ER_status_by_IHC=er,PR_status_by_IHC=pr,HER2_IHC_status=ihc,HER2_FISH_status=fish,HER2_IHC_score=score,
            ER_negative=er=='Negative',PR_negative=pr=='Negative',HER2_negative=not positive and negative,
            primary_TNBC_eligible=main,strict_IHC_negative_eligible=strict,included_in_analysis=included,
            exclusion_reason=reason,present_in_source_expression_matrix=True))
    output=project/'results/bulk';write_tsv(output/'TCGA_TNBC_sample_definition_manifest.tsv',rows,list(rows[0]))
    selected=[r for r in rows if r['included_in_analysis']]
    if len(selected)!=len({r['patient_id'] for r in selected}):raise RuntimeError('Duplicate patient entered final cohort')
    values=dict(definition='complete_documented_expression_universe_primary01_clinical_ER_PR_negative_resolved_HER2_negative',
        expression_source='UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2',expression_source_version='2017-10-13',expression_source_sha256=SOURCE_SHA256,
        expression_matrix_samples_n=len(samples),primary_tumour_columns_n=sum(r['is_primary_tumour'] for r in rows),
        unique_primary_patients_n=len(by_patient),primary_patients_with_multiple_samples_n=sum(len(v)>1 for v in by_patient.values()),
        clinical_primary_patients_matched_n=sum(p in patients for p in by_patient),
        clinical_source_commit=CLINICAL_COMMIT,clinical_source_sha256=CLINICAL_SHA256,
        ER_PR_both_negative_primary_n=sum(r['selected_primary_sample'] and r['ER_negative'] and r['PR_negative'] for r in rows),
        primary_TNBC_eligible_n=sum(r['primary_TNBC_eligible'] for r in rows),strict_IHC_negative_n=sum(r['strict_IHC_negative_eligible'] for r in rows),
        TNBC_matched_to_expression_matrix_n=len(selected),strict_IHC_sensitivity_mode=strict_ihc,PAM50_used_for_inclusion=False,
        duplicate_primary_sample_rule='lexicographically first source barcode among type01 samples for each patient; outcome/expression independent')
    write_tsv(output/'TCGA_TNBC_definition_audit.tsv',[dict(metric=k,value=v) for k,v in values.items()],['metric','value'])
    (output/'TCGA_TNBC_selected_samples.txt').write_text('\n'.join(r['sample_id'] for r in selected)+'\n')
    print(json.dumps(dict(eligible_n=len(selected),primary_TNBC_n=values['primary_TNBC_eligible_n'],strict_IHC_n=values['strict_IHC_negative_n'],full_expression_columns=len(samples),selection_uses_outcome=False)))
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--expression-matrix',type=Path);p.add_argument('--source-receipt',type=Path);p.add_argument('--clinical-file',type=Path)
    p.add_argument('--strict-ihc',action='store_true');a=p.parse_args()
    reconcile_clinical(a.project.resolve(),a.strict_ihc,a.expression_matrix,a.source_receipt,a.clinical_file)


if __name__=='__main__':main()
