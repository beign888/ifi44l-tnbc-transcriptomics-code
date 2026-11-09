#!/usr/bin/env python3
"""Run the version16 extensions after the frozen core public-data workflow."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, datetime

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--plan-dir',type=Path)
    parser.add_argument('--stage',choices=['context','scanb','all'],default='all')
    args=parser.parse_args();p=args.project.resolve()
    plans=args.plan_dir.resolve() if args.plan_dir else (p/'protocols/v16' if (p/'protocols/v16').is_dir() else p/'audit_v16')
    jobs=[];rscript=shutil.which('Rscript')
    if not rscript:raise SystemExit('Rscript is required; see the recorded R environment.')
    if args.stage in ['context','all']:
        plan=plans/'ANALYSIS_PLAN_V16.md'
        amendment=plans/'ANALYSIS_PLAN_V16_AMENDMENT1.md'
        if sha(plan)!='8d24b3a4971f981ad4651d7f332e229da8ec2e25b9695af67e3b8961d947c97e':raise SystemExit('Context protocol hash mismatch')
        if sha(amendment)!='3dbc05693df636aef345ff557fa3b91fb5c89a29b26167ea2f3438a5a46d7ad2':raise SystemExit('Context clarification hash mismatch')
        jobs.extend([
          ('30_conditional_context',[rscript,str(p/'scripts/30_IFI44L_IFN_conditional_context.R'),'--project',str(p),'--analysis-plan',str(plan),'--expected-plan-sha256',sha(plan),'--analysis-plan-amendment',str(amendment),'--expected-amendment-sha256',sha(amendment)]),
          ('32_single_cell_context',[sys.executable,str(p/'scripts/32_single_cell_within_lineage_IFN_context.py'),'--project',str(p),'--plan-dir',str(plans)]),
          ('31_conditional_figure',[sys.executable,str(p/'scripts/31_make_IFI44L_IFN_context_figure.py'),'--project',str(p)])])
    if args.stage in ['scanb','all']:
        plan=plans/'GSE96058_ANALYSIS_PLAN.md'
        amendment=plans/'GSE96058_PLATFORM_AMENDMENT.md'
        if sha(plan)!='f174e513f1c3091a19d0ab3233512048c64e0b02f564a46101fd324ab550a7ae':raise SystemExit('SCAN-B protocol hash mismatch')
        if sha(amendment)!='06b06ed2d8879e2b1aecc2df094eca9baabcd317c1686d55c0801d2029f8db22':raise SystemExit('SCAN-B platform amendment hash mismatch')
        scanb_input=p/'source_metadata/GSE96058/GSE96058_TNBC143_frozen_all_genes.tsv.gz'
        if not scanb_input.is_file():
            jobs.append(('34_SCANB_offline_restore',[sys.executable,str(p/'scripts/34_prepare_GSE96058_inputs.py'),'--project',str(p),'--restore-parts']))
        jobs.extend([
          ('35_SCANB_validation',[rscript,str(p/'scripts/35_GSE96058_validation.R'),'--project',str(p),'--analysis-plan',str(plan),'--expected-plan-sha256',sha(plan),'--analysis-plan-amendment',str(amendment),'--expected-amendment-sha256',sha(amendment)]),
          ('36_SCANB_figure',[sys.executable,str(p/'scripts/36_make_GSE96058_validation_figure.py'),'--project',str(p)])])
    for _,command in jobs:
        if not Path(command[1]).is_file():raise SystemExit('Missing extension script: '+command[1])
    logs=p/'logs/v16_extensions';logs.mkdir(parents=True,exist_ok=True)
    receipt=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),project=str(p),stage=args.stage,python=sys.executable,protocol_directory=str(plans),jobs=[],status='RUNNING')
    target=logs/'execution_receipt.json'
    target.write_text(json.dumps(receipt,indent=2)+'\n')
    for name,command in jobs:
        log=logs/(name+'.log');print('Running '+name,flush=True)
        with log.open('w') as f:result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,cwd=p)
        receipt['jobs'].append(dict(name=name,command=command,script_sha256=sha(Path(command[1])),log=log.name,log_sha256=sha(log),exit_code=result.returncode))
        target.write_text(json.dumps(receipt,indent=2)+'\n')
        if result.returncode:
            receipt['status']='FAILED';target.write_text(json.dumps(receipt,indent=2)+'\n');raise SystemExit('Failed '+name+'; see '+str(log))
    receipt.update(status='PASS_ACTUAL_EXTENSION_EXECUTION',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),limits='Core analyses are prerequisites and are not rerun by this extension-only wrapper; input/model/visual audits remain separate.')
    target.write_text(json.dumps(receipt,indent=2)+'\n');print(receipt['status'],flush=True)

if __name__=='__main__':main()
