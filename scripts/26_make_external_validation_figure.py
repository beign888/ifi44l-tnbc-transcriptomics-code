#!/usr/bin/env python3
"""Plot the full prespecified GSE21653 validation, including imprecise results."""
from pathlib import Path
import argparse, importlib.util, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--project',type=Path,required=True)
    args=parser.parse_args()
    project=args.project.resolve()
    helper_path=Path(__file__).with_name('15_make_bulk_validation_figures_v11.py')
    spec=importlib.util.spec_from_file_location('plot_helpers',helper_path)
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    h.set_style()
    folder=project/'results/validation_GSE21653'
    output=project/'figures/v15_supplementary';output.mkdir(parents=True,exist_ok=True)
    modules=pd.read_csv(folder/'GSE21653_IFI44L_immune_signature_correlations.tsv',sep='\t')
    cox=pd.read_csv(folder/'GSE21653_continuous_IFI44L_Cox.tsv',sep='\t')
    clinical=pd.read_csv(folder/'GSE21653_clinical_IFI44L.tsv',sep='\t')
    survival=clinical.loc[np.isfinite(clinical.dfs_days)&(clinical.dfs_days>0)&clinical.dfs_event.isin([0,1])].copy()
    fig=plt.figure(figsize=(7.2,6.6),constrained_layout=True)
    grid=fig.add_gridspec(2,2,height_ratios=[1,1.15],width_ratios=[1,1])
    a=fig.add_subplot(grid[0,0]);b=fig.add_subplot(grid[0,1])
    nested=grid[1,:].subgridspec(2,1,height_ratios=[4,1],hspace=.03)
    c=fig.add_subplot(nested[0]);risk=fig.add_subplot(nested[1])
    ordered=modules.sort_values('spearman_rho')
    display={'IFN_ALPHA_RESPONSE':'IFN-alpha response','IFN_GAMMA_RESPONSE':'IFN-gamma response',
      'T_CELL_INFLAMED':'T-cell inflamed','T_CELL_EXHAUSTION':'T-cell exhaustion',
      'NK_CELL':'NK cell','CD8_T_CELL':'CD8 T cell','TREG':'Treg',
      'M2_MACROPHAGE':'M2-macrophage module'}
    labels=[display.get(s,s.replace('_',' ').lower()) for s in ordered.signature]
    y=np.arange(len(ordered))
    colours=[h.ORANGE if s=='EPITHELIAL_CONTENT' else h.TEAL for s in ordered.signature]
    a.hlines(y,0,ordered.spearman_rho,color=colours,lw=1.5)
    a.scatter(ordered.spearman_rho,y,color=colours,s=19,zorder=3)
    a.axvline(0,color=h.GREY,lw=.7)
    a.set_yticks(y,labels);a.set_xlim(-.5,1)
    a.set_xlabel('Spearman rho with IFI44L')
    a.set_title(f'GSE21653 expression modules (n={int(modules.n.iloc[0])})',loc='left',fontweight='bold')
    h.panel_label(a,'a',x=-.25)
    order=['Unadjusted','Age-adjusted','Clinicopathological sensitivity']
    selected=cox.set_index('model').loc[order].reset_index()
    labels=['Unadjusted','Age-adjusted','Age, pT, pN and grade']
    y=np.arange(3)[::-1]
    for yi,(_,row),label in zip(y,selected.iterrows(),labels):
        b.hlines(yi,row.ci95_low,row.ci95_high,color=h.BLUE,lw=1.2)
        b.scatter(row.hazard_ratio_per_IFI44L_SD,yi,color=h.BLUE,marker='s',s=25,zorder=3)
    b.axvline(1,color=h.GREY,lw=.8,ls='--')
    b.set_yticks(y,[f'{label}\nn={int(row.n)}, events={int(row.events)}' for label,(_,row) in zip(labels,selected.iterrows())])
    b.set_xscale('log')
    b.set_xlim(min(.6,float(selected.ci95_low.min())*.8),max(1.4,float(selected.ci95_high.max())*1.2))
    ticks=[v for v in [.2,.3,.4,.6,.8,1.,1.2,1.6,2.,3.,4.] if b.get_xlim()[0]<=v<=b.get_xlim()[1]]
    b.set_xticks(ticks,[f'{v:g}' for v in ticks])
    b.set_ylim(-.7,2.5)
    b.set_xlabel('Hazard ratio per IFI44L s.d. (95% CI)')
    b.set_title('Continuous DFS models',loc='left',fontweight='bold')
    h.panel_label(b,'b',x=-.28)
    result=h.draw_km(c,risk,survival,'dfs_days','dfs_event','IFI44L_group','GSE21653 disease-free survival')
    h.panel_label(c,'c',x=-.08)
    base=output/'Supplementary_Figure_S1_GSE21653_validation'
    h.save_bundle(fig,base)
    for index,name in [(1,'KM_confidence'),(2,'KM_censoring'),(3,'KM_risk_counts')]:
        pd.DataFrame(result[index]).to_csv(output/f'Supplementary_Figure_S1_{name}.tsv',sep='\t',index=False)
    audit={'figure':'Supplementary Figure S1','module_patients':int(modules.n.iloc[0]),'DFS_positive_time_patients':len(survival),'DFS_events':int(survival.dfs_event.sum()),'logrank_p':result[0],'models':selected.to_dict(orient='records'),'interpretation':'Observational validation; small event count limits adjusted estimates; zero-time sensitivity remains in complete model tables.'}
    (output/'Supplementary_Figure_S1_build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'figure':str(base),'patients':len(survival),'events':int(survival.dfs_event.sum())}))


if __name__=='__main__':main()
