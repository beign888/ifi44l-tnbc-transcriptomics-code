#!/usr/bin/env python3
"""Plot saved script30 results only; no model fitting or result modification."""
from __future__ import annotations
import argparse,hashlib,json,datetime,platform
from pathlib import Path
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
import numpy as np
import pandas as pd

BLUE='#2B6CA3';TEAL='#23877D';ORANGE='#D28A35';GREY='#747474';BLACK='#222222';LIGHT='#D8D8D8'
COHORTS=['TCGA','GSE58812','GSE21653']
PRIMARY=[('TCGA','OS'),('TCGA','PFI'),('GSE58812','MFS'),('GSE58812','OS'),('GSE21653','DFS')]
SECONDARY=[('TCGA','OS'),('TCGA','PFI'),('GSE21653','DFS')]
GENES=['IFI44L','IFI44','MX1','ISG15']
BASE='Supplementary_Figure_S2_IFI44L_IFN_conditional_context'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def set_style():
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],
      'font.size':7.6,'axes.labelsize':7.6,'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,
      'legend.fontsize':7.2,'axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False,
      'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def unique_row(frame,**keys):
    selected=frame
    for key,value in keys.items():selected=selected[selected[key]==value]
    if len(selected)!=1:raise ValueError(f'Expected one source row: {keys}, found {len(selected)}')
    return selected.iloc[0]
def draw_panel_letter(ax,label):ax.text(-.12,1.10,label,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def fmt_q(value):return f'{value:.3f}' if value>=.001 else f'{value:.1e}'
def save_tsv(frame,path):frame.to_csv(path,sep='\t',index=False,float_format='%.17g',lineterminator='\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--results-dir',type=Path);parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args();project=args.project.resolve();source=(args.results_dir or project/'results/conditional_context_v16').resolve()
    out=(args.output_dir or project/'figures/v16_supplementary').resolve();out.mkdir(parents=True,exist_ok=True)
    protected={p.name:sha(p) for p in source.iterdir() if p.is_file()}
    inputs={}
    def read(name):
        path=source/name;inputs[name]=sha(path);return pd.read_csv(path,sep='\t')
    metadata={c:read(c+'_full_expression_cohort_inputs.tsv') for c in COHORTS}
    coefficients=read('Cox_all_coefficients.tsv');summary=read('Cox_model_summary.tsv');ph=read('Cox_PH_diagnostics.tsv')
    primary_lrt=read('primary5_conditional_LRT.tsv');secondary_lrt=read('secondary3_conditional_LRT.tsv')
    diagnostic=read('predictor_correlation_diagnostics.tsv');redundancy=read('fixed_expression_redundancy10.tsv')
    scatter=[]
    for panel,cohort in zip('abc',COHORTS):
        frame=metadata[cohort]
        if not frame.sample_id.is_unique:raise ValueError('Scatter observations must be unique cohort samples')
        d=unique_row(diagnostic,cohort=cohort,endpoint='full_expression_cohort',variable1='IFI44L_z',variable2='IFNa19_z')
        if len(frame)!=d['n']:raise ValueError('Saved scatter/diagnostic n mismatch')
        points=frame[['sample_id','IFI44L_z','IFNa19_z']].copy();points.insert(0,'cohort',cohort);points.insert(0,'panel',panel)
        points['displayed_spearman_rho']=d.spearman_rho;points['score_definition']='original IFI44L-excluded19-gene IFNa score; standardized in full expression cohort'
        scatter.append(points)
    scatter=pd.concat(scatter,ignore_index=True)

    def forest_sources(family,order,models,panel):
        result=[];lrt_table=primary_lrt if family=='primary' else secondary_lrt
        for row_order,(cohort,endpoint) in enumerate(order):
            test=unique_row(lrt_table,cohort=cohort,endpoint=endpoint)
            for model in models:
                c=unique_row(coefficients,cohort=cohort,endpoint=endpoint,family=family,model=model,term='IFI44L_z')
                s=unique_row(summary,cohort=cohort,endpoint=endpoint,family=family,model=model)
                p=unique_row(ph,cohort=cohort,endpoint=endpoint,family=family,model=model,term='GLOBAL')
                result.append(dict(panel=panel,row_order=row_order,cohort=cohort,endpoint=endpoint,family=family,model=model,
                    n=int(c['n']),events=int(c.events),coefficient_count=int(s.coefficient_count),events_per_coefficient=s.events_per_coefficient,
                    hazard_ratio=c.hazard_ratio,ci95_low=c.ci95_low,ci95_high=c.ci95_high,nominal_Wald_p=c.nominal_p_value,
                    conditional_LRT_p=test.p_value,conditional_LRT_fdr_bh=test.fdr_bh,conditional_LRT_BH_n=int(test.BH_family_size),
                    nominal_global_PH_p=p.nominal_p_value,nominal_global_PH_flag=bool(p.nominal_p_value<.05),
                    effect_scale='IFI44L HR per full-expression-cohort SD; horizontal bars95% CI'))
        return pd.DataFrame(result)
    primary=forest_sources('primary',PRIMARY,['age_plus_IFI44L','age_plus_IFNa_plus_IFI44L'],'d')
    secondary=forest_sources('secondary',SECONDARY,['clinical_base_plus_IFNa_plus_IFI44L'],'f')
    limits=pd.concat([primary,secondary]);xlow=float(limits.ci95_low.min()*.72);xhigh=float(max(limits.ci95_high.max()*1.20,1.5))
    if not (xlow>0 and np.isfinite(xhigh)):raise ValueError('Invalid/non-finite forest CI: cannot silently truncate')
    pairs=[(a,b) for i,a in enumerate(GENES) for b in GENES[i+1:]]+[(a,'IFNa16_excluding_IFI44L_IFI44_MX1_ISG15') for a in GENES]
    heat_rows=[]
    for i,(a,b) in enumerate(pairs):
        for cohort in COHORTS:
            r=unique_row(redundancy,cohort=cohort,variable1=a,variable2=b)
            heat_rows.append(dict(panel='e',row_order=i,column_order=COHORTS.index(cohort),cohort=cohort,variable1=a,variable2=b,
                n=int(r['n']),spearman_rho=r.spearman_rho,p_value=r.p_value,fdr_bh=r.fdr_bh,BH_family_size=int(r.BH_family_size),
                label=a+'–'+('IFNa16' if b.startswith('IFNa16') else b)))
    heat=pd.DataFrame(heat_rows)
    if len(heat)!=30 or len(primary)!=10 or len(secondary)!=3:raise ValueError('Planned panels are incomplete')
    for label,frame in [('scatter',scatter),('primary_forest',primary),('redundancy_heatmap',heat),('secondary_forest',secondary)]:
        save_tsv(frame,out/f'{BASE}_{label}_source.tsv')

    set_style();fig=plt.figure(figsize=(9.1,9.5))
    outer=fig.add_gridspec(3,1,left=.075,right=.975,bottom=.115,top=.90,height_ratios=[1.10,1.86,1.02],hspace=.63)
    top=outer[0].subgridspec(1,3,wspace=.40)
    xvals=scatter.IFI44L_z.to_numpy();yvals=scatter.IFNa19_z.to_numpy()
    xpad=.10*(xvals.max()-xvals.min());ypad=.10*(yvals.max()-yvals.min())
    for col,(panel,cohort,colour) in enumerate(zip('abc',COHORTS,[BLUE,TEAL,ORANGE])):
        ax=fig.add_subplot(top[0,col]);data=scatter[scatter.cohort==cohort]
        ax.scatter(data.IFI44L_z,data.IFNa19_z,s=12,color=colour,alpha=.67,edgecolors='none')
        ax.axhline(0,lw=.45,color=LIGHT,zorder=0);ax.axvline(0,lw=.45,color=LIGHT,zorder=0)
        ax.set_xlim(xvals.min()-xpad,xvals.max()+xpad);ax.set_ylim(yvals.min()-ypad,yvals.max()+ypad)
        ax.set_xlabel('IFI44L expression (SD)');ax.set_ylabel('Original IFNa19 score (SD)')
        ax.set_title(cohort,fontweight='bold',pad=8)
        ax.text(.04,.96,f"n = {len(data)}; Spearman ρ = {data.displayed_spearman_rho.iloc[0]:.3f}",transform=ax.transAxes,ha='left',va='top',fontsize=7.4,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.86,pad=2))
        draw_panel_letter(ax,panel)

    middle=outer[1].subgridspec(1,2,width_ratios=[1.30,1],wspace=.62)
    forest_grid=middle[0,0].subgridspec(1,3,width_ratios=[1.0,1.3,1.22],wspace=.06)
    label_ax=fig.add_subplot(forest_grid[0,0]);axd=fig.add_subplot(forest_grid[0,1]);value_ax=fig.add_subplot(forest_grid[0,2])
    for ax in [label_ax,value_ax]:ax.set_axis_off();ax.set_ylim(-.65,4.65)
    axd.set_ylim(-.65,4.65);axd.set_xscale('log');axd.set_xlim(xlow,xhigh)
    axd.axvline(1,color=GREY,lw=.65,linestyle='--');axd.set_yticks([])
    ticks=[t for t in [.05,.1,.2,.5,1,2] if xlow<=t<=xhigh]
    axd.xaxis.set_major_locator(FixedLocator(ticks));axd.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'));axd.xaxis.set_minor_locator(NullLocator())
    axd.set_xlabel('IFI44L HR (log scale)')
    label_ax.text(0,1.13,'d',transform=label_ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
    label_ax.text(0,1.06,'Primary association models',transform=label_ax.transAxes,fontweight='bold',fontsize=8.5,va='bottom')
    value_ax.text(0,1.01,'HR [95% CI]',transform=value_ax.transAxes,fontsize=7.4,va='bottom')
    legend=[Line2D([0],[0],marker='o',color=GREY,lw=.8,markersize=4,label='Age + IFI44L'),
      Line2D([0],[0],marker='s',color=BLUE,lw=.8,markersize=4,label='Age + IFNa19 + IFI44L')]
    label_ax.legend(handles=legend,loc='lower left',bbox_to_anchor=(0,.985),ncol=2,columnspacing=.9,handlelength=1.4,handletextpad=.4,borderaxespad=0)
    for index,(cohort,endpoint) in enumerate(PRIMARY):
        group=primary[(primary.cohort==cohort)&(primary.endpoint==endpoint)];y=4-index
        flag=bool(group.nominal_global_PH_flag.any());first=group.iloc[0]
        label_ax.text(0,y,f'{cohort} {endpoint}'+('†' if flag else '')+f'\nn={first.n}; events={first.events}\nLRT q={fmt_q(first.conditional_LRT_fdr_bh)}',va='center',fontsize=7.3,linespacing=1.25)
        for model,dy,colour,marker in [('age_plus_IFI44L',.15,GREY,'o'),('age_plus_IFNa_plus_IFI44L',-.15,BLUE,'s')]:
            r=group[group.model==model].iloc[0]
            axd.errorbar(r.hazard_ratio,y+dy,xerr=[[r.hazard_ratio-r.ci95_low],[r.ci95_high-r.hazard_ratio]],fmt=marker,
                color=colour,markersize=4,elinewidth=.9,capsize=2,capthick=.8)
            value_ax.text(0,y+dy,f'{r.hazard_ratio:.2f} [{r.ci95_low:.2f}–{r.ci95_high:.2f}]',va='center',fontsize=7.3,color=colour)

    axe=fig.add_subplot(middle[0,1]);matrix=heat.pivot(index='row_order',columns='cohort',values='spearman_rho')[COHORTS].to_numpy()
    cmap=LinearSegmentedColormap.from_list('ifi_context',['#F4F7FA','#A9C8E2',BLUE])
    image=axe.imshow(matrix,vmin=0,vmax=1,cmap=cmap,aspect='auto')
    axe.set_yticks(range(10));axe.set_yticklabels([heat[heat.row_order==i].label.iloc[0] for i in range(10)],fontsize=7.3)
    axe.set_xticks(range(3));axe.set_xticklabels([f'{c}\nn={int(heat[heat.cohort==c].n.iloc[0])}' for c in COHORTS],fontsize=7.2)
    axe.tick_params(length=0);axe.axhline(5.5,color='white',lw=2)
    for row in range(10):
        for col in range(3):axe.text(col,row,f'{matrix[row,col]:.3f}',ha='center',va='center',fontsize=7.4,color='white' if matrix[row,col]>.78 else BLACK)
    axe.set_title('Fixed expression comparisons',fontweight='bold',pad=14)
    draw_panel_letter(axe,'e')
    colourbar=fig.colorbar(image,ax=axe,orientation='horizontal',fraction=.045,pad=.13,aspect=22,ticks=[0,.5,1])
    colourbar.set_label('Spearman ρ (fixed 10 tests per cohort)',fontsize=7.2,labelpad=3);colourbar.ax.tick_params(labelsize=7.2,length=2)

    bottom=outer[2].subgridspec(1,3,width_ratios=[1.7,1.9,1.1],wspace=.08)
    labelf=fig.add_subplot(bottom[0,0]);axf=fig.add_subplot(bottom[0,1]);valuesf=fig.add_subplot(bottom[0,2])
    for ax in [labelf,valuesf]:ax.set_axis_off();ax.set_ylim(-.5,2.5)
    axf.set_ylim(-.5,2.5);axf.set_xscale('log');axf.set_xlim(xlow,xhigh);axf.set_yticks([])
    axf.axvline(1,color=GREY,lw=.65,linestyle='--');axf.xaxis.set_major_locator(FixedLocator(ticks));axf.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'));axf.xaxis.set_minor_locator(NullLocator())
    axf.set_xlabel('Conditional IFI44L HR (95% CI; log scale)')
    labelf.text(0,1.14,'f',transform=labelf.transAxes,fontweight='bold',fontsize=10,va='bottom')
    labelf.text(0,1.04,'Secondary: age + IFNa19 + IFI44L + clinical context',transform=labelf.transAxes,fontweight='bold',fontsize=8.5,va='bottom')
    valuesf.text(0,1.01,'HR [95% CI]; LRT q',transform=valuesf.transAxes,fontsize=7.4,va='bottom')
    for index,(cohort,endpoint) in enumerate(SECONDARY):
        r=secondary[(secondary.cohort==cohort)&(secondary.endpoint==endpoint)].iloc[0];y=2-index
        note='stage-stratified' if cohort=='TCGA' else 'pT, pN and grade adjusted'
        flag='†' if r.nominal_global_PH_flag else ''
        labelf.text(0,y,f'{cohort} {endpoint}{flag} · {note}\nn={r.n}; events={r.events}; coefficients={r.coefficient_count}',va='center',fontsize=7.5,linespacing=1.3)
        axf.errorbar(r.hazard_ratio,y,xerr=[[r.hazard_ratio-r.ci95_low],[r.ci95_high-r.hazard_ratio]],fmt='s',color=TEAL,markersize=4.3,elinewidth=1,capsize=2,capthick=.8)
        valuesf.text(0,y,f'{r.hazard_ratio:.2f} [{r.ci95_low:.2f}–{r.ci95_high:.2f}]\nq={fmt_q(r.conditional_LRT_fdr_bh)}',va='center',fontsize=7.4,color=TEAL,linespacing=1.3)

    fig.text(.075,.985,'Supplementary Figure S2 | IFI44L in the broader interferon-expression context',fontweight='bold',fontsize=10,va='top')
    fig.text(.075,.953,'Independent cohort samples; existing IFNa19 scores for Cox; fixed comparator-excluded IFNa16 score for panel e.',fontsize=7.5,va='top')
    primary_flags=primary[(primary.model=='age_plus_IFNa_plus_IFI44L')&primary.nominal_global_PH_flag]
    secondary_flags=secondary[secondary.nominal_global_PH_flag]
    flag_text='; '.join([f'{r.cohort} {r.endpoint} '+('primary' if r.family=='primary' else 'secondary')+f' joint model (P={r.nominal_global_PH_p:.3f})' for r in pd.concat([primary_flags,secondary_flags]).itertuples()])
    fig.text(.075,.068,'IFNa16 excludes IFI44L, IFI44, MX1 and ISG15. All six gene pairs and four gene–module correlations are retained.',fontsize=7.3,va='top')
    fig.text(.075,.047,'† Nominal global PH flag: '+flag_text+'.',fontsize=7.3,va='top')
    fig.text(.075,.026,'Conditional comparisons are exploratory association checks, not validated prediction. Sparse events and correlated predictors limit interpretation.',fontsize=7.3,va='top')
    fig.canvas.draw()
    font_sizes=[t.get_fontsize() for t in fig.findobj(mpl.text.Text) if t.get_text()]
    if min(font_sizes)<7:raise ValueError('A nonempty plot label is smaller than7pt')
    for ax in [axd,axf]:
        if ax.get_xlim()[0]>=limits.ci95_low.min() or ax.get_xlim()[1]<=limits.ci95_high.max():raise ValueError('Forest CI clipping guard failed')
    for suffix in ['pdf','svg','png']:
        fig.savefig(out/f'{BASE}.{suffix}',dpi=300,bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
    changed=[name for name,digest in protected.items() if sha(source/name)!=digest]
    if changed:raise RuntimeError('Script30 result changed while plotting: '+str(changed))
    receipt=dict(status='FIGURE_GENERATED_FROM_SAVED_RESULTS_REQUIRES_VISUAL_QA',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        data_source=str(source),input_sha256=inputs,script_sha256=sha(Path(__file__).resolve()),script30_result_files_checked_unchanged=len(protected),result_modifications=changed,
        source_row_counts=dict(scatter=len(scatter),primary_forest=len(primary),redundancy_heatmap=len(heat),secondary_forest=len(secondary)),
        figure_size_inches=[9.1,9.5],minimum_nonempty_text_font_pt=min(font_sizes),pdf_fonttype=42,svg_fonttype='none',PNG_dpi=300,
        forest_common_log_xlim=[xlow,xhigh],full_CI_range=[float(limits.ci95_low.min()),float(limits.ci95_high.max())],no_forest_CI_truncated=True,
        nominal_global_PH_flags=pd.concat([primary_flags,secondary_flags])[['panel','cohort','endpoint','family','nominal_global_PH_p']].to_dict('records'),
        no_models_fitted=True,correlation_values_read_from_script30=True,
        software=dict(python=platform.python_version(),matplotlib=mpl.__version__,pandas=pd.__version__,numpy=np.__version__),
        output_files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.glob(BASE+'*') if p.is_file() and p.name!=BASE+'_build_receipt.json'})
    (out/f'{BASE}_build_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','source_row_counts','minimum_nonempty_text_font_pt','forest_common_log_xlim','nominal_global_PH_flags']}))

if __name__=='__main__':main()
