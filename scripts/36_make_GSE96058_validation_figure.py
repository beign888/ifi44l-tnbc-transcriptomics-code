#!/usr/bin/env python3
"""Render saved SCAN-B results only, retaining all tests and exact plot sources."""
from pathlib import Path
import argparse,datetime,hashlib,json,platform
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
import numpy as np
import pandas as pd

BASE='Supplementary_Figure_S4_GSE96058_validation'
BLUE='#2B6CA3';TEAL='#23877D';ORANGE='#D28A35';GREY='#747474';BLACK='#222222'
MODULES=[('IFN_ALPHA_RESPONSE','IFN-alpha response'),('IFN_GAMMA_RESPONSE','IFN-gamma response'),
 ('ANTIGEN_PRESENTATION','Antigen presentation'),('CHECKPOINT_AXIS','Checkpoint axis'),
 ('T_CELL_INFLAMED','T-cell inflamed'),('CYTOLYTIC_ACTIVITY','Cytolytic activity'),('CD8_T_CELL','CD8 T cell'),
 ('NK_CELL','NK cell'),('T_CELL_EXHAUSTION','T-cell exhaustion'),('TREG','Regulatory T cell'),
 ('M2_MACROPHAGE','M2 macrophage'),('DENDRITIC_CELL','Dendritic cell'),('EPITHELIAL_CONTENT','Epithelial content')]
MODEL_ORDER=['M1','M2','M3','M5','M7','M9']
MODEL_LABELS={'M1':'Unadjusted IFI44L','M2':'Age + IFI44L','M3':'Clinical + IFI44L',
 'M5':'Age + IFNa19 + IFI44L','M7':'Clinical + IFNa19 + IFI44L','M9':'Age + IFNa19 + IFI44L; platform-stratified'}

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def one(d,**keys):
    for k,v in keys.items():d=d[d[k]==v]
    if len(d)!=1:raise ValueError(f'Nonunique/missing source key {keys}: {len(d)} rows')
    return d.iloc[0]

def save_tsv(d,p):d.to_csv(p,sep='\t',index=False,float_format='%.17g',lineterminator='\n',na_rep='NA')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--results-dir',type=Path);ap.add_argument('--output-dir',type=Path)
    a=ap.parse_args();root=a.project.resolve();src=(a.results_dir or root/'results/validation_GSE96058_v16').resolve();out=(a.output_dir or root/'figures/v16_supplementary').resolve();out.mkdir(parents=True,exist_ok=True)
    protected={p.name:sha(p) for p in src.iterdir() if p.is_file()};inputs={}
    def read(name):
        p=src/name;inputs[name]=sha(p);return pd.read_csv(p,sep='\t')
    clinical=read('GSE96058_clinical_IFI44L.tsv');cor=read('GSE96058_IFI44L_immune_signature_correlations.tsv')
    coverage=read('GSE96058_immune_signature_definitions.tsv');models=read('GSE96058_continuous_IFI44L_Cox.tsv')
    coeff=read('GSE96058_all_Cox_coefficients.tsv');lrt=read('GSE96058_conditional_LRT.tsv');ph=read('GSE96058_all_Cox_PH_diagnostics.tsv')
    inputs['35_GSE96058_validation.execution.R']=sha(src/'35_GSE96058_validation.execution.R')
    if len(clinical)!=143 or not clinical.sample_id.is_unique or len(models)!=9 or len(lrt)!=3 or len(cor)!=13:raise ValueError('Incomplete frozen nine-model/three-LRT/13-module SCAN-B results')
    r=one(cor,signature='IFN_ALPHA_RESPONSE')
    scatter=clinical[['sample_id','patient_source_id','instrument','IFI44L_expression','IFI44L_z','IFNa19_score','IFNa19_z']].copy()
    scatter.insert(0,'panel','a');scatter['displayed_spearman_rho']=r.spearman_rho;scatter['n_full_expression_cohort']=int(r['n'])
    scatter['score_definition']='existing IFI44L-excluded19-gene IFN-alpha ssGSEA score; full143-case sample-SD standardization'
    if scatter[['IFI44L_z','IFNa19_z']].isna().any().any():raise ValueError('Cannot silently omit missing scatter observations')
    module_rows=[]
    for i,(signature,label) in enumerate(MODULES):
        r=one(cor,signature=signature);cov=coverage[coverage.signature==signature]
        module_rows.append(dict(panel='b',row_order=i,signature=signature,label=label,spearman_rho=r.spearman_rho,p_value=r.p_value,
          fdr_bh=r.fdr_bh,BH_family_size=int(r.BH_family_size),n=int(r['n']),fixed_gene_count=len(cov),
          present_gene_count=int(cov.present.sum()),variable_gene_count=int(cov.variable_in_cohort.sum())))
    modules=pd.DataFrame(module_rows)
    forest=[]
    for i,mid in enumerate(MODEL_ORDER):
        r=one(models,model_id=mid);c=one(coeff,model_id=mid,term='IFI44L_z');q=lrt[lrt.full_model_id==mid]
        if len(q)>1:raise ValueError('Ambiguous conditional LRT source')
        forest.append(dict(panel='c',row_order=i,model_id=mid,label=MODEL_LABELS[mid],n=int(r['n']),events=int(r.events),
          coefficient_count=int(r.coefficient_count),events_per_coefficient=r.events_per_coefficient,
          hazard_ratio=c.hazard_ratio,ci95_low=c.ci95_low,ci95_high=c.ci95_high,nominal_Wald_p=c.p_value,
          conditional_LRT_p=float(q.p_value.iloc[0]) if len(q) else np.nan,
          conditional_LRT_fdr_bh=float(q.fdr_bh.iloc[0]) if len(q) else np.nan,conditional_LRT_BH_n=int(q.BH_family_size.iloc[0]) if len(q) else np.nan,
          nominal_IFI44L_PH_p=r.proportional_hazards_p,nominal_global_PH_p=r.global_proportional_hazards_p,
          nominal_PH_flag=bool(min(r.proportional_hazards_p,r.global_proportional_hazards_p)<.05),
          scale='IFI44L HR per full143-case sample SD;95% CI; OS days; Efron ties'))
    forest=pd.DataFrame(forest);all9=models.copy();all9.insert(0,'plotted_in_panel_c',all9.model_id.isin(MODEL_ORDER))
    for name,data in [('scatter',scatter),('modules',modules),('Cox',forest),('all9_models',all9)]:save_tsv(data,out/f'{BASE}_{name}_source.tsv')
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'font.size':7.6,
      'axes.labelsize':7.6,'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'legend.fontsize':7.2,
      'axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
    fig=plt.figure(figsize=(9.1,8.5));outer=fig.add_gridspec(2,1,left=.075,right=.975,bottom=.15,top=.89,height_ratios=[1.03,1.22],hspace=.48)
    top=outer[0].subgridspec(1,2,width_ratios=[1,1.62],wspace=.35);axa=fig.add_subplot(top[0,0])
    for inst,col in [('HiSeq 2000',BLUE),('NextSeq 500',ORANGE)]:
        data=scatter[scatter.instrument==inst]
        axa.scatter(data.IFI44L_z,data.IFNa19_z,s=13,color=col,edgecolors='none',alpha=.7,label=f'{inst} (n={len(data)})')
    if set(scatter.instrument)!= {'HiSeq 2000','NextSeq 500'}:raise ValueError('Unrecognized platform would be silently dropped from scatter')
    axa.axhline(0,color='#D8D8D8',lw=.45,zorder=0);axa.axvline(0,color='#D8D8D8',lw=.45,zorder=0)
    axa.margins(.10);axa.set_xlabel('IFI44L expression (SD)');axa.set_ylabel('Original IFNa19 score (SD)')
    axa.set_title('SCAN-B expression context',fontweight='bold',pad=11);axa.text(-.16,1.10,'a',transform=axa.transAxes,fontsize=10,fontweight='bold',va='bottom')
    axa.text(.04,.96,f'n = {len(scatter)}; Spearman ρ = {scatter.displayed_spearman_rho.iloc[0]:.3f}',transform=axa.transAxes,va='top',fontsize=7.4,
      bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
    axa.legend(loc='lower right',handlelength=1,handletextpad=.3,borderaxespad=.4)
    bg=top[0,1].subgridspec(1,3,width_ratios=[1.55,1.65,.46],wspace=.06)
    bl=fig.add_subplot(bg[0,0]);axb=fig.add_subplot(bg[0,1]);bv=fig.add_subplot(bg[0,2])
    for axis in [bl,axb,bv]:axis.set_ylim(-.6,12.6)
    for axis in [bl,bv]:axis.set_axis_off()
    axb.set_xlim(-1,1);axb.axvline(0,color=GREY,lw=.65,ls='--');axb.set_yticks([]);axb.set_xticks([-1,-.5,0,.5,1]);axb.set_xlabel('Spearman ρ')
    bl.text(0,1.10,'b',transform=bl.transAxes,fontsize=10,fontweight='bold',va='bottom');bl.text(0,1.035,'All 13 fixed modules',transform=bl.transAxes,fontsize=8.5,fontweight='bold',va='bottom');bv.text(0,1.015,'ρ',transform=bv.transAxes,fontsize=7.6,va='bottom')
    for _,row in modules.iterrows():
        y=12-row.row_order;colour=BLUE if row.signature=='IFN_ALPHA_RESPONSE' else TEAL
        bl.text(0,y,row.label,va='center',fontsize=7.3);axb.plot([0,row.spearman_rho],[y,y],color=colour,lw=1.1);axb.scatter([row.spearman_rho],[y],color=colour,s=15,zorder=3)
        bv.text(0,y,f'{row.spearman_rho:.3f}',va='center',fontsize=7.3,color=colour)
    cg=outer[1].subgridspec(1,3,width_ratios=[1.64,1.7,1.02],wspace=.08)
    cl=fig.add_subplot(cg[0,0]);axc=fig.add_subplot(cg[0,1]);cv=fig.add_subplot(cg[0,2])
    for axis in [cl,axc,cv]:axis.set_ylim(-.55,5.55)
    for axis in [cl,cv]:axis.set_axis_off()
    low=float(forest.ci95_low.min()*.8);high=float(forest.ci95_high.max()*1.2)
    if not(0<low<high<np.inf):raise ValueError('Invalid full forest CI range')
    axc.set_xscale('log');axc.set_xlim(low,high);axc.set_yticks([]);axc.axvline(1,color=GREY,lw=.65,ls='--')
    ticks=[x for x in [.1,.2,.25,.5,1,2,3,5] if low<=x<=high];axc.xaxis.set_major_locator(FixedLocator(ticks));axc.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'));axc.xaxis.set_minor_locator(NullLocator());axc.set_xlabel('IFI44L HR per SD (log scale)')
    cl.text(0,1.10,'c',transform=cl.transAxes,fontsize=10,fontweight='bold',va='bottom');cl.text(0,1.025,'OS associations: all six models containing IFI44L',transform=cl.transAxes,fontsize=8.5,fontweight='bold',va='bottom');cv.text(0,1.025,'HR [95% CI]; conditional LRT q',transform=cv.transAxes,fontsize=7.4,va='bottom')
    for _,row in forest.iterrows():
        y=5-row.row_order;conditional=np.isfinite(row.conditional_LRT_fdr_bh);colour=ORANGE if row.model_id=='M9' else BLUE if conditional else GREY
        label=row.label.replace('; platform-stratified','\nplatform-stratified')
        label=f'{row.model_id}  '+label+('†' if row.nominal_PH_flag else '')+f'\nn={row.n}; events={row.events}; k={row.coefficient_count}'
        cl.text(0,y,label,va='center',fontsize=7.3,linespacing=1.15)
        axc.errorbar(row.hazard_ratio,y,xerr=[[row.hazard_ratio-row.ci95_low],[row.ci95_high-row.hazard_ratio]],fmt='s' if conditional else 'o',color=colour,markersize=4.2,elinewidth=1,capsize=2,capthick=.8)
        value=f'{row.hazard_ratio:.3f} [{row.ci95_low:.3f}–{row.ci95_high:.3f}]'
        if conditional:value+=f'\nq={row.conditional_LRT_fdr_bh:.3f}'
        cv.text(0,y,value,va='center',fontsize=7.3,color=colour,linespacing=1.15)
    fig.text(.075,.985,'Supplementary Figure S4 | SCAN-B RNA-seq context and OS associations',fontsize=10,fontweight='bold',va='top')
    fig.text(.075,.953,'Clinical-receptor-defined GSE96058: 143 primary cases, 26 OS events; explicit technical replicates excluded.',fontsize=7.5,va='top')
    fig.text(.075,.097,'Clinical = age, tumour size, node status and NHG. k = fitted coefficients; the clinical conditional model has 26 events / 6 coefficients.',fontsize=7.3,va='top')
    fig.text(.075,.073,'Conditional LRT q-values use a separate three-test BH family. All 13 module tests, nine models and complete PH diagnostics are retained.',fontsize=7.3,va='top')
    fig.text(.075,.049,'Source log2(FPKM + 0.1) is unchanged; IFNa19 excludes IFI44L. These exploratory associations do not establish incremental prediction.',fontsize=7.3,va='top')
    flags=ph[ph.p_value<.05]
    if len(flags):fig.text(.075,.025,'† Nominal term/global PH flag; see the complete PH source table. Non-significant tests do not verify the PH assumption.',fontsize=7.2,va='top')
    fig.canvas.draw();fonts=[x.get_fontsize() for x in fig.findobj(mpl.text.Text) if x.get_text()]
    if min(fonts)<7:raise ValueError('Nonempty figure label smaller than7pt')
    if not(axc.get_xlim()[0]<forest.ci95_low.min()<forest.ci95_high.max()<axc.get_xlim()[1]):raise ValueError('Forest CIs would be clipped')
    for ext in ['pdf','svg','png']:fig.savefig(out/f'{BASE}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
    changed=[name for name,h in protected.items() if sha(src/name)!=h]
    if changed:raise RuntimeError('Analysis output changed during plotting: '+str(changed))
    receipt=dict(status='SAVED_RESULT_FIGURE_GENERATED_REQUIRES_ACTUAL_VISUAL_QA',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),script_sha256=sha(Path(__file__)),input_sha256=inputs,
      source_result_files_protected=len(protected),source_result_modifications=changed,source_row_counts=dict(scatter=len(scatter),module_tests=len(modules),plotted_IFI44L_models=len(forest),all_model_source=len(all9)),
      minimum_nonempty_font_pt=min(fonts),figure_inches=[9.1,8.5],PDF_fonttype=42,SVG_fonttype='none',PNG_dpi=300,
      forest_axis_scale='log',full_IFI44L_CI_range=[float(forest.ci95_low.min()),float(forest.ci95_high.max())],forest_xlim=[low,high],all_CIs_unclipped=True,
      conditional_LRT_family_size=3,all_module_tests_retained=True,no_model_fitting=True,no_correlation_recomputation=True,PH_flags=flags.to_dict('records'),
      software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=mpl.__version__),
      output_files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.glob(BASE+'*') if p.is_file() and p.name!=BASE+'_build_receipt.json'})
    (out/f'{BASE}_build_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','source_row_counts','minimum_nonempty_font_pt','forest_xlim','source_result_modifications']}))

if __name__=='__main__':main()
