#!/usr/bin/env Rscript
# Frozen independent-cohort validation; no cutpoint or covariate selection by P.
suppressPackageStartupMessages({
  library(data.table);library(AnnotationDbi);library(hgu133plus2.db)
  library(limma);library(GSVA);library(survival);library(jsonlite);library(digest)
})
args<-commandArgs(trailingOnly=TRUE)
value<-function(flag,default=NULL){i<-match(flag,args);if(is.na(i))return(default);if(i==length(args))stop('Missing ',flag);args[[i+1L]]}
project<-normalizePath(value('--project','.'),mustWork=TRUE)
work<-file.path(project,'source_metadata','GSE21653')
out<-value('--output-dir',file.path(project,'results','validation_GSE21653'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(work,'GSE21653_ANALYSIS_DESIGN.json')))
input_file<-value('--input',file.path(work,'GSE21653_TNBC87_frozen_all_probes.tsv.gz'))
provenance<-fromJSON(file.path(work,'GSE21653_frozen_input_provenance.json'))
is_full<-!grepl('frozen_all_probes',basename(input_file),fixed=TRUE)
expected_input_sha<-if(is_full)provenance$full_source_sha256 else provenance$subset_sha256
stopifnot(digest(input_file,algo='sha256',file=TRUE)==expected_input_sha,
  digest(file.path(work,'GSE21653_all266_metadata.tsv'),algo='sha256',file=TRUE)==provenance$all266_metadata_sha256,
  digest(file.path(work,'baseline_module_definitions.tsv'),algo='sha256',file=TRUE)==
    '8bc3ff195ccaef81d53e6f4543d44f584e16ade8a3b0e04fb8231c3068562bb3')
metadata<-fread(file.path(work,'GSE21653_all266_metadata.tsv'),na.strings=c('NA',''))
clinical<-metadata[eligible_receptor_TNBC==TRUE]
stopifnot(nrow(clinical)==87L,!anyDuplicated(clinical$sample_id),
          all(clinical[['er ihc']]=='0'),all(clinical[['pr ihc']]=='0'),all(clinical$erbb2=='0'),
          all(clinical$data_processing=='RMA using a quantile normalization'))

# A full GEO series matrix or an exact all-probe subset can be used. The source
# values are already RMA-processed; applying log2(x+1) here would be a second log.
input<-fread(input_file,skip=if(is_full) 'ID_REF' else 0,nrows=54613,check.names=FALSE)
stopifnot(nrow(input)==54613L,names(input)[1]=='ID_REF',!anyDuplicated(input[[1]]))
stopifnot(all(clinical$sample_id%in%names(input)))
probe_expr<-as.matrix(input[,clinical$sample_id,with=FALSE]);storage.mode(probe_expr)<-'double'
rownames(probe_expr)<-input[[1]]
stopifnot(all(is.finite(probe_expr)),'204439_at'%in%rownames(probe_expr))
annotation<-as.data.table(AnnotationDbi::select(hgu133plus2.db,keys=rownames(probe_expr),
                 keytype='PROBEID',columns=c('PROBEID','SYMBOL','ENTREZID')))
annotation<-unique(annotation[!is.na(SYMBOL)&nzchar(SYMBOL)])
fwrite(annotation,file.path(out,'GSE21653_probe_annotation.tsv'),sep='\t')
fixed_annotation<-annotation[SYMBOL=='IFI44L']
stopifnot(nrow(fixed_annotation)==1L,fixed_annotation$PROBEID=='204439_at',fixed_annotation$ENTREZID=='10964')
mapped<-probe_expr[match(annotation$PROBEID,rownames(probe_expr)),,drop=FALSE]
rownames(mapped)<-annotation$SYMBOL
gene_expr<-avereps(mapped,ID=rownames(mapped))
fwrite(data.table(gene=rownames(gene_expr),as.data.table(gene_expr)),
  file.path(out,'GSE21653_normalized_gene_expression.tsv.gz'),sep='\t',compress='gzip')
ifi<-as.numeric(probe_expr['204439_at',clinical$sample_id])
stopifnot(identical(ifi,as.numeric(gene_expr['IFI44L',clinical$sample_id])))
number<-function(x)suppressWarnings(as.numeric(as.character(x)))
clinical[,`:=`(patient_source_id=source_case_id,age_years=number(`age at diagnosis`),
  dfs_months=number(`dfs time (months)`),dfs_event=number(`dfs evt`),
  pT_raw=pt,pN_raw=pn,grade_raw=`sbr grade`,IFI44L_expression=ifi,
  IFI44L_z=as.numeric(scale(ifi)),age10=number(`age at diagnosis`)/10)]
clinical[,dfs_days:=dfs_months*365.25/12]
clinical[,`:=`(IFI44L_group=ifelse(IFI44L_expression>=median(IFI44L_expression),'High','Low'),
  pT_group=factor(ifelse(pT_raw=='pT1','T1',ifelse(pT_raw%in%c('pT2','pT3'),'T2/3',NA)),levels=c('T1','T2/3')),
  pN_group=factor(ifelse(pN_raw=='0','N0',ifelse(pN_raw=='1','N-positive',NA)),levels=c('N0','N-positive')),
  grade_group=factor(ifelse(grade_raw%in%c('1','2'),'G1/2',ifelse(grade_raw=='3','G3',NA)),levels=c('G1/2','G3')))]
clinical[,primary_DFS_eligible:=is.finite(dfs_months)&dfs_months>0&dfs_event%in%c(0,1)]
clinical[,zero_time_sensitivity_eligible:=is.finite(dfs_months)&dfs_months>=0&dfs_event%in%c(0,1)]
clinical[,primary_exclusion_reason:=ifelse(primary_DFS_eligible,'included',
  ifelse(is.finite(dfs_months)&dfs_months==0&dfs_event%in%c(0,1),'zero recorded time; retained in zero-time sensitivity',
         'missing or invalid recorded DFS time/event'))]
fwrite(clinical,file.path(out,'GSE21653_clinical_IFI44L.tsv'),sep='\t',na='NA')
fwrite(clinical[primary_DFS_eligible==FALSE],file.path(out,'GSE21653_DFS_primary_exclusions.tsv'),sep='\t',na='NA')
fwrite(clinical[is.finite(dfs_months)&dfs_months==0],file.path(out,'GSE21653_zero_time_records.tsv'),sep='\t',na='NA')

warning_rows<-list();coef_rows<-list();ph_rows<-list();definitions<-list()
fit_one<-function(label,rhs,include_zero=FALSE){
  eligible<-if(include_zero)clinical$zero_time_sensitivity_eligible else clinical$primary_DFS_eligible
  d<-as.data.frame(clinical[which(eligible)])
  needed<-c('dfs_months','dfs_event',all.vars(as.formula(paste('~',rhs))))
  d<-d[complete.cases(d[,needed,drop=FALSE]),,drop=FALSE]
  warnings<-character()
  fit<-withCallingHandlers(coxph(as.formula(paste('Surv(dfs_months,dfs_event)~',rhs)),
                    data=d,ties='efron',x=TRUE,y=TRUE,model=TRUE),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  summary_fit<-summary(fit);ph<-cox.zph(fit)
  terms<-rownames(summary_fit$coefficients)
  coef_rows[[label]]<<-data.table(model=label,term=terms,n=fit$n,events=fit$nevent,
    coefficient=summary_fit$coefficients[,'coef'],hazard_ratio=summary_fit$conf.int[,'exp(coef)'],
    ci95_low=summary_fit$conf.int[,'lower .95'],ci95_high=summary_fit$conf.int[,'upper .95'],
    p_value=summary_fit$coefficients[,'Pr(>|z|)'])
  ph_rows[[label]]<<-data.table(model=label,term=rownames(ph$table),as.data.table(ph$table))
  if(length(warnings))warning_rows[[label]]<<-data.table(model=label,warning=warnings)
  definitions[[label]]<<-data.table(model=label,formula=paste('Surv(dfs_months,dfs_event)~',rhs),
    n=fit$n,events=fit$nevent,coefficient_count=length(coef(fit)),events_per_coefficient=fit$nevent/length(coef(fit)),
    time_rule=if(include_zero)'finite DFS >=0; original zero times unchanged' else 'finite DFS >0',
    IFI44L_standardization='sample SD across all 87 eligible patients; fixed probe 204439_at; original RMA values',
    missingness='complete cases for variables in each model; no imputation',ties='efron',
    categorical_coding='pT1 vs pT2/3; pN0 vs pN1; SBR grades1/2 vs3',
    warning_count=length(warnings),interpretation=if(grepl('Clinicopathological',label))
      'exploratory unstable clinical sensitivity; 5 coefficients and limited events; not independent prediction' else 'association only')
  row<-match('IFI44L_z',terms)
  data.table(endpoint='Disease-free survival',model=label,n=fit$n,events=fit$nevent,
    hazard_ratio_per_IFI44L_SD=summary_fit$conf.int[row,'exp(coef)'],
    ci95_low=summary_fit$conf.int[row,'lower .95'],ci95_high=summary_fit$conf.int[row,'upper .95'],
    p_value=summary_fit$coefficients[row,'Pr(>|z|)'],proportional_hazards_p=ph$table['IFI44L_z','p'],
    global_proportional_hazards_p=ph$table['GLOBAL','p'],warning_count=length(warnings),
    coefficients=length(coef(fit)),events_per_coefficient=fit$nevent/length(coef(fit)))
}
cox_results<-rbindlist(list(
  fit_one('Unadjusted','IFI44L_z'),
  fit_one('Age-adjusted','IFI44L_z+age10'),
  fit_one('Clinicopathological sensitivity','IFI44L_z+age10+pT_group+pN_group+grade_group'),
  fit_one('Age-adjusted including zero-time events (sensitivity)','IFI44L_z+age10',TRUE)))
cox_results[,fdr_bh:=p.adjust(p_value,'BH')]
fwrite(cox_results,file.path(out,'GSE21653_continuous_IFI44L_Cox.tsv'),sep='\t')
fwrite(rbindlist(coef_rows),file.path(out,'GSE21653_all_Cox_coefficients.tsv'),sep='\t')
fwrite(rbindlist(ph_rows),file.path(out,'GSE21653_all_Cox_PH_diagnostics.tsv'),sep='\t')
fwrite(rbindlist(definitions),file.path(out,'GSE21653_Cox_model_definitions.tsv'),sep='\t')
warnings<-if(length(warning_rows))rbindlist(warning_rows) else data.table(model=character(),warning=character())
fwrite(warnings,file.path(out,'GSE21653_Cox_warnings.tsv'),sep='\t')
frequency<-rbindlist(lapply(c('pT_group','pN_group','grade_group'),function(v){
 d<-clinical[primary_DFS_eligible&complete.cases(age10,pT_group,pN_group,grade_group)]
 d[,.(n=.N,events=sum(dfs_event)),by=c(v)][,variable:=v]
}),fill=TRUE)
fwrite(frequency,file.path(out,'GSE21653_clinical_model_group_counts.tsv'),sep='\t',na='NA')

definition<-fread(file.path(work,'baseline_module_definitions.tsv'))
signatures<-lapply(split(as.character(definition$gene),definition$signature),function(x)setdiff(unique(x),'IFI44L'))
members<-rbindlist(lapply(names(signatures),function(nm)data.table(signature=nm,gene=signatures[[nm]],present=signatures[[nm]]%in%rownames(gene_expr))))
fwrite(members,file.path(out,'GSE21653_immune_signature_definitions.tsv'),sep='\t')
available<-lapply(signatures,intersect,y=rownames(gene_expr))
stopifnot(length(available)==13L,all(lengths(available)>=2L),!any(vapply(available,function(x)'IFI44L'%in%x,logical(1))))
scores<-gsva(ssgseaParam(gene_expr,available,minSize=2,maxSize=500,normalize=TRUE),verbose=FALSE)
fwrite(as.data.table(t(scores),keep.rownames='sample_id'),file.path(out,'GSE21653_immune_signature_ssGSEA_scores.tsv'),sep='\t')
correlations<-rbindlist(lapply(rownames(scores),function(nm){
 ct<-suppressWarnings(cor.test(ifi,as.numeric(scores[nm,clinical$sample_id]),method='spearman',exact=FALSE))
 data.table(signature=nm,spearman_rho=unname(ct$estimate),p_value=ct$p.value,n=length(ifi))
}))
correlations[,fdr_bh:=p.adjust(p_value,'BH')]
fwrite(correlations,file.path(out,'GSE21653_IFI44L_immune_signature_correlations.tsv'),sep='\t')

km_data<-as.data.frame(clinical[primary_DFS_eligible==TRUE])
km_data$IFI44L_group<-factor(km_data$IFI44L_group,levels=c('Low','High'))
km<-survfit(Surv(dfs_days,dfs_event)~IFI44L_group,data=km_data,conf.type='log-log')
s<-summary(km,censored=TRUE)
fwrite(data.table(group=as.character(s$strata),time_days=s$time,n_risk=s$n.risk,n_event=s$n.event,
 n_censor=s$n.censor,survival=s$surv,ci95_low=s$lower,ci95_high=s$upper),
 file.path(out,'GSE21653_KM_curve.tsv'),sep='\t')
logrank<-survdiff(Surv(dfs_days,dfs_event)~IFI44L_group,data=km_data)
fwrite(data.table(endpoint='DFS',n=nrow(km_data),events=sum(km_data$dfs_event),
 chi_square=logrank$chisq,p_value=pchisq(logrank$chisq,df=1,lower.tail=FALSE),
 grouping='High >= median across87; Low < median; graphical sensitivity'),
 file.path(out,'GSE21653_KM_logrank.tsv'),sep='\t')
expression_audit<-data.table(dataset='GSE21653',source='deposited RMA with quantile normalization',
 additional_log_transform=FALSE,input_probe_rows=nrow(probe_expr),TNBC_patients=ncol(probe_expr),
 mapped_genes=nrow(gene_expr),fixed_IFI44L_probe='204439_at',IFI44L_EntrezID='10964',
 min_expression=min(probe_expr),max_expression=max(probe_expr),
 q01=unname(quantile(probe_expr,.01)),q50=unname(quantile(probe_expr,.5)),q99=unname(quantile(probe_expr,.99)),
 module_count=length(available),IFI44L_in_any_module=FALSE,
 inference_unit='source patient/case, not probe; unique identifiers do not verify identity across repositories')
fwrite(expression_audit,file.path(out,'GSE21653_audit.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'GSE21653_R_sessionInfo.txt'))
message('GSE21653 frozen validation completed; all models and all 13 module correlations retained.')
