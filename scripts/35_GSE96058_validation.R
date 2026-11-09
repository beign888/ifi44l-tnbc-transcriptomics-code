#!/usr/bin/env Rscript
# Source-defined SCAN-B clinical TNBC; specifications frozen before expression review.
suppressPackageStartupMessages({
  library(data.table);library(GSVA);library(survival);library(jsonlite);library(digest)
})
options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE)
value<-function(flag,default=NULL){i<-match(flag,args);if(is.na(i))return(default);if(i==length(args))stop('Missing ',flag);args[[i+1L]]}
project<-normalizePath(value('--project','.'),mustWork=TRUE)
work<-file.path(project,'source_metadata','GSE96058')
out<-value('--output-dir',file.path(project,'results','validation_GSE96058_v16'))
default_protocol<-file.path(project,'protocols','v16','GSE96058_ANALYSIS_PLAN.md')
if(!file.exists(default_protocol))default_protocol<-file.path(project,'audit_v16','GSE96058_ANALYSIS_PLAN.md')
plan<-value('--analysis-plan',default_protocol)
expected_plan<-value('--expected-plan-sha256')
if(is.null(expected_plan))stop('Supply frozen --expected-plan-sha256 before fitting')
plan_sha<-digest(file=plan,algo='sha256')
if(!identical(plan_sha,expected_plan))stop('Frozen analysis-plan SHA256 mismatch; no models fitted')
default_amendment<-file.path(project,'protocols','v16','GSE96058_PLATFORM_AMENDMENT.md')
if(!file.exists(default_amendment))default_amendment<-file.path(project,'audit_v16','GSE96058_PLATFORM_AMENDMENT.md')
amendment<-value('--analysis-plan-amendment',default_amendment)
expected_amendment<-value('--expected-amendment-sha256')
if(is.null(expected_amendment))stop('Supply frozen --expected-amendment-sha256 before fitting')
amendment_sha<-digest(file=amendment,algo='sha256')
if(!identical(amendment_sha,expected_amendment))stop('Frozen platform-amendment SHA256 mismatch; no models fitted')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
provenance_file<-file.path(work,'GSE96058_frozen_input_provenance.json')
provenance<-fromJSON(provenance_file)
stopifnot(provenance$eligible_clinical_TNBC==143L,provenance$gene_rows==30865L,
          provenance$analysis_plan_sha256==plan_sha)
input_file<-value('--input',file.path(work,provenance$subset_file))
input_sha<-digest(file=input_file,algo='sha256')
is_full<-identical(input_sha,provenance$full_source_sha256)
if(!is_full&&!identical(input_sha,provenance$subset_sha256))stop('Expression input is neither verified full source nor frozen subset')
input_paths<-c(expression=input_file,metadata=file.path(work,'GSE96058_all3409_metadata.tsv'),
  linkage=file.path(work,'GSE96058_expression_linkage.tsv'),modules=file.path(work,'baseline_module_definitions.tsv'),
  protocol=plan,amendment=amendment,provenance=provenance_file)
input_manifest<-data.table(role=names(input_paths),path=unname(input_paths),bytes=file.info(input_paths)$size,
  sha256=vapply(input_paths,function(p)digest(file=p,algo='sha256'),character(1)))
stopifnot(input_manifest[role=='metadata',sha256]==provenance$metadata_sha256,
          input_manifest[role=='linkage',sha256]==provenance$linkage_sha256,
          input_manifest[role=='modules',sha256]==provenance$module_definitions_sha256)
fwrite(input_manifest,file.path(out,'GSE96058_input_manifest.tsv'),sep='\t')
file.copy(plan,file.path(out,'GSE96058_ANALYSIS_PLAN.execution.md'),overwrite=TRUE)
file.copy(amendment,file.path(out,'GSE96058_PLATFORM_AMENDMENT.execution.md'),overwrite=TRUE)
script_arg<-grep('^--file=',commandArgs(),value=TRUE)
if(length(script_arg)==1L){script<-sub('^--file=','',script_arg);file.copy(script,file.path(out,'35_GSE96058_validation.execution.R'),overwrite=TRUE)}
write_json(list(status='STARTED_BEFORE_MODELS',recorded_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
  plan_sha256=plan_sha,amendment_sha256=amendment_sha,input_sha256=input_sha,input_route=if(is_full)'full official CSV' else 'exact143-column frozen subset',
  prespecified_models=9L,prespecified_LRT_family=3L,IFI44L_expression_inspected_before_protocol=FALSE),
  file.path(out,'GSE96058_STARTUP_RECEIPT.json'),auto_unbox=TRUE,pretty=TRUE)

warnings_rows<-list();errors_rows<-list();coef_rows<-list();PH_rows<-list();definition_rows<-list()
summary_rows<-list();design_rows<-list();VIF_rows<-list();scaling_rows<-list();eligibility_rows<-list();LRT_rows<-list()
warn<-function(model,step){function(w){
  warnings_rows[[length(warnings_rows)+1L]]<<-data.table(model=model,step=step,warning=conditionMessage(w));invokeRestart('muffleWarning')
}}
metadata<-fread(input_paths['metadata'],check.names=FALSE,na.strings=c('NA',''))
linkage<-fread(input_paths['linkage'],check.names=FALSE)
clinical<-metadata[eligible_receptor_TNBC==TRUE]
stopifnot(nrow(metadata)==3409L,nrow(clinical)==143L,nrow(linkage)==143L,
          !anyDuplicated(clinical$sample_id),!anyDuplicated(clinical$patient_source_id),
          all(clinical[['er status']]=='0'),all(clinical[['pgr status']]=='0'),all(clinical[['her2 status']]=='0'),
          !any(clinical$is_technical_replicate))
clinical<-clinical[match(linkage$sample_id,sample_id)]
stopifnot(identical(as.character(clinical$title),as.character(linkage$expression_title)))
if(is_full){
  # fread selects only the fixed143 columns while streaming/parsing the full CSV;
  # it never constructs a3409-sample numeric matrix.
  selected<-c(1L,as.integer(linkage$source_column_index_0based)+1L)
  raw<-fread(input_file,select=selected,check.names=FALSE,na.strings=c('NA',''))
  stopifnot(identical(names(raw)[-1L],as.character(linkage$expression_title)))
  setnames(raw,c('gene',as.character(linkage$sample_id)))
}else raw<-fread(input_file,check.names=FALSE,na.strings=c('NA',''))
stopifnot(nrow(raw)==30865L,names(raw)[1]=='gene',!anyDuplicated(raw$gene),
          identical(names(raw)[-1L],as.character(clinical$sample_id)))
gene_expr<-as.matrix(raw[,-1L,with=FALSE]);storage.mode(gene_expr)<-'double';rownames(gene_expr)<-raw$gene
stopifnot(all(is.finite(gene_expr)),'IFI44L'%in%rownames(gene_expr))
rm(raw);gc(verbose=FALSE)
variable<-apply(gene_expr,1L,sd)>0
module_def<-fread(input_paths['modules'],check.names=FALSE)
signatures<-lapply(split(as.character(module_def$gene),module_def$signature),function(g)setdiff(unique(g),'IFI44L'))
stopifnot(length(signatures)==13L,length(signatures[['IFN_ALPHA_RESPONSE']])==19L)
coverage<-rbindlist(lapply(names(signatures),function(nm)data.table(signature=nm,gene=signatures[[nm]],
  present=signatures[[nm]]%in%rownames(gene_expr),variable_in_cohort=signatures[[nm]]%in%rownames(gene_expr)[variable])))
available<-lapply(signatures,intersect,y=rownames(gene_expr))
stopifnot(all(lengths(available)>=2L),!any(vapply(available,function(g)'IFI44L'%in%g,logical(1))))
# Same GSVA ssGSEA implementation as the existing external cohorts. All source
# background rows are supplied; GSVA's own method-specific filters are logged.
scores<-withCallingHandlers(gsva(ssgseaParam(gene_expr,available,minSize=2,maxSize=500,normalize=TRUE),verbose=FALSE),warning=warn('fixed13modules','GSVA'))
scores<-scores[,clinical$sample_id,drop=FALSE]
standardize<-function(x,label){s<-sd(x);if(any(!is.finite(x))||!is.finite(s)||s<=0)stop('Invalid full-cohort scaling: ',label)
  scaling_rows[[length(scaling_rows)+1L]]<<-data.table(predictor=label,n=length(x),mean=mean(x),sample_SD=s,
    scope='all143 eligible expression cases before endpoint/clinical exclusions');(x-mean(x))/s}
number<-function(x)suppressWarnings(as.numeric(as.character(x)))
clinical[,`:=`(IFI44L_expression=as.numeric(gene_expr['IFI44L',]),
  IFNa19_score=as.numeric(scores['IFN_ALPHA_RESPONSE',]),age_years=number(`age at diagnosis`),
  age10=number(`age at diagnosis`)/10,size_mm=number(`tumor size`),size10=number(`tumor size`)/10,
  os_days=number(`overall survival days`),os_event=number(`overall survival event`),
  node_group=factor(`lymph node status`,levels=c('NodeNegative','NodePositive')),
  grade_group=factor(ifelse(nhg%in%c('G1','G2'),'G1/2',ifelse(nhg=='G3','G3',NA)),levels=c('G1/2','G3')),
  instrument=factor(`instrument model`,levels=c('HiSeq 2000','NextSeq 500')))]
clinical[,IFI44L_z:=standardize(IFI44L_expression,'IFI44L')]
clinical[,IFNa19_z:=standardize(IFNa19_score,'original_IFNa19_score')]
clinical[,positive_OS_eligible:=is.finite(os_days)&os_days>0&os_event%in%c(0,1)]
fwrite(clinical,file.path(out,'GSE96058_clinical_IFI44L.tsv'),sep='\t',na='NA')
fwrite(coverage,file.path(out,'GSE96058_immune_signature_definitions.tsv'),sep='\t',na='NA')
fwrite(as.data.table(t(scores),keep.rownames='sample_id'),file.path(out,'GSE96058_immune_signature_ssGSEA_scores.tsv'),sep='\t')
correlations<-rbindlist(lapply(rownames(scores),function(nm){
  ct<-suppressWarnings(cor.test(clinical$IFI44L_expression,as.numeric(scores[nm,]),method='spearman',exact=FALSE))
  data.table(signature=nm,spearman_rho=unname(ct$estimate),p_value=ct$p.value,n=nrow(clinical))}))
correlations[,`:=`(fdr_bh=p.adjust(p_value,'BH',n=13L),BH_family_size=13L)]
fwrite(correlations,file.path(out,'GSE96058_IFI44L_immune_signature_correlations.tsv'),sep='\t')
fwrite(data.table(gene=rownames(gene_expr),as.data.table(gene_expr)),file.path(out,'GSE96058_normalized_gene_expression.tsv.gz'),sep='\t',compress='gzip')
fwrite(rbindlist(scaling_rows),file.path(out,'GSE96058_expression_cohort_standardization.tsv'),sep='\t')

make_cases<-function(family,needed){
  d<-copy(clinical)
  reasons<-vapply(seq_len(nrow(d)),function(i){
    why<-character()
    if(!is.finite(d$os_days[i])||d$os_days[i]<=0)why<-c(why,'missing/nonpositive OS days')
    if(is.na(d$os_event[i])||!d$os_event[i]%in%c(0,1))why<-c(why,'missing/invalid OS event')
    for(name in needed){x<-d[[name]][i];if(is.na(x)||(is.numeric(x)&&!is.finite(x)))why<-c(why,paste0('missing/nonfinite ',name))}
    if(length(why))paste(why,collapse=';') else 'included'
  },character(1))
  d[,`:=`(case_family=family,eligible_common_complete_case=reasons=='included',exclusion_reason=reasons)]
  eligibility_rows[[length(eligibility_rows)+1L]]<<-d
  as.data.frame(d[eligible_common_complete_case==TRUE])
}
age_cases<-make_cases('age_IFNa_common',c('age10','IFNa19_z','IFI44L_z'))
clinical_cases<-make_cases('clinical_IFNa_common',c('age10','size10','node_group','grade_group','IFNa19_z','IFI44L_z'))
platform_cases<-make_cases('platform_stratified_age_IFNa_common',c('age10','IFNa19_z','IFI44L_z','instrument'))
unadjusted_cases<-make_cases('unadjusted',c('IFI44L_z'))
stopifnot(nrow(age_cases)==143L,sum(age_cases$os_event)==26L,
          nrow(clinical_cases)==137L,sum(clinical_cases$os_event)==26L,nrow(platform_cases)==143L)

diagnostics<-function(d,rhs,label){
  X<-model.matrix(as.formula(paste('~',rhs)),data=d);X<-X[,colnames(X)!='(Intercept)',drop=FALSE]
  valid<-all(is.finite(apply(X,2L,sd))&apply(X,2L,sd)>0)
  if(valid){Z<-scale(X);sing<-svd(Z,nu=0,nv=0)$d;tol<-max(dim(Z))*.Machine$double.eps*max(sing);rank<-sum(sing>tol);condition<-max(sing)/min(sing)}
  else{Z<-NULL;sing<-numeric();tol<-NA_real_;rank<-NA_integer_;condition<-Inf}
  design_rows[[length(design_rows)+1L]]<<-data.table(model=label,n=nrow(d),covariate_columns=ncol(X),
    standardized_design_rank=rank,full_rank=valid&&!is.na(rank)&&rank==ncol(X),standardized_design_condition_number=condition,
    rank_tolerance=tol,singular_values=paste(format(sing,digits=17),collapse=';'),
    definition='sample-SD-standardized complete-case covariate design; intercept/strata excluded')
  for(j in seq_len(ncol(X))){
    vif<-NA_real_
    if(valid){if(ncol(X)==1L)vif<-1 else{res<-qr.resid(qr(cbind(Intercept=1,Z[,-j,drop=FALSE])),Z[,j]);frac<-sum(res^2)/sum(Z[,j]^2);vif<-if(frac>0)1/frac else Inf}}
    VIF_rows[[length(VIF_rows)+1L]]<<-data.table(model=label,term=colnames(X)[j],VIF=vif)
  }
  valid&&!is.na(rank)&&rank==ncol(X)
}
fits<-list()
fit_one<-function(id,label,rhs,d,family,stratify=FALSE){
  rhs_fit<-if(stratify)paste0(rhs,'+strata(instrument)') else rhs
  formula<-as.formula(paste('Surv(os_days,os_event)~',rhs_fit));before<-length(warnings_rows)
  fit<-tryCatch({
    if(!diagnostics(d,rhs,label))stop('Rank deficient covariate design; no variable dropped')
    withCallingHandlers(coxph(formula,data=d,ties='efron',singular.ok=FALSE,x=TRUE,y=TRUE,model=TRUE),warning=warn(label,'coxph'))
  },error=function(e){errors_rows[[length(errors_rows)+1L]]<<-data.table(model=label,step='coxph',error=conditionMessage(e));NULL})
  if(is.null(fit)){
    summary_rows[[length(summary_rows)+1L]]<<-data.table(model_id=id,model=label,status='not_estimable',endpoint='Overall survival',case_family=family,n=nrow(d),events=sum(d$os_event));return(invisible(NULL))
  }
  fits[[id]]<<-fit;s<-summary(fit)
  coef_rows[[length(coef_rows)+1L]]<<-data.table(model_id=id,model=label,term=rownames(s$coefficients),
    coefficient=s$coefficients[,'coef'],SE=s$coefficients[,'se(coef)'],hazard_ratio=s$conf.int[,'exp(coef)'],
    ci95_low=s$conf.int[,'lower .95'],ci95_high=s$conf.int[,'upper .95'],p_value=s$coefficients[,'Pr(>|z|)'],n=fit$n,events=fit$nevent)
  ph<-tryCatch(withCallingHandlers(cox.zph(fit,transform='km'),warning=warn(label,'cox.zph')),
    error=function(e){errors_rows[[length(errors_rows)+1L]]<<-data.table(model=label,step='cox.zph',error=conditionMessage(e));NULL})
  if(!is.null(ph))PH_rows[[length(PH_rows)+1L]]<<-data.table(model_id=id,model=label,term=rownames(ph$table),
    chisq=ph$table[,'chisq'],df=ph$table[,'df'],p_value=ph$table[,'p'])
  has_ifi<-'IFI44L_z'%in%rownames(s$coefficients)
  summary_rows[[length(summary_rows)+1L]]<<-data.table(model_id=id,model=label,status='fitted',endpoint='Overall survival',case_family=family,
    n=fit$n,events=fit$nevent,coefficient_count=length(coef(fit)),events_per_coefficient=fit$nevent/length(coef(fit)),
    log_likelihood=as.numeric(logLik(fit)),model_df=attr(logLik(fit),'df'),
    hazard_ratio_per_IFI44L_SD=if(has_ifi)s$conf.int['IFI44L_z','exp(coef)'] else NA_real_,
    ci95_low=if(has_ifi)s$conf.int['IFI44L_z','lower .95'] else NA_real_,ci95_high=if(has_ifi)s$conf.int['IFI44L_z','upper .95'] else NA_real_,
    p_value=if(has_ifi)s$coefficients['IFI44L_z','Pr(>|z|)'] else NA_real_,
    proportional_hazards_p=if(has_ifi&&!is.null(ph))ph$table['IFI44L_z','p'] else NA_real_,
    global_proportional_hazards_p=if(!is.null(ph))ph$table['GLOBAL','p'] else NA_real_,warning_count=length(warnings_rows)-before)
  definition_rows[[length(definition_rows)+1L]]<<-data.table(model_id=id,model=label,formula=paste(deparse(formula),collapse=' '),case_family=family,
    n=fit$n,events=fit$nevent,time_source='overall survival days, strictly positive',event_source='overall survival event, source0/1',
    predictor_scale='IFI44L/IFNa standardized over all143; age per10years; size per10mm; source log2(FPKM+0.1) unchanged',
    categorical_reference='explicit source lymph node status NodeNegative vs NodePositive only; all other status labels missing/excluded; NHG G1/2 vsG3; instrument stratifies baseline hazard when specified',
    missingness='complete cases; common cases across paired IFNa nested models; no imputation',ties='Efron',
    PH='all-term/global Schoenfeld tests with km transform; non-significance does not verify PH',
    role='frozen exploratory revision extension;26 events/5-6clinical coefficients; conditional association, not validated incremental prediction')
  invisible(fit)
}
fit_one('M1','Unadjusted','IFI44L_z',unadjusted_cases,'unadjusted')
fit_one('M2','Age-adjusted','age10+IFI44L_z',age_cases,'age_IFNa_common')
fit_one('M3','Clinicopathological sensitivity','age10+size10+node_group+grade_group+IFI44L_z',clinical_cases,'clinical_IFNa_common')
fit_one('M4','Age+IFNa','age10+IFNa19_z',age_cases,'age_IFNa_common')
fit_one('M5','Age+IFNa+IFI44L','age10+IFNa19_z+IFI44L_z',age_cases,'age_IFNa_common')
fit_one('M6','Clinical+IFNa','age10+size10+node_group+grade_group+IFNa19_z',clinical_cases,'clinical_IFNa_common')
fit_one('M7','Clinical+IFNa+IFI44L','age10+size10+node_group+grade_group+IFNa19_z+IFI44L_z',clinical_cases,'clinical_IFNa_common')
fit_one('M8','Platform-stratified age+IFNa','age10+IFNa19_z',platform_cases,'platform_stratified_age_IFNa_common',TRUE)
fit_one('M9','Platform-stratified age+IFNa+IFI44L','age10+IFNa19_z+IFI44L_z',platform_cases,'platform_stratified_age_IFNa_common',TRUE)
for(pair in list(c('M4','M5'),c('M6','M7'),c('M8','M9'))){
  b<-fits[[pair[1]]];f<-fits[[pair[2]]]
  if(is.null(b)||is.null(f)){LRT_rows[[length(LRT_rows)+1L]]<-data.table(base_model_id=pair[1],full_model_id=pair[2],status='not_estimable',p_value=NA_real_);next}
  stopifnot(identical(rownames(b$model),rownames(f$model)),b$n==f$n,b$nevent==f$nevent)
  ll0<-as.numeric(logLik(b));ll1<-as.numeric(logLik(f));df<-attr(logLik(f),'df')-attr(logLik(b),'df');stat<-2*(ll1-ll0)
  stopifnot(df==1L,stat>= -1e-8)
  LRT_rows[[length(LRT_rows)+1L]]<-data.table(base_model_id=pair[1],full_model_id=pair[2],status='estimable',n=f$n,events=f$nevent,
    base_log_likelihood=ll0,full_log_likelihood=ll1,raw_LRT_statistic=stat,LRT_statistic=max(0,stat),df=df,
    p_value=pchisq(max(0,stat),df=df,lower.tail=FALSE),role='additional cohort exploratory conditional association; not independent predictive utility')
}
lrt<-rbindlist(LRT_rows,fill=TRUE);stopifnot(nrow(lrt)==3L)
lrt[,`:=`(fdr_bh=p.adjust(p_value,'BH',n=3L),BH_family_size=3L,family='additional_GSE96058_exploratory3')]
eligibility<-rbindlist(eligibility_rows,fill=TRUE)
empty_warnings<-data.table(model=character(),step=character(),warning=character())
empty_errors<-data.table(model=character(),step=character(),error=character())
platform_counts<-rbindlist(lapply(list(all143=as.data.frame(clinical),clinical137=clinical_cases),function(d){as.data.table(d)[,.(n=.N,events=sum(os_event)),by=instrument]}),idcol='case_set')
group_counts<-rbindlist(lapply(c('node_group','grade_group'),function(v){as.data.table(clinical_cases)[,.(n=.N,events=sum(os_event)),by=v][,variable:=v]}),fill=TRUE)
predictors<-model.matrix(~age10+size10+node_group+grade_group+IFNa19_z+IFI44L_z,data=clinical_cases)[,-1L,drop=FALSE]
predictor_correlation<-rbindlist(lapply(combn(colnames(predictors),2L,simplify=FALSE),function(p)data.table(variable1=p[1],variable2=p[2],n=nrow(predictors),
  pearson_r=cor(predictors[,p[1]],predictors[,p[2]]),spearman_rho=cor(predictors[,p[1]],predictors[,p[2]],method='spearman'))))
tables<-list(continuous_IFI44L_Cox=rbindlist(summary_rows,fill=TRUE),all_Cox_coefficients=rbindlist(coef_rows,fill=TRUE),
  all_Cox_PH_diagnostics=rbindlist(PH_rows,fill=TRUE),Cox_model_definitions=rbindlist(definition_rows,fill=TRUE),
  Cox_design_diagnostics=rbindlist(design_rows,fill=TRUE),Cox_VIF_diagnostics=rbindlist(VIF_rows,fill=TRUE),
  conditional_LRT=lrt,complete_case_and_exclusion_records=eligibility,excluded_records=eligibility[eligible_common_complete_case==FALSE],
  platform_case_event_counts=platform_counts,clinical_model_group_counts=group_counts,predictor_correlation_diagnostics=predictor_correlation,
  Cox_warnings=if(length(warnings_rows))rbindlist(warnings_rows) else empty_warnings,Cox_errors=if(length(errors_rows))rbindlist(errors_rows) else empty_errors)
for(nm in names(tables))fwrite(tables[[nm]],file.path(out,paste0('GSE96058_',nm,'.tsv')),sep='\t',na='NA')
fwrite(data.table(dataset='GSE96058',full_source_gene_rows=nrow(gene_expr),variable_genes=sum(variable),expression_cases=ncol(gene_expr),
  source_scale='log2(FPKM+0.1); no further log or count normalization',additional_log_transform=FALSE,
  min_expression=min(gene_expr),max_expression=max(gene_expr),module_count=length(signatures),IFI44L_in_any_module=FALSE,
  primary_OS_n=nrow(age_cases),primary_OS_events=sum(age_cases$os_event),clinical_n=nrow(clinical_cases),clinical_events=sum(clinical_cases$os_event)),
  file.path(out,'GSE96058_audit.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'GSE96058_R_sessionInfo.txt'))
receipt<-list(status='ACTUAL_MODELS_EXECUTED_REQUIRES_INDEPENDENT_REVIEW',plan_sha256=plan_sha,amendment_sha256=amendment_sha,
  expression_input_route=if(is_full)'full official CSV' else 'frozen exact143 subset',expression_input_sha256=input_sha,
  expression_cases=143L,source_gene_rows=30865L,planned_models=9L,fitted_models=sum(tables$continuous_IFI44L_Cox$status=='fitted'),
  planned_LRT=3L,module_tests=13L,warnings_n=nrow(tables$Cox_warnings),errors_n=nrow(tables$Cox_errors),
  limitations=c('26OS events; clinical models5-6coefficients.','Source diagnostic receptor calls retained; not retrospectively assigned a modern uniform1% cutoff.',
  'Platform strata permit different baseline hazards; zero-event stratum has no independent event information.',
  'No KM cutpoint, independent clinical-utility claim, treatment-response analysis or significance-based specification selection.'))
write_json(receipt,file.path(out,'GSE96058_EXECUTION_RECEIPT.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
files<-list.files(out,full.names=TRUE);files<-files[!file.info(files)$isdir&basename(files)!='OUTPUT_MANIFEST_SHA256.tsv']
fwrite(data.table(file=basename(files),bytes=file.info(files)$size,sha256=vapply(files,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'OUTPUT_MANIFEST_SHA256.tsv'),sep='\t')
print(lrt);message('GSE96058 frozen9-model/3-LRT/13-module execution completed; independent verification remains required.')
