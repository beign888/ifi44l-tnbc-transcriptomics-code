#!/usr/bin/env Rscript
# v16 targeted conditional-association and fixed-expression redundancy analysis.
# Existing v15 results and the public v1.3.0 code snapshot are never overwritten.
suppressPackageStartupMessages({
  library(data.table);library(survival);library(GSVA);library(limma)
  library(jsonlite);library(digest)
})
options(stringsAsFactors=FALSE)
set.seed(20260831)
args<-commandArgs(trailingOnly=TRUE)
value<-function(flag,default=NULL){i<-match(flag,args);if(is.na(i))return(default);if(i==length(args))stop('Missing ',flag);args[[i+1L]]}
project<-normalizePath(value('--project','.'),mustWork=TRUE)
out<-value('--output-dir',file.path(project,'results','conditional_context_v16'))
plan<-value('--analysis-plan',file.path(project,'audit_v16','ANALYSIS_PLAN_V16.md'))
expected_plan<-value('--expected-plan-sha256')
if(is.null(expected_plan))stop('Supply the frozen --expected-plan-sha256 before fitting any model')
plan_sha<-digest(file=plan,algo='sha256')
if(!identical(plan_sha,expected_plan))stop('Frozen protocol SHA-256 mismatch; no models fitted')
amendment<-value('--analysis-plan-amendment',file.path(project,'audit_v16','ANALYSIS_PLAN_V16_AMENDMENT1.md'))
expected_amendment<-value('--expected-amendment-sha256')
if(is.null(expected_amendment))stop('Supply --expected-amendment-sha256 before fitting any model')
amendment_sha<-digest(file=amendment,algo='sha256')
if(!identical(amendment_sha,expected_amendment))stop('Frozen clarification SHA-256 mismatch; no models fitted')
dir.create(out,recursive=TRUE,showWarnings=FALSE)

input_records<-list();warnings_rows<-list();errors_rows<-list()
summary_rows<-list();coef_rows<-list();PH_rows<-list();definition_rows<-list()
design_rows<-list();VIF_rows<-list();predictor_rows<-list();eligibility_rows<-list()
LRT_rows<-list();scaling_rows<-list();score_rows<-list();coverage_rows<-list();redundancy_rows<-list()
add_warning<-function(cohort,endpoint,family,model,step){function(w){
  warnings_rows[[length(warnings_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=model,step=step,warning=conditionMessage(w))
  invokeRestart('muffleWarning')
}}
read_input<-function(relative){
  path<-file.path(project,relative)
  if(!file.exists(path))stop('Missing input: ',relative)
  input_records[[relative]]<<-data.table(path=relative,bytes=file.info(path)$size,sha256=digest(file=path,algo='sha256'))
  fread(path,check.names=FALSE,na.strings=c('NA',''))
}
standardize<-function(x,label,cohort){
  if(any(!is.finite(x)))stop('Non-finite full-expression-cohort scaling input: ',cohort,'/',label)
  s<-sd(x);if(!is.finite(s)||s<=0)stop('Zero/invalid sample SD: ',cohort,'/',label)
  scaling_rows[[length(scaling_rows)+1L]]<<-data.table(cohort=cohort,predictor=label,n=length(x),mean=mean(x),sample_SD=s,
    scope='full eligible expression cohort before any endpoint/clinical complete-case exclusion')
  (x-mean(x))/s
}
definition<-read_input('results/bulk/immune_signature_definitions.tsv')
IFNa19<-unique(as.character(definition[signature=='IFN_ALPHA_RESPONSE',gene]))
fixed_genes<-c('IFI44L','IFI44','MX1','ISG15')
IFNa16<-setdiff(IFNa19,fixed_genes)
if(length(IFNa19)!=19L||length(IFNa16)!=16L||'IFI44L'%in%IFNa19)stop('Expected fixed 19-gene IFI44L-excluded IFNa and its fixed 16-gene comparator-excluded subset')

configs<-list(
 TCGA=list(clinical='results/bulk/TCGA_TNBC_survival_merged.tsv',scores='results/bulk/immune_signature_ssGSEA_scores_wide.tsv',
   expression='source_metadata/TCGA_documented_primary/TCGA_BRCA_documented_primary_subset.tsv.gz',expression_field='IFI44L_expression',
   scale='supplied log2(RSEM normalized count + 1); no additional log',n=157L),
 GSE58812=list(clinical='results/validation_GSE58812/GSE58812_clinical_IFI44L.tsv',scores='results/validation_GSE58812/GSE58812_immune_signature_ssGSEA_scores.tsv',
   expression='results/validation_GSE58812/GSE58812_normalized_gene_expression.tsv.gz',expression_field='IFI44L_log2',
   scale='already log2(MAS5 intensity + 1), mapped/averaged by symbol; no further log',n=107L),
 GSE21653=list(clinical='results/validation_GSE21653/GSE21653_clinical_IFI44L.tsv',scores='results/validation_GSE21653/GSE21653_immune_signature_ssGSEA_scores.tsv',
   expression='results/validation_GSE21653/GSE21653_normalized_gene_expression.tsv.gz',expression_field='IFI44L_expression',
   scale='deposited RMA/quantile-normalized expression, mapped/averaged by symbol; no additional log',n=87L))
data_by_cohort<-list()
for(cohort in names(configs)){
  cfg<-configs[[cohort]];d<-read_input(cfg$clinical);s<-read_input(cfg$scores)
  if(nrow(d)!=cfg$n||anyDuplicated(d$sample_id))stop('Clinical expression-cohort identity mismatch: ',cohort)
  if(cohort=='TCGA'){
    score<-s[signature=='IFN_ALPHA_RESPONSE'];if(nrow(score)!=1L)stop('Non-unique TCGA IFNa score')
    ids<-setdiff(names(score),'signature');vals<-as.numeric(score[,..ids]);names(vals)<-ids
  }else{
    if(anyDuplicated(s$sample_id))stop('Duplicated score sample identifiers')
    vals<-s[['IFN_ALPHA_RESPONSE']];names(vals)<-s$sample_id
  }
  if(!setequal(names(vals),d$sample_id))stop('Clinical/module sample universe mismatch: ',cohort)
  d[,cohort:=cohort];d[,age10:=as.numeric(age_years)/10]
  d[,IFI44L_source_expression:=as.numeric(get(cfg$expression_field))]
  iz<-standardize(d$IFI44L_source_expression,'IFI44L',cohort)
  if(!isTRUE(all.equal(iz,as.numeric(d$IFI44L_z),tolerance=1e-10)))stop('Stored IFI44L scaling differs from full cohort: ',cohort)
  d[,IFI44L_z:=iz];d[,IFNa19_score:=as.numeric(vals[sample_id])]
  nz<-standardize(d$IFNa19_score,'original_IFNa19_score',cohort)
  d[,IFNa19_z:=nz]
  if(cohort=='TCGA')d[,stage_group:=factor(stage_group,levels=c('I','II','III'))]
  if(cohort=='GSE21653')d[,`:=`(pT_group=factor(pT_group,levels=c('T1','T2/3')),
    pN_group=factor(pN_group,levels=c('N0','N-positive')),grade_group=factor(grade_group,levels=c('G1/2','G3')))]
  fwrite(d,file.path(out,paste0(cohort,'_full_expression_cohort_inputs.tsv')),sep='\t',na='NA')
  data_by_cohort[[cohort]]<-d

  # The old 19-gene scores are read as supplied for Cox. Only the fixed 16-gene
  # expression-comparator score is newly calculated, using the full background.
  raw<-read_input(cfg$expression);gene<-as.character(raw[[1]])
  if(!setequal(names(raw)[-1L],d$sample_id))stop('Expression columns differ from clinical cohort: ',cohort)
  ids<-d$sample_id;expression<-as.matrix(raw[,..ids]);storage.mode(expression)<-'double';rownames(expression)<-gene
  expression<-avereps(expression,ID=rownames(expression))
  expression<-expression[rowSums(is.finite(expression))==ncol(expression),,drop=FALSE]
  expression<-expression[apply(expression,1L,sd)>0,,drop=FALSE]
  if(!all(fixed_genes%in%rownames(expression)))stop('A fixed comparator is unavailable/nonvariable: ',cohort)
  if(!isTRUE(all.equal(as.numeric(expression['IFI44L',]),d$IFI44L_source_expression,tolerance=1e-10)))stop('IFI44L gene matrix and clinical metadata differ: ',cohort)
  available<-intersect(IFNa16,rownames(expression));if(length(available)<2L)stop('Insufficient common16-gene IFNa coverage')
  coverage_rows[[cohort]]<-data.table(cohort=cohort,gene=IFNa19,in_original_IFNa19=TRUE,
    removed_comparator=IFNa19%in%fixed_genes,in_fixed_IFNa16=IFNa19%in%IFNa16,
    available_variable_gene=IFNa19%in%rownames(expression),used_in_IFNa16=IFNa19%in%available,
    expression_scale=cfg$scale,full_variable_background_genes=nrow(expression))
  m16<-withCallingHandlers(gsva(ssgseaParam(expression,list(IFNa16_excluding_IFI44L_IFI44_MX1_ISG15=available),
    minSize=2,maxSize=500,normalize=TRUE),verbose=FALSE),warning=add_warning(cohort,'expression','fixed_redundancy','IFNa16','GSVA'))
  m16<-as.numeric(m16[1L,ids]);z16<-standardize(m16,'fixed_IFNa16_score',cohort)
  score_rows[[cohort]]<-data.table(cohort=cohort,sample_id=ids,IFI44L=as.numeric(expression['IFI44L',]),
    IFI44=as.numeric(expression['IFI44',]),MX1=as.numeric(expression['MX1',]),ISG15=as.numeric(expression['ISG15',]),
    IFNa19_score=d$IFNa19_score,IFNa19_z=d$IFNa19_z,IFNa16_score=m16,IFNa16_z=z16)
  for(pair in combn(fixed_genes,2L,simplify=FALSE)){
    test<-suppressWarnings(cor.test(expression[pair[1],],expression[pair[2],],method='spearman',exact=FALSE))
    redundancy_rows[[length(redundancy_rows)+1L]]<-data.table(cohort=cohort,test_type='fixed_gene_pair',variable1=pair[1],variable2=pair[2],
      n=ncol(expression),spearman_rho=unname(test$estimate),p_value=test$p.value)
  }
  for(gene in fixed_genes){
    test<-suppressWarnings(cor.test(expression[gene,],m16,method='spearman',exact=FALSE))
    redundancy_rows[[length(redundancy_rows)+1L]]<-data.table(cohort=cohort,test_type='fixed_gene_vs_common_IFNa16',variable1=gene,
      variable2='IFNa16_excluding_IFI44L_IFI44_MX1_ISG15',n=ncol(expression),spearman_rho=unname(test$estimate),p_value=test$p.value)
  }
  predictor_rows[[length(predictor_rows)+1L]]<-data.table(cohort=cohort,endpoint='full_expression_cohort',family='diagnostic',variable1='IFI44L_z',variable2='IFNa19_z',
    n=nrow(d),pearson_r=cor(d$IFI44L_z,d$IFNa19_z),spearman_rho=cor(d$IFI44L_z,d$IFNa19_z,method='spearman'),role='descriptive collinearity diagnostic; no significance screening')
  rm(raw,expression);gc(verbose=FALSE)
}

design_diagnostics<-function(data,rhs,cohort,endpoint,family,model){
  X<-model.matrix(as.formula(paste('~',rhs)),data=data);X<-X[,colnames(X)!='(Intercept)',drop=FALSE]
  sds<-apply(X,2L,sd);valid<-all(is.finite(sds)&sds>0)
  if(valid){
    Z<-scale(X);singular_values<-svd(Z,nu=0,nv=0)$d
    tol<-max(dim(Z))*.Machine$double.eps*max(singular_values)
    rank<-sum(singular_values>tol);condition<-if(min(singular_values)>0)max(singular_values)/min(singular_values) else Inf
  }else{Z<-NULL;singular_values<-numeric();tol<-NA_real_;rank<-NA_integer_;condition<-Inf}
  design_rows[[length(design_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=model,n=nrow(data),
    covariate_columns=ncol(X),standardized_design_rank=rank,full_rank=valid&&!is.na(rank)&&rank==ncol(X),
    standardized_design_condition_number=condition,rank_tolerance=tol,
    singular_values=paste(format(singular_values,digits=17),collapse=';'),
    diagnostic_definition='center/sample-SD-standardize each covariate-design column within common complete cases; exclude intercept and stage strata')
  for(j in seq_len(ncol(X))){
    vif<-NA_real_
    if(valid){
      if(ncol(X)==1L)vif<-1 else{
        y<-Z[,j];other<-cbind(Intercept=1,Z[,-j,drop=FALSE]);residual<-qr.resid(qr(other),y)
        variance_fraction<-sum(residual^2)/sum((y-mean(y))^2);vif<-if(variance_fraction>0)1/variance_fraction else Inf
      }
    }
    VIF_rows[[length(VIF_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=model,term=colnames(X)[j],VIF=vif,
      role='ordinary standardized-design variance-inflation diagnostic; not a causal or predictive-performance metric')
  }
  list(X=X,rank_ok=valid&&!is.na(rank)&&rank==ncol(X),rank=rank,condition=condition)
}

run_scenario<-function(cohort,endpoint,time_col,event_col,time_unit,family,clinical_rhs='age10',strata_term=NULL){
  d<-copy(data_by_cohort[[cohort]])
  d[,`:=`(analysis_time=as.numeric(get(time_col)),analysis_event=as.numeric(get(event_col)))]
  full_design<-paste(clinical_rhs,'IFNa19_z','IFI44L_z',sep='+')
  needed<-unique(all.vars(as.formula(paste('~',full_design))))
  if(!is.null(strata_term))needed<-unique(c(needed,strata_term))
  reasons<-vapply(seq_len(nrow(d)),function(i){
    why<-character()
    if(!is.finite(d$analysis_time[i])||d$analysis_time[i]<=0)why<-c(why,'missing/nonpositive endpoint time')
    if(is.na(d$analysis_event[i])||!d$analysis_event[i]%in%c(0,1))why<-c(why,'missing/invalid endpoint event')
    for(name in needed){x<-d[[name]][i];if(is.na(x)||(is.numeric(x)&&!is.finite(x)))why<-c(why,paste0('missing/nonfinite ',name))}
    if(length(why))paste(unique(why),collapse=';') else 'included'
  },character(1))
  eligible<-reasons=='included';dat<-as.data.frame(d[eligible])
  inclusion<-copy(d);inclusion[,`:=`(endpoint=endpoint,family=family,time_unit=time_unit,eligible_common_complete_case=eligible,exclusion_reason=reasons)]
  eligibility_rows[[length(eligibility_rows)+1L]]<<-inclusion
  if(family=='primary'){
    models<-list(age=clinical_rhs,age_plus_IFI44L=paste(clinical_rhs,'IFI44L_z',sep='+'),
      age_plus_IFNa=paste(clinical_rhs,'IFNa19_z',sep='+'),age_plus_IFNa_plus_IFI44L=full_design)
    base_name<-'age_plus_IFNa';full_name<-'age_plus_IFNa_plus_IFI44L'
  }else{
    models<-list(clinical_base=clinical_rhs,clinical_base_plus_IFNa=paste(clinical_rhs,'IFNa19_z',sep='+'),
      clinical_base_plus_IFNa_plus_IFI44L=full_design)
    base_name<-'clinical_base_plus_IFNa';full_name<-'clinical_base_plus_IFNa_plus_IFI44L'
  }
  fits<-list()
  for(label in names(models)){
    rhs<-models[[label]];di<-design_diagnostics(dat,rhs,cohort,endpoint,family,label)
    rhs_fit<-if(is.null(strata_term))rhs else paste0(rhs,'+strata(',strata_term,')')
    formula<-as.formula(paste('Surv(analysis_time,analysis_event)~',rhs_fit))
    before<-length(warnings_rows);model_error<-NULL
    fit<-tryCatch({
      if(!di$rank_ok)stop('Standardized common-case design is rank deficient; no covariate dropped')
      withCallingHandlers(coxph(formula,data=dat,ties='efron',singular.ok=FALSE,x=TRUE,y=TRUE,model=TRUE),
        warning=add_warning(cohort,endpoint,family,label,'coxph'))
    },error=function(e){model_error<<-conditionMessage(e);NULL})
    if(is.null(fit)){
      errors_rows[[length(errors_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,step='coxph',error=model_error)
      summary_rows[[length(summary_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,status='not_estimable',
        n=nrow(dat),events=sum(dat$analysis_event),coefficient_count=ncol(di$X),events_per_coefficient=sum(dat$analysis_event)/ncol(di$X),log_likelihood=NA_real_)
      next
    }
    fits[[label]]<-fit;s<-summary(fit)
    summary_rows[[length(summary_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,status='fitted',
      n=fit$n,events=fit$nevent,coefficient_count=length(coef(fit)),events_per_coefficient=fit$nevent/length(coef(fit)),
      log_likelihood=unname(as.numeric(logLik(fit))),model_df=attr(logLik(fit),'df'),warning_count=length(warnings_rows)-before)
    coef_rows[[length(coef_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,
      term=rownames(s$coefficients),coefficient=s$coefficients[,'coef'],SE=s$coefficients[,'se(coef)'],
      hazard_ratio=s$conf.int[,'exp(coef)'],ci95_low=s$conf.int[,'lower .95'],ci95_high=s$conf.int[,'upper .95'],
      nominal_p_value=s$coefficients[,'Pr(>|z|)'],n=fit$n,events=fit$nevent)
    ph<-tryCatch(withCallingHandlers(cox.zph(fit,transform='km'),warning=add_warning(cohort,endpoint,family,label,'cox.zph')),
      error=function(e){errors_rows[[length(errors_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,step='cox.zph',error=conditionMessage(e));NULL})
    if(!is.null(ph))PH_rows[[length(PH_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,
      term=rownames(ph$table),chisq=ph$table[,'chisq'],df=ph$table[,'df'],nominal_p_value=ph$table[,'p'])
    summary_rows[[length(summary_rows)]]$warning_count<<-length(warnings_rows)-before
    definition_rows[[length(definition_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,model=label,
      formula=paste(deparse(formula),collapse=' '),time_source=time_col,event_source=event_col,time_unit=time_unit,
      n=fit$n,events=fit$nevent,coefficient_count=length(coef(fit)),
      scale='IFI44L and original19-gene IFNa standardized by full-expression-cohort sample SD; age/10',
      complete_cases='identical full-model cases used for every nested model within this scenario; endpoint time strictly>0; no imputation',
      categorical_reference='TCGA strataI/II/III; GSE21653 pT1, pN0 and grade1/2 when used',
      ties='Efron',PH='term/global Schoenfeld tests with km transform; non-significant global test does not certify all terms',
      interpretation='conditional association and internal nested partial-likelihood comparison only; sparse events and near-collinearity can limit precision; no validated prediction/causal claim')
  }
  Xfull<-model.matrix(as.formula(paste('~',full_design)),data=dat)
  Xfull<-Xfull[,colnames(Xfull)!='(Intercept)',drop=FALSE]
  if(ncol(Xfull)>1L)for(pair in combn(colnames(Xfull),2L,simplify=FALSE))predictor_rows[[length(predictor_rows)+1L]]<<-data.table(
    cohort=cohort,endpoint=endpoint,family=family,variable1=pair[1],variable2=pair[2],n=nrow(dat),
    pearson_r=cor(Xfull[,pair[1]],Xfull[,pair[2]]),spearman_rho=cor(Xfull[,pair[1]],Xfull[,pair[2]],method='spearman'),
    role='descriptive collinearity diagnostic; no significance screening')
  base<-fits[[base_name]];full<-fits[[full_name]]
  if(!is.null(base)&&!is.null(full)){
    ll0<-unname(as.numeric(logLik(base)));ll1<-unname(as.numeric(logLik(full)))
    df<-attr(logLik(full),'df')-attr(logLik(base),'df');raw_stat<-2*(ll1-ll0)
    if(df!=1L||base$n!=full$n||base$nevent!=full$nevent||raw_stat< -1e-8)stop('Invalid nested one-df likelihood-ratio comparison')
    statistic<-max(0,raw_stat)
    LRT_rows[[length(LRT_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,status='estimable',
      base_model=base_name,full_model=full_name,n=full$n,events=full$nevent,base_log_likelihood=ll0,full_log_likelihood=ll1,
      raw_LRT_statistic=raw_stat,LRT_statistic=statistic,df=df,p_value=pchisq(statistic,df=df,lower.tail=FALSE),
      interpretation='within-cohort conditional association check; not external predictive validation')
  }else LRT_rows[[length(LRT_rows)+1L]]<<-data.table(cohort=cohort,endpoint=endpoint,family=family,status='not_estimable',
    base_model=base_name,full_model=full_name,n=nrow(dat),events=sum(dat$analysis_event),p_value=NA_real_)
}

run_scenario('TCGA','OS','OS_days','OS_event','days','primary')
run_scenario('TCGA','PFI','PFI_days','PFI_event','days','primary')
run_scenario('GSE58812','MFS','mfs_days','metastasis_event','days','primary')
run_scenario('GSE58812','OS','os_days','death_event','days','primary')
run_scenario('GSE21653','DFS','dfs_months','dfs_event','months','primary')
run_scenario('TCGA','OS','OS_days','OS_event','days','secondary',strata_term='stage_group')
run_scenario('TCGA','PFI','PFI_days','PFI_event','days','secondary',strata_term='stage_group')
run_scenario('GSE21653','DFS','dfs_months','dfs_event','months','secondary',clinical_rhs='age10+pT_group+pN_group+grade_group')

LRT<-rbindlist(LRT_rows,fill=TRUE)
if(sum(LRT$family=='primary')!=5L||sum(LRT$family=='secondary')!=3L)stop('Planned LRT family count changed')
LRT[,BH_family_size:=ifelse(family=='primary',5L,3L)]
LRT[,fdr_bh:=NA_real_]
for(fam in c('primary','secondary')){
  idx<-which(LRT$family==fam&is.finite(LRT$p_value))
  LRT$fdr_bh[idx]<-p.adjust(LRT$p_value[idx],method='BH',n=if(fam=='primary')5L else 3L)
}
redundancy<-rbindlist(redundancy_rows)
if(any(redundancy[,.N,by=cohort]$N!=10L))stop('Fixed10-expression-test family changed')
redundancy[,`:=`(fdr_bh=p.adjust(p_value,method='BH',n=10L),BH_family_size=10L),by=cohort]
summary<-rbindlist(summary_rows,fill=TRUE);coefficients<-rbindlist(coef_rows,fill=TRUE)
eligibility<-rbindlist(eligibility_rows,fill=TRUE)
empty_warnings<-data.table(cohort=character(),endpoint=character(),family=character(),model=character(),step=character(),warning=character())
empty_errors<-data.table(cohort=character(),endpoint=character(),family=character(),model=character(),step=character(),error=character())
files<-list(Cox_model_summary=summary,Cox_all_coefficients=coefficients,Cox_PH_diagnostics=rbindlist(PH_rows,fill=TRUE),
  Cox_model_definitions=rbindlist(definition_rows,fill=TRUE),Cox_design_diagnostics=rbindlist(design_rows,fill=TRUE),
  Cox_VIF_diagnostics=rbindlist(VIF_rows,fill=TRUE),predictor_correlation_diagnostics=rbindlist(predictor_rows,fill=TRUE),
  complete_case_and_exclusion_records=eligibility,excluded_records=eligibility[eligible_common_complete_case==FALSE],
  expression_cohort_standardization=rbindlist(scaling_rows),IFNa16_gene_coverage=rbindlist(coverage_rows),
  fixed_gene_and_module_scores=rbindlist(score_rows),fixed_expression_redundancy10=redundancy,
  primary5_conditional_LRT=LRT[family=='primary'],secondary3_conditional_LRT=LRT[family=='secondary'],
  warnings=if(length(warnings_rows))rbindlist(warnings_rows) else empty_warnings,
  errors=if(length(errors_rows))rbindlist(errors_rows) else empty_errors,input_manifest=rbindlist(input_records))
for(name in names(files))fwrite(files[[name]],file.path(out,paste0(name,'.tsv')),sep='\t',na='NA')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
receipt<-list(status='ACTUAL_MODELS_EXECUTED_NOT_INDEPENDENTLY_VALIDATED',protocol_file=plan,protocol_sha256=plan_sha,
  amendment_file=amendment,amendment_sha256=amendment_sha,
  protocol_role='local revision decision frozen before these models; not public preregistration',full_expression_cohort_n=list(TCGA=157L,GSE58812=107L,GSE21653=87L),
  fitted_models=sum(summary$status=='fitted'),planned_models=29L,primary_LRT_n=5L,secondary_LRT_n=3L,
  fixed_expression_comparators=c('IFI44','MX1','ISG15'),redundancy_primary_method='Spearman; ten fixed tests per cohort with BH10',
  Cox_IFNa_module='existing original19-gene IFI44L-excluded score; not rescored',
  redundancy_IFNa_module='fixed common16-gene list after simultaneously excluding IFI44L and the three comparators; ssGSEA full variable-gene background; normalize=TRUE',
  old_outputs_overwritten=FALSE,public_v1_3_0_code_modified=FALSE,warnings_n=nrow(files$warnings),errors_n=nrow(files$errors),
  limits=c('High predictor correlation and sparse events can limit conditional precision.','A nonsignificant test is not proof of no information.','LRT improvement is not external validated predictive performance.','No C-index/AUC, optimized score/cutpoint, interaction or result-selected subgroup was calculated.'))
write_json(receipt,file.path(out,'EXECUTION_RECEIPT.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
output_files<-list.files(out,full.names=TRUE);output_files<-output_files[file.info(output_files)$isdir==FALSE]
output_files<-output_files[basename(output_files)!='OUTPUT_MANIFEST_SHA256.tsv']
fwrite(data.table(file=basename(output_files),bytes=file.info(output_files)$size,
  sha256=vapply(output_files,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'OUTPUT_MANIFEST_SHA256.tsv'),sep='\t')
print(LRT);message('v16 execution completed; results require independent source/numerical review before interpretation.')
