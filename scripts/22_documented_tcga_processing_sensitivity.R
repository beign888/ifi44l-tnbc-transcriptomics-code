#!/usr/bin/env Rscript
# Processing sensitivity only: the same clinically selected patients, not an
# independent patient cohort. No outcome modelling or GSEA is performed here.
suppressPackageStartupMessages({library(data.table);library(limma);library(GSVA)})

args <- commandArgs(trailingOnly=TRUE)
value <- function(flag, default=NULL) {
  index <- match(flag,args)
  if (is.na(index)) return(default)
  if (index==length(args)) stop('Missing value after ',flag)
  args[[index+1L]]
}
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub('^--file=','',script_arg[[1L]]))) else getwd()
work <- normalizePath(value('--work-dir',script_dir),mustWork=TRUE)
out <- value('--output-dir',work)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
default_input <- file.path(work,'documented_HiSeqV2_matched151_expression.tsv.gz')
if(!file.exists(default_input)) default_input <- file.path(work,'TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz')
matrix_file <- value('--input',value('--matrix',default_input))
baseline_project <- value('--baseline-project',value('--project',NULL))
baseline_files <- if(is.null(baseline_project)) {
  list(selected=file.path(work,'baseline_TCGA151_selected_samples.txt'),
       modules=file.path(work,'baseline_module_definitions.tsv'),
       ifi=file.path(work,'baseline_TCGA_IFI44L.tsv'),
       scores=file.path(work,'baseline_TCGA_ssGSEA_scores_wide.tsv'))
} else {
  bulk <- file.path(normalizePath(baseline_project,mustWork=TRUE),'results','bulk')
  list(selected=file.path(bulk,'TCGA_TNBC_selected_samples.txt'),
       modules=file.path(bulk,'immune_signature_definitions.tsv'),
       ifi=file.path(bulk,'tnbc_sample_metadata.tsv'),
       scores=file.path(bulk,'immune_signature_ssGSEA_scores_wide.tsv'))
}
stopifnot(file.exists(file.path(work,'PROCESSING_SENSITIVITY_DESIGN.json')))
requested <- unique(readLines(baseline_files$selected))
stopifnot(length(requested)==151L, all(substr(requested,14,15)=='01'))
header <- names(fread(matrix_file,nrows=0,check.names=FALSE))
available <- header[-1L]
mapping <- rbindlist(lapply(requested,function(sample) {
  exact <- which(available==sample)
  candidates <- if(length(exact)) exact else which(substr(available,1,15)==sample)
  data.table(sample_id=sample,source_sample_id=if(length(candidates)==1L) available[candidates] else NA_character_,
             mapping_rule=if(length(exact)==1L) 'exact barcode' else if(length(candidates)==1L) 'unique first-15 primary-tumour barcode' else 'missing or ambiguous; excluded',
             candidate_n=length(candidates),included=length(candidates)==1L)
}))
fwrite(mapping,file.path(out,'HiSeqV2_sample_mapping.tsv'),sep='\t')
fwrite(mapping[included==FALSE],file.path(out,'HiSeqV2_missing_or_ambiguous_samples.tsv'),sep='\t')
included <- mapping[included==TRUE]
fwrite(data.table(source_sample_id=setdiff(available,included$source_sample_id)),file.path(out,'HiSeqV2_source_samples_not_selected.tsv'),sep='\t')
if(nrow(included)<50L) stop('Insufficient matched samples for the bounded processing check')
stopifnot(!anyDuplicated(substr(included$sample_id,1,12)))

message('Reading documented RSEM matrix for ',nrow(included),' frozen clinical patients')
input <- fread(matrix_file,select=c(header[1L],included$source_sample_id),check.names=FALSE)
expr <- as.matrix(input[,-1L,with=FALSE]);storage.mode(expr)<-'double'
rownames(expr)<-as.character(input[[1L]])
colnames(expr)<-included$sample_id[match(colnames(expr),included$source_sample_id)]
expr<-expr[,included$sample_id,drop=FALSE]
gene_rows_input <- nrow(expr);duplicate_symbol_rows<-sum(duplicated(rownames(expr)))
expr<-avereps(expr,ID=rownames(expr))
nonfinite_rows<-sum(rowSums(is.finite(expr))!=ncol(expr))
expr<-expr[rowSums(is.finite(expr))==ncol(expr),,drop=FALSE]
invariant_rows<-sum(apply(expr,1,sd)==0)
expr<-expr[apply(expr,1,sd)>0,,drop=FALSE]
stopifnot('IFI44L'%in%rownames(expr))

definition<-fread(baseline_files$modules)
signatures<-lapply(split(as.character(definition$gene),definition$signature),unique)
signatures<-lapply(signatures,function(genes)setdiff(genes,'IFI44L'))
members<-rbindlist(lapply(names(signatures),function(name)
  data.table(signature=name,gene=signatures[[name]],present_in_documented_matrix=signatures[[name]]%in%rownames(expr))))
fwrite(members,file.path(out,'HiSeqV2_module_gene_coverage.tsv'),sep='\t')
signatures_present<-lapply(signatures,function(genes)intersect(genes,rownames(expr)))
if(any(lengths(signatures_present)<2L)) stop('An original module has fewer than two mapped genes')
stopifnot(length(signatures_present)==13L)

# This is the identical ssGSEA specification used in primary script 01.
scores<-gsva(ssgseaParam(expr,signatures_present,minSize=2,maxSize=500,normalize=TRUE),verbose=FALSE)
fwrite(data.table(signature=rownames(scores),as.data.table(scores)),file.path(out,'HiSeqV2_ssGSEA_scores_wide.tsv'),sep='\t')
ifi<-as.numeric(expr['IFI44L',]);names(ifi)<-colnames(expr)
fwrite(data.table(sample_id=names(ifi),IFI44L_log2_RSEM_normalized_count_plus1=ifi),file.path(out,'HiSeqV2_IFI44L_matched151.tsv'),sep='\t')
correlations<-function(x,score_matrix) {
  output<-rbindlist(lapply(rownames(score_matrix),function(name) {
    result<-suppressWarnings(cor.test(x,score_matrix[name,],method='spearman',exact=FALSE))
    data.table(signature=name,n=length(x),spearman_rho=unname(result$estimate),p_value=result$p.value)
  }))
  output[,fdr_bh:=p.adjust(p_value,'BH')]
  output
}
documented<-correlations(ifi,scores)
fwrite(documented,file.path(out,'HiSeqV2_IFI44L_module_correlations.tsv'),sep='\t')
old_metadata<-fread(baseline_files$ifi)
old_ifi<-old_metadata$IFI44L_expression[match(included$sample_id,old_metadata$sample_id)]
stopifnot(all(is.finite(old_ifi)))
old_wide<-fread(baseline_files$scores,check.names=FALSE)
old_scores<-as.matrix(old_wide[,-1L,with=FALSE]);rownames(old_scores)<-old_wide[[1L]]
old_scores<-old_scores[,included$sample_id,drop=FALSE]
archived<-correlations(old_ifi,old_scores)
fwrite(archived,file.path(out,'Archived_matched151_IFI44L_module_correlations.tsv'),sep='\t')
combined<-merge(archived,documented,by='signature',suffixes=c('_archived','_documented_RSEM'))
combined[,same_direction:=sign(spearman_rho_archived)==sign(spearman_rho_documented_RSEM)]
fwrite(combined,file.path(out,'Processing_sensitivity_module_comparison.tsv'),sep='\t')
rank_test<-suppressWarnings(cor.test(old_ifi,ifi,method='spearman',exact=FALSE))
fwrite(data.table(comparison='IFI44L ranks: archived units versus documented log2(RSEM normalized count + 1)',
                  n=length(ifi),spearman_rho=unname(rank_test$estimate),p_value=rank_test$p.value),
       file.path(out,'IFI44L_cross_processing_rank_concordance.tsv'),sep='\t')
audit<-data.table(metric=c('requested_clinical_patients','matched_primary_tumour_samples','missing_or_ambiguous_samples',
                          'input_gene_rows','duplicated_symbol_rows','nonfinite_rows_removed','invariant_rows_removed',
                          'retained_variable_gene_rows','module_count','IFI44L_in_any_module','additional_expression_transform',
                          'documented_expression_units','prognostic_models_run','enrichment_models_run','independent_patient_validation'),
                  value=as.character(c(length(requested),nrow(included),sum(!mapping$included),gene_rows_input,
                                       duplicate_symbol_rows,nonfinite_rows,invariant_rows,nrow(expr),length(signatures_present),
                                       any(vapply(signatures,function(genes)'IFI44L'%in%genes,logical(1))),
                                       'none; use supplied documented log2 matrix','log2(RSEM normalized count + 1)',
                                       FALSE,FALSE,FALSE)))
fwrite(audit,file.path(out,'PROCESSING_SENSITIVITY_MODEL_AUDIT.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'PROCESSING_SENSITIVITY_R_sessionInfo.txt'))
message('Bounded processing sensitivity complete; all ',length(signatures_present),' original modules retained.')
print(combined[signature%in%c('IFN_ALPHA_RESPONSE','EPITHELIAL_CONTENT')])
print(data.table(IFI44L_cross_processing_rho=unname(rank_test$estimate),n=length(ifi)))
