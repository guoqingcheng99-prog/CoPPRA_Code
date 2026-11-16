source('CoPPRA_config.R')
library(matrixStats)
out <- output_root
 q<-readRDS(file.path(out,'QC','data','quantification.rds'))
 counts<-fread(file.path(out,'QC','tables','QC_stepwise_peptide_protein_counts.csv'))
 detection<-fread(file.path(out,'QC','tables','QC_primary_biological_detection.csv'))
 center<-fread(file.path(out,'QC','tables','QC_centering_validation.csv'))
 protein_idx<-which(grepl('Proteins supported',counts$Stage));missing_idx<-which(grepl('Peptides detected in >=20%',counts$Stage,fixed=TRUE))
 stopifnot(protein_idx<missing_idx,all(diff(counts$Unique_peptides)<=0),all(diff(counts$Protein_accessions)<=0),
  nrow(q$runs)==45,nrow(q$sample_map)==15,all(q$runs[, .N,by=Experiment]$N==3),
  max(abs(center$Median_normalized_common))<1e-10,
  setequal(detection[Retained==TRUE,Sequence],q$primary_sequences),
  all(detection[Retained==TRUE,R]>=2),all(detection[Retained==TRUE,U]>=2),
  all(detection[Retained==TRUE,RU]>=2),all(detection[Retained==TRUE,Ctrl]>=2),
  counts[grepl('Exclude references',Stage),Runs]==36)
 x<-fread(file.path(out,'Limma','tables','peptide_limma_matrix_used.csv'))
 lim<-fread(file.path(out,'Limma','tables','all_peptide_differential_results.csv'))
 m<-as.matrix(x[,-1]);rownames(m)<-x$Sequence
 observed<-q$biological_log2[x$Sequence,colnames(m),drop=FALSE]
 missing<-!is.finite(observed)
 med<-rowMedians(observed,na.rm=TRUE)
 stopifnot(ncol(m)==12,all(is.finite(m)),!any(grepl('Carrier',colnames(m))),
  setequal(x$Sequence,q$primary_sequences),all(abs(m[!missing]-observed[!missing])<1e-10),
  all(abs(m[missing]-med[row(m)[missing]])<1e-10))
 for(g in c('Ru','Ul','RU')){
  expected<-lim[[paste0(g,'_vs_Control_pvalue')]]<0.05 & abs(lim[[paste0(g,'_vs_Control_limma_logFC')]])>1
  stopifnot(identical(lim[[paste0(g,'_sig')]],expected),
    max(abs(lim[[paste0(g,'_log2FC')]]-lim[[paste0(g,'_vs_Control_limma_logFC')]]))<1e-12)
 }
 prot<-fread(file.path(out,'Limma','tables','significant_proteins.csv'))
 ct<-fread(file.path(out,'Limma','tables','conjunctional_targeting_proteins.csv'))
 stopifnot(setequal(ct$ProteinID,prot[RU_sig & !Ru_sig & !Ul_sig,ProteinID]))
 checks <- data.table(Analysis=analysis_label,Raw_files=45,Biological_samples=12,Primary_peptides=nrow(m),Imputed_values=sum(missing),
  Normalization_max_abs_error=max(abs(center$Median_normalized_common)),QC_filter_order='PASS',
  Observed_value_preservation='PASS',Imputation='PASS',Limma_call_consistency='PASS',CT_definition='PASS')
fwrite(checks,file.path(output_root,'validation_summary.csv'))
print(checks)
