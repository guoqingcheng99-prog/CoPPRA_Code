# Joint sample-label permutations preserve peptide correlation and n=3/group.
source('CoPPRA_config.R')
library(limma)
library(matrixStats)
library(ggplot2)
paths <- stage_paths('Permutation')
tab <- fread(file.path(output_root,'Limma','tables','peptide_limma_matrix_used.csv'))
y <- as.matrix(tab[, -1]);rownames(y)<-tab$Sequence
meta <- fread(file.path(output_root,'Limma','tables','differential_sample_metadata.csv'))
meta <- meta[match(colnames(y),Experiment)]
g <- factor(meta$Condition,levels=c('Control','Ru','Ul','RU'))
design <- model.matrix(~0+g);colnames(design)<-levels(g)
stopifnot(all(table(g)==3),ncol(y)==12,all(is.finite(y)))
sumsq <- rowSums(y*y)
calculate <- function(d,details=FALSE){
 means <- y %*% d / 3
 s2 <- pmax(0,(sumsq - 3*rowSums(means*means))/8)
 sv <- limma::squeezeVar(s2,df=8)
 df <- pmin(8+sv$df.prior,nrow(y)*8)
 effects <- means[,2:4,drop=FALSE] - means[,1]
 p <- 2*pt(-abs(effects)/sqrt(sv$var.post*2/3),df=df)
 colnames(p)<-colnames(effects)<-c('R','U','RU')
 counts <- colSums(p<0.05 & abs(effects)>1)
 if(details) return(list(effects=effects,p=p,counts=counts))
 counts
}
observed <- calculate(design,TRUE)
saved <- fread(file.path(output_root,'Limma','tables','all_peptide_differential_results.csv'))
saved <- saved[match(rownames(y),Sequence)]
for(i in seq_along(c('Ru','Ul','RU'))){
 nm <- c('Ru','Ul','RU')[i]
 stopifnot(max(abs(observed$effects[,i]-saved[[paste0(nm,'_vs_Control_limma_logFC')]]))<1e-8,
 max(abs(observed$p[,i]-saved[[paste0(nm,'_vs_Control_pvalue')]]))<1e-8)
}
set.seed(20261029)
B <- as.integer(Sys.getenv('COPPRA_PERMUTATIONS','10000'))
stopifnot(length(B) == 1L, !is.na(B), B > 0L)
permutations <- replicate(B,sample.int(12L))
saveRDS(permutations,file.path(paths$data,'joint_sample_label_permutations.rds'))
chunks <- split(seq_len(B),ceiling(seq_len(B)/100))
res <- parallel::mclapply(seq_along(chunks),function(k){
 ans<-vapply(chunks[[k]],function(i) calculate(design[permutations[,i],,drop=FALSE]),numeric(3))
 message(analysis_label,' permutation chunk ',k,'/',length(chunks))
 ans
},mc.cores=if (.Platform$OS.type == "windows") 1L else 2L,mc.preschedule=FALSE,mc.set.seed=FALSE)
null <- t(do.call(cbind,res));colnames(null)<-c('R','U','RU')
fwrite(cbind(data.table(Permutation=seq_len(B)),as.data.table(null)),file.path(paths$tables,'null_responsive_peptide_counts.csv'))
summary <- rbindlist(lapply(seq_len(3),function(i) data.table(Condition=colnames(null)[i],
 Observed=observed$counts[i],Null_mean=mean(null[,i]),Null_median=median(null[,i]),
 Null_95_low=quantile(null[,i],.025),Null_95_high=quantile(null[,i],.975),
 Empirical_P=(1+sum(null[,i]>=observed$counts[i]))/(B+1),Permutations=B)))
fwrite(summary,file.path(paths$tables,'global_permutation_summary.csv'))
long<-melt(as.data.table(null),measure.vars=c('R','U','RU'),variable.name='Condition',value.name='Responsive_peptides')
p<-ggplot(long,aes(Responsive_peptides,fill=Condition))+geom_histogram(bins=40,color='white')+
 geom_vline(data=summary,aes(xintercept=Observed),linetype=2,linewidth=.8)+
 facet_wrap(~Condition,scales='free_x')+scale_fill_manual(values=condition_colors)+
 labs(title='Joint sample-label permutation test',x='Responsive peptide count',y='Permutations')+
 theme_bw(base_size=14)+theme(legend.position='none')
save_plot(p,'global_permutation_null',12,4,folder=paths$plots)
writeLines(c(sprintf('%d joint permutations across 12 biological samples; three per group.', B),
 'Primary observed-data eligibility and peptide-median imputation fixed before permutation.',
 'Group means, residual variances and empirical-Bayes prior refitted for every permutation.',
 'Fast balanced-design calculation verified against saved limma effects and P values (<1e-8).',
 'Count rule: nominal P<0.05 and |limma log2FC|>1; empirical P=(1+null>=observed)/(B+1).'),
 file.path(paths$root,'method.txt'))
save_session(paths$root)
