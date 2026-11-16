source('CoPPRA_config.R')
paths<-stage_paths('QC')
x<-fread(file.path(paths$tables,'QC_stepwise_peptide_protein_counts.csv'))
tab<-x[,.(Step=Step+1L,Stage,Peptides=format(Unique_peptides,big.mark=','),
 `Protein groups`=format(Protein_groups,big.mark=','),Accessions=format(Protein_accessions,big.mark=','))]
tab[,Stage:=vapply(Stage,function(z)paste(strwrap(z,width=50),collapse='\n'),character(1))]
g<-gridExtra::tableGrob(tab,rows=NULL,theme=gridExtra::ttheme_minimal(base_size=10,
 core=list(fg_params=list(hjust=0,x=.03),bg_params=list(fill=rep(c('white','#EEF4FA'),length.out=nrow(tab)),col='grey85')),
 colhead=list(fg_params=list(col='white',fontface='bold'),bg_params=list(fill='#1F77B4',col=NA))))
plot<-gridExtra::arrangeGrob(grid::textGrob(paste(analysis_label,'- stepwise peptide and protein retention'),gp=grid::gpar(fontsize=16,fontface='bold')),g,
 grid::textGrob('Peptides = distinct unmodified sequences. Proteins = inferred groups and member accessions.\nCounts reflect observed identifications; imputation adds no peptides or proteins.',gp=grid::gpar(fontsize=9)),ncol=1,heights=c(.07,.86,.07))
save_plot(plot,'QC_stepwise_counts_table',13,10,folder=paths$plots)
