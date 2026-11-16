source("CoPPRA_config.R")
source("CoPPRA_heatmap_12_samples.R")
run_coppra_heatmaps(file.path(output_root,"QC","data","quantification.rds"),file.path(output_root,"Heatmap"),analysis_label)
