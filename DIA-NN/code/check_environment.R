source("CoPPRA_config.R")
packages <- c("data.table", "knitr", "jsonlite", "R.utils", "matrixStats", "ggplot2", "pheatmap", "gridExtra", "PCAtools", "ggrepel", "dplyr", "tidyr", "limma", "ComplexUpset", "PEIMAN2", "fgsea", "openxlsx", "readr", "AnnotationDbi", "clusterProfiler", "org.Hs.eg.db", "GO.db", "topGO", "enrichplot", "ComplexHeatmap", "circlize", "ragg", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Missing R packages: ",paste(missing,collapse=", "))
versions <- data.table(Package=packages, Version=vapply(packages,function(x) as.character(packageVersion(x)), character(1)))
methods_dir <- file.path(output_root, "Methods")
dir.create(methods_dir, recursive=TRUE, showWarnings=FALSE)
fwrite(versions,file.path(methods_dir,"package_versions.csv"))
if (as.character(packageVersion("limma")) != "3.68.4") {
  stop("This workflow was validated with limma 3.68.4. Install that version in ",
       "code/R_library or set COPPRA_R_LIBRARY. See ../limma_source.json.")
}
print(versions)
