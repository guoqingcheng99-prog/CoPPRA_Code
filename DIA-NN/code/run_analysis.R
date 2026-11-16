args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
setwd(dirname(normalizePath(script)))
if (!nzchar(Sys.getenv("COPPRA_PACKAGE"))) Sys.setenv(COPPRA_PACKAGE = dirname(getwd()))
selection <- if (length(args)) args[1] else "all"
stages <- c(
  "1.1" = "1.1 CoPPRA_QC.rmd",
  "1.2" = "1.2 CoPPRA_QCplot.rmd",
  "qc-table" = "01_qc_table.R",
  "2.1" = "2.1 CoPPRA_PCA.rmd",
  "2.2" = "2.2 CoPPRA_PCAplot.rmd",
  "heatmaps" = "02_heatmaps.R",
  "3.1" = "3.1 CoPPRA_Limma.rmd",
  "3.2" = "3.2 CoPPRA_Limmaplot.rmd",
  "5.1" = "5.1 CoPPRA_PEIMAN_PSEA.rmd",
  "5.2" = "5.2 CoPPRA_PEIMAN_SEA.rmd",
  "7.1" = "7.1 CoPPRA_TopGO.rmd",
  "7.2" = "7.2 CoPPRA_TopGOplot.rmd",
  "permutation" = "08_global_permutation.R",
  "validate" = "10_validate.R"
)
if (selection %in% c("--help", "-h")) {
  cat("Usage: Rscript run_analysis.R [all|STAGE]\n",
      "Stages: ", paste(names(stages), collapse = ", "), "\n", sep = "")
  quit(status = 0L)
}
if (length(args) > 1L || !(selection %in% c("all", names(stages)))) {
  stop("Unknown stage or extra arguments. Use --help to list stages.")
}
source("CoPPRA_config.R")
files <- if (selection == "all") stages else stages[selection]
logdir <- file.path(output_root, "logs")
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
for (f in files) {
  message("START ", f, " ", Sys.time())
  dest <- if (grepl("rmd$", f)) file.path(logdir,
    paste0(tools::file_path_sans_ext(f), ".R")) else f
  if (grepl("rmd$", f)) knitr::purl(f, output = dest, documentation = 0, quiet = TRUE)
  grDevices::pdf(NULL)
  e <- new.env(parent = globalenv())
  stage_warnings <- character()
  tryCatch(withCallingHandlers(source(dest, local = e, echo = FALSE), warning = function(w) {
    stage_warnings <<- c(stage_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), finally = {
    writeLines(unique(stage_warnings), file.path(logdir,
      paste0(tools::file_path_sans_ext(basename(f)), "_warnings.txt")))
    grDevices::graphics.off()
  })
  cat(f, as.character(Sys.time()), sep = "\t",
      file = file.path(logdir, "completed_stages.tsv"), append = TRUE)
  cat("\n", file = file.path(logdir, "completed_stages.tsv"), append = TRUE)
  rm(e)
  gc()
  message("DONE ", f, " ", Sys.time())
}
