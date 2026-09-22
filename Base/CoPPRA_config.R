# Shared settings. Keep this file beside all eight R Markdown files.
# Inputs: evidence_new.txt and proteinGroups_base.txt in the Base folder.
project_dir <- path.expand(Sys.getenv(
  "COPPRA_PROJECT_DIR", unset = "~/Desktop/Jafari_Lab/CoPPRA/Base"
))
output_root <- path.expand(Sys.getenv("COPPRA_BASE_OUTPUT",
  unset = Sys.getenv("COPPRA_OUTPUT_ROOT", unset = file.path(project_dir, "Results"))))
evidence_file <- file.path(project_dir, "evidence_new.txt")
protein_groups_file <- path.expand(Sys.getenv("COPPRA_BASE_PROTEIN_GROUPS",
  unset = file.path(project_dir, "proteinGroups_base.txt")))

library(data.table)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", error = FALSE, dpi = 600, fig.retina = 1)
report_input <- knitr::current_input()
if (!length(report_input) || !nzchar(report_input)) report_input <- "CoPPRA"
knitr::opts_chunk$set(fig.path = file.path(output_root, "reports", "figures",
  paste0(tools::file_path_sans_ext(basename(report_input)), "-")))

stage_paths <- function(stage) {
  root <- file.path(output_root, stage)
  paths <- list(root = root, tables = file.path(root, "tables"),
                plots = file.path(root, "plots"), data = file.path(root, "data"))
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

write_tables <- function(tables, folder) {
  for (name in names(tables)) {
    fwrite(tables[[name]], file.path(folder, paste0(name, ".csv")))
  }
}

matrix_table <- function(x) {
  cbind(data.table(Sequence = rownames(x)), as.data.table(x))
}

to_matrix <- function(data, column, value, column_order) {
  wide <- dcast(data, as.formula(paste("Sequence ~", column)),
                value.var = value, fill = NA_real_)
  # Preserve metadata columns even if a run has no quantified evidence rows.
  for (name in setdiff(column_order, names(wide))) wide[, (name) := NA_real_]
  mat <- as.matrix(wide[, ..column_order])
  rownames(mat) <- wide$Sequence
  mat
}

save_plot <- function(plot, name, width, height, folder = plot_folder, dpi = 600) {
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(folder, paste0(name, ".", ext)), plot,
                    width = width, height = height, dpi = dpi, bg = "white")
  }
}

save_session <- function(folder) {
  writeLines(capture.output(sessionInfo()), file.path(folder, "sessionInfo.txt"))
}

condition_colors <- c(R = "#D62728", U = "#2CA02C", RU = "#FF7F0E",
                       Ctrl = "#1F77B4", Carrier = "#9467BD")
