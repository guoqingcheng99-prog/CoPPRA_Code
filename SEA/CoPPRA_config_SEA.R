# Shared settings. Keep this file beside all SEA R Markdown files.
# Inputs: evidence_SEA_v2.txt and proteinGroups_SEA_v2.txt in the SEA folder.
project_dir_SEA <- path.expand(Sys.getenv(
  "COPPRA_SEA_DIR", unset = "~/Desktop/Jafari_Lab/CoPPRA/SEA"
))
output_root_SEA <- path.expand(Sys.getenv("COPPRA_SEA_OUTPUT", unset = file.path(project_dir_SEA, "Results_v2")))
evidence_file_SEA <- file.path(project_dir_SEA, "evidence_SEA_v2.txt")
protein_groups_file_SEA <- path.expand(Sys.getenv("COPPRA_SEA_PROTEIN_GROUPS",
  unset = file.path(project_dir_SEA, "proteinGroups_SEA_v2.txt")))

library(data.table)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", error = FALSE, dpi = 600, fig.retina = 1)
report_input_SEA <- knitr::current_input()
if (!length(report_input_SEA) || !nzchar(report_input_SEA)) report_input_SEA <- "CoPPRA_SEA"
knitr::opts_chunk$set(fig.path = file.path(output_root_SEA, "reports", "figures",
  paste0(tools::file_path_sans_ext(basename(report_input_SEA)), "-")))

stage_paths_SEA <- function(stage) {
  root <- file.path(output_root_SEA, stage)
  paths_SEA <- list(root = root, tables = file.path(root, "tables"),
                plots = file.path(root, "plots"), data = file.path(root, "data"))
  invisible(lapply(paths_SEA, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths_SEA
}

write_tables_SEA <- function(tables, folder) {
  for (name in names(tables)) {
    fwrite(tables[[name]], file.path(folder, paste0(name, "_SEA.csv")))
  }
}

matrix_table_SEA <- function(x) {
  cbind(data.table(Sequence = rownames(x)), as.data.table(x))
}

to_matrix_SEA <- function(data, column, value, column_order) {
  wide <- dcast(data, as.formula(paste("Sequence ~", column)),
                value.var = value, fill = NA_real_)
  # Preserve metadata columns even if a run has no quantified evidence rows.
  for (name in setdiff(column_order, names(wide))) wide[, (name) := NA_real_]
  mat <- as.matrix(wide[, ..column_order])
  rownames(mat) <- wide$Sequence
  mat
}

save_plot_SEA <- function(plot, name, width, height, folder = plot_folder_SEA, dpi = 600) {
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(folder, paste0(name, "_SEA.", ext)), plot,
                    width = width, height = height, dpi = dpi, bg = "white")
  }
}

save_session_SEA <- function(folder) {
  writeLines(capture.output(sessionInfo()), file.path(folder, "sessionInfo_SEA.txt"))
}

condition_colors_SEA <- c(R = "#D62728", U = "#2CA02C", RU = "#FF7F0E",
                       Ctrl = "#1F77B4", Carrier = "#9467BD")

# Keep saved list fields readable; tag the objects loaded into an R session.
read_data_SEA <- function(file) {
  data <- readRDS(file)
  setNames(data, paste0(names(data), "_SEA"))
}

# Display the SEA dataset as PTM in figure titles, axes and legends.
# File paths, saved data and internal dataset keys continue to use SEA.
display_label_SEA <- function(x) {
  gsub("\\bSEA\\b", "PTM", as.character(x), perl = TRUE)
}
