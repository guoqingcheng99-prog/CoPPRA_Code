# Run core stages 1.1-4.2 with: Rscript /path/to/Base/Code/run_Base.R
args_BASE <- commandArgs(trailingOnly = FALSE)
script_BASE <- sub('^--file=', '', args_BASE[grepl('^--file=', args_BASE)])
code_dir_BASE <- if (length(script_BASE)) dirname(normalizePath(script_BASE)) else getwd()
Sys.setenv(COPPRA_PROJECT_DIR = dirname(code_dir_BASE))
source(file.path(code_dir_BASE, "CoPPRA_config.R"), local = TRUE)
if (!rmarkdown::pandoc_available()) {
  pandoc_BASE <- Sys.glob('/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/*/pandoc')
  arch_BASE <- if (grepl('aarch64|arm64', R.version$arch)) 'aarch64' else 'x86_64'
  pandoc_BASE <- pandoc_BASE[grepl(paste0('/', arch_BASE, '/'), pandoc_BASE)]
  if (length(pandoc_BASE)) Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_BASE[1]))
}
stopifnot(rmarkdown::pandoc_available())
reports_BASE <- file.path(output_root, 'reports')
logs_BASE <- file.path(output_root, 'logs')
dir.create(reports_BASE, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_BASE, recursive = TRUE, showWarnings = FALSE)
files_BASE <- sort(list.files(code_dir_BASE, pattern = '^[1-4]\\.[12].*\\.rmd$', full.names = TRUE))
# An optional argument selects a starting stage when resuming a stopped run.
start_BASE <- commandArgs(trailingOnly = TRUE)
if (length(start_BASE)) files_BASE <- files_BASE[substr(basename(files_BASE), 1, 3) >= start_BASE[1]]
status_file_BASE <- file.path(logs_BASE, 'completed_stages_BASE.csv')
status_BASE <- if (file.exists(status_file_BASE)) read.csv(status_file_BASE) else
  data.frame(File = character(), Started = character(), Finished = character())
render_dir_BASE <- file.path(logs_BASE, 'render')
dir.create(render_dir_BASE, showWarnings = FALSE)
setwd(code_dir_BASE)
for (file_BASE in files_BASE) {
  started_BASE <- as.character(Sys.time())
  message('\nSTART ', basename(file_BASE), ' at ', started_BASE)
  rmarkdown::render(file_BASE, output_dir = reports_BASE,
                    knit_root_dir = code_dir_BASE, intermediates_dir = render_dir_BASE, envir = new.env(parent = globalenv()),
                    clean = TRUE, quiet = FALSE)
  status_BASE <- status_BASE[status_BASE$File != basename(file_BASE), , drop = FALSE]
  status_BASE <- rbind(status_BASE, data.frame(File = basename(file_BASE),
    Started = started_BASE, Finished = as.character(Sys.time())))
  write.csv(status_BASE, status_file_BASE, row.names = FALSE)
  gc()
  message('DONE ', basename(file_BASE))
}
