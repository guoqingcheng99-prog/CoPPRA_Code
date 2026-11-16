# Reproducible heatmaps from the existing CoPPRA QC handoff.
# This script never writes to QC, PCA, or Limma outputs.
run_coppra_heatmaps <- function(q_file, output_dir, analysis_label) {
  for (pkg in c("data.table", "matrixStats", "ComplexHeatmap", "circlize", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  dirs <- list(plots = file.path(plot_root,"Heatmap"), tables=file.path(output_dir,"tables"), data=file.path(output_dir,"data"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  q <- readRDS(q_file)
  source_hash <- unname(tools::md5sum(q_file))
  sample_map <- as.data.frame(q$sample_map)
  sample_map$Condition <- as.character(sample_map$Condition)
  sample_map <- sample_map[sample_map$Condition %in% c("Ctrl", "R", "U", "RU"), ]
  sample_map <- sample_map[order(match(sample_map$Condition, c("Ctrl", "R", "U", "RU")),
                                 sample_map$Bio_rep), ]
  samples <- as.character(sample_map$Experiment)
  stopifnot(length(samples) == 12L, !anyDuplicated(samples),
            all(table(sample_map$Condition) == 3L),
            all(samples %in% colnames(q$biological_log2)))
  # The source is already log2 transformed and raw-file median normalized.
  source <- q$biological_log2[, samples, drop = FALSE]
  stopifnot(!anyDuplicated(rownames(source)))
  observed_n <- rowSums(is.finite(source))
  keep <- rownames(source) %in% q$primary_sequences
  before <- source[keep, , drop = FALSE]
  missing <- !is.finite(before)
  medians <- matrixStats::rowMedians(before, na.rm = TRUE)
  x <- before
  idx <- which(missing, arr.ind = TRUE)
  x[idx] <- medians[idx[, 1L]]
  stopifnot(all(is.finite(x)), identical(x[!missing], before[!missing]))
  variances <- matrixStats::rowVars(x)
  correlation <- cor(x, method = "pearson")
  distance <- pmax(0, 1 - correlation)
  dim(distance) <- dim(correlation)
  dimnames(distance) <- dimnames(correlation)
  diag(distance) <- 0
  sample_tree <- hclust(as.dist(distance), method = "average")
  # A secondary view only: select by unsupervised variance, not significance.
  variable_rows <- which(is.finite(variances) & variances > 0)
  selected <- head(variable_rows[order(-variances[variable_rows],
                                       rownames(x)[variable_rows])], 1000L)
  z <- sweep(x[selected, , drop = FALSE], 1L, rowMeans(x[selected, , drop = FALSE]), "-")
  z <- sweep(z, 1L, sqrt(variances[selected]), "/")
  row_distance <- 1 - cor(t(z))
  row_distance[row_distance < 0] <- 0
  row_dist <- as.dist(row_distance)
  peptide_tree <- hclust(row_dist, method = "average")
  matrix_table <- function(m, key = "Sequence") {
    ans <- data.table::as.data.table(m)
    ans <- cbind(data.table::data.table(ID = rownames(m)), ans)
    data.table::setnames(ans, "ID", key)
    ans
  }
  write_table <- function(x, name) data.table::fwrite(x, file.path(dirs$tables, paste0(name, ".csv")))
  write_table(sample_map, "sample_metadata")
  write_table(matrix_table(x), "peptide_normalized_log2_imputed_matrix")
  write_table(matrix_table(correlation, "Sample"), "sample_Pearson_correlation")
  write_table(matrix_table(z), "top1000_variable_peptides_row_zscore")
  write_table(data.table::data.table(
    Sequence = rownames(source), Observed_samples = observed_n, Included = keep,
    Exclusion = ifelse(keep, "", "Fails at least 2/3 observed biological replicates in every condition")),
    "peptide_inclusion")
  write_table(data.table::data.table(
    Sample = samples, Condition = sample_map$Condition,
    Observed_values = colSums(!missing), Imputed_values = colSums(missing),
    Imputed_percent = 100 * colMeans(missing)), "imputation_by_sample")
  write_table(data.table::data.table(
    Order = seq_along(sample_tree$order), Sample = samples[sample_tree$order]),
    "sample_cluster_order")
  write_table(data.table::data.table(
    Sequence = rownames(z), Variance_before_row_scaling = variances[selected]),
    "top1000_peptide_selection")
  all_row_order <- order(-rowMeans(x), rownames(x))
  write_table(data.table::data.table(Order = seq_len(nrow(x)), Sequence = rownames(x)[all_row_order]),
    "all_peptide_display_order")
  metrics <- data.table::data.table(
    Metric = c("Biological samples", "QC peptide rows", "Peptide rows retained",
               "Excluded rows: fewer than 2 observations", "Values imputed", "Values imputed (%)",
               "Minimum off-diagonal Pearson r", "Maximum off-diagonal Pearson r"),
    Value = c(ncol(x), nrow(source), nrow(x), sum(!keep), sum(missing), 100 * mean(missing),
              range(correlation[lower.tri(correlation)])))
  write_table(metrics, "analysis_summary")
  colors <- c(Ctrl = "#1F77B4", R = "#D62728", U = "#2CA02C", RU = "#FF7F0E")
  annotation <- function() ComplexHeatmap::HeatmapAnnotation(
    Condition = factor(sample_map$Condition, levels = names(colors)),
    col = list(Condition = colors), show_legend = FALSE,
    annotation_name_gp = grid::gpar(fontsize = 9), simple_anno_size = grid::unit(3, "mm"))
  n_label <- format(nrow(x), big.mark = ",")
  footer <- function(note) grid::grid.text(note, x = grid::unit(0.02, "npc"),
    y = grid::unit(0.025, "npc"), just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 8.5, col = "#444444"))
  draw_heatmap <- function(ht, note) {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
      padding = grid::unit(c(18, 6, 6, 6), "mm"))
    footer(note)
  }
  render <- function(stem, make_plot, width = 9, height = 8.5) {
    grDevices::pdf(file.path(dirs$plots, paste0(stem, ".pdf")), width = width,
                   height = height, useDingbats = FALSE)
    tryCatch(make_plot(), finally = grDevices::dev.off())
    for (dpi in c(600L, 160L)) {
      suffix <- if (dpi == 600L) ".png" else "_preview.png"
      ragg::agg_png(file.path(dirs$plots, paste0(stem, suffix)), width = width,
                    height = height, units = "in", res = dpi, background = "white")
      tryCatch(make_plot(), finally = grDevices::dev.off())
    }
  }
  corr_low <- floor(min(correlation) * 10) / 10
  render("sample_correlation_12", function() {
    ht <- ComplexHeatmap::Heatmap(correlation, name = "Pearson r",
      col = circlize::colorRamp2(c(corr_low, (corr_low + 1) / 2, 1),
                               c("#F3F6FA", "#88B5D6", "#164B7B")),
      cluster_rows = sample_tree, cluster_columns = sample_tree,
      row_dend_reorder = FALSE, column_dend_reorder = FALSE,
      top_annotation = annotation(), column_names_rot = 45,
      row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
      column_title = paste0(analysis_label, " | 12-sample correlation\n", n_label,
                            " normalized, log2-transformed, imputed peptides"),
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      rect_gp = grid::gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(at = seq(corr_low, 1, by = 0.05)),
      cell_fun = function(j, i, xx, yy, w, h, fill) {
        grid::grid.text(sprintf("%.3f", correlation[i, j]), xx, yy,
          gp = grid::gpar(fontsize = 8.5, col = if (correlation[i, j] > 0.93) "white" else "#182C40"))
      })
    draw_heatmap(ht, paste0("Average-linkage clustering; distance = 1 - Pearson r. References excluded.\n",
      "Missing values: peptide median across the 12 samples (", sprintf("%.2f", 100 * mean(missing)), "% imputed)."))
  })
  log2_limit <- ceiling(max(abs(range(x))))
  render("all_peptide_abundance_12", function() {
    ht <- ComplexHeatmap::Heatmap(x, name = "Normalized\nlog2 intensity",
      col = circlize::colorRamp2(c(-log2_limit, 0, log2_limit), c("#2166AC", "#F7F7F7", "#B2182B")),
      cluster_rows = FALSE, row_order = all_row_order,
      cluster_columns = sample_tree, column_dend_reorder = FALSE,
      show_row_names = FALSE, top_annotation = annotation(), column_names_rot = 45,
      column_names_gp = grid::gpar(fontsize = 11),
      column_title = paste0(analysis_label, " | All retained peptide abundances\n", n_label, " peptides x 12 biological samples"),
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      row_title = paste0(n_label, " peptides (ordered by mean abundance)"),
      use_raster = TRUE, raster_device = "agg_png", raster_quality = 3,
      raster_by_magick = FALSE, raster_resize_mat = FALSE)
    draw_heatmap(ht, "Samples: average-linkage clustering, distance = 1 - Pearson r, using all retained peptides.\nRows ordered by decreasing mean abundance; no row clustering or row scaling. References excluded.")
  }, width = 9, height = 9)
  render("top1000_peptide_abundance_12", function() {
    ht <- ComplexHeatmap::Heatmap(z, name = "Peptide\nrow z-score",
      col = circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B")),
      cluster_rows = peptide_tree, cluster_columns = sample_tree,
      row_dend_reorder = FALSE, column_dend_reorder = FALSE,
      show_row_names = FALSE, top_annotation = annotation(), column_names_rot = 45,
      column_names_gp = grid::gpar(fontsize = 11),
      column_title = paste0(analysis_label, " | Most variable peptide abundances\n", nrow(z),
                            " peptides x 12 biological samples"),
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      row_title = paste0(nrow(z), " most variable peptides"),
      heatmap_legend_param = list(at = c(-3, 0, 3), labels = c("<= -3", "0", ">= 3")),
      use_raster = TRUE, raster_device = "agg_png", raster_quality = 3,
      raster_by_magick = FALSE, raster_resize_mat = FALSE)
    draw_heatmap(ht, "Rows: selected by variance after imputation, row z-scored, then clustered (1 - Pearson r; average linkage).\nSample tree uses all retained peptides, as in the correlation heatmap. Color saturates at +/-3 SD.")
  }, width = 9, height = 9)
  saveRDS(list(matrix_before_imputation = before, matrix_imputed = x, missing_mask = missing,
               sample_metadata = sample_map, correlation = correlation, sample_tree = sample_tree,
               top_peptide_zscores = z, top_peptide_tree = peptide_tree,
               source_path = normalizePath(q_file), source_md5 = source_hash,
               source_method = q$method), file.path(dirs$data, "heatmap_data.rds"))
  writeLines(c(
    paste0("# ", analysis_label, ": heatmaps of all 12 biological samples"), "",
    paste0("Input: `", normalizePath(q_file), "`"),
    paste0("Input MD5: `", source_hash, "`"), "",
    "Samples: Ctrl1-3, R1-3, U1-3, RU1-3. Three carrier/reference samples are excluded.",
    "The saved QC input was already log2 transformed and median centered per raw file,",
    "collapsed by unmodified peptide Sequence, then median aggregated across technical runs.",
    "No second log transformation or normalization was applied.", "",
    paste0("Retained ", n_label, " peptides observed in at least 2/3 biological replicates in each of the four conditions; excluded ", sum(!keep), " rows."),
    paste0("Filled ", sum(missing), " missing values (", sprintf("%.4f", 100 * mean(missing)),
           "%) with each peptide's observed median across those 12 samples, matching the existing PCA approach."),
    "Observed values are unchanged. No significance filtering was used.", "",
    "- sample_correlation_12: Pearson correlation of all retained imputed peptides. Both axes use the same average-linkage tree on 1 - r.",
    "- all_peptide_abundance_12: all retained normalized log2 peptide values. Samples use that tree; rows are ordered by decreasing mean abundance, not clustered.",
    "- top1000_peptide_abundance_12: 1,000 highest-variance peptides after imputation, row z-scored for display. Rows use average linkage on 1 - r; sample tree still uses all retained peptides. Colors saturate at +/-3 SD.", "",
    "PNG exports are 600 dpi; PDF exports and smaller previews are included.",
    "Full numerical matrices, inclusion decisions, imputation counts, display order, and trees are saved in tables/ and data/.",
    "Peptide rows in overview figures are rasterized and compressed at page scale; use the matrix exports for individual values.",
    "Clustering is descriptive. Median imputation can affect correlations and clustering.",
    "Existing QC, PCA, and differential-analysis files were not modified."),
    file.path(output_dir, "README.md"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  stopifnot(identical(unname(tools::md5sum(q_file)), source_hash))
  message(analysis_label, ": completed ", nrow(x), " peptides x ", ncol(x), " samples.")
  invisible(list(matrix = x, sample_tree = sample_tree, summary = metrics))
}
