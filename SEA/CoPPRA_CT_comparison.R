# Shared CT tables and schematic Venn style for Base, SEA and PSEA.
ct_venn_plot <- function(sets) {
  stopifnot(length(sets) %in% c(2L, 3L), !is.null(names(sets)))
  labels <- names(sets)
  display_labels <- ifelse(labels == "SEA", "PTM", labels)
  palette <- c(Base = "#2878B5", SEA = "#E88931", PSEA = "#43A047")
  stopifnot(all(labels %in% names(palette)))
  angle <- seq(0, 2 * pi, length.out = 361)
  if (length(sets) == 2L) {
    centers <- data.table::data.table(x = c(-0.65, 0.65), y = c(0, 0), Dataset = labels)
    positions <- data.table::data.table(x = c(-1.1, 0, 1.1), y = 0,
      Count = c(length(setdiff(sets[[1]], sets[[2]])), length(intersect(sets[[1]], sets[[2]])),
                length(setdiff(sets[[2]], sets[[1]]))))
    titles <- data.table::data.table(x = c(-1, 1), y = 1.5,
      Label = paste0(display_labels, "\nCT = ", lengths(sets)))
    ylim <- c(-1.45, 1.85)
  } else {
    centers <- data.table::data.table(x = c(-0.65, 0.65, 0), y = c(0.4, 0.4, -0.7), Dataset = labels)
    a <- sets[[1]]; b <- sets[[2]]; c <- sets[[3]]
    regions <- list(setdiff(a, union(b, c)), setdiff(b, union(a, c)), setdiff(c, union(a, b)),
      setdiff(intersect(a, b), c), setdiff(intersect(a, c), b), setdiff(intersect(b, c), a),
      Reduce(intersect, sets))
    positions <- data.table::data.table(x = c(-1.1, 1.1, 0, 0, -0.65, 0.65, 0),
      y = c(0.7, 0.7, -1.45, 1, -0.65, -0.65, 0), Count = lengths(regions))
    titles <- data.table::data.table(x = c(-1.1, 1.1, 0), y = c(1.95, 1.95, -2.15),
      Label = paste0(display_labels, "\nCT = ", lengths(sets)))
    ylim <- c(-2.5, 2.35)
  }
  circles <- data.table::rbindlist(lapply(seq_len(nrow(centers)), function(i)
    data.table::data.table(x = centers$x[i] + 1.2 * cos(angle),
      y = centers$y[i] + 1.2 * sin(angle), Dataset = centers$Dataset[i])))
  ggplot2::ggplot(circles, ggplot2::aes(x, y, group = Dataset, fill = Dataset, colour = Dataset)) +
    ggplot2::geom_polygon(alpha = 0.25, linewidth = 0.9) +
    ggplot2::geom_text(data = positions, ggplot2::aes(x, y, label = Count), inherit.aes = FALSE, size = 6) +
    ggplot2::geom_text(data = titles, ggplot2::aes(x, y, label = Label), inherit.aes = FALSE, size = 5) +
    ggplot2::scale_fill_manual(values = palette) + ggplot2::scale_colour_manual(values = palette) +
    ggplot2::coord_fixed(xlim = c(-2.1, 2.1), ylim = ylim, clip = "off") +
    ggplot2::labs(title = "CT protein overlap") + ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin = ggplot2::margin(15, 15, 15, 15))
}

compare_ct_results <- function(files, paths, prefix) {
  if (!all(file.exists(files))) stop("Run the corresponding 3.1 limma analyses first. Missing: ",
    paste(files[!file.exists(files)], collapse = ", "))
  results <- lapply(files, readRDS)
  stopifnot(all(vapply(results, function(x) identical(x$pval_threshold, results[[1]]$pval_threshold) &&
    identical(x$log2fc_threshold, results[[1]]$log2fc_threshold), logical(1))))
  sets <- lapply(results, function(x) {
    p <- data.table::as.data.table(x$protein_significance)
    stopifnot(!anyDuplicated(p$ProteinID), !anyNA(p[, .(ProteinID, Ru_sig, Ul_sig, RU_sig)]))
    sort(p[RU_sig & !Ru_sig & !Ul_sig, ProteinID])
  })
  labels <- names(sets)
  all_ids <- sort(unique(as.character(unlist(sets, use.names = FALSE))))
  membership <- data.table::data.table(ProteinID = all_ids)
  for (label in labels) membership[, (label) := ProteinID %in% sets[[label]]]
  membership[, Intersection := vapply(seq_len(.N), function(i)
    paste(labels[vapply(sets, function(ids) ProteinID[i] %in% ids, logical(1))], collapse = " & "), character(1))]
  common <- membership[ProteinID %in% Reduce(intersect, sets)]
  counts <- data.table::data.table(Set = c(labels, "Common to all", "Union"),
    CT_proteins = c(lengths(sets), nrow(common), nrow(membership)))
  regions <- lapply(seq_len(2^length(sets) - 1L), function(mask) {
    included <- as.logical(intToBits(mask)[seq_along(sets)])
    selected <- Reduce(intersect, sets[included])
    exclusive <- setdiff(selected, unlist(sets[!included], use.names = FALSE))
    data.table::data.table(Intersection = paste(labels[included], collapse = " & "), CT_proteins = length(exclusive))
  })
  intersections <- data.table::rbindlist(regions)
  stopifnot(sum(intersections$CT_proteins) == nrow(membership))
  for (folder in paths) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  tables <- list(counts = counts, intersections = intersections, membership = membership, common_CT = common)
  for (name in names(tables)) data.table::fwrite(tables[[name]],
    file.path(paths$tables, paste0(prefix, "_", name, ".csv")))
  for (label in labels) data.table::fwrite(membership[ProteinID %in% sets[[label]]],
    file.path(paths$tables, paste0(prefix, "_", label, "_CT.csv")))
  data.table::fwrite(data.table::data.table(Dataset = names(files), Result_file = unname(files),
    MD5 = unname(tools::md5sum(files))), file.path(paths$tables, paste0(prefix, "_inputs.csv")))
  plot <- ct_venn_plot(sets)
  for (ext in c("png", "pdf")) ggplot2::ggsave(file.path(paths$plots, paste0(prefix, "_Venn.", ext)),
    plot, width = 8, height = if (length(sets) == 2L) 6 else 8, units = "in", dpi = 600, bg = "white")
  saveRDS(list(sets = sets, counts = counts, intersections = intersections, membership = membership,
    common = common, files = files), file.path(paths$data, paste0(prefix, ".rds")))
  list(plot = plot, counts = counts, intersections = intersections, common = common, membership = membership)
}
