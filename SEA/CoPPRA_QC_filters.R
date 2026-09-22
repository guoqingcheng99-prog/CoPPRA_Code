# Shared CoPPRA QC filters. Requires data.table; does not change normalization.
# A missing proteinGroups table produces an explicitly incomplete coverage audit.
# The calling QC Rmd must not save a downstream handoff from an incomplete audit.
coppra_qc <- function(evidence, sample_map, protein_groups_file, audit_dir,
                      suffix = "", evidence_q_max = 0.01, protein_q_max = 0.01,
                      min_bio_replicates = 2L, min_protein_peptides = 2L,
                      filter_peptide_recurrence = FALSE,
                      filter_protein_peptide_count = TRUE,
                      min_raw_file_fraction = 0.20) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  ev <- data.table::copy(data.table::as.data.table(evidence))
  data.table::setnames(ev, make.names(names(ev)))
  required <- c("Sequence", "Raw.file", "Experiment", "Intensity", "Q.value",
                "Proteins", "Leading.proteins", "Leading.razor.proteins",
                "Protein.group.IDs", "Peptide.ID", "Decoy", "Potential.contaminant")
  missing <- setdiff(required, names(ev))
  if (length(missing)) stop("Missing evidence columns: ", paste(missing, collapse = ", "))
  if (any(!ev$Experiment %in% sample_map$Experiment)) stop("Unknown evidence Experiment label.")
  stopifnot(min_bio_replicates >= 1L, min_protein_peptides >= 1L,
            evidence_q_max >= 0, evidence_q_max <= 1,
            protein_q_max >= 0, protein_q_max <= 1,
            min_raw_file_fraction >= 0, min_raw_file_fraction <= 1)
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  emit <- function(x, name) data.table::fwrite(x, file.path(audit_dir, paste0(name, suffix, ".csv")))
  split_ids <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    unique(trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)))
  }
  nonempty_count <- function(x) length(unique(x[!is.na(x) & nzchar(as.character(x))]))
  flagged <- function(x) !is.na(x) & trimws(as.character(x)) == "+"
  ev[, Protein.group.IDs := as.character(Protein.group.IDs)]
  ev[, QC_row := .I]
  ev[, QC_condition := as.character(sample_map$Condition[match(Experiment, sample_map$Experiment)])]
  complete <- length(protein_groups_file) == 1L && !is.na(protein_groups_file) && file.exists(protein_groups_file)
  trace <- list()
  record <- function(d, stage, status = "applied") {
    trace[[length(trace) + 1L]] <<- data.table::data.table(
      Stage = stage, Status = status, Evidence_rows = nrow(d),
      Peptide_sequences = nonempty_count(d$Sequence),
      Protein_accessions = length(split_ids(unique(d$Proteins))),
      Candidate_protein_group_IDs = length(split_ids(d$Protein.group.IDs)),
      Leading_protein_combinations = nonempty_count(d$Leading.proteins),
      Assigned_groups = if ("QC_group" %in% names(d) && complete) nonempty_count(d$QC_group) else NA_integer_,
      Razor_representatives = nonempty_count(d$Leading.razor.proteins))
  }
  before <- data.table::copy(ev)
  record(ev, "01 Imported evidence", "baseline")
  ev <- ev[!flagged(Decoy)]
  record(ev, "02 Remove evidence decoys")
  ev <- ev[!flagged(Potential.contaminant)]
  record(ev, "03 Remove evidence contaminants")
  ev <- ev[!is.na(Sequence) & Sequence != "" & is.finite(Intensity) & Intensity > 0]
  legacy <- data.table::copy(ev)
  record(ev, "04 Valid sequence and positive finite intensity (previous QC)")
  # This is an explicit evidence/sample q-value cutoff, separate from global
  # peptide/protein identification FDR and separate from an individual PEP.
  ev <- ev[is.finite(Q.value) & Q.value >= 0 & Q.value <= evidence_q_max]
  record(ev, "05 Evidence Q-value <= threshold")

  if (complete) {
    pg <- data.table::fread(protein_groups_file, na.strings = c("", "NA"))
    data.table::setnames(pg, make.names(names(pg)))
    needed <- c("id", "Protein.IDs", "Only.identified.by.site", "Q.value", "Potential.contaminant",
                "Peptide.IDs", "Peptide.is.razor", "Peptide.sequences")
    missing <- setdiff(needed, names(pg))
    if (length(missing)) stop("Missing proteinGroups columns: ", paste(missing, collapse = ", "))
    decoy_cols <- intersect(c("Reverse", "Decoy"), names(pg))
    if (!length(decoy_cols)) stop("proteinGroups needs a Reverse or Decoy flag.")
    pg[, id := as.character(id)]
    if (anyNA(pg$id) || anyDuplicated(pg$id)) stop("Invalid/duplicate protein-group IDs.")
    if (length(setdiff(split_ids(before$Protein.group.IDs), pg$id)))
      stop("Some evidence protein-group IDs are missing from proteinGroups; use the same run.")
    pg[, QC_decoy := Reduce(`|`, lapply(.SD, flagged)), .SDcols = decoy_cols]
    pg[, QC_contaminant := flagged(Potential.contaminant)]
    pg[, QC_site_only := flagged(Only.identified.by.site)]
    pg[, QC_bad_q := !is.finite(Q.value) | Q.value < 0 | Q.value > protein_q_max]
    peptide_ids <- strsplit(as.character(pg$Peptide.IDs), ";", fixed = TRUE)
    peptide_razor <- strsplit(as.character(pg$Peptide.is.razor), ";", fixed = TRUE)
    peptide_seqs <- strsplit(as.character(pg$Peptide.sequences), ";", fixed = TRUE)
    if (any(lengths(peptide_ids) != lengths(peptide_razor)) ||
        any(lengths(peptide_ids) != lengths(peptide_seqs)))
      stop("Unpaired peptide IDs, razor flags or peptide sequences in proteinGroups.")
    pg_peptides <- data.table::data.table(
      QC_group = rep(pg$id, lengths(peptide_ids)),
      Peptide.ID = unlist(peptide_ids, use.names = FALSE),
      Sequence = unlist(peptide_seqs, use.names = FALSE),
      QC_razor = tolower(trimws(unlist(peptide_razor, use.names = FALSE))))
    if (any(!pg_peptides$QC_razor %in% c("true", "false"))) stop("Unknown proteinGroups razor flag.")
    ev[, Peptide.ID := as.character(Peptide.ID)]
    # Numeric IDs are run-local: verify sequence identity, not just ID overlap.
    sequence_check <- unique(before[, .(Peptide.ID = as.character(Peptide.ID), Sequence)])
    known_sequences <- unique(pg_peptides[, .(Peptide.ID, Sequence)])
    if (nrow(sequence_check[!known_sequences, on = .(Peptide.ID, Sequence)]))
      stop("Evidence peptide IDs/sequences disagree with proteinGroups; use the same run.")

    # Resolve the existing razor assignment against its candidate group IDs.
    # Never reassign a peptide to another group merely because its chosen group fails QC.
    members <- lapply(pg$Protein.IDs, split_ids)
    names(members) <- pg$id
    tuples <- unique(before[, .(Leading.razor.proteins, Protein.group.IDs)])
    tuples[, QC_group := vapply(seq_len(.N), function(i) {
      razor <- split_ids(Leading.razor.proteins[i])
      candidate <- split_ids(Protein.group.IDs[i])
      if (!length(razor) || !length(candidate)) return(NA_character_)
      matched <- candidate[vapply(candidate, function(g) any(razor %in% members[[g]]), logical(1))]
      if (length(matched) == 1L) matched else NA_character_
    }, character(1))]
    # Detect a mismatched table instead of silently discarding all evidence.
    one <- !grepl(";", tuples$Protein.group.IDs) & !is.na(tuples$Leading.razor.proteins)
    if (any(one & is.na(tuples$QC_group)))
      stop("Razor accessions disagree with proteinGroups membership; verify the paired input files.")
    before <- merge(before, tuples, by = c("Leading.razor.proteins", "Protein.group.IDs"), all.x = TRUE, sort = FALSE)
    legacy <- merge(legacy, tuples, by = c("Leading.razor.proteins", "Protein.group.IDs"), all.x = TRUE, sort = FALSE)
    ev <- merge(ev, tuples, by = c("Leading.razor.proteins", "Protein.group.IDs"), all.x = TRUE, sort = FALSE)
    ev <- ev[!is.na(QC_group)]
    record(ev, "06 Resolve unique-plus-razor group assignment")
    eligible <- pg_peptides[QC_razor == "true", .(QC_group, Peptide.ID, Sequence)]
    ev <- ev[eligible, on = .(QC_group, Peptide.ID, Sequence), nomatch = 0]
    record(ev, "06b Require proteinGroups peptide-is-razor flag")
    pg_index <- match(ev$QC_group, pg$id)
    ev <- ev[!pg$QC_decoy[pg_index]]
    record(ev, "07 Remove reverse protein groups")
    ev <- ev[!pg$QC_contaminant[match(QC_group, pg$id)]]
    record(ev, "08 Remove contaminant protein groups")
    ev <- ev[!pg$QC_site_only[match(QC_group, pg$id)]]
    record(ev, "09 Remove only-identified-by-site groups")
    ev <- ev[!pg$QC_bad_q[match(QC_group, pg$id)]]
    record(ev, "10 Protein-group Q-value <= threshold")
    emit(pg[, .(id, Protein.IDs, QC_decoy, QC_contaminant, QC_site_only, Q.value, QC_bad_q)], "QC_protein_group_flags")
  } else {
    # Provisional razor representatives are not asserted to be verified protein groups.
    ev[, QC_group := as.character(Leading.razor.proteins)]
    ev <- ev[!is.na(QC_group) & nzchar(QC_group) & !grepl(";", QC_group)]
    record(ev, "06 Provisional single razor representative", "provisional; proteinGroups missing")
    for (stage in c("07 Remove reverse protein groups", "08 Remove contaminant protein groups",
                    "09 Remove only-identified-by-site groups", "10 Protein-group Q-value <= threshold"))
      record(ev, stage, "NOT APPLIED: proteinGroups missing")
  }
  # Sequence pooling cannot safely retain inconsistent group assignments.
  assignment <- ev[, .(N_assignments = data.table::uniqueN(QC_group)), by = Sequence]
  ev <- ev[Sequence %in% assignment[N_assignments == 1L, Sequence]]
  record(ev, "11 Require consistent assignment per pooled Sequence")
  recurrence <- unique(ev[QC_condition != "Carrier", .(Sequence, Experiment, QC_condition)])[
    , .(Biological_replicates = .N), by = .(Sequence, QC_condition)]
  pass_sequences <- unique(recurrence[Biological_replicates >= min_bio_replicates, Sequence])
  recurrence[, Pass_condition := Biological_replicates >= min_bio_replicates]
  emit(recurrence, "QC_peptide_recurrence")
  if (filter_peptide_recurrence) ev <- ev[Sequence %in% pass_sequences]
  record(ev, sprintf("12 >=%d biological replicates in >=1 non-carrier condition", min_bio_replicates),
         if (filter_peptide_recurrence) "applied" else "DISABLED: no biological-replicate requirement")

  # Preserve SEA's verified unique-plus-razor mapping before counting support.
  if (complete) {
    ev[, Original.Proteins := Proteins]
    ev[, Proteins := pg$Protein.IDs[match(QC_group, pg$id)]]
    ev[, Leading.proteins := Proteins]
    ev[, Protein.group.IDs := QC_group]
    if ("Gene.names" %in% names(pg)) ev[, Gene.Names := pg$Gene.names[match(QC_group, pg$id)]]
    if ("Protein.names" %in% names(pg)) ev[, Protein.Names := pg$Protein.names[match(QC_group, pg$id)]]
  }
  overview <- function(d, stage) data.table::rbindlist(lapply(c("All_samples", "Non_carrier"), function(scope) {
    z <- if (scope == "Non_carrier") d[QC_condition != "Carrier"] else d
    input <- if (scope == "Non_carrier") before[QC_condition != "Carrier"] else before
    n_raw <- data.table::uniqueN(input$Raw.file)
    n_peptides <- nonempty_count(z$Sequence)
    n_detected <- nrow(unique(z[, .(Sequence, Raw.file)]))
    data.table::data.table(Stage = stage, Scope = scope, Evidence_rows = nrow(z),
      Peptide_sequences = n_peptides,
      Protein_accessions = length(split_ids(unique(z$Proteins))),
      Candidate_protein_group_IDs = length(split_ids(z$Protein.group.IDs)),
      Razor_representatives = nonempty_count(z$Leading.razor.proteins),
      Verified_assigned_groups = if (complete && "QC_group" %in% names(z)) nonempty_count(z$QC_group) else NA_integer_,
      Total_raw_files = n_raw, Raw_files_with_evidence = data.table::uniqueN(z$Raw.file),
      Peptide_raw_file_detections = n_detected,
      Peptide_completeness = if (n_peptides > 0L) n_detected / (n_peptides * n_raw) else NA_real_,
      Peptide_missingness = if (n_peptides > 0L) 1 - n_detected / (n_peptides * n_raw) else NA_real_)
  }))
  coverage_before_added <- overview(ev, "Before_protein_and_detection_filters")

  # Count each accession's distinct pooled sequences before the detection filter.
  protein_map <- unique(ev[, .(Sequence, Proteins)])
  protein_map <- unique(protein_map[, .(ProteinID = split_ids(Proteins)), by = Sequence])
  support <- protein_map[, .(Distinct_peptide_sequences = data.table::uniqueN(Sequence)), by = ProteinID]
  support[, Pass_peptide_count := Distinct_peptide_sequences >= min_protein_peptides]
  emit(support, "QC_protein_peptide_support")
  retained_proteins <- support[Pass_peptide_count == TRUE, ProteinID]
  ev[, Original.QC.group := QC_group]
  if (filter_protein_peptide_count) {
    annotations <- unique(ev[, .(Proteins, Leading.proteins)])
    retain_ids <- function(x) vapply(x, function(value) {
      ids <- split_ids(value)
      paste(ids[ids %in% retained_proteins], collapse = ";")
    }, character(1))
    annotations[, `:=`(Filtered_proteins = retain_ids(Proteins),
      Filtered_leading_proteins = retain_ids(Leading.proteins),
      Had_protein_annotation = !is.na(Proteins) & nzchar(trimws(Proteins)))]
    ev[annotations, on = .(Proteins, Leading.proteins), `:=`(
      Proteins = i.Filtered_proteins, Leading.proteins = i.Filtered_leading_proteins,
      Had_protein_annotation = i.Had_protein_annotation)]
    ev <- ev[!Had_protein_annotation | nzchar(Proteins)]
    ev[, Had_protein_annotation := NULL]
  }
  ev[, Protein_group_retained := TRUE]
  record(ev, sprintf("13 Retain protein accessions with >=%d distinct peptides", min_protein_peptides),
    if (filter_protein_peptide_count) "applied; remove rows mapped only to excluded accessions" else "DISABLED")
  coverage_after_protein <- overview(ev, "After_protein_filter")

  # The denominator is all input raw files, including carriers and files that
  # have no surviving observations after quality filtering.
  total_raw_files <- data.table::uniqueN(before$Raw.file)
  minimum_raw_files <- as.integer(ceiling(min_raw_file_fraction * total_raw_files))
  raw_detection <- ev[, .(Raw_files_detected = data.table::uniqueN(Raw.file)), by = Sequence]
  raw_detection[, `:=`(Total_raw_files = total_raw_files, Minimum_raw_files = minimum_raw_files,
    Detection_fraction = Raw_files_detected / total_raw_files,
    Pass_raw_file_fraction = Raw_files_detected >= minimum_raw_files)]
  emit(raw_detection, "QC_peptide_raw_file_detection")
  ev <- ev[Sequence %in% raw_detection[Pass_raw_file_fraction == TRUE, Sequence]]
  record(ev, sprintf("14 Peptide observed in >=%d/%d raw files (%.0f%%)",
    minimum_raw_files, total_raw_files, 100 * min_raw_file_fraction))
  data.table::setorder(ev, QC_row)

  stages <- data.table::rbindlist(trace)
  stages[, Evidence_rows_removed := c(0L, head(Evidence_rows, -1L) - tail(Evidence_rows, -1L))]
  stages[, Peptide_sequences_removed := c(0L, head(Peptide_sequences, -1L) - tail(Peptide_sequences, -1L))]
  stages[, Candidate_protein_group_IDs_removed := c(0L, head(Candidate_protein_group_IDs, -1L) - tail(Candidate_protein_group_IDs, -1L))]
  emit(stages, "QC_filter_stages")
  coverage <- data.table::rbindlist(list(overview(before, "Before_QC"), overview(legacy, "Previous_QC"),
    coverage_before_added, coverage_after_protein,
    overview(ev, if (complete) "After_complete_QC" else "After_available_QC_INCOMPLETE")))
  emit(coverage, "QC_coverage_before_after")
  scoped <- function(d, stage, key, universe) {
    levels <- unique(before[[key]])
    counts <- d[, .(Evidence_rows = .N, Peptide_sequences = nonempty_count(Sequence),
      Razor_representatives = nonempty_count(Leading.razor.proteins),
      Verified_assigned_groups = if (complete && "QC_group" %in% names(d)) nonempty_count(QC_group) else NA_integer_), by = key]
    out <- merge(data.table::data.table(value = levels), counts, by.x = "value", by.y = key, all.x = TRUE, sort = FALSE)
    for (j in c("Evidence_rows", "Peptide_sequences", "Razor_representatives")) data.table::set(out, which(is.na(out[[j]])), j, 0L)
    if (complete) out[is.na(Verified_assigned_groups), Verified_assigned_groups := 0L]
    data.table::setnames(out, "value", key)
    out[, `:=`(Stage = stage, Completeness_denominator_sequences = universe,
      Peptide_completeness = if (universe > 0L) Peptide_sequences / universe else NA_real_)]
    out
  }
  for (key in c("Raw.file", "Experiment", "QC_condition")) {
    out <- data.table::rbindlist(list(scoped(before, "Before_QC", key, nonempty_count(before$Sequence)),
      scoped(legacy, "Previous_QC", key, nonempty_count(legacy$Sequence)),
      scoped(ev, if (complete) "After_complete_QC" else "After_available_QC_INCOMPLETE", key, nonempty_count(ev$Sequence))))
    emit(out, paste0("QC_before_after_by_", key))
  }
  status <- data.table::data.table(QC_complete = complete,
    Protein_groups_file = protein_groups_file,
    Evidence_q_max = evidence_q_max, Protein_group_q_max = protein_q_max,
    Peptide_recurrence_filter_enabled = filter_peptide_recurrence,
    Protein_accession_peptide_count_filter_enabled = filter_protein_peptide_count,
    Minimum_raw_file_fraction = min_raw_file_fraction,
    Total_raw_files = total_raw_files, Minimum_raw_files = minimum_raw_files,
    Protein_peptide_count_scope = "each accession in the existing SEA mapping, before raw-file detection filtering",
    Protein_peptide_count_action = "exclude accessions and rows mapped only to excluded accessions",
    Filter_order = "protein peptide support, then peptide raw-file detection; no second protein filter",
    Minimum_biological_replicates = min_bio_replicates, Minimum_distinct_peptides = min_protein_peptides,
    Missing_filters = if (complete) "" else "protein-group Reverse/Decoy, contaminant, Only identified by site, protein-group Q-value; group assignment unverified")
  emit(status, "QC_filter_status")
  emit(ev[, .(QC_row, Sequence, Raw.file, Experiment, Original.QC.group, QC_group,
              Protein_group_retained, Leading.razor.proteins, Q.value)], "QC_retained_evidence")
  list(evidence = ev, complete = complete, stages = stages, coverage = coverage, status = status,
       n_input = nrow(before), n_previous_qc = nrow(legacy))
}
