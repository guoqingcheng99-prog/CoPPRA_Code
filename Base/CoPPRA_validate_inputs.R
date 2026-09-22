# Verify the paired MaxQuant files without changing evidence rows or mappings.
validate_coppra_input_pair <- function(evidence, protein_groups_file) {
  if (!file.exists(protein_groups_file)) stop("Missing protein-group input: ", protein_groups_file)
  pg <- data.table::fread(protein_groups_file, na.strings = c("", "NA"),
    select = c("id", "Protein IDs", "Peptide IDs", "Peptide sequences"))
  data.table::setnames(pg, make.names(names(pg)))
  pg[, id := as.character(id)]
  if (anyNA(pg$id) || anyDuplicated(pg$id)) stop("Invalid protein-group IDs in paired input.")
  peptide_ids <- strsplit(as.character(pg$Peptide.IDs), ";", fixed = TRUE)
  sequences <- strsplit(as.character(pg$Peptide.sequences), ";", fixed = TRUE)
  if (any(lengths(peptide_ids) != lengths(sequences))) stop("Unpaired peptide IDs/sequences in proteinGroups.")
  known <- unique(data.table::data.table(
    GroupID = rep(pg$id, lengths(peptide_ids)),
    Peptide.ID = unlist(peptide_ids, use.names = FALSE),
    Sequence = unlist(sequences, use.names = FALSE)))
  observed <- unique(evidence[, .(Protein.group.IDs, Peptide.ID = as.character(Peptide.ID), Sequence)])
  observed <- observed[!is.na(Protein.group.IDs) & nzchar(Protein.group.IDs),
    .(GroupID = trimws(unlist(strsplit(as.character(Protein.group.IDs), ";", fixed = TRUE)))),
    by = .(Peptide.ID, Sequence)]
  if (length(setdiff(observed$GroupID, pg$id)) ||
      nrow(observed[!known, on = .(GroupID, Peptide.ID, Sequence)])) {
    stop("Evidence group IDs and peptide IDs/sequences do not match proteinGroups; use the paired run.")
  }
  invisible(TRUE)
}
