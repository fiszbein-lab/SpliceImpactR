rank_test_reference <- function(strand = "+") {
  definitions <- list(
    A_INC_distant = matrix(c(101L, 400L, 1001L, 1600L), ncol = 2L, byrow = TRUE),
    E_EXC = matrix(c(101L, 199L, 1001L, 1300L), ncol = 2L, byrow = TRUE),
    Z_INC_close = matrix(c(101L, 400L, 1001L, 1300L), ncol = 2L, byrow = TRUE)
  )
  annotations <- sequences <- list()
  for (id in names(definitions)) {
    spans <- definitions[[id]]
    exon_seq <- c(paste0("ATG", strrep("GCT", (spans[1L, 2L] - spans[1L, 1L] - 2L) / 3L)),
                  paste0(strrep("GCT", 100L), if (id == "A_INC_distant") strrep("CCC", 100L) else ""))
    if (strand == "-") {
      spans <- cbind(2001L - spans[, 2L], 2001L - spans[, 1L])
    }
    lengths <- spans[, 2L] - spans[, 1L] + 1L
    offset <- c(0L, head(cumsum(lengths), -1L))
    ex <- data.table::data.table(type = "exon", transcript_id = id,
      exon_id = paste0(id, "_", 1:2), exon_number = 1:2,
      start = spans[, 1L], end = spans[, 2L], cds_gen_start = spans[, 1L],
      cds_gen_stop = spans[, 2L], cds_rel_start = offset + 1L, cds_rel_stop = cumsum(lengths),
      cds_has = TRUE, cds_len = lengths, feature_length = lengths,
      absolute_exon_position = 1:2, coding_exon_position = 1:2,
      absolute_exon_class = c("first", "last"), coding_exon_class = c("first", "last"),
      start_frame = offset %% 3L, stop_frame = (cumsum(lengths) - 1L) %% 3L)
    tr <- data.table::data.table(type = "transcript", transcript_id = id,
      start = min(ex$start), end = max(ex$end))
    ann <- data.table::rbindlist(list(tr, ex), fill = TRUE)
    ann[, `:=`(gene_id = "G", gene_name = "G", chr = "chr1", strand = strand,
                transcript_name = id, transcript_type = "protein_coding",
                transcript_support_level = "1", protein_id = paste0("P_", id))]
    annotations[[id]] <- ann
    dna <- paste(exon_seq, collapse = "")
    sequences[[id]] <- data.table::data.table(transcript_id = id, protein_id = paste0("P_", id),
      transcript_seq = dna, protein_seq = as.character(Biostrings::translate(Biostrings::DNAString(dna))))
  }
  ann <- data.table::rbindlist(annotations, fill = TRUE)[, row_uid := .I]
  events <- data.table::data.table(event_id = "E", event_type = "A5SS", form = c("INC", "EXC"),
    gene_id = "G", chr = "chr1", strand = strand,
    inc = if (strand == "+") c("101-400", "101-199") else c("1601-1900", "1802-1900"),
    exc = if (strand == "+") c("", "200-400") else c("", "1601-1801"),
    delta_psi = c(0.3, -0.3), p.value = 0.001, padj = 0.01,
    n_samples = 4L, n_control = 2L, n_case = 2L)
  list(events = events, annotations = ann, sequences = data.table::rbindlist(sequences))
}

# Plus-strand transcripts from exon matrices. `cds` is one coding span for all
# transcripts or a list of spans by transcript ID. Coding bases come from one
# shared genome, so identical positions give identical sequence.
rank_structure_reference <- function(transcripts, cds = NULL) {
  genome <- paste(rep(c("ACGTTGCA", "GATTACAC", "CCGGTTAA"), length.out = 150L), collapse = "")
  ann <- data.table::rbindlist(lapply(names(transcripts), function(id) {
    e <- transcripts[[id]]
    ex <- data.table::data.table(type = "exon", transcript_id = id, start = e[, 1L],
                                 end = e[, 2L], exon_id = paste0(id, "_", seq_len(nrow(e))))
    span <- if (is.list(cds)) cds[[id]] else cds
    if (!is.null(span)) {
      ex[, `:=`(cds_gen_start = pmax(start, span[1L]), cds_gen_stop = pmin(end, span[2L]))]
      ex[cds_gen_start > cds_gen_stop, `:=`(cds_gen_start = NA_real_, cds_gen_stop = NA_real_)]
    }
    tr <- data.table::data.table(type = "transcript", transcript_id = id,
                                 start = min(e[, 1L]), end = max(e[, 2L]))
    data.table::rbindlist(list(tr, ex), fill = TRUE)
  }), fill = TRUE)
  ann[, `:=`(gene_id = "G", gene_name = "G", chr = "chr1", strand = "+",
             transcript_name = transcript_id, transcript_type = "protein_coding",
             transcript_support_level = "1", protein_id = paste0("P_", transcript_id),
             row_uid = .I)]
  coding <- if (is.null(cds)) NULL else ann[type == "exon" & !is.na(cds_gen_start)]
  sequences <- data.table::data.table(transcript_id = names(transcripts),
    protein_id = paste0("P_", names(transcripts)), transcript_seq = NA_character_, protein_seq = "M")
  if (!is.null(coding)) {
    sequences[, transcript_seq := vapply(transcript_id, function(id) {
      x <- coding[transcript_id == id][order(start)]
      if (!nrow(x)) return(NA_character_)
      paste(substring(genome, x$cds_gen_start, x$cds_gen_stop), collapse = "")
    }, character(1))]
  }
  list(annotations = ann, sequences = sequences)
}

rank_events <- function(event_id, event_type, form, inc, exc, delta_psi) {
  data.table::data.table(event_id = event_id, event_type = event_type, form = form,
    gene_id = "G", chr = "chr1", strand = "+", inc = inc, exc = exc,
    delta_psi = delta_psi, p.value = 0.001, padj = 0.01,
    n_samples = 4L, n_control = 2L, n_case = 2L)
}
