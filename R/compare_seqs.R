#' @title Internal sequence alignment and exon utility functions
#' @description
#' A collection of internal helpers for alignment, exon parsing, and
#' length summarization used across SpliceImpactR.
#' These are not exported.
#' @name internal_helpers
#' @noRd
#' @keywords internal

#' Compute percent identity of an alignment
#' @param aln A \link[pwalign]{pairwiseAlignment} object.
#' @importFrom pwalign alignedPattern alignedSubject
#' @return Numeric percent identity (0–100) or NA.
.aln_pid <- function(aln) {
  ap <- as.character(pwalign::alignedPattern(aln))
  as <- as.character(pwalign::alignedSubject(aln))
  if (!nzchar(ap) || !nzchar(as)) return(NA_real_)
  a <- strsplit(ap, "", fixed = TRUE)[[1]]
  b <- strsplit(as, "", fixed = TRUE)[[1]]
  keep <- (a != "-" & b != "-")
  if (!any(keep)) return(NA_real_)
  100 * sum(a[keep] == b[keep]) / sum(keep)
}

#' Align two DNA sequences and compute identity, score, and width.
#' @param a,b Character strings of DNA sequences.
#' @return List with fields \code{pid}, \code{score}, \code{width}.
#' @importFrom pwalign pairwiseAlignment
#' @importFrom Biostrings DNAString
#' @noRd
#' @keywords internal
.align_dna <- function(a, b, alignmentMat) {
  # Coerce to character safely
  a <- if (is(a, "DNAString")) as.character(a) else as.character(a)
  b <- if (is(b, "DNAString")) as.character(b) else as.character(b)

  # Guard against NA / empty
  if (is.null(a) || is.null(b) ||
      anyNA(c(a, b)) ||
      !nzchar(a) || !nzchar(b)) {
    return(list(pid = NA_real_, score = NA_real_, width = NA_integer_))
  }

  # Now rebuild DNAString objects after cleaning
  a <- DNAString(toupper(a))
  b <- DNAString(toupper(b))

  # Fast path: identical
  if (identical(a, b)) {
    n <- length(a)  # works on DNAString
    return(list(pid = 100, score = n, width = n))
  }

  aln <- pwalign::pairwiseAlignment(
    a, b, type = "global",
    substitutionMatrix = alignmentMat,
    gapOpening = 10, gapExtension = 0.5
  )

  list(
    pid   = .aln_pid(aln),
    score = pwalign::score(aln),
    width = as.integer(width(pwalign::aligned(aln))[1])
  )
}

#' Align two protein sequences and compute identity, score, and width.
#' @param a,b Character strings of amino acid sequences.
#' @return List with fields \code{pid}, \code{score}, \code{width}.
#'
#' @importFrom Biostrings AAString
#' @noRd
#' @keywords internal
.align_aa <- function(a, b, alignmentMat) {
  # Coerce to character safely
  a <- if (is(a, "AAString")) as.character(a) else as.character(a)
  b <- if (is(b, "AAString")) as.character(b) else as.character(b)

  if (is.null(a) || is.null(b) ||
      anyNA(c(a, b)) ||
      !nzchar(a) || !nzchar(b)) {
    return(list(pid = NA_real_, score = NA_real_, width = NA_integer_))
  }

  a <- AAString(toupper(a))
  b <- AAString(toupper(b))

  if (identical(a, b)) {
    n <- length(a)
    return(list(pid = 100, score = n, width = n))
  }

  aln <- pwalign::pairwiseAlignment(
    a, b, type = "global",
    substitutionMatrix = alignmentMat,
    gapOpening = 10, gapExtension = 0.5
  )

  list(
    pid   = .aln_pid(aln),
    score = pwalign::score(aln),
    width = as.integer(width(pwalign::aligned(aln))[1])
  )
}

#' Sum exon CDS and feature lengths for each event row
#' @param H Event table containing exon IDs.
#' @param annotations Annotation data.table with exon_id, cds_len, and feature_length.
#' @param col_exons Column name in \code{H} containing exon identifiers.
#' @param out_prefix Character vector specifying which outputs to compute ("cds", "exon").
#' @param exon_delim Regex for splitting exon ID lists.
#' @noRd
#' @keywords internal
#' @return data.table with summed exon and CDS lengths per event.
.sum_exon_lengths <- function(H,
                              annotations,
                              col_exons,
                              out_prefix = c("cds","exon"),
                              exon_delim = "[,;|[:space:]]+") {

  L <- as.data.table(annotations)[, .(exon_id = as.character(exon_id),
                                      cds_len = as.integer(cds_len),
                                      feature_length = as.integer(feature_length))]
  setkey(L, exon_id)

  DT <- as.data.table(H)[, .(row_id = .I,
                             event_type = event_type_inc,
                             exons = get(col_exons))]

  # explode to long
  LONG <- DT[, {
    s <- as.character(exons); s[is.na(s)] <- ""
    ids <- unlist(strsplit(s, exon_delim))
    ids <- ids[nzchar(ids)]
    .(exon_id = ids)
  }, by = .(row_id, event_type)]

  if (!nrow(LONG)) {
    # no exon ids: return NA vectors aligned to H
    res <- data.table(row_id = seq_len(nrow(H)))
    if ("cds" %in% out_prefix)  res[, cds_len := NA_integer_]
    if ("exon" %in% out_prefix) res[, exon_len := NA_integer_]
    return(res[order(row_id)])
  }

  # join once, then sum
  J <- L[LONG, on = "exon_id", nomatch = 0L]
  SUM <- J[, .(
    cds_len   = sum(cds_len,         na.rm = TRUE),
    exon_len  = sum(feature_length,  na.rm = TRUE)
  ), by = row_id]

  # ensure every input row has an entry
  ALL <- data.table(row_id = seq_len(nrow(H)))[SUM, on = "row_id"]
  if ("cds" %notin% out_prefix)  ALL[, cds_len := NULL]
  if ("exon" %notin% out_prefix) ALL[, exon_len := NULL]
  ALL[order(row_id)]
}


#' Negated %in% operator
#' @noRd
#' @keywords internal
`%notin%` <- function(x, y) !(x %in% y)

#' Classify pair by protein-coding status
#' @param p1,p2 Peptide or protein IDs (character).
#' @return "protein_coding", "onePC", or "noPC".
#' @noRd
#' @keywords internal
.pc_class <- function(p1, p2) {
  has1 <- !is.na(p1) && nzchar(p1)
  has2 <- !is.na(p2) && nzchar(p2)
  if (has1 && has2) "protein_coding" else if (has1 || has2) "onePC" else "noPC"
}

#' Parse exon identifier strings
#' @param x Character vector of delimited exon IDs.
#' @return Character vector of parsed exon IDs.
#' @noRd
#' @keywords internal
.parse_exon_ids <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  y <- as.character(x); y[is.na(y)] <- ""
  out <- unlist(strsplit(y, "[,;|[:space:]]+"))
  out <- trimws(out)
  out[nzchar(out)]
}

#' Safe version of nchar() returning NA for missing or empty strings
#' @param z Character vector.
#' @return Integer vector of string lengths or NA..
#' @return Character vector of parsed exon IDs.
#' @noRd
#' @keywords internal
.safe_nchar <- function(z) ifelse(is.na(z) | !nzchar(z), NA_integer_, nchar(z))




#' Compare isoform nucleotide and protein sequences by pairwise alignment
#'
#' Performs transcript- and protein-level global alignments between
#' included and excluded isoform sequences to quantify sequence similarity,
#' coding differences, and exon coverage differences.
#'
#' This function wraps internal helpers for alignment
#' (\code{.align_dna()}, \code{.align_aa()}),
#' exon length summarization (\code{.sum_exon_lengths()}),
#' and protein-coding classification (\code{.pc_class()}).
#'
#' @param hits A \link[data.table]{data.table} or data.frame containing
#'   isoform pairs. Must include columns:
#'   \code{transcript_seq_inc}, \code{transcript_seq_exc},
#'   \code{protein_seq_inc}, and \code{protein_seq_exc}. From prior analysis.
#' @param annotations output from get_annotations (annotations)
#' @param include_sequences Logical; if \code{TRUE}, retains raw sequences in the output.
#'   Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, prints processing messages.
#'
#' @return A \link[data.table]{data.table} with added columns describing
#'   sequence length, coding length, and alignment similarity metrics, including:
#'   \itemize{
#'     \item \code{pc_class}: protein-coding status ("protein_coding", "onePC", "noPC")
#'     \item \code{prot_len_inc/exc}: protein sequence lengths
#'     \item \code{tx_len_inc/exc}: transcript sequence lengths
#'     \item \code{exon_cds_len_*}, \code{exon_len_*}: summed exon/CDS lengths
#'     \item \code{dna_pid}, \code{dna_score}, \code{prot_pid}, \code{prot_score}: alignment metrics
#'   }
#'
#'
#' @examples
#' sample_frame <- data.frame(path = c(check_extdata_dir('rawData/control_S5/'),
#'                                     check_extdata_dir('rawData/control_S6/'),
#'                                     check_extdata_dir('rawData/control_S7/'),
#'                                     check_extdata_dir('rawData/control_S8/'),
#'                                     check_extdata_dir('rawData/case_S1/'),
#'                                     check_extdata_dir('rawData/case_S2/'),
#'                                     check_extdata_dir('rawData/case_S3/'),
#'                                     check_extdata_dir('rawData/case_S4/')),
#'                            sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
#'                            condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
#'                            stringsAsFactors = FALSE)
#' hit_index <- get_hitindex(sample_frame)
#' res <- get_differential_inclusion(hit_index)
#' annots <- get_annotation(load = "test")
#' matched <- get_matched_events_chunked(res, annots$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annots$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' aligned <- compare_sequences_alignment(pairs, annots$annotations)
#'
#' @export
#' @importFrom data.table as.data.table setcolorder
#' @importFrom pwalign nucleotideSubstitutionMatrix
#' @importFrom utils data
compare_sequences_alignment <- function(hits, annotations, include_sequences = FALSE, verbose = TRUE) {
  data("BLOSUM62", package="Biostrings", envir=environment())
  NUC44 <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  valid_aa <- colnames(BLOSUM62)
  valid_nt <- colnames(NUC44)

  ## Fix aa and nt sequences in rare occurence of nonBLOSUM62/NUC44 chars
  hits$protein_seq_inc <- vapply(
    hits$protein_seq_inc,
    function(s) {
      if (is.na(s)) return(NA_character_)
      paste0(strsplit(s, "")[[1]][ strsplit(s, "")[[1]] %chin% valid_aa ], collapse = "")
    },
    FUN.VALUE = character(1)
  )
  hits$protein_seq_exc <- vapply(
    hits$protein_seq_exc,
    function(s) {
      if (is.na(s)) return(NA_character_)
      paste0(strsplit(s, "")[[1]][ strsplit(s, "")[[1]] %chin% valid_aa ], collapse = "")
    },
    FUN.VALUE = character(1)
  )

  hits$transcript_seq_inc <- vapply(
    hits$transcript_seq_inc,
    function(s) {
      if (is.na(s)) return(NA_character_)
      paste0(strsplit(s, "")[[1]][ strsplit(s, "")[[1]] %chin% valid_nt ], collapse = "")
    },
    FUN.VALUE = character(1)
  )

  hits$transcript_seq_exc <- vapply(
    hits$transcript_seq_exc,
    function(s) {
      if (is.na(s)) return(NA_character_)
      paste0(strsplit(s, "")[[1]][ strsplit(s, "")[[1]] %chin% valid_nt ], collapse = "")
    },
    FUN.VALUE = character(1)
  )

  DT <- as.data.table(hits)

  if (verbose == TRUE) print(paste0("[INFO] Processing ", dim(DT)[1], " transcript and protein sequence alignments, this may take a little..."))
  inc_sum <- .sum_exon_lengths(DT, annotations, col_exons = "exons_inc",
                               out_prefix = c("cds","exon"),
                               exon_delim = "[,;|[:space:]]+")

  exc_sum <- .sum_exon_lengths(DT, annotations, col_exons = "exons_exc",
                               out_prefix = c("cds","exon"),
                               exon_delim = "[,;|[:space:]]+")

  # lengths & classes
  prot_len_inc <- .safe_nchar(DT$protein_seq_inc)
  prot_len_exc <- .safe_nchar(DT$protein_seq_exc)
  tx_len_inc   <- .safe_nchar(DT$transcript_seq_inc)
  tx_len_exc   <- .safe_nchar(DT$transcript_seq_exc)
  pc_class     <- mapply(.pc_class, DT$protein_seq_inc, DT$protein_seq_exc)

  DT[, transcript_seq_inc := vapply(transcript_seq_inc, as.character, character(1))]
  DT[, transcript_seq_exc := vapply(transcript_seq_exc, as.character, character(1))]
  DT[, protein_seq_inc    := vapply(protein_seq_inc,    as.character, character(1))]
  DT[, protein_seq_exc    := vapply(protein_seq_exc,    as.character, character(1))]

  # handle rows with missing sequences
  dna_res <- vector("list", nrow(DT))
  aa_res  <- vector("list", nrow(DT))

  valid_dna <- !is.na(DT$transcript_seq_inc) & !is.na(DT$transcript_seq_exc)
  valid_aa  <- !is.na(DT$protein_seq_inc)    & !is.na(DT$protein_seq_exc)

  # fill non-valid rows with default NA result list
  for (i in which(!valid_dna)) {
    dna_res[[i]] <- list(pid = NA_real_, score = NA_real_, width = NA_integer_)
  }

  for (i in which(!valid_aa)) {
    aa_res[[i]] <- list(pid = NA_real_, score = NA_real_, width = NA_integer_)
  }

  # align only valid rows
  dna_res[valid_dna] <- mapply(
    function(a, b) .align_dna(a, b, NUC44),
    DT$transcript_seq_inc[valid_dna],
    DT$transcript_seq_exc[valid_dna],
    SIMPLIFY = FALSE
  )

  aa_res[valid_aa] <- mapply(
    function(a, b) .align_aa(a, b, BLOSUM62),
    DT$protein_seq_inc[valid_aa],
    DT$protein_seq_exc[valid_aa],
    SIMPLIFY = FALSE
  )

  DT[, `:=`(
    pc_class   = as.character(pc_class),

    prot_len_inc      = prot_len_inc,
    prot_len_exc      = prot_len_exc,
    prot_len_diff = prot_len_inc - prot_len_exc,
    prot_len_diff_abs = abs(prot_len_inc - prot_len_exc),

    tx_len_inc        = tx_len_inc,
    tx_len_exc        = tx_len_exc,
    tx_len_diff   = tx_len_inc - tx_len_exc,
    tx_len_diff_abs   = abs(tx_len_inc - tx_len_exc),

    # coding lengths (CDS)
    exon_cds_len_inc       = inc_sum$cds_len,
    exon_cds_len_exc       = exc_sum$cds_len,
    exon_cds_len_diff      = inc_sum$cds_len - exc_sum$cds_len,
    exon_cds_len_diff_abs  = abs(inc_sum$cds_len - exc_sum$cds_len),

    # total exon lengths (feature length)
    exon_len_inc      = inc_sum$exon_len,
    exon_len_exc      = exc_sum$exon_len,
    exon_len_diff     = inc_sum$exon_len - exc_sum$exon_len,
    exon_len_diff_abs = abs(inc_sum$exon_le - exc_sum$exon_len),

    dna_pid    = vapply(dna_res, `[[`, numeric(1), "pid"),
    dna_score  = vapply(dna_res, `[[`, numeric(1), "score"),
    dna_width  = vapply(dna_res, `[[`, integer(1), "width"),

    prot_pid   = vapply(aa_res,  `[[`, numeric(1), "pid"),
    prot_score = vapply(aa_res,  `[[`, numeric(1), "score"),
    prot_width = vapply(aa_res,  `[[`, integer(1), "width")

  )]

  if (!isTRUE(include_sequences)) {
    keep <- setdiff(names(DT),
                    c("transcript_seq_inc","transcript_seq_exc",
                      "protein_seq_inc","protein_seq_exc"))
    DT <- DT[, ..keep]
  }

  metric_cols <- c("dna_pid","dna_score","dna_width",
                   "prot_pid","prot_score","prot_width")
  data.table::setcolorder(DT, c(setdiff(names(DT), metric_cols), metric_cols))
  DT[]
}








