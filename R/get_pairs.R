#' Pair inclusion and exclusion forms of splicing events
#'
#' Builds paired tables of inclusion/exclusion forms for splicing events from
#' rMATS-like or HITindex-like inputs. In rMATS mode, events are paired when
#' both INC and EXC forms exist for a given event ID. In HITindex mode, all
#' positive and negative deltaPSI rows within each event are cross-joined.
#'
#' @param x A data.frame, data.table, or `SpliceImpactResult` containing
#'   splicing event information.
#' @param source Character string specifying input structure:
#'   \describe{
#'     \item{\code{"paired"}}{(rMATS-like) requires exactly one INC and one EXC
#'       per event ID.}
#'     \item{\code{"multi"}}{(HITindex-like) pairs all positive and negative
#'       \code{delta_psi} values within each event.}
#'   }
#' @param return_class Character. Output mode: `"data.table"`, `"S4"`, or
#'   `"auto"` (default). In `auto`, S4 input returns updated S4 output.
#'
#' @return A \link[data.table]{data.table} (or updated `SpliceImpactResult`
#' when `return_class` resolves to S4) where each row represents an
#' inclusion-exclusion pair of the same event. S4 output sets
#' `metadata$matching` to `"legacy"` and removes `matching_diagnostics` left by
#' [get_ranked_pairs()].
#' @details
#' In \code{source="paired"} mode, only events with exactly one INC and one EXC
#' row are retained. When delta PSI is supplied, the positive form is the case
#' and the negative form is the control; without it, INC is the case.
#' In \code{source="multi"} mode, all positive deltaPSI rows are
#' joined with all negative deltaPSI rows (cartesian join) within each event.
#' Missing sample-count diagnostics are retained as `NA`. With matched transcript
#' input, both modes return the same columns, including when no pairs remain.
#' A `transcript_type` column on the input (added by both matchers) becomes
#' `transcript_type_case` and `transcript_type_control`, for example to show
#' which selected transcripts are annotated as NMD.
#'
#' An event with several rising or falling sites, such as an AFE with three
#' first exons, gives one comparison per rising x falling pair of sites, so a
#' site can appear in several rows. `n_event_comparisons` gives the number of
#' comparisons its event defines (1 for a two-form event). After
#' [keep_sig_pairs()], `site_significant_case` and `site_significant_control`
#' show whether each form or site met the thresholds itself; the other sites
#' of a passing event are kept as comparators. Counts of pairs, including the
#' domain-enrichment foreground, are therefore counts of comparisons.
#' Forms without a matched transcript ID cannot produce a pair. An unmatched
#' form is not interpreted as an absent or noncoding isoform.
#'
#' @importFrom data.table as.data.table setkeyv setnames setorderv
#' @export
#'
#' @examples
#' ex <- load_example_data("sample_frame")
#' sample_frame <- ex$sample_frame
#' hit_index <- get_hitindex(sample_frame)
#' res <- get_differential_inclusion(hit_index)
#' annots <- load_example_data("annotation_df")$annotation_df
#' matched <- get_matched_events_chunked(res, annots$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annots$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' print(pairs)
get_pairs <- function(x,
                      source = c("paired","multi"),
                      return_class = c("auto", "data.table", "S4")) {

  source  <- match.arg(source)
  return_class <- match.arg(return_class)
  if (methods::is(x, "SpliceImpactResult")) {
    .spi_obj <- x
    DT <- as.data.table(as_dt_from_s4(x, "matched"))
    if (!ncol(DT)) DT <- as.data.table(as_dt_from_s4(x, "res_di"))
    if (!ncol(DT)) DT <- as.data.table(as_dt_from_s4(x, "di_events"))
  } else {
    .spi_in <- .resolve_splice_input(x, what = "di_events")
    .spi_obj <- .spi_in$obj
    DT <- data.table::copy(as.data.table(.spi_in$dt))
  }

  required <- if (source == "paired") c("event_id", "form") else c(
    "event_id", "gene_id", "transcript_id", "chr", "strand", "event_type",
    "form", "exons", "protein_id", "inc", "exc", "delta_psi",
    "p.value", "padj", "transcript_seq", "protein_seq"
  )
  missing <- setdiff(required, names(DT))
  if (length(missing)) {
    stop("get_pairs(source='", source, "') missing required columns: ",
         paste(missing, collapse = ", "),
         ". Run annotation matching + sequence attachment before pairing.")
  }
  if (anyNA(DT$event_id) || any(!nzchar(trimws(as.character(DT$event_id))))) {
    stop("event_id must be non-missing and non-empty.")
  }

  # Imported DI may not have sample counts; unavailable is not zero.
  for (column in c("n_samples", "n_control", "n_case")) {
    if (!column %in% names(DT)) DT[, (column) := NA_integer_]
  }
  shared <- intersect(c("gene_id", "chr", "strand", "event_type"), names(DT))
  if (length(shared) && nrow(DT)) {
    # An event with more than one distinct metadata combination is inconsistent.
    combos <- unique(DT, by = c("event_id", shared))$event_id
    conflicting <- unique(combos[duplicated(combos)])
    if (length(conflicting)) {
      conflicting <- unique(DT$event_id)[unique(DT$event_id) %chin% conflicting]
      stop("An event_id must identify one gene, chromosome, strand and event type. Conflicting IDs: ",
           paste(utils::head(conflicting, 5L), collapse = ", "))
    }
  }

  if (source == "paired") {
    # A paired event has exactly one INC and one EXC row.
    eligible <- DT[, .(keep = .N == 2L && sum(form == "INC", na.rm = TRUE) == 1L &&
                         sum(form == "EXC", na.rm = TRUE) == 1L), by = event_id]
    DT <- DT[event_id %chin% eligible[keep == TRUE, event_id]]
  }

  if ("delta_psi" %in% names(DT)) {
    if (!is.numeric(DT$delta_psi)) stop("delta_psi must be numeric.")
    signs <- DT[, .(keep = any(delta_psi > 0, na.rm = TRUE) &&
                     any(delta_psi < 0, na.rm = TRUE)), by = event_id]
    DT <- DT[event_id %chin% signs[keep == TRUE, event_id]]
    POS <- data.table::copy(DT[!is.na(delta_psi) & delta_psi > 0])
    NEG <- data.table::copy(DT[!is.na(delta_psi) & delta_psi < 0])
  } else {
    # Preserve form-based pairing for callers without differential statistics.
    POS <- data.table::copy(DT[form == "INC"])
    NEG <- data.table::copy(DT[form == "EXC"])
  }

  # A multi-site event defines one comparison per rising x falling site pair,
  # whether or not each form found a transcript.
  n_comparisons <- merge(POS[, .(n_pos = .N), by = event_id], NEG[, .(n_neg = .N), by = event_id],
                         by = "event_id")[, .(event_id, n = n_pos * n_neg)]

  # A failed annotation match is not evidence for an absent/noncoding isoform.
  if ("transcript_id" %in% names(DT)) {
    POS <- POS[!is.na(transcript_id) & nzchar(transcript_id)]
    NEG <- NEG[!is.na(transcript_id) & nzchar(transcript_id)]
  }
  left_cols <- setdiff(names(POS), "event_id")
  right_cols <- setdiff(names(NEG), "event_id")
  setnames(POS, left_cols, paste0(left_cols, "_case"))
  setnames(NEG, right_cols, paste0(right_cols, "_control"))
  out <- merge(POS, NEG, by = "event_id", allow.cartesian = TRUE)

  if ("delta_psi" %in% names(DT)) {
    out[, `:=`(ordA = -abs(delta_psi_case), ordB = -abs(delta_psi_control))]
    data.table::setorderv(out, c("event_id", "ordA", "ordB"))
    out[, c("ordA", "ordB") := NULL]
  } else {
    data.table::setorderv(out, "event_id")
  }

  # Both modes use the established multi-mode result schema when given matched
  # transcripts, including a typed zero-row result when there are no pairs.
  form_cols <- c("form", "exons", "protein_id", "inc", "exc", "delta_psi",
                 "p.value", "padj", "n_samples", "n_control", "n_case",
                 "transcript_seq", "protein_seq", "transcript_type", "site_significant")
  paired_cols <- unlist(lapply(form_cols, function(column) paste0(column, c("_control", "_case"))))
  cols_old <- c("event_id", "gene_id_control", "transcript_id_control", "transcript_id_case",
                "chr_control", "strand_control", "event_type_control", paired_cols)
  cols_new <- c("event_id", "gene_id", "transcript_id_control", "transcript_id_case",
                "chr", "strand", "event_type", paired_cols)
  present <- cols_old %in% names(out)
  cols_old <- cols_old[present]
  cols_new <- cols_new[present]
  out <- out[, ..cols_old]
  data.table::setnames(out, cols_old, cols_new)
  out[, n_event_comparisons := n_comparisons$n[match(event_id, n_comparisons$event_id)]]
  out <- .return_splice_output(out[], obj = .spi_obj, what = "paired_hits", return_class = return_class)
  # These pairs replace any from get_ranked_pairs(), so their provenance does too.
  if (methods::is(out, "SpliceImpactResult")) {
    out@metadata$matching <- "legacy"
    out@metadata$matching_diagnostics <- NULL
  }
  out
}
