#' Parse closed event intervals for the alternate matcher (internal)
#' @keywords internal
#' @noRd
.si_rank_intervals <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(IRanges::IRanges())
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1L]])
  if (any(!grepl("^[0-9]+-[0-9]+$", parts))) {
    stop("Event intervals must be one-based closed start-end spans separated by semicolons.")
  }
  bounds <- do.call(rbind, strsplit(parts, "-", fixed = TRUE))
  lo <- as.numeric(bounds[, 1L]); hi <- as.numeric(bounds[, 2L])
  if (any(!is.finite(lo) | !is.finite(hi) | lo < 1 | hi < lo |
          hi > .Machine$integer.max)) stop("Invalid event interval bounds.")
  out <- IRanges::IRanges(start = lo, end = hi)
  out[order(IRanges::start(out), IRanges::end(out))]
}

#' Stop on an input problem, naming the events that have it (internal)
#' @keywords internal
#' @noRd
.si_rank_input_error <- function(problem, ids) {
  ids <- unique(as.character(ids))
  more <- if (length(ids) > 5L) sprintf(" and %d more", length(ids) - 5L) else ""
  stop(problem, " Affected: ", paste(utils::head(ids, 5L), collapse = ", "), more, ".", call. = FALSE)
}

#' Event-side boundary positions for relaxed structural matching (internal)
#'
#' Outer ends of flanking exons, the shared A3SS/A5SS boundary and terminal
#' TSS/PAS ends are not part of an event. In relaxed matching, the intervals
#' that end there are mapped to the exon containing their event-side boundary.
#' @return Named positions keyed by interval index; empty when nothing relaxes.
#' @keywords internal
#' @noRd
.si_rank_relaxed_anchors <- function(event_type, strand, lo, hi) {
  ni <- length(lo)
  if (event_type %in% c("SE", "MXE", "RI", "A5SS", "A3SS") && ni >= 2L) {
    return(stats::setNames(c(hi[1L], lo[ni]), c(1L, ni)))
  }
  if (ni != 1L) return(stats::setNames(numeric(), character()))
  if (event_type %in% c("A5SS", "A3SS")) {
    variable_end <- (event_type == "A5SS" && strand == "+") ||
      (event_type == "A3SS" && strand == "-")
    return(stats::setNames(if (variable_end) hi else lo, 1L))
  }
  if (event_type %in% c("AFE", "HFE", "ALE", "HLE")) {
    splice_at_end <- xor(event_type %in% c("AFE", "HFE"), strand == "-")
    return(stats::setNames(if (splice_at_end) hi else lo, 1L))
  }
  stats::setNames(numeric(), character())
}

#' Check segment coverage and event-specific exon structure (internal)
#'
#' Exact matching requires every required interval inside one exon. Relaxed
#' matching keeps every event-defining boundary exact but lets the outer ends
#' named by [.si_rank_relaxed_anchors()] differ.
#' @keywords internal
#' @noRd
.si_rank_compatible <- function(ex, inc, exc, event_type, form, strand, relaxed = FALSE) {
  fail <- function(reason) list(eligible = FALSE, reason = reason,
                                context = NA_character_, exons = "", match = NA_character_)
  if (!length(inc)) return(fail("no_required_interval"))
  if (!nrow(ex)) return(fail("no_exons"))
  er <- IRanges::IRanges(ex$start, ex$end)
  if (length(exc) && length(IRanges::findOverlaps(er, exc))) {
    return(fail("forbidden_interval_overlap"))
  }
  ni <- length(inc); ne <- nrow(ex)
  lo <- IRanges::start(inc); hi <- IRanges::end(inc)
  hits <- IRanges::findOverlaps(inc, er, type = "within")
  mapping <- split(S4Vectors::subjectHits(hits), S4Vectors::queryHits(hits))
  ix <- vapply(seq_along(inc), function(i) {
    j <- mapping[[as.character(i)]]
    if (is.null(j)) NA_integer_ else j[1L]
  }, integer(1))
  if (relaxed) {
    anchor <- .si_rank_relaxed_anchors(event_type, strand, lo, hi)
    for (k in seq_along(anchor)) {
      j <- which(ex$start <= anchor[[k]] & ex$end >= anchor[[k]])
      ix[as.integer(names(anchor)[k])] <- if (length(j)) j[1L] else NA_integer_
    }
  }
  if (anyNA(ix)) return(fail("required_interval_not_covered"))
  context <- "interval_only"
  # Generic DI imports can preserve event type while using SITE labels.
  # Infer a branch only when the supplied interval topology identifies it.
  if (form == "SITE" && event_type == "SE" && ni == 2L && length(exc)) form <- "EXC"
  if (form == "SITE" && event_type == "RI" &&
      (ni == 1L || (ni == 3L && all(lo[-1L] == hi[-ni] + 1L)))) form <- "INC"

  # In genomic order, exact inner boundaries and adjacent annotated exons
  # encode the same junction chain on either strand.
  chain_ok <- function() {
    if (ni < 2L || any(diff(ix) != 1L)) return(FALSE)
    all(ex$end[ix[-ni]] == hi[-ni] & ex$start[ix[-1L]] == lo[-1L] &
          lo[-1L] > hi[-ni] + 1L)
  }
  if (event_type %in% c("SE", "MXE")) {
    expected <- if (event_type == "SE" && form == "EXC") 2L else 3L
    if (ni == expected) {
      if (!chain_ok()) return(fail("splice_junction_chain_mismatch"))
      context <- "junction_chain"
    } else if (ni == 1L && !(event_type == "SE" && form == "EXC")) {
      if (ix == 1L || ix == ne || ex$start[ix] != lo || ex$end[ix] != hi) {
        return(fail("internal_exon_boundary_mismatch"))
      }
      context <- "partial_event"
    } else return(fail("insufficient_event_context"))
  } else if (event_type == "RI") {
    if (form == "INC") {
      contiguous <- ni == 1L || all(lo[-1L] == hi[-ni] + 1L)
      if (!contiguous || length(unique(ix)) != 1L) {
        return(fail("retained_intron_not_in_one_exon"))
      }
      context <- if (ni == 3L) "retained_intron" else "partial_event"
    } else {
      if (ni != 2L || !chain_ok()) return(fail("retained_intron_spliced_chain_mismatch"))
      context <- "junction_chain"
    }
  } else if (event_type %in% c("A5SS", "A3SS")) {
    # The rMATS converters supply the variable exon and its partner, so the
    # junction is checked; a lone variable exon verifies only the splice site.
    if (!ni %in% c(1L, 2L)) return(fail("insufficient_event_context"))
    variable_end <- (event_type == "A5SS" && strand == "+") ||
      (event_type == "A3SS" && strand == "-")
    anchor <- if (ni == 1L || variable_end) 1L else ni
    j <- ix[anchor]
    boundary_ok <- if (variable_end) ex$end[j] == hi[anchor] && j < ne else
      ex$start[j] == lo[anchor] && j > 1L
    if (!boundary_ok) return(fail("alternative_splice_boundary_mismatch"))
    if (ni == 2L && !chain_ok()) return(fail("splice_junction_chain_mismatch"))
    context <- if (ni == 2L) "junction_chain" else "splice_site_only"
  } else if (event_type %in% c("AFE", "HFE", "ALE", "HLE")) {
    if (ni != 1L) return(fail("insufficient_event_context"))
    first <- event_type %in% c("AFE", "HFE")
    splice_at_end <- xor(first, strand == "-")
    terminal <- if (splice_at_end) 1L else ne
    # Relaxed: only the splice-site end must match; the TSS/PAS end may differ.
    ends_ok <- if (relaxed) {
      ne > 1L && (if (splice_at_end) ex$end[ix] == hi else ex$start[ix] == lo)
    } else {
      ex$start[ix] == lo && ex$end[ix] == hi
    }
    if (ix != terminal || !ends_ok) {
      return(fail("terminal_exon_boundary_mismatch"))
    }
    context <- "terminal_exon"
  }
  list(eligible = TRUE, reason = "eligible", context = context,
       exons = paste(ex$exon_id[ix], collapse = ";"),
       match = if (relaxed) "relaxed" else "exact")
}

#' Extract a coding sequence whose length agrees with annotation (internal)
#' @keywords internal
#' @noRd
.si_rank_cds <- function(ex, sequence, strand) {
  empty <- list(ranges = IRanges::IRanges(), sequence = NA_character_,
                status = "no_cds_annotation")
  if (!all(c("cds_gen_start", "cds_gen_stop") %in% names(ex))) return(empty)
  valid <- !is.na(ex$cds_gen_start) & !is.na(ex$cds_gen_stop) &
    ex$cds_gen_start >= ex$start & ex$cds_gen_stop <= ex$end &
    ex$cds_gen_stop >= ex$cds_gen_start
  cds <- ex[which(valid)]
  if (!nrow(cds)) return(empty)
  if (strand == "-") cds <- cds[rev(seq_len(nrow(cds)))]
  rr <- IRanges::IRanges(cds$cds_gen_start, cds$cds_gen_stop)
  result <- list(ranges = rr, sequence = NA_character_, status = "missing_sequence")
  if (is.na(sequence) || !nzchar(sequence)) return(result)
  sequence <- chartr("U", "T", toupper(sequence))
  if (!grepl("^[ACGTRYSWKMBDHVN]+$", sequence)) {
    result$status <- "invalid_nucleotide_sequence"
    return(result)
  }
  n_cds <- sum(IRanges::width(rr))
  if (nchar(sequence) == n_cds) {
    result$sequence <- sequence
    result$status <- "cds_length_agrees"
  } else if (nchar(sequence) == n_cds + 3L &&
             substr(sequence, n_cds + 1L, n_cds + 3L) %in% c("TAA", "TAG", "TGA")) {
    result$sequence <- substr(sequence, 1L, n_cds)
    result$status <- "terminal_stop_removed"
  } else if (nchar(sequence) == sum(ex$end - ex$start + 1L)) {
    ordered <- if (strand == "-") ex[rev(seq_len(nrow(ex)))] else ex
    offset <- c(0L, utils::head(cumsum(ordered$end - ordered$start + 1L), -1L))
    k <- match(cds$exon_id, ordered$exon_id)
    a <- offset[k] + if (strand == "+") cds$cds_gen_start - cds$start + 1L else
      cds$end - cds$cds_gen_stop + 1L
    result$sequence <- paste(substring(sequence, a, a + IRanges::width(rr) - 1L), collapse = "")
    result$status <- "cds_extracted_from_transcript"
  } else result$status <- "sequence_length_mismatch"
  result
}

#' Genomic span of a set of ranges (internal)
#' @keywords internal
#' @noRd
.si_rank_span <- function(r) {
  if (!length(r)) return(IRanges::IRanges())
  IRanges::IRanges(min(IRanges::start(r)), max(IRanges::end(r)))
}

#' Frame-agnostic structural similarity of two transcripts (internal)
#'
#' Compares exonic genomic positions inside `window`, minus `mask`. Positions
#' in one genome determine sequence, so this measures how much of the coding
#' region the two transcripts share without depending on where each annotated
#' reading frame ends: a frameshifted isoform keeps its downstream exons.
#' @return Base-pair Jaccard `similarity`, the shared fraction of each
#'   transcript's compared positions, and the number of compared positions.
#' @keywords internal
#' @noRd
.si_rank_positions <- function(a, b, window, mask = IRanges::IRanges()) {
  pa <- IRanges::intersect(IRanges::setdiff(IRanges::reduce(a), mask), window)
  pb <- IRanges::intersect(IRanges::setdiff(IRanges::reduce(b), mask), window)
  shared <- sum(IRanges::width(IRanges::intersect(pa, pb)))
  total <- sum(IRanges::width(IRanges::union(pa, pb)))
  na <- sum(IRanges::width(pa)); nb <- sum(IRanges::width(pb))
  list(similarity = if (total) shared / total else NA_real_,
       coverage_case = if (na) shared / na else NA_real_,
       coverage_control = if (nb) shared / nb else NA_real_,
       compared_nt = as.integer(total))
}

#' Fraction of a region covered by exons (internal)
#' @keywords internal
#' @noRd
.si_rank_coverage <- function(exons, region) {
  if (!length(region)) return(NA_real_)
  sum(IRanges::width(IRanges::intersect(region, exons))) / sum(IRanges::width(region))
}

#' Legacy-style approximate eligibility (internal)
#'
#' The legacy matcher's gates: every required interval overlaps one exon by at
#' least `min_overlap` of its width, and no exon overlaps a forbidden interval
#' by that much. Event boundaries are not checked.
#' @keywords internal
#' @noRd
.si_rank_approximate <- function(exons, inc, exc, min_overlap = 0.05) {
  if (!length(inc) || !length(exons)) return(FALSE)
  best <- function(q) {
    h <- IRanges::findOverlaps(q, exons)
    out <- numeric(length(q))
    if (length(h)) {
      w <- IRanges::width(IRanges::pintersect(q[S4Vectors::queryHits(h)], exons[S4Vectors::subjectHits(h)]))
      top <- tapply(w, S4Vectors::queryHits(h), max)
      out[as.integer(names(top))] <- top
    }
    out / IRanges::width(q)
  }
  eps <- 1e-9
  all(best(inc) + eps >= min_overlap) && (!length(exc) || all(best(exc) + eps < min_overlap))
}

#' Jointly rank compatible transcript pairs by coding-region structure
#'
#' An opt-in alternative to independent transcript selection. Candidates must
#' represent the supplied AS form before coding status, transcript support and
#' coding-region similarity are considered. The legacy matcher is unchanged as
#' the default in [get_splicing_impact()].
#'
#' @param events Canonical differential event table, or a [SpliceImpactResult]
#'   with differential events (its significance-filtered `res_di` events when
#'   present, otherwise `di_events`). Required columns: `event_id`, `event_type`, `form`,
#'   `gene_id`, `chr`, `strand`, `inc`, `exc`. Coordinates are one-based, closed.
#'   In multi mode, positive and negative `delta_psi` forms are paired within
#'   each event. Filter significance before calling this function if desired.
#' @param annotations Processed transcript/exon table from [get_annotation()].
#'   An optional `transcript_support_level` column sets TSL tiers (unknown when
#'   absent). An optional `transcript_type` column is reported as each selected
#'   transcript's biotype; it is not used for ranking. Exons require `exon_id`,
#'   `start`, `end`. Exon
#'   `cds_gen_start`/`cds_gen_stop` define which transcripts are coding and the
#'   coding window used for similarity.
#' @param sequences Sequence table with `transcript_id`, `transcript_seq`,
#'   `protein_seq` and optional `protein_id`. CDS or full transcript sequences
#'   are recognized by agreement with annotation lengths. A terminal stop
#'   codon may be present. Conflicting duplicate transcript records are errors.
#'   Sequences are attached to selected pairs for downstream analysis; similarity
#'   itself is computed from annotated coordinates.
#' @param source Pairing mode, as in [get_pairs()].
#' @param max_candidates Maximum candidate transcripts per form in one
#'   comparison (default 100, so at most 10,000 pairs are scored). A form with
#'   more keeps the best supported, in the ranking's own order: annotated CDS,
#'   TSL tier, exact structure, numeric TSL, then transcript ID. The number
#'   dropped is reported in `candidates_dropped_case` and
#'   `candidates_dropped_control`, and a message counts the affected comparisons.
#' @param fallback Logical. When no structurally matching pair exists, select an
#'   `approximate` pair with the legacy matcher's overlap rules (default `TRUE`).
#'   Set `FALSE` to leave such comparisons unresolved.
#' @param verbose Emit progress messages.
#'
#' @details
#' Ranking is lexicographic: number of coding transcripts (descending),
#' worst TSL tier (1--3, then 4--5, then unknown), evidence (both CDS sequences
#' usable, CDS annotated without a usable sequence, then a member without CDS),
#' context similarity, whole-window similarity (`orf_similarity`), number of
#' members that are not exact structural matches, worst numeric TSL, summed
#' TSL, and transcript IDs. TSL 1--3 is a package policy tier; within it
#' similarity precedes exact TSL. This is not a transcript-expression
#' likelihood model.
#'
#' A transcript is coding when it has an annotated CDS; its biotype is not
#' used. Transcripts annotated as nonsense-mediated decay therefore count as
#' coding: when an event shifts the frame, the transcript that carries it
#' faithfully is often annotated that way, and requiring the `protein_coding`
#' biotype would favour a frame-restoring isoform.
#'
#' The event mask is the symmetric difference between required intervals of the
#' two forms, plus their explicit forbidden intervals. Context similarity is
#' the base-pair Jaccard overlap of the two transcripts' exonic positions in a
#' coding window, outside the mask: the span coded by either transcript when
#' both are coding, otherwise one window shared by every pair of the comparison
#' (the span coded by any candidate, or the candidates' full span without CDS).
#' Positions in one genome determine sequence, so this compares the coding
#' region without depending on where each annotated reading frame ends: a
#' frameshifted isoform keeps credit for its downstream exons and is not
#' outranked by a frame-restoring isoform. `orf_similarity` is the same score
#' including the event. Protein identity is deliberately not a criterion.
#'
#' Forbidden segments cannot overlap any exon, and event-defining boundaries
#' must match exactly: SE/MXE and spliced RI check adjacent exons and junction
#' boundaries, retained RI pieces must lie within one exon, A3SS/A5SS check the
#' variable boundary and adjacent exon existence (an omitted partner interval
#' cannot be checked), and terminal events check the first/last exon.
#' A transcript that also covers every required segment is an `exact` match.
#' The outer ends of flanking exons, the shared A3SS/A5SS boundary and terminal
#' TSS/PAS ends are not part of an event, so a transcript that differs only
#' there is a `relaxed` match (`structural_match`). Relaxed matches must still
#' contain their form's alternative region and lose ties to exact matches.
#' For terminal events they are used only when no transcript matches exactly,
#' so sites sharing a splice site stay distinguishable. Partial and generic
#' interval-only definitions are labelled in diagnostics.
#'
#' With `fallback = TRUE`, a comparison without a structural pair is retried
#' with the legacy matcher's gates: each required interval overlaps an exon by
#' at least 5%, and no exon overlaps a forbidden interval by 5%. A pair must use
#' two transcripts and differ at the event in the expected direction:
#' `event_agreement` adds, for each form, how much more of that form's
#' event-specific bases its own transcript covers than its partner (0 to 2), and
#' must be positive. Approximate pairs are ranked by `event_agreement` first,
#' then by the criteria above, and are labelled `matching_tier = "approximate"`
#' with `approximate` structural matches; their event boundaries are unverified.
#'
#' Input errors (missing metadata, strands other than `+`/`-`, unknown forms,
#' malformed intervals, duplicate definitions) stop the run before any matching
#' and name the affected events.
#'
#' @return A list with `pairs` (one selected pair per positive/negative form
#'   combination, compatible with [compare_sequence_frame()]; members' biotypes
#'   are in `transcript_type_case` and `transcript_type_control`, `NA` without a
#'   `transcript_type` column, and `n_event_comparisons` counts the comparisons
#'   its event defines, resolved or not), `matched` (selected
#'   form rows, qualified by `matching_group`), `candidates` (eligibility,
#'   structural match, approximate eligibility and rejection reasons), `rankings`
#'   (all scored combinations, by tier), `events` (all input rows and candidate
#'   counts), `unmatched` (form pairs unresolved by every tier, with the
#'   structural `reason`, the `fallback` outcome and any candidates dropped by
#'   `max_candidates`),
#'   `pair_rejections` (self-pairs, combinations unable to distinguish the forms,
#'   and approximate combinations without an event difference), and `settings`.
#'   Ties and context-score margins are descriptive diagnostics, not statistical
#'   confidence or probabilities.
#'
#' @seealso [get_matched_events_chunked()], [get_pairs()], [get_splicing_impact()]
#' @examples
#' ex <- load_example_data(c("sample_frame", "annotation_df"))
#' hit_index <- get_hitindex(ex$sample_frame)
#' res <- get_differential_inclusion(hit_index)
#' res_di <- keep_sig_pairs(res)
#' ranked <- get_ranked_pairs(res_di, ex$annotation_df$annotations,
#'                            ex$annotation_df$sequences)
#' # One case/control transcript pair per resolved comparison
#' head(ranked$pairs[, c("event_id", "transcript_id_case",
#'                       "transcript_id_control", "matching_tier")])
#' table(ranked$pairs$matching_tier)
#' # Comparisons no tier resolved, and why
#' table(ranked$unmatched$reason)
#' @export
get_ranked_pairs <- function(events, annotations, sequences,
                             source = c("multi", "paired"),
                             max_candidates = 100L,
                             fallback = TRUE,
                             verbose = TRUE) {
  source <- match.arg(source)
  if (!is.numeric(max_candidates) || length(max_candidates) != 1L ||
      is.na(max_candidates) || max_candidates < 1) {
    stop("max_candidates must be a positive scalar.")
  }
  max_candidates <- as.integer(min(max_candidates, .Machine$integer.max))
  if (!is.logical(fallback) || length(fallback) != 1L || is.na(fallback)) {
    stop("fallback must be TRUE or FALSE.")
  }
  # As in legacy matching, S4 input uses significant events when available.
  input <- .resolve_splice_input(events, what = "res_di")
  if (!ncol(data.table::as.data.table(input$dt))) input <- .resolve_splice_input(events, what = "di_events")
  ev <- data.table::copy(data.table::as.data.table(input$dt))
  reserved <- intersect(c("event_row", "n_candidates", "matching_status"), names(ev))
  if (length(reserved)) ev[, (reserved) := NULL]
  required <- c("event_id", "event_type", "form", "gene_id", "chr", "strand", "inc", "exc")
  if (!all(required %in% names(ev))) stop("Missing canonical event columns: ", paste(setdiff(required, names(ev)), collapse = ", "))
  for (nm in required) data.table::set(ev, j = nm, value = as.character(ev[[nm]]))
  # Input errors stop before any matching and name the events to fix.
  for (nm in setdiff(required, c("inc", "exc"))) {
    bad <- is.na(ev[[nm]]) | !nzchar(ev[[nm]])
    if (any(bad)) {
      .si_rank_input_error(paste0("Event metadata cannot be missing: ", nm, "."),
                           if (nm == "event_id") paste("row", which(bad)) else ev$event_id[bad])
    }
  }
  bad <- !ev$strand %in% c("+", "-")
  if (any(bad)) {
    .si_rank_input_error(sprintf("Ranked matching requires explicit + or - strands (found %s).",
                                 paste(unique(ev$strand[bad]), collapse = ", ")), ev$event_id[bad])
  }
  bad <- !ev$form %in% c("INC", "EXC", "SITE")
  if (any(bad)) {
    .si_rank_input_error(sprintf("Unknown event form (found %s); use INC, EXC or SITE.",
                                 paste(unique(ev$form[bad]), collapse = ", ")), ev$event_id[bad])
  }
  ev[, event_row := .I]
  parse_all <- function(x) lapply(x, function(z) tryCatch(.si_rank_intervals(z), error = function(e) NULL))
  incs <- parse_all(ev$inc)
  excs <- parse_all(ev$exc)
  bad <- vapply(incs, is.null, logical(1)) | vapply(excs, is.null, logical(1))
  if (any(bad)) {
    .si_rank_input_error(paste("Event intervals must be one-based closed start-end spans",
                               "(start <= end) separated by semicolons."), ev$event_id[bad])
  }
  keys <- .mk_di_key(ev$event_id, ev$form, ev$inc, ev$exc)
  bad <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  if (any(bad)) {
    .si_rank_input_error("Each event/form/interval definition must identify one differential row.",
                         ev$event_id[bad])
  }
  for (nm in c("p.value", "padj")) if (!nm %in% names(ev)) ev[, (nm) := NA_real_]
  for (nm in c("n_samples", "n_control", "n_case")) if (!nm %in% names(ev)) ev[, (nm) := NA_integer_]
  # Use the established event grouping/direction contract, with row IDs as
  # temporary transcript IDs. No annotation winner is selected at this step.
  skeleton <- data.table::copy(ev)
  skeleton[, `:=`(transcript_id = as.character(event_row), exons = "", protein_id = NA_character_,
                  transcript_seq = NA_character_, protein_seq = NA_character_)]
  groups <- get_pairs(skeleton, source = source, return_class = "data.table")
  if (!"delta_psi" %in% names(ev)) ev[, delta_psi := NA_real_]
  groups[, `:=`(case_row = as.integer(transcript_id_case), control_row = as.integer(transcript_id_control))]

  ann <- data.table::as.data.table(annotations)
  needed <- c("type", "transcript_id", "gene_id", "chr", "strand", "start", "end", "exon_id")
  if (!all(needed %in% names(ann))) stop("Missing annotation columns: ", paste(setdiff(needed, names(ann)), collapse = ", "))
  ann <- data.table::copy(ann[gene_id %chin% ev$gene_id & type %chin% c("transcript", "exon")])
  if (!"transcript_support_level" %in% names(ann)) ann[, transcript_support_level := NA_character_]
  # The biotype is reported as a label only; it does not affect ranking.
  if (!"transcript_type" %in% names(ann)) ann[, transcript_type := NA_character_]
  meta_cols <- c("transcript_id", "gene_id", "chr", "strand", "transcript_type", "transcript_support_level")
  tx <- unique(ann[type == "transcript", ..meta_cols])
  if (anyNA(tx$transcript_id) || anyDuplicated(tx$transcript_id)) stop("Conflicting or missing transcript metadata.")
  tx[, tsl_rank := .si_tsl_rank(transcript_support_level)]
  tx[, tsl_tier := data.table::fifelse(tsl_rank <= 3L, 0L,
                    data.table::fifelse(tsl_rank <= 5L, 1L, 2L))]
  data.table::setorder(tx, transcript_id)
  data.table::setindexv(tx, c("gene_id", "chr", "strand"))
  exon_cols <- intersect(c("transcript_id", "exon_id", "gene_id", "chr", "strand",
                           "start", "end", "cds_gen_start", "cds_gen_stop"), names(ann))
  ex <- unique(ann[type == "exon", ..exon_cols])
  if (nrow(ex) && any(is.na(ex$start) | is.na(ex$end) | ex$start < 1 | ex$end < ex$start)) {
    stop("Invalid annotation exon bounds.")
  }
  data.table::setorder(ex, transcript_id, start, end, exon_id)
  ex_by_tx <- split(ex, ex$transcript_id, drop = TRUE)
  for (id in names(ex_by_tx)) {
    e <- ex_by_tx[[id]]
    m <- tx[list(id), on = "transcript_id"]
    if (!nrow(m) || is.na(m$transcript_id) || any(e$gene_id != m$gene_id | e$chr != m$chr | e$strand != m$strand)) {
      stop("Exon/transcript locus metadata disagree for ", id)
    }
    if (anyDuplicated(e$exon_id) || (nrow(e) > 1L && any(e$start[-1L] <= utils::head(e$end, -1L)))) {
      stop("Conflicting or overlapping annotated exons for ", id)
    }
  }
  sq <- data.table::as.data.table(sequences)
  if (!all(c("transcript_id", "transcript_seq", "protein_seq") %in% names(sq))) {
    stop("Sequence table requires transcript_id, transcript_seq and protein_seq.")
  }
  sq_cols <- c("transcript_id", "protein_id", "transcript_seq", "protein_seq")
  present <- intersect(sq_cols, names(sq))
  sq <- unique(sq[transcript_id %chin% tx$transcript_id, ..present])
  if (!"protein_id" %in% names(sq)) sq[, protein_id := NA_character_]
  if (anyDuplicated(sq$transcript_id)) stop("Conflicting sequence records for a transcript ID.")
  sq <- sq[match(tx$transcript_id, sq$transcript_id)]
  sq[, transcript_id := tx$transcript_id]
  cds <- stats::setNames(lapply(seq_len(nrow(tx)), function(i) {
    e <- ex_by_tx[[tx$transcript_id[i]]]
    if (is.null(e)) e <- ex[0]
    .si_rank_cds(e, sq$transcript_seq[i], tx$strand[i])
  }), tx$transcript_id)
  # Coding means an annotated CDS, whatever the biotype: NMD-annotated
  # transcripts often carry an event's frameshift faithfully.
  tx[, `:=`(coding = vapply(cds, function(z) length(z$ranges) > 0L, logical(1), USE.NAMES = FALSE))]

  exons_of <- function(id) {
    exons <- ex_by_tx[[id]]
    if (is.null(exons)) ex[0] else exons
  }
  ranges_of <- function(id) {
    exons <- exons_of(id)
    IRanges::IRanges(exons$start, exons$end)
  }
  # Outer flank extents never define an event, so non-terminal forms always
  # admit relaxed matches. Terminal forms use them only without an exact
  # match, so sites that share a splice site stay distinguishable.
  terminal_types <- c("AFE", "HFE", "ALE", "HLE")
  relaxed_ok <- logical(nrow(ev))
  fits <- function(id, row, allow_relaxed = relaxed_ok[row]) {
    args <- list(exons_of(id), incs[[row]], excs[[row]], ev$event_type[row], ev$form[row], ev$strand[row])
    check <- do.call(.si_rank_compatible, args)
    if (!check$eligible && allow_relaxed) check <- do.call(.si_rank_compatible, c(args, relaxed = TRUE))
    check
  }
  candidate_parts <- vector("list", nrow(ev))
  for (i in seq_len(nrow(ev))) {
    e <- ev[i]
    eligible_tx <- tx[list(e$gene_id, e$chr, e$strand),
                      on = c("gene_id", "chr", "strand"), which = TRUE, nomatch = 0L]
    checks <- lapply(tx$transcript_id[eligible_tx], fits, row = i, allow_relaxed = FALSE)
    relaxed_ok[i] <- !e$event_type %in% terminal_types ||
      !any(vapply(checks, function(check) check$eligible, logical(1)))
    if (relaxed_ok[i]) {
      checks <- lapply(seq_along(eligible_tx), function(k) {
        if (checks[[k]]$eligible) return(checks[[k]])
        .si_rank_compatible(exons_of(tx$transcript_id[eligible_tx[k]]), incs[[i]], excs[[i]],
                            e$event_type, e$form, e$strand, relaxed = TRUE)
      })
    }
    candidate_parts[[i]] <- data.table::rbindlist(lapply(seq_along(eligible_tx), function(k) {
      j <- eligible_tx[k]
      check <- checks[[k]]
      data.table::data.table(event_row = i, event_id = e$event_id, form = e$form,
        transcript_id = tx$transcript_id[j], eligible = check$eligible,
        reason = check$reason, structural_context = check$context,
        structural_match = check$match, exons = check$exons,
        approximate_ok = .si_rank_approximate(ranges_of(tx$transcript_id[j]), incs[[i]], excs[[i]]),
        coding = tx$coding[j], transcript_type = tx$transcript_type[j],
        tsl_rank = tx$tsl_rank[j], tsl_tier = tx$tsl_tier[j],
        sequence_status = cds[[tx$transcript_id[j]]]$status,
        protein_available = !is.na(sq$protein_seq[j]) && nzchar(sq$protein_seq[j]))
    }), fill = TRUE)
  }
  candidates <- data.table::rbindlist(candidate_parts, fill = TRUE)
  if (!ncol(candidates)) candidates <- data.table::data.table(event_row = integer(), event_id = character(),
    form = character(), transcript_id = character(), eligible = logical(), reason = character(),
    structural_context = character(), structural_match = character(), exons = character(),
    approximate_ok = logical(), coding = logical(), transcript_type = character(),
    tsl_rank = integer(), tsl_tier = integer(),
    sequence_status = character(), protein_available = logical())
  good <- candidates[eligible == TRUE]
  data.table::setkey(good, event_row)
  approx <- candidates[approximate_ok == TRUE]
  approx[eligible == FALSE, structural_match := "approximate"]
  data.table::setkey(approx, event_row)
  counts <- good[, .(n_candidates = .N), by = event_row]
  ev <- counts[ev, on = "event_row"]
  ev[is.na(n_candidates), n_candidates := 0L]
  ev[, matching_status := ifelse(n_candidates > 0L, "compatible_candidates", "no_compatible_transcript")]
  rank_parts <- selected <- selected_forms <- unresolved <- rejected <- list()
  # Exact outer boundaries only break ties left by coding-region similarity.
  rank_cols <- c("coding_count", "tsl_tier", "evidence_rank", "context_similarity",
                 "orf_similarity", "n_inexact", "tsl_worst", "tsl_sum",
                 "transcript_id_case", "transcript_id_control")
  rank_order <- c(-1L, 1L, 1L, -1L, -1L, 1L, 1L, 1L, 1L, 1L)
  evidence_labels <- c("cds_sequences", "cds_annotation", "exon_annotation", "no_outside_event_context")
  exon_ids_for <- function(id, row) {
    x <- exons_of(id)
    h <- IRanges::findOverlaps(incs[[row]], IRanges::IRanges(x$start, x$end))
    paste(unique(x$exon_id[S4Vectors::subjectHits(h)[order(S4Vectors::queryHits(h))]]), collapse = ";")
  }
  # Score one tier's candidate combinations for a form comparison. Returns the
  # ordered ranking, or NULL with the reason no combination survived.
  score_tier <- function(ca, cb, row_a, row_b, label, event_id, tier) {
    mask <- IRanges::union(IRanges::setdiff(incs[[row_a]], incs[[row_b]]),
                           IRanges::setdiff(incs[[row_b]], incs[[row_a]]))
    mask <- IRanges::union(mask, IRanges::union(excs[[row_a]], excs[[row_b]]))
    # Bases that only one form requires define each form's side of the event.
    own_a <- IRanges::setdiff(incs[[row_a]], incs[[row_b]])
    own_b <- IRanges::setdiff(incs[[row_b]], incs[[row_a]])
    ids <- unique(c(ca$transcript_id, cb$transcript_id))
    rng <- stats::setNames(lapply(ids, ranges_of), ids)
    cov_a <- vapply(rng, .si_rank_coverage, numeric(1), region = own_a)
    cov_b <- vapply(rng, .si_rank_coverage, numeric(1), region = own_b)
    if (tier == "approximate" && !ev$event_type[row_a] %in% terminal_types) {
      # "relaxed" promises the form's alternative region; without it, a member
      # chosen by the fallback is only an approximation.
      if (length(own_a)) ca[structural_match == "relaxed" & cov_a[transcript_id] < 1, structural_match := "approximate"]
      if (length(own_b)) cb[structural_match == "relaxed" & cov_b[transcript_id] < 1, structural_match := "approximate"]
    }
    # A pair with a non-coding member is scored on one window shared by the
    # comparison, so its pairs stay comparable.
    spans <- Filter(length, lapply(ids, function(id) .si_rank_span(cds[[id]]$ranges)))
    event_window <- if (length(spans)) IRanges::reduce(do.call(c, spans)) else
      IRanges::reduce(do.call(c, unname(lapply(rng, .si_rank_span))))
    window_type <- if (length(spans)) "event_cds_span" else "event_exon_span"

    combinations <- expand.grid(a = seq_len(nrow(ca)), b = seq_len(nrow(cb)))
    ta_all <- ca$transcript_id[combinations$a]; tb_all <- cb$transcript_id[combinations$b]
    agreement <- (if (length(own_a)) cov_a[ta_all] - cov_a[tb_all] else 0) +
      (if (length(own_b)) cov_b[tb_all] - cov_b[ta_all] else 0)
    reason <- ifelse(ta_all == tb_all, "same_transcript", "")
    if (tier == "structural") {
      # Cross-compatibility uses the partner form's own exact/relaxed rule.
      cross_a <- vapply(ca$transcript_id, function(id) fits(id, row_b)$eligible, logical(1))
      cross_b <- vapply(cb$transcript_id, function(id) fits(id, row_a)$eligible, logical(1))
      reason[!nzchar(reason) & cross_a[combinations$a] & cross_b[combinations$b]] <- "both_forms_compatible"
      if (!ev$event_type[row_a] %in% terminal_types) {
        # Relaxed outer boundaries never excuse a partial alternative region:
        # each transcript must contain the bases specific to its own form.
        partial <- (length(own_a) > 0L & cov_a[ta_all] < 1) | (length(own_b) > 0L & cov_b[tb_all] < 1)
        reason[!nzchar(reason) & partial] <- "event_region_not_covered"
      }
    } else {
      reason[!nzchar(reason) & !(agreement > 1e-12)] <- "no_event_difference"
    }
    bad <- which(nzchar(reason))
    if (length(bad)) rejected[[length(rejected) + 1L]] <<- data.table::data.table(
      matching_group = label, event_id = event_id, matching_tier = tier,
      transcript_id_case = ta_all[bad], transcript_id_control = tb_all[bad], reason = reason[bad])
    keep <- which(!nzchar(reason))
    if (!length(keep)) {
      why <- if (any(reason == "both_forms_compatible")) "no_form_distinguishing_pair" else
        if (any(reason == "event_region_not_covered")) "event_region_not_covered" else
          if (any(reason == "no_event_difference")) "no_event_difference" else "no_distinct_transcript_pair"
      return(list(ranking = NULL, reason = why))
    }
    ranking <- data.table::rbindlist(lapply(keep, function(k) {
      a <- ca[combinations$a[k]]; b <- cb[combinations$b[k]]
      ta <- a$transcript_id; tb <- b$transcript_id
      ra <- cds[[ta]]$ranges; rb <- cds[[tb]]$ranges
      both_coding <- length(ra) > 0L && length(rb) > 0L
      window <- if (both_coding) IRanges::reduce(c(.si_rank_span(ra), .si_rank_span(rb))) else event_window
      ctx <- .si_rank_positions(rng[[ta]], rng[[tb]], window, mask)
      whole <- .si_rank_positions(rng[[ta]], rng[[tb]], window)
      evidence_rank <- if (!both_coding) 2L else
        if (!is.na(cds[[ta]]$sequence) && !is.na(cds[[tb]]$sequence)) 0L else 1L
      if (is.na(ctx$similarity)) evidence_rank <- 3L
      data.table::data.table(matching_group = label, event_id = event_id,
        case_row = row_a, control_row = row_b, transcript_id_case = ta, transcript_id_control = tb,
        matching_tier = tier, event_agreement = agreement[k],
        coding_count = as.integer(a$coding) + as.integer(b$coding),
        tsl_tier = max(a$tsl_tier, b$tsl_tier), evidence_rank = evidence_rank,
        similarity_source = evidence_labels[evidence_rank + 1L],
        context_window = if (both_coding) "pair_cds_span" else window_type,
        context_similarity = ctx$similarity, orf_similarity = whole$similarity,
        context_coverage_case = ctx$coverage_case, context_coverage_control = ctx$coverage_control,
        context_compared_nt = ctx$compared_nt,
        n_inexact = as.integer(a$structural_match != "exact") + as.integer(b$structural_match != "exact"),
        tsl_worst = max(a$tsl_rank, b$tsl_rank), tsl_sum = a$tsl_rank + b$tsl_rank)
    }), fill = TRUE)
    # Approximate pairs are ordered by how faithfully they differ at the event.
    cols <- if (tier == "approximate") c("event_agreement", rank_cols) else rank_cols
    data.table::setorderv(ranking, cols, if (tier == "approximate") c(-1L, rank_order) else rank_order,
                          na.last = TRUE)
    ranking[, pair_rank := seq_len(.N)]
    tie_cols <- utils::head(cols, -2L)
    best <- ranking[1L]
    tied <- vapply(seq_len(nrow(ranking)), function(i) {
      all(vapply(tie_cols, function(nm) isTRUE(all.equal(ranking[[nm]][i], best[[nm]][1L], tolerance = 1e-12)), logical(1)))
    }, logical(1))
    ranking[, `:=`(selected = pair_rank == 1L, tied_best = tied)]
    list(ranking = ranking, reason = NA_character_)
  }

  # A form with more than `max_candidates` candidates keeps the best supported,
  # in the pair ranking's own order: annotated CDS, TSL tier, exact structure,
  # numeric TSL, then transcript ID. The number dropped is reported.
  cap_candidates <- function(cand) {
    n <- nrow(cand)
    if (n <= max_candidates) return(list(kept = cand, dropped = 0L))
    ord <- order(!cand$coding, cand$tsl_tier, cand$structural_match != "exact", cand$tsl_rank,
                 cand$transcript_id, method = "radix")
    list(kept = cand[sort(ord[seq_len(max_candidates)])], dropped = n - max_candidates)
  }

  for (g in seq_len(nrow(groups))) {
    row_a <- groups$case_row[g]; row_b <- groups$control_row[g]
    event_id <- groups$event_id[g]
    label <- .mk_pair_key(event_id, ev$inc[row_a], ev$inc[row_b], ev$exc[row_a], ev$exc[row_b])
    if (verbose && (g == 1L || g %% 25L == 0L)) message("[RANK] Form pair ", g, "/", nrow(groups))
    cap_a <- cap_candidates(good[list(row_a), nomatch = 0L])
    cap_b <- cap_candidates(good[list(row_b), nomatch = 0L])
    ca <- cap_a$kept; cb <- cap_b$kept
    dropped <- dropped_any <- c(cap_a$dropped, cap_b$dropped)
    result <- if (nrow(ca) && nrow(cb)) score_tier(ca, cb, row_a, row_b, label, event_id, "structural") else
      list(ranking = NULL, reason = "missing_form_candidate")
    tier <- "structural"
    fallback_status <- "not_run"
    if (is.null(result$ranking) && fallback) {
      cap_fa <- cap_candidates(approx[list(row_a), nomatch = 0L])
      cap_fb <- cap_candidates(approx[list(row_b), nomatch = 0L])
      fa <- cap_fa$kept; fb <- cap_fb$kept
      dropped_any <- pmax(dropped_any, c(cap_fa$dropped, cap_fb$dropped))
      if (nrow(fa) && nrow(fb)) {
        approximate <- score_tier(fa, fb, row_a, row_b, label, event_id, "approximate")
        fallback_status <- if (is.null(approximate$ranking)) approximate$reason else "selected"
        if (!is.null(approximate$ranking)) {
          result <- approximate
          ca <- fa; cb <- fb; tier <- "approximate"
          dropped <- c(cap_fa$dropped, cap_fb$dropped)
        }
      } else fallback_status <- "missing_form_candidate"
    }
    if (is.null(result$ranking)) {
      unresolved[[length(unresolved) + 1L]] <- data.table::data.table(matching_group = label,
        event_id = event_id, case_row = row_a, control_row = row_b, reason = result$reason,
        fallback = fallback_status, candidates_dropped_case = dropped_any[1L],
        candidates_dropped_control = dropped_any[2L])
      next
    }
    ranking <- result$ranking
    rank_parts[[length(rank_parts) + 1L]] <- ranking
    best <- ranking[1L]
    form_rows <- data.table::copy(ev[c(row_a, row_b)])
    choices <- c(best$transcript_id_case, best$transcript_id_control)
    form_rows[, `:=`(transcript_id = choices,
                     transcript_type = tx$transcript_type[match(choices, tx$transcript_id)])]
    chosen_exons <- c(ca[match(choices[1L], transcript_id), exons], cb[match(choices[2L], transcript_id), exons])
    # Approximate members have no verified structure; list their overlapping exons.
    for (k in which(is.na(chosen_exons) | !nzchar(chosen_exons))) {
      chosen_exons[k] <- exon_ids_for(choices[k], c(row_a, row_b)[k])
    }
    form_rows[, exons := chosen_exons]
    form_rows[, delta_psi := c(1, -1)]
    # get_pairs handles the established output schema; restore measured PSI
    # afterwards, including NA for a form-based call without DI statistics.
    with_seq <- attach_sequences(form_rows, sq, return_class = "data.table")
    pair <- get_pairs(with_seq, source = "multi", return_class = "data.table")
    pair[, `:=`(delta_psi_case = ev$delta_psi[row_a], delta_psi_control = ev$delta_psi[row_b],
                n_event_comparisons = groups$n_event_comparisons[g],
                matching_group = label, matching_method = "orf", matching_tier = tier,
                n_candidate_pairs = nrow(ranking), n_tied_best = sum(ranking$tied_best),
                candidates_dropped_case = dropped[1L], candidates_dropped_control = dropped[2L],
                event_agreement = best$event_agreement,
                similarity_source = best$similarity_source, context_window = best$context_window,
                context_similarity = best$context_similarity, orf_similarity = best$orf_similarity,
                coding_count = best$coding_count, tsl_tier = best$tsl_tier, tsl_worst = best$tsl_worst,
                structural_context_case = ca[match(best$transcript_id_case, transcript_id), structural_context],
                structural_context_control = cb[match(best$transcript_id_control, transcript_id), structural_context],
                structural_match_case = ca[match(best$transcript_id_case, transcript_id), structural_match],
                structural_match_control = cb[match(best$transcript_id_control, transcript_id), structural_match])]
    comparable <- ranking[-1L][coding_count == best$coding_count & tsl_tier == best$tsl_tier &
                                evidence_rank == best$evidence_rank]
    margin <- if (!nrow(comparable) || is.na(best$context_similarity)) NA_real_ else
      best$context_similarity - max(comparable$context_similarity, na.rm = TRUE)
    pair[, context_score_margin := if (is.finite(margin)) margin else NA_real_]
    pair[, matching_pair_id := paste(label, transcript_id_case, transcript_id_control, sep = "|")]
    with_seq[, `:=`(matching_group = label, delta_psi = ev$delta_psi[c(row_a, row_b)])]
    selected[[length(selected) + 1L]] <- pair
    selected_forms[[length(selected_forms) + 1L]] <- with_seq
  }
  pairs <- data.table::rbindlist(selected, fill = TRUE)
  if (!ncol(pairs)) {
    empty <- data.table::copy(skeleton[0])[, transcript_type := character()]
    if (!"delta_psi" %in% names(empty)) empty[, delta_psi := numeric()]
    pairs <- get_pairs(empty, source = "multi", return_class = "data.table")
    pairs[, `:=`(matching_group = character(), matching_method = character(), matching_tier = character(),
                  n_candidate_pairs = integer(), n_tied_best = integer(),
                  candidates_dropped_case = integer(), candidates_dropped_control = integer(),
                  event_agreement = numeric(),
                  similarity_source = character(), context_window = character(),
                  context_similarity = numeric(), orf_similarity = numeric(), coding_count = integer(),
                  tsl_tier = integer(), tsl_worst = integer(),
                  structural_context_case = character(), structural_context_control = character(),
                  structural_match_case = character(), structural_match_control = character(),
                  context_score_margin = numeric(), matching_pair_id = character())]
  }
  unmatched_rows <- setdiff(ev$event_row, unique(c(groups$case_row, groups$control_row)))
  if (length(unmatched_rows)) {
    ev[event_row %in% unmatched_rows, matching_status := "no_opposite_direction_form"]
  }
  matched <- data.table::rbindlist(selected_forms, fill = TRUE)
  if (!ncol(matched)) matched <- data.table::copy(skeleton[0])[, `:=`(transcript_type = character(),
                                                                       matching_group = character())]
  unmatched <- data.table::rbindlist(unresolved, fill = TRUE)
  if (!ncol(unmatched)) unmatched <- data.table::data.table(matching_group = character(),
    event_id = character(), case_row = integer(), control_row = integer(), reason = character(),
    fallback = character(), candidates_dropped_case = integer(), candidates_dropped_control = integer())
  n_capped <- sum(pairs$candidates_dropped_case > 0L | pairs$candidates_dropped_control > 0L) +
    sum(unmatched$candidates_dropped_case > 0L | unmatched$candidates_dropped_control > 0L)
  if (n_capped) {
    message(sprintf(paste("get_ranked_pairs: %d comparison(s) had a form with more than max_candidates = %d",
                          "candidates; the best supported were kept (see candidates_dropped_case/_control)."),
                    n_capped, max_candidates))
  }
  rejections <- data.table::rbindlist(rejected, fill = TRUE)
  if (!ncol(rejections)) rejections <- data.table::data.table(matching_group = character(),
    event_id = character(), matching_tier = character(), transcript_id_case = character(),
    transcript_id_control = character(), reason = character())
  rankings <- data.table::rbindlist(rank_parts, fill = TRUE)
  if (!ncol(rankings)) rankings <- data.table::data.table(matching_group = character(),
    event_id = character(), case_row = integer(), control_row = integer(),
    transcript_id_case = character(), transcript_id_control = character(),
    matching_tier = character(), event_agreement = numeric(), coding_count = integer(),
    tsl_tier = integer(), evidence_rank = integer(), similarity_source = character(),
    context_window = character(), context_similarity = numeric(), orf_similarity = numeric(),
    context_coverage_case = numeric(), context_coverage_control = numeric(),
    context_compared_nt = integer(), n_inexact = integer(), tsl_worst = integer(),
    tsl_sum = integer(), pair_rank = integer(), selected = logical(), tied_best = logical())
  list(pairs = pairs, matched = matched, candidates = candidates, rankings = rankings,
       pair_rejections = rejections, events = ev, unmatched = unmatched,
       settings = list(method = "orf", protocol_version = 5L, source = source,
                       coding = "annotated CDS; biotype not used",
                       outer_boundaries = "relaxed; terminal TSS/PAS only without an exact match",
                       similarity = "exonic base-pair Jaccard in the coding window, outside the event mask",
                       fallback = fallback, fallback_min_overlap = 0.05,
                       max_candidates = max_candidates,
                       candidate_cut = "annotated CDS, TSL tier, exact structure, TSL, transcript ID",
                       tsl_tiers = list(high = seq_len(3L), lower = 4:5, unknown = 6L)))
}
