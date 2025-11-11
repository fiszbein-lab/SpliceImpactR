#' Build exon and transcript tables from a parsed annotation
#'
#' Constructs standardized exon and transcript data tables from a
#' flattened annotation object such as one produced by
#' [rtracklayer::import()] on a GTF file.
#'
#' @param ann A `data.frame` or `data.table` containing parsed gene
#'   annotation records with at least the following columns:
#'   \code{type}, \code{chr}, \code{start}, \code{end},
#'   \code{strand}, \code{gene_id}, \code{gene_name},
#'   \code{transcript_id}, \code{transcript_name},
#'   \code{transcript_type}, \code{protein_id},
#'   \code{exon_id}, and \code{exon_number}.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{`exons`}{A `data.table` of exon records with added
#'     `classification` (`"first"`, `"last"`, or `"internal"`)
#'     determined per transcript.}
#'   \item{`transcripts`}{A `data.table` of transcript records.}
#'   \item{`protein_tx`}{Character vector of transcript IDs
#'     annotated as protein-coding (non-NA `protein_id`).}
#' }
#'
#' @details
#' This function is used internally to prepare lightweight annotation
#' summaries for downstream event mapping and HIT-index aggregation.
#' It does not query external resources and assumes the input
#' annotation already includes both exon and transcript entries.
#'
#' @importFrom data.table as.data.table setkey fifelse
#' @keywords internal
build_from_annotations <- function(ann) {
  ann <- data.table::as.data.table(ann)

  tx <- ann[type == "transcript",
            .(rownum = row_uid,
              chr, start, stop = end, strand,
              gene_id, gene_name,
              transcript_id, transcript_name, transcript_type,
              protein_id)]
  data.table::setkey(tx, transcript_id)

  ex <- ann[type == "exon",
            .(rownum = row_uid,
              chr, start, stop = end, strand,
              gene_id, gene_name,
              transcript_id, exon_id, exon_number)]

  mx <- ex[!is.na(exon_number), .(max_ex = max(exon_number)), by = transcript_id]
  ex <- mx[ex, on = "transcript_id"]
  ex[, classification :=
       data.table::fifelse(!is.na(exon_number) & exon_number == 1L, "first",
               data.table::fifelse(!is.na(exon_number) & exon_number == max_ex, "last", "internal"))]
  ex[is.na(exon_number) | is.na(max_ex), classification := NA_character_]
  ex[, max_ex := NULL]

  protein_tx <- unique(na.omit(tx$transcript_id[!is.na(tx$protein_id)]))
  list(exons = ex, transcripts = tx, protein_tx = protein_tx)
}

#' Expand semicolon-delimited coordinate strings into long-format tables
#'
#' Parses inclusion or exclusion coordinate strings of the form
#' `"start-end;start-end;..."` into per-exon start/stop rows.
#'
#' @param ev A `data.frame` or `data.table` containing at least
#'   `inc`, `exc`, `event_type`, `gene_id`, `chr`, `strand`, and `form` columns.
#' @param which Character. One of `"inc"` or `"exc"` specifying which coordinate
#'   column to explode (default `"inc"`).
#'
#' @return A `data.table` with one row per coordinate interval containing:
#'   `event_row`, `event_type`, `gene_id`, `chr`, `strand`, `form`, `event_id`,
#'   `start`, `stop`, and optionally `inc_idx` (ordinal index for inclusion parts).
#'
#' @details
#' This helper is used internally to vectorize exon coordinate parsing for
#' downstream construction of `GRanges` or other genomic features.
#'
#' @importFrom data.table as.data.table rbindlist %chin%
#' @keywords internal
explode_coords <- function(ev, which = c("inc","exc")) {
  which <- match.arg(which)
  ev <- data.table::as.data.table(ev)
  col  <- which
  # make sure ev has an integer row id
  if (!(".__row__" %chin% names(ev))) ev[, .__row__ := .I]
  # split strings; empty -> no rows
  split_vec <- strsplit(
    ifelse(is.na(ev[[col]]) | !nzchar(ev[[col]]), "", ev[[col]]),
    ";", fixed = TRUE)

  # build long rows
  out_list <- vector("list", nrow(ev))
  for (i in seq_len(nrow(ev))) {
    parts <- split_vec[[i]]
    if (length(parts) == 0L || (length(parts)==1L && parts[1]== "")) next
    xy <- matrix(as.integer(unlist(strsplit(parts, "-", fixed = TRUE))), ncol = 2, byrow = TRUE)
    out_list[[i]] <- data.table(
      event_row = ev$.__row__[i],
      event_type= ev$event_type[i],
      gene_id   = ev$gene_id[i] %||% ev$geneR[i] %||% ev$gene[i],
      chr       = ev$chr[i],
      strand    = ev$strand[i],
      form      = ev$form[i],
      event_id = ev$event_id[i],
      start     = xy[,1],
      stop      = xy[,2],
      inc_idx   = if (which=="inc") seq_len(nrow(xy)) else NA_integer_
    )
  }
  if (length(Filter(NROW, out_list)) == 0L) {
    out <- data.table(
      event_row = integer(),
      event_type=character(),
      gene_id=character(),
      chr=character(),
      strand=character(),
      form=character(),
      event_id=character(),
      start=integer(),
      stop=integer(),
      inc_idx=integer(),
      key=character()
    )
  } else {
    out <- data.table::rbindlist(out_list, fill=TRUE)
    out[, key := paste(chr, strand, gene_id, sep="|")]
  }
  out[]
}


#' Match alternative splicing events to annotated exons and transcripts
#'
#' Performs vectorized overlap mapping between event coordinates (e.g. inclusion
#' or exclusion intervals) and exons from an annotation resource. Returns the
#' best-matching transcript and exon set per event based on coverage, exon
#' classification, and protein-coding preference.
#'
#' @param events A `data.frame` or `data.table` containing splicing events with
#'   columns `chr`, `strand`, `gene_id`, and coordinate fields `inc` and `exc`
#'   (semicolon-delimited `"start-end"` strings).
#' @param annotations A gene annotation table (e.g. from
#'   [rtracklayer::import()] on a GTF file) with exon and transcript rows, passed
#'   to [build_from_annotations()].
#' @param minOverlap Minimum fractional overlap (0-1) required between an
#'   inclusion segment and an annotated exon to count as a hit. Default `0.05`.
#'
#' @return A `data.table` containing one row per input event with the following
#'   appended columns:
#'   \describe{
#'     \item{`transcript_id`}{The best matching transcript ID.}
#'     \item{`exons`}{Semicolon delimited list of exon IDs covered by the event.}
#'     \item{`inc_exons_by_idx`}{Semicolon delimited exon IDs per inclusion
#'       index (maintaining order).}
#'     \item{`inc_rows_by_idx`}{Semicolon delimited exon row indices matching
#'       inclusion order.}
#'   }
#'
#' @details
#' The function uses [GenomicRanges::findOverlaps()] to match event inclusion
#' intervals to exons. Candidate transcripts are filtered to ensure sufficient
#' coverage and absence of overlaps with exclusion coordinates.
#'
#' The algorithm prioritizes exon classification (`first`, `internal`, `last`)
#' consistent with the event type (AFE, ALE, SE, etc.), followed by reciprocal
#' overlap fractions, intersection width, and protein-coding status.
#'
#' @importFrom data.table as.data.table setkey setorder rbindlist uniqueN
#' @importFrom data.table setnames %chin%
#' @importFrom GenomicRanges GRanges findOverlaps pintersect
#' @importFrom IRanges IRanges width
#' @importFrom S4Vectors queryHits subjectHits
#' @keywords internal
match_events_to_annotations_vec <- function(events,
                                            annotations,
                                            minOverlap = 0.05) {

  eps <- 1e-9  # numerical guard
  want_cols <- c("event_id","event_type","form", "gene_id","chr","strand",
                 "inc","exc","delta_psi","p.value","padj",
                 "n_samples","n_control","n_case")

  AA <- build_from_annotations(annotations)
  EX <- AA$exons
  TX <- AA$transcripts
  protein_tx <- AA$protein_tx
  # 1) events id + meta
  ev <- data.table::as.data.table(events)
  if (!(".__row__" %chin% names(ev))) ev[, .__row__ := .I]
  meta_cols <- intersect(want_cols, names(ev))
  ev_meta   <- ev[, c(".__row__", meta_cols), with = FALSE]
  setnames(ev_meta, ".__row__", "event_row")
  setkey(ev_meta, event_row)

  # 2) long coords
  inc_long <- explode_coords(ev, "inc")
  exc_long <- explode_coords(ev, "exc")

  if (!nrow(inc_long) && !nrow(exc_long)) {
    return(ev_meta[, `:=`(
      transcript_id = NA_character_,
      exons = "",
      inc_exons_by_idx = "",
      inc_rows_by_idx = ""
    )][order(event_row)])
  }

  if (!nrow(inc_long)) {
    inc_long <- data.table(
      event_row = ev$.__row__,
      event_type=ev$event_type,
      gene_id=ev$gene_id,
      chr=ev$chr,
      strand=ev$strand,
      form=ev$form,
      event_id=ev$event_id,
      start=integer(),
      stop=integer(),
      inc_idx=integer(),
      key=paste(ev$chr, ev$strand, ev$gene_id, sep="|"),
      inc_coord=character()
    )
  } else if (!("inc_coord" %in% names(inc_long))) {
    inc_long[, inc_coord := sprintf("%d%d", start, stop)]
  }

  # n_inc per event
  ninc_tbl <- inc_long[, .N, by = event_row]
  data.table::setnames(ninc_tbl, "N", "n_inc")
  ninc_map <- setNames(ninc_tbl$n_inc, ninc_tbl$event_row)

  # 3) prefilter by (chr,strand,gene)
  valid_keys <- unique(paste(EX$chr, EX$strand, EX$gene_id, sep="|"))
  inc_long   <- inc_long[key %chin% valid_keys]
  if (nrow(exc_long)) exc_long <- exc_long[key %chin% valid_keys]
  if (!nrow(inc_long)) return(ev_meta[, `:=`
                                      (transcript_id = NA_character_,
                                        exons = "",
                                        inc_exons_by_idx = "",
                                        inc_rows_by_idx = "")][order(event_row)])

  # 4) restrict EX; build GRanges once
  keys_needed <- unique(inc_long$key)
  EX_sub <- EX[paste(chr, strand, gene_id, sep="|") %chin% keys_needed]
  gr_inc <- GenomicRanges::GRanges(inc_long$chr,
                                   IRanges::IRanges(inc_long$start, inc_long$stop),
                                   inc_long$strand)
  gr_ex  <- GenomicRanges::GRanges(EX_sub$chr,
                                   IRanges::IRanges(EX_sub$start,   EX_sub$stop),
                                   EX_sub$strand)


  # 5) global INC overlaps (gate on INC side)
  H  <- GenomicRanges::findOverlaps(gr_inc, gr_ex, ignore.strand = FALSE)
  if (!length(H)) return(ev_meta[, `:=`
                                 (transcript_id = NA_character_,
                                   exons = "",
                                   inc_exons_by_idx = "",
                                   inc_rows_by_idx = "")][order(event_row)])

  q  <- gr_inc[queryHits(H)]
  s <- gr_ex[subjectHits(H)]
  wI <- IRanges::width(GenomicRanges::pintersect(q, s))
  rq <- as.numeric(wI / IRanges::width(q))
  rs <- as.numeric(wI / IRanges::width(s))
  keep <- which(rq + eps >= minOverlap)
  if (!length(keep)) return(ev_meta[, `:=`(transcript_id = NA_character_, exons = "", inc_exons_by_idx = "", inc_rows_by_idx = "")][order(event_row)])

  # 6) candidates (inc_hits)
  inc_hits <- data.table(
    event_row      = inc_long$event_row[queryHits(H)[keep]],
    event_type     = inc_long$event_type[queryHits(H)[keep]],
    inc_idx        = inc_long$inc_idx[queryHits(H)[keep]],
    gene_id        = inc_long$gene_id[queryHits(H)[keep]],
    chr            = inc_long$chr[queryHits(H)[keep]],
    strand         = inc_long$strand[queryHits(H)[keep]],
    form           = inc_long$form[queryHits(H)[keep]],
    exon_row       = EX_sub$rownum[subjectHits(H)[keep]],
    exon_id        = EX_sub$exon_id[subjectHits(H)[keep]],
    classification = EX_sub$classification[subjectHits(H)[keep]],
    transcript_id  = EX_sub$transcript_id[subjectHits(H)[keep]],
    recip_q        = rq[keep],
    recip_s        = rs[keep],
    int_w          = wI[keep]
  )
  inc_hits[, protein_link := transcript_id %chin% protein_tx]



  # 7) coverage gate
  cov_by_tx <- inc_hits[, .(covered = uniqueN(inc_idx)), by = .(event_row, transcript_id)]
  setkey(cov_by_tx, event_row); setkey(ninc_tbl, event_row)
  cov_by_tx <- cov_by_tx[ninc_tbl]
  good_cov  <- cov_by_tx[covered == n_inc]
  if (!nrow(good_cov)) return(ev_meta[, `:=`
                                      (transcript_id = NA_character_,
                                        exons = "",
                                        inc_exons_by_idx = "",
                                        inc_rows_by_idx = "")][order(event_row)])

  data.table::setkey(good_cov, event_row, transcript_id)
  data.table::setkey(inc_hits, event_row, transcript_id)
  inc_hits <- inc_hits[good_cov, nomatch = 0L]

  # 8) EXC veto (EXC side)
  if (nrow(exc_long)) {
    gr_exc <- GenomicRanges::GRanges(exc_long$chr,
                                     IRanges::IRanges(exc_long$start, exc_long$stop),
                                     exc_long$strand)
    HX <- GenomicRanges::findOverlaps(gr_exc, gr_ex, ignore.strand = FALSE)
    if (length(HX)) {
      qx  <- gr_exc[queryHits(HX)]
      sx  <- gr_ex [subjectHits(HX)]
      wIx <- GenomicRanges::width(GenomicRanges::pintersect(qx, sx))
      rqx <- as.numeric(wIx / GenomicRanges::width(qx))
      keep_bad <- which(rqx + eps >= minOverlap)
      if (length(keep_bad)) {
        bad_pairs <- data.table(
          event_row     = exc_long$event_row[queryHits(HX)[keep_bad]],
          transcript_id = EX_sub$transcript_id[subjectHits(HX)[keep_bad]]
        )
        setkey(bad_pairs, event_row, transcript_id)
        inc_hits <- inc_hits[!bad_pairs]
      }
    }
    if (!nrow(inc_hits)) return(ev_meta[, `:=`(transcript_id = NA_character_, exons = "", inc_exons_by_idx = "", inc_rows_by_idx = "")][order(event_row)])
  }

  # 9) scoring & winner
  class_pref_map <- list(
    SE="internal", MXE="internal", RI="internal",
    A5SS=c("first","internal","last"),
    A3SS=c("last","internal","first"),
    AFE="first", HFE="first", ALE="last", HLE="last"
  )
  inc_hits[, class_pref_hit := mapply(function(et, cls) cls %chin% (class_pref_map[[et]] %||% character()),
                                      event_type, classification)]

  data.table::setorder(inc_hits,
                       event_row, transcript_id, inc_idx,
                       -class_pref_hit, -recip_q, -recip_s, -int_w)
  inc_hits <- inc_hits[, .SD[1L], by = .(event_row, transcript_id, inc_idx)]

  # recompute coverage for require_all_inc
  cov_by_tx <- inc_hits[, .(covered = uniqueN(inc_idx)), by = .(event_row, transcript_id)]
  good_cov  <- cov_by_tx[covered == ninc_map[as.character(event_row)]]
  inc_hits  <- inc_hits[good_cov, on = .(event_row, transcript_id), nomatch = 0L]

  grp_has_pro <- inc_hits[, .(has_pro = any(protein_link)), by = .(event_row, transcript_id)]
  data.table::setkey(grp_has_pro, event_row, transcript_id)
  data.table::setkey(inc_hits,     event_row, transcript_id)
  inc_hits <- inc_hits[grp_has_pro][ has_pro == FALSE | protein_link == TRUE ]
  inc_hits[, has_pro := NULL]

  data.table::setorder(inc_hits, event_row, transcript_id, -class_pref_hit, -recip_q, -int_w)
  best_per_tx <- inc_hits[, .SD[1L], by = .(event_row, transcript_id)]

  supp <- inc_hits[, .(
    inc_exons_all = paste(unique(exon_id),  collapse=";"),
    inc_rows_all  = paste(unique(exon_row), collapse=";"),
    n_inc_covered = uniqueN(inc_idx)
  ), by = .(event_row, transcript_id)]

  data.table::setkey(supp, event_row, transcript_id)
  data.table::setkey(best_per_tx, event_row, transcript_id)
  best_per_tx <- best_per_tx[supp]

  data.table::setorder(best_per_tx, event_row, -protein_link, -class_pref_hit, -recip_q, -int_w)
  final <- best_per_tx[, .SD[1L], by = event_row]

  # 10) per-INC mapping for the chosen transcript
  data.table::setkey(final, event_row, transcript_id)
  data.table::setkey(inc_hits, event_row, transcript_id)
  winner_hits <- inc_hits[final, nomatch = 0L]

  setkey(inc_long, event_row, inc_idx)
  winner_hits <- inc_long[winner_hits, on = .(event_row, inc_idx)]

  per_inc <- winner_hits[
    , .(exon_ids  = paste(unique(exon_id),  collapse = ","),
        exon_rows = paste(unique(exon_row), collapse = ",")),
    by = .(event_row, inc_idx)
  ]

  all_slots <- data.table(event_row = unique(per_inc$event_row))
  all_slots <- all_slots[, .(event_row, inc_idx = seq_len(ninc_map[as.character(event_row)])), by = event_row]
  per_inc   <- all_slots[per_inc, on = .(event_row, inc_idx)]
  per_inc[is.na(exon_ids),  exon_ids := ""]
  per_inc[is.na(exon_rows), exon_rows := ""]

  inc_exons_by_idx <- per_inc[, .(inc_exons_by_idx = paste(exon_ids,  collapse = ";")), by = event_row]
  inc_rows_by_idx  <- per_inc[, .(inc_rows_by_idx  = paste(exon_rows, collapse = ";")), by = event_row]
  data.table::setkey(inc_exons_by_idx, event_row)
  data.table::setkey(inc_rows_by_idx,  event_row)

  # 11) build the simple table (winner transcript + exons) and attach meta + per-INC lists
  winner_exons_long <- winner_hits[
    , .(event_row, inc_idx, exon_id, recip_q, int_w, transcript_id)
  ]
  if (nrow(winner_exons_long)) {
    data.table::setorder(winner_exons_long, event_row, inc_idx, -recip_q, -int_w)
    simple_core <- winner_exons_long[
      , .(transcript_id = transcript_id[1L],
          exons = paste(unique(exon_id), collapse = ";")),
      by = event_row
    ]
  } else {
    simple_core <- data.table::data.table(event_row = integer(), transcript_id = character(), exons = character())
  }

  data.table::setkey(simple_core, event_row)
  out <- ev_meta[simple_core]                 # keep all events
  out <- inc_exons_by_idx[out]
  out <- inc_rows_by_idx[out]
  out[is.na(exons), `:=`(exons = "", transcript_id = NA_character_,
                         inc_exons_by_idx = "", inc_rows_by_idx = "")]
  data.table::setorder(out, event_row)
  out
}

#' Match splicing events to transcript annotations in chunks
#'
#' This is a wrapper around match_events_to_annotations_vec that
#' processes large event tables in manageable chunks to reduce memory usage.
#' It sequentially runs matching per chunk and concatenates the results.
#'
#' @param events A data.frame or data.table of event definitions. Must include
#'   coordinates compatible with match_events_to_annotations_vec
#' @param annotations A data.frame of transcript/exon annotation rows.
#'   Typically generated by get_annotations
#' @param chunk_size Integer, number of event rows to process per chunk
#'   (default = 50,000).
#' @param minOverlap double ranging from 0 to 1 (default 0.05), required minimum
#' overlap to consider a match
#'
#' @return A data.table with matched transcripts and exons for all events.
#'   The output order matches the original event order.
#'
#' @details
#' This function is intended for large-scale event matching across many
#' splicing events, where running the full table at once may exceed memory
#' limits. It can later be parallelized using \pkg{future.apply} or similar.
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
#'
#' @importFrom data.table as.data.table rbindlist setorder %chin% data.table
#' @export
get_matched_events_chunked <- function(events,
                                 annotations,
                                 chunk_size = 50000,
                                 minOverlap = 0.05
) {
  ev <- data.table::as.data.table(events)
  if (!(".__row__" %chin% names(ev))) ev[, .__row__ := .I]

  # indices per chunk
  n <- nrow(ev)
  if (n == 0L) return(data.table::data.table())
  chunk_id <- ceiling(seq_len(n) / chunk_size)
  idx_list <- split(seq_len(n), chunk_id)

  # run chunk-by-chunk
  out_list <- vector("list", length(idx_list))
  for (k in seq_along(idx_list)) {
    ii <- idx_list[[k]]
    cat(sprintf("Chunk %d/%d: rows %d..%d\n", k, length(idx_list), min(ii), max(ii)))
    out_list[[k]] <- match_events_to_annotations_vec(
      events      = ev[ii],
      annotations = annotations,
      minOverlap = minOverlap
    )
  }

  # stitch and restore event order
  ans <- data.table::rbindlist(out_list, fill = TRUE)
  if ("event_row" %chin% names(ans)) data.table::setorder(ans, event_row)
  ans
}


#' Attach transcript and protein sequences to an event or annotation table
#'
#' Merges sequence data (transcript and protein sequences) onto an input table
#' of splicing events or transcript annotations, by matching on
#' \code{transcript_id}. When multiple sequences share the same
#' \code{transcript_id}, the function keeps the entry with a non-missing
#' \code{protein_id} and the longest \code{protein_seq}.
#'
#' @param x A data.frame or data.table containing a \code{transcript_id} column.
#' @param sequences A data.frame or data.table with at least the columns:
#'   \code{transcript_id}, \code{protein_id}, \code{transcript_seq},
#'   and \code{protein_seq}.
#'
#' @return A \link[data.table]{data.table} with the same rows as \code{x} and
#'   appended sequence columns (\code{transcript_seq}, \code{protein_seq}, etc.).
#' @details The join is left-sided: all rows from \code{x} are preserved.
#'   Duplicate \code{transcript_id}s in \code{sequences} are resolved internally
#'   based on protein presence and sequence length.
#'
#' @importFrom data.table as.data.table setorder setkey setcolorder fifelse
#' @export
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
attach_sequences <- function(x, sequences) {
  x  <- data.table::as.data.table(x)
  seq <- data.table::as.data.table(sequences)

  # Keep only what we need and normalize types
  seq <- seq[, .(transcript_id = as.character(transcript_id),
                 protein_id    = as.character(protein_id),
                 transcript_seq= as.character(transcript_seq),
                 protein_seq   = as.character(protein_seq))]

  # If the sequence table has duplicate transcript_ids, keep the "best" one:
  # prioritize rows with a non-NA protein_id, and longer protein_seq (if ties)
  if (any(duplicated(seq$transcript_id))) {
    seq[, has_prot := as.integer(!is.na(protein_id) & nzchar(protein_id))]
    seq[, prot_len := data.table::fifelse(is.na(protein_seq), 0L, nchar(protein_seq))]
    data.table::setorder(seq, transcript_id, -has_prot, -prot_len)
    seq <- seq[!duplicated(transcript_id)]
    seq[, c("has_prot","prot_len") := NULL]
  }

  if (!("transcript_id" %in% names(x)))
    stop("`x` must contain a column named `transcript_id`.")

  x[, transcript_id := as.character(transcript_id)]
  data.table::setkey(seq, transcript_id)

  out <- seq[x, on = "transcript_id"]   # left join; preserves x’s row order
  data.table::setcolorder(out, c(names(x), setdiff(names(out), names(x))))
  out[]
}




