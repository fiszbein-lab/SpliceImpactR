for (pkg in c("methods", "S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package `", pkg, "` is required for summarized_experiment_utils.R")
  }
}

#' SpliceImpact result container (S4)
#'
#' @slot raw_events `SummarizedExperiment` of sample/form-level rows.
#' @slot di_events `GRanges` of differential inclusion rows.
#' @slot paired_hits `GRanges` of paired event-level rows.
#' @slot segments `GRangesList` of per-event segment parts (`inc_case`, etc.).
#' @slot metadata list with flags and optional provenance.
setClass(
  "SpliceImpactResult",
  slots = c(
    raw_events = "SummarizedExperiment",
    di_events = "GRanges",
    paired_hits = "GRanges",
    segments = "GRangesList",
    metadata = "list"
  )
)

setValidity("SpliceImpactResult", function(object) {
  msg <- character()

  has_hits <- length(object@paired_hits) > 0L
  has_segments <- length(object@segments) > 0L

  if (has_hits) {
    ph_key <- as.character(S4Vectors::mcols(object@paired_hits)$pair_key)
    if (length(ph_key) == 0L) {
      msg <- c(msg, "paired_hits must contain mcols(pair_key) when non-empty.")
    } else {
      if (anyDuplicated(ph_key)) msg <- c(msg, "paired_hits pair_key must be unique.")
      if (has_segments) {
        seg_names <- names(object@segments)
        if (is.null(seg_names)) msg <- c(msg, "segments must be named when non-empty.")
        else if (!setequal(seg_names, ph_key)) msg <- c(msg, "segments names must match paired_hits$pair_key.")
      }
    }
  } else if (has_segments && is.null(names(object@segments))) {
    msg <- c(msg, "segments must be named when non-empty.")
  }

  di_key <- as.character(S4Vectors::mcols(object@di_events)$di_key)
  if (length(di_key) && anyDuplicated(di_key)) msg <- c(msg, "di_events di_key must be unique.")

  rrd <- SummarizedExperiment::rowData(object@raw_events)
  if ("raw_key" %in% names(rrd)) {
    raw_key <- as.character(rrd$raw_key)
    if (anyDuplicated(raw_key)) msg <- c(msg, "raw_events raw_key must be unique.")
  }

  if (length(msg)) msg else TRUE
})

.parse_span <- function(x) {
  if (is.na(x) || !nzchar(x)) return(IRanges::IRanges())
  p <- strsplit(x, "-", fixed = TRUE)[[1]]
  if (length(p) != 2) return(IRanges::IRanges())
  s <- suppressWarnings(as.integer(p[1]))
  e <- suppressWarnings(as.integer(p[2]))
  if (is.na(s) || is.na(e)) return(IRanges::IRanges())
  IRanges::IRanges(start = min(s, e), end = max(s, e))
}

.parse_spans <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(IRanges::IRanges())
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(IRanges::IRanges())

  starts <- integer(0)
  ends <- integer(0)
  for (p in parts) {
    rr <- .parse_span(p)
    if (length(rr)) {
      starts <- c(starts, IRanges::start(rr))
      ends <- c(ends, IRanges::end(rr))
    }
  }
  if (!length(starts)) return(IRanges::IRanges())
  IRanges::IRanges(start = starts, end = ends)
}

.norm_spans <- function(x) {
  rr <- .parse_spans(x)
  if (!length(rr)) return("")
  ord <- order(IRanges::start(rr), IRanges::end(rr))
  rr <- rr[ord]
  paste0(IRanges::start(rr), "-", IRanges::end(rr), collapse = ";")
}

.mk_pair_key <- function(event_id, inc_case, inc_control, exc_case, exc_control) {
  paste(
    as.character(event_id),
    vapply(inc_case, .norm_spans, character(1)),
    vapply(inc_control, .norm_spans, character(1)),
    vapply(exc_case, .norm_spans, character(1)),
    vapply(exc_control, .norm_spans, character(1)),
    sep = "|"
  )
}

.mk_di_key <- function(event_id, form, inc, exc) {
  paste(
    as.character(event_id),
    as.character(form),
    vapply(inc, .norm_spans, character(1)),
    vapply(exc, .norm_spans, character(1)),
    sep = "|"
  )
}

.mk_raw_key <- function(event_id, form, sample, inc, exc) {
  paste(
    as.character(event_id),
    as.character(form),
    as.character(sample),
    vapply(inc, .norm_spans, character(1)),
    vapply(exc, .norm_spans, character(1)),
    sep = "|"
  )
}

.gr_from_span_col <- function(dt, span_col, chr_col = "chr", strand_col = "strand") {
  n <- nrow(dt)
  starts <- rep.int(1L, n)
  ends <- rep.int(0L, n)
  for (i in seq_len(n)) {
    rr <- .parse_spans(dt[[span_col]][i])
    if (length(rr)) {
      starts[i] <- min(IRanges::start(rr))
      ends[i] <- max(IRanges::end(rr))
    }
  }
  ir <- IRanges::IRanges(start = starts, end = ends)

  chr <- as.character(dt[[chr_col]])
  chr[is.na(chr) | !nzchar(chr)] <- "unknown"
  st <- as.character(dt[[strand_col]])
  st[is.na(st) | !(st %in% c("+", "-", "*"))] <- "*"

  GenomicRanges::GRanges(seqnames = chr, ranges = ir, strand = st)
}

.empty_gr <- function() GenomicRanges::GRanges()
.empty_grl <- function() GenomicRanges::GRangesList()
.empty_se <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      psi = matrix(numeric(0), nrow = 0, ncol = 1, dimnames = list(character(), "value")),
      inclusion_reads = matrix(numeric(0), nrow = 0, ncol = 1, dimnames = list(character(), "value")),
      exclusion_reads = matrix(numeric(0), nrow = 0, ncol = 1, dimnames = list(character(), "value"))
    ),
    rowRanges = .empty_gr()
  )
}

#' Convert differential table to GRanges
#' @keywords internal
as_granges_res <- function(res_dt) {
  dt <- data.table::as.data.table(res_dt)
  need <- c("event_id", "form", "inc", "exc", "chr", "strand")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("as_granges_res: missing required columns: ", paste(miss, collapse = ", "))
  gr <- .gr_from_span_col(dt, span_col = "inc", chr_col = "chr", strand_col = "strand")
  dt[, di_key := .mk_di_key(event_id, form, inc, exc)]
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(dt)
  gr
}

#' Convert paired hits to GRanges
#' @keywords internal
as_granges_hits <- function(hits_dt) {
  dt <- data.table::as.data.table(hits_dt)
  need <- c("event_id", "chr", "strand", "inc_case", "inc_control", "exc_case", "exc_control")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("as_granges_hits: missing required columns: ", paste(miss, collapse = ", "))

  dt[, anchor_span := ifelse(nzchar(inc_case), inc_case, inc_control)]
  dt[, pair_key := .mk_pair_key(event_id, inc_case, inc_control, exc_case, exc_control)]
  gr <- .gr_from_span_col(dt, span_col = "anchor_span", chr_col = "chr", strand_col = "strand")
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(dt[, !"anchor_span"])
  gr
}

#' Convert paired hit segments to GRangesList
#' @keywords internal
as_segments_grl <- function(hits_dt) {
  dt <- data.table::as.data.table(hits_dt)
  need <- c("event_id", "chr", "strand", "inc_case", "inc_control", "exc_case", "exc_control")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("as_segments_grl: missing required columns: ", paste(miss, collapse = ", "))
  dt[, pair_key := .mk_pair_key(event_id, inc_case, inc_control, exc_case, exc_control)]

  mk_one <- function(i, col) {
    rr <- .parse_spans(dt[[col]][i])
    if (!length(rr)) return(GenomicRanges::GRanges())
    st <- as.character(dt$strand[i]); if (is.na(st) || !(st %in% c("+", "-", "*"))) st <- "*"
    chr <- as.character(dt$chr[i]); if (is.na(chr) || !nzchar(chr)) chr <- "unknown"
    GenomicRanges::GRanges(
      seqnames = chr,
      ranges = rr,
      strand = st,
      segment = col,
      part_idx = seq_along(rr),
      raw_span = as.character(dt[[col]][i])
    )
  }

  lst <- lapply(seq_len(nrow(dt)), function(i) c(
    mk_one(i, "inc_case"),
    mk_one(i, "inc_control"),
    mk_one(i, "exc_case"),
    mk_one(i, "exc_control")
  ))
  out <- GenomicRanges::GRangesList(lst)
  names(out) <- dt$pair_key
  out
}

#' Convert sample-level table to SummarizedExperiment
#' @keywords internal
as_se_raw_events <- function(data_dt) {
  dt <- data.table::as.data.table(data_dt)
  need <- c("event_id", "form", "sample", "chr", "strand", "psi", "inclusion_reads", "exclusion_reads")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("as_se_raw_events: missing required columns: ", paste(miss, collapse = ", "))

  rr <- if ("inc" %in% names(dt)) dt$inc else rep("", nrow(dt))
  ex <- if ("exc" %in% names(dt)) dt$exc else rep("", nrow(dt))
  gr <- .gr_from_span_col(
    data.table::data.table(chr = dt$chr, strand = dt$strand, span = rr),
    span_col = "span", chr_col = "chr", strand_col = "strand"
  )

  row_id <- paste(dt$event_id, dt$form, dt$sample, seq_len(nrow(dt)), sep = "|")
  dt[, raw_key := .mk_raw_key(event_id, form, sample, rr, ex)]

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      psi = matrix(as.numeric(dt$psi), ncol = 1, dimnames = list(row_id, "value")),
      inclusion_reads = matrix(as.numeric(dt$inclusion_reads), ncol = 1, dimnames = list(row_id, "value")),
      exclusion_reads = matrix(as.numeric(dt$exclusion_reads), ncol = 1, dimnames = list(row_id, "value"))
    ),
    rowRanges = gr
  )
  rd <- S4Vectors::DataFrame(dt)
  rownames(rd) <- row_id
  SummarizedExperiment::rowData(se) <- rd
  se
}

#' Build S4 SpliceImpact container
#'
#' Accepts any subset of `data`, `res`, and `hits_final`. Missing pieces are
#' stored as empty valid slots and can be added later with [add_splice_part()].
#'
#' @param data Optional sample-level input table.
#' @param res Optional differential inclusion result table.
#' @param hits_final Optional paired/final hit table.
#' @param metadata Optional list.
#' @return `SpliceImpactResult`
#' @export
as_splice_impact_result <- function(data = NULL, res = NULL, hits_final = NULL, metadata = list()) {
  if (!is.null(data) && !is.data.frame(data)) stop("as_splice_impact_result: `data` must be a data.frame/data.table.")
  if (!is.null(res) && !is.data.frame(res)) stop("as_splice_impact_result: `res` must be a data.frame/data.table.")
  if (!is.null(hits_final) && !is.data.frame(hits_final)) stop("as_splice_impact_result: `hits_final` must be a data.frame/data.table.")

  md <- c(metadata, list(
    has_raw = !is.null(data),
    has_di = !is.null(res),
    has_hits = !is.null(hits_final)
  ))

  obj <- methods::new(
    "SpliceImpactResult",
    raw_events = if (!is.null(data)) as_se_raw_events(data) else .empty_se(),
    di_events = if (!is.null(res)) as_granges_res(res) else .empty_gr(),
    paired_hits = if (!is.null(hits_final)) as_granges_hits(hits_final) else .empty_gr(),
    segments = if (!is.null(hits_final)) as_segments_grl(hits_final) else .empty_grl(),
    metadata = md
  )
  methods::validObject(obj)
  obj
}

#' Add one part to an existing SpliceImpactResult
#'
#' @param obj `SpliceImpactResult`
#' @param data Optional raw sample-level table.
#' @param res Optional differential result table.
#' @param hits_final Optional paired/final hits table.
#' @return Updated `SpliceImpactResult`
#' @export
add_splice_part <- function(obj, data = NULL, res = NULL, hits_final = NULL) {
  if (!methods::is(obj, "SpliceImpactResult")) stop("add_splice_part: `obj` must be a SpliceImpactResult.")
  n_parts <- sum(!vapply(list(data, res, hits_final), is.null, logical(1)))
  if (n_parts == 0L) stop("add_splice_part: provide one of `data`, `res`, or `hits_final`.")
  if (n_parts > 1L) stop("add_splice_part: provide only one part at a time.")

  out <- obj
  if (!is.null(data)) {
    out@raw_events <- as_se_raw_events(data)
    out@metadata$has_raw <- TRUE
  } else if (!is.null(res)) {
    out@di_events <- as_granges_res(res)
    out@metadata$has_di <- TRUE
  } else {
    out@paired_hits <- as_granges_hits(hits_final)
    out@segments <- as_segments_grl(hits_final)
    out@metadata$has_hits <- TRUE
  }
  methods::validObject(out)
  out
}

#' Convert S4 slots back to data.table
#'
#' @param x `SpliceImpactResult`
#' @param slot One of `raw_events`, `di_events`, `paired_hits`.
#' @param keep_internal_keys Keep internal key columns (`raw_key`, `di_key`, `pair_key`).
#' @return `data.table`
#' @export
as_dt_from_s4 <- function(x, slot = c("raw_events", "di_events", "paired_hits"), keep_internal_keys = FALSE) {
  slot <- match.arg(slot)
  stopifnot(methods::is(x, "SpliceImpactResult"))

  if (slot == "raw_events") {
    out <- data.table::as.data.table(as.data.frame(SummarizedExperiment::rowData(x@raw_events)))
    if (!keep_internal_keys && "raw_key" %in% names(out)) out[, raw_key := NULL]
    return(out)
  }
  if (slot == "di_events") {
    out <- data.table::as.data.table(as.data.frame(S4Vectors::mcols(x@di_events)))
    if (!keep_internal_keys && "di_key" %in% names(out)) out[, di_key := NULL]
    return(out)
  }
  out <- data.table::as.data.table(as.data.frame(S4Vectors::mcols(x@paired_hits)))
  if (!keep_internal_keys && "pair_key" %in% names(out)) out[, pair_key := NULL]
  out
}

#' Coerce input to data.table for internal pipelines
#' @keywords internal
coerce_to_dt <- function(x, what = c("raw_events", "di_events", "paired_hits")) {
  what <- match.arg(what)
  if (methods::is(x, "SpliceImpactResult")) return(as_dt_from_s4(x, slot = what))
  data.table::as.data.table(x)
}


#' S4 slot and key schema for SpliceImpactResult
#'
#' @return Named list describing slots, core key columns, and assay names.
#' @export
spliceimpact_s4_schema <- function() {
  list(
    slots = list(
      raw_events = "SummarizedExperiment: sample/form rows with rowRanges + rowData",
      di_events = "GRanges: differential inclusion rows with mcols",
      paired_hits = "GRanges: paired case/control rows with mcols",
      segments = "GRangesList: per-pair genomic segments (inc_case/inc_control/exc_case/exc_control)",
      metadata = "list: provenance and pipeline metadata"
    ),
    keys = list(
      raw_events = "raw_key (rowData)",
      di_events = "di_key (mcols)",
      paired_hits = "pair_key (mcols)",
      segments = "names(segments) == pair_key"
    ),
    assays = c("psi", "inclusion_reads", "exclusion_reads")
  )
}

#' Detailed guide for using SpliceImpactResult
#'
#' Prints practical guidance on slots, assays, key columns, and common access
#' patterns for conversion between S4 and `data.table`.
#'
#' @param as_markdown Logical; if `TRUE`, returns guide text instead of printing.
#' @return Invisible character guide text (or visible text if `as_markdown = TRUE`).
#' @export
spliceimpact_s4_guide <- function(as_markdown = FALSE) {
  schema <- spliceimpact_s4_schema()

  txt <- paste(
    "# SpliceImpactResult Guide",
    "",
    "## 1) Slots",
    paste0("- `raw_events`: ", schema$slots$raw_events),
    paste0("- `di_events`: ", schema$slots$di_events),
    paste0("- `paired_hits`: ", schema$slots$paired_hits),
    paste0("- `segments`: ", schema$slots$segments),
    paste0("- `metadata`: ", schema$slots$metadata),
    "",
    "## 2) Assays and row metadata",
    "- Assays in `raw_events`: `psi`, `inclusion_reads`, `exclusion_reads`.",
    "- Per-row columns live in `SummarizedExperiment::rowData(obj@raw_events)`.",
    "- Genomic spans for raw rows are in `SummarizedExperiment::rowRanges(obj@raw_events)`.",
    "",
    "## 3) Keys",
    "- `raw_key`: unique row identity for sample/form-level rows.",
    "- `di_key`: unique row identity for differential rows.",
    "- `pair_key`: unique row identity for paired case/control rows.",
    "- `names(obj@segments)` use the same `pair_key` values.",
    "",
    "## 4) Access patterns",
    "- Differential GRanges metadata: `S4Vectors::mcols(obj@di_events)`.",
    "- Paired-hits metadata: `S4Vectors::mcols(obj@paired_hits)`.",
    "- Segment ranges for one pair: `obj@segments[[pair_key]]`.",
    "",
    "## 5) Convert back to data.table",
    "- Raw rows: `as_dt_from_s4(obj, 'raw_events')`",
    "- Differential rows: `as_dt_from_s4(obj, 'di_events')`",
    "- Paired rows: `as_dt_from_s4(obj, 'paired_hits')`",
    "",
    "## 6) Build and add incrementally",
    "- Build all at once: `as_splice_impact_result(data, res, hits_final)`.",
    "- Add later: `obj <- add_splice_part(obj, hits_final = hits_final)`.",
    sep = "\n"
  )

  if (isTRUE(as_markdown)) return(txt)
  cat(txt, "\n", sep = "")
  invisible(txt)
}

#' List predefined paired-hit column subsets
#'
#' @return Named list of predefined column vectors for paired-hit accessors.
#' @export
spliceimpact_hit_colsets <- function() {
  list(
    core = c(
      "event_id", "event_type", "gene_id", "chr", "strand",
      "transcript_id_control", "transcript_id_case",
      "protein_id_control", "protein_id_case",
      "inc_control", "inc_case", "exc_control", "exc_case",
      "delta_psi_control", "delta_psi_case",
      "padj_control", "padj_case",
      "prot_pid", "frame_call", "summary_classification",
      "diff_n", "n_ppi"
    ),
    domain = c(
      "event_id", "event_type", "gene_id", "chr", "strand",
      "transcript_id_control", "transcript_id_case",
      "exons_control", "exons_case",
      "domains_exons_case", "domains_exons_control",
      "case_only_domains", "control_only_domains",
      "case_only_domains_list", "control_only_domains_list", "either_domains_list",
      "case_only_n", "control_only_n", "diff_n"
    ),
    ppi = c(
      "event_id", "event_type", "gene_id", "chr", "strand",
      "transcript_id_control", "transcript_id_case",
      "protein_id_control", "protein_id_case",
      "case_ppi", "control_ppi",
      "n_case_ppi", "n_control_ppi", "n_ppi",
      "case_pfam_changed", "control_pfam_changed",
      "case_elm_changed", "control_elm_changed"
    ),
    sequence = c(
      "event_id", "event_type", "gene_id", "chr", "strand",
      "transcript_id_control", "transcript_id_case",
      "protein_id_control", "protein_id_case",
      "transcript_seq_control", "transcript_seq_case",
      "protein_seq_control", "protein_seq_case",
      "prot_len_control", "prot_len_case", "prot_len_diff", "prot_len_diff_abs",
      "tx_len_control", "tx_len_case", "tx_len_diff", "tx_len_diff_abs",
      "dna_pid", "dna_score", "dna_width",
      "prot_pid", "prot_score", "prot_width",
      "frame_call", "rescue", "summary_classification"
    )
  )
}

#' Access paired-hits as compact data.table subsets
#'
#' Extracts `paired_hits` columns from a [SpliceImpactResult] (or a paired-hits
#' `data.table`) using predefined subset groups such as `core`, `domain`, `ppi`,
#' and `sequence`.
#'
#' @param x A [SpliceImpactResult] object or a paired-hits `data.frame`/`data.table`.
#' @param col_subset Character vector of subset names. Any of `"core"`,
#'   `"domain"`, `"ppi"`, `"sequence"`, or `"all"`.
#' @param cols Optional explicit column vector. If supplied, `col_subset` is ignored.
#' @param drop_missing Logical; if `TRUE`, silently drops requested columns that
#'   are absent. If `FALSE`, errors on missing columns.
#' @param keep_internal_keys Passed through when `x` is S4. Default `FALSE`.
#'
#' @return A `data.table` containing the selected columns.
#' @export
get_hits_final_view <- function(
    x,
    col_subset = c("core"),
    cols = NULL,
    drop_missing = TRUE,
    keep_internal_keys = FALSE
) {
  if (methods::is(x, "SpliceImpactResult")) {
    dt <- as_dt_from_s4(x, slot = "paired_hits", keep_internal_keys = keep_internal_keys)
  } else {
    dt <- data.table::as.data.table(x)
  }

  if (!is.null(cols)) {
    want <- unique(as.character(cols))
  } else {
    col_subset <- unique(as.character(col_subset))
    if (!length(col_subset)) col_subset <- "core"

    if ("all" %in% col_subset) {
      want <- names(dt)
    } else {
      presets <- spliceimpact_hit_colsets()
      bad <- setdiff(col_subset, names(presets))
      if (length(bad)) {
        stop(
          "get_hits_final_view: unknown col_subset value(s): ",
          paste(bad, collapse = ", "),
          ". Use one or more of: ",
          paste(c(names(presets), "all"), collapse = ", ")
        )
      }
      want <- unique(unlist(presets[col_subset], use.names = FALSE))
    }
  }

  missing <- setdiff(want, names(dt))
  if (length(missing) && !isTRUE(drop_missing)) {
    stop(
      "get_hits_final_view: requested columns missing: ",
      paste(missing, collapse = ", ")
    )
  }

  keep <- intersect(want, names(dt))
  dt[, ..keep]
}

#' Convenience accessor for core paired-hit columns
#'
#' @inheritParams get_hits_final_view
#' @return `data.table` with the `core` subset.
#' @export
get_hits_core <- function(x, drop_missing = TRUE, keep_internal_keys = FALSE) {
  get_hits_final_view(
    x = x,
    col_subset = "core",
    drop_missing = drop_missing,
    keep_internal_keys = keep_internal_keys
  )
}

#' Convenience accessor for domain-focused paired-hit columns
#'
#' @inheritParams get_hits_final_view
#' @return `data.table` with the `domain` subset.
#' @export
get_hits_domain <- function(x, drop_missing = TRUE, keep_internal_keys = FALSE) {
  get_hits_final_view(
    x = x,
    col_subset = "domain",
    drop_missing = drop_missing,
    keep_internal_keys = keep_internal_keys
  )
}

#' Convenience accessor for PPI-focused paired-hit columns
#'
#' @inheritParams get_hits_final_view
#' @return `data.table` with the `ppi` subset.
#' @export
get_hits_ppi <- function(x, drop_missing = TRUE, keep_internal_keys = FALSE) {
  get_hits_final_view(
    x = x,
    col_subset = "ppi",
    drop_missing = drop_missing,
    keep_internal_keys = keep_internal_keys
  )
}

#' Convenience accessor for sequence/frame-focused paired-hit columns
#'
#' @inheritParams get_hits_final_view
#' @return `data.table` with the `sequence` subset.
#' @export
get_hits_sequence <- function(x, drop_missing = TRUE, keep_internal_keys = FALSE) {
  get_hits_final_view(
    x = x,
    col_subset = "sequence",
    drop_missing = drop_missing,
    keep_internal_keys = keep_internal_keys
  )
}
