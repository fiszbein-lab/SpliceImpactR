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
#' inclusion-exclusion pair of the same event.
#' @details
#' In \code{source="paired"} mode, only events with exactly one INC and one EXC
#' row are retained. In \code{source="multi"} mode, all positive deltaPSI rows are
#' joined with all negative deltaPSI rows (cartesian join) within each event.
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
    if (!nrow(DT)) DT <- as.data.table(as_dt_from_s4(x, "res_di"))
    if (!nrow(DT)) DT <- as.data.table(as_dt_from_s4(x, "di_events"))
  } else {
    .spi_in <- .resolve_splice_input(x, what = "di_events")
    .spi_obj <- .spi_in$obj
    DT <- as.data.table(.spi_in$dt)
  }

  # --- guards ---
  miss <- setdiff("event_id", names(DT))
  if (length(miss)) stop("Missing key columns: ", paste(miss, collapse=", "))

  setkeyv(DT, "event_id")

  if (source == "paired") {
    need_paired <- c("event_id", "form")
    miss_paired <- setdiff(need_paired, names(DT))
    if (length(miss_paired)) {
      stop("get_pairs(source='paired') missing required columns: ",
           paste(miss_paired, collapse = ", "))
    }

    # keep keys that appear **exactly twice** (one INC + one EXC)
    cnt <- DT[, .N, by = "event_id"]
    keep_keys <- cnt[N == 2L, "event_id"]
    if (!nrow(keep_keys)) return(.return_splice_output(DT[0], obj = .spi_obj, what = "paired_hits", return_class = return_class))

    DT2 <- DT[keep_keys, on = event_id]

    INC <- DT2[form == "INC"]
    EXC <- DT2[form == "EXC"]

    # restrict to keys present in BOTH forms
    common_keys <- intersect(
      do.call(paste, c(INC[, "event_id"], sep = "\r")),
      do.call(paste, c(EXC[, "event_id"], sep = "\r"))
    )
    if (!length(common_keys)) return(.return_splice_output(DT[0], obj = .spi_obj, what = "paired_hits", return_class = return_class))

    INC <- INC[do.call(paste, c(.SD, sep="\r")) %chin% common_keys, .SDcols = event_id]
    EXC <- EXC[do.call(paste, c(.SD, sep="\r")) %chin% common_keys, .SDcols = event_id]

    # drop form col (we’re going to suffix all non-keys)
    INC[, (form) := NULL]
    EXC[, (form) := NULL]

    # suffix and join
    inc_cols <- setdiff(names(INC), "event_id")
    exc_cols <- setdiff(names(EXC), "event_id")
    setnames(INC, inc_cols, paste0(inc_cols, "_case"))
    setnames(EXC, exc_cols, paste0(exc_cols, "_control"))

    setkeyv(INC, "event_id")
    setkeyv(EXC, "event_id")
    out <- EXC[INC, on = event_id, nomatch = 0L]
    data.table::setorderv(out, event_id)
    return(.return_splice_output(out[], obj = .spi_obj, what = "paired_hits", return_class = return_class))
  }

  # ---------- HITindex mode: pair ALL positive with ALL negative within each key ----------
  need_multi <- c(
    "event_id", "gene_id", "transcript_id", "chr", "strand", "event_type",
    "form", "exons", "protein_id", "inc", "exc", "delta_psi",
    "p.value", "padj", "n_samples", "n_control", "n_case",
    "transcript_seq", "protein_seq"
  )
  miss_multi <- setdiff(need_multi, names(DT))
  if (length(miss_multi)) {
    stop("get_pairs(source='multi') missing required columns: ",
         paste(miss_multi, collapse = ", "),
         ". Run annotation matching + sequence attachment before pairing.")
  }

  # keep only keys that have at least one + and one – deltaPSI
  sign_tbl   <- DT[, .(has_pos = any(delta_psi > 0, na.rm=TRUE),
                       has_neg = any(delta_psi < 0, na.rm=TRUE)), by = "event_id"]
  valid_keys <- sign_tbl[has_pos & has_neg, "event_id"]
  if (!nrow(valid_keys)) return(.return_splice_output(DT[0], obj = .spi_obj, what = "paired_hits", return_class = return_class))

  DT2 <- DT[valid_keys, on = 'event_id']

  POS <- DT2[delta_psi > 0]
  NEG <- DT2[delta_psi < 0]

  # suffix and pair (cartesian join)
  left_cols  <- setdiff(names(POS), "event_id")
  right_cols <- setdiff(names(NEG), "event_id")
  setnames(POS, left_cols,  paste0(left_cols,  "_case"))
  setnames(NEG, right_cols, paste0(right_cols, "_control"))

  setkeyv(POS, "event_id")
  setkeyv(NEG, "event_id")
  out <- merge(POS, NEG, by = 'event_id', allow.cartesian = TRUE)

  # order like combine_inc_exc(): by key (then strongest |deltaPSI| on each side if available)
  dA <- paste0('delta_psi', "_case")
  dB <- paste0('delta_psi', "_control")
  out[, ordA := -abs(get(dA))]
  out[, ordB := -abs(get(dB))]
  data.table::setorderv(out, c("event_id", "ordA", "ordB"))
  out[, c("ordA","ordB") := NULL]
  col_names_ord <- c("form", "exons", "protein_id", "inc", "exc", "delta_psi", "p.value", "padj", "n_samples", "n_control", "n_case", "transcript_seq", "protein_seq")
  cn_adjust <- unlist(lapply(col_names_ord, function(x) paste0(x, c("_control","_case"))))
  
  cols_old <- c("event_id", "gene_id_control", "transcript_id_control", "transcript_id_case",
                "chr_control", "strand_control", "event_type_control", cn_adjust)
  
  cols_new <- c("event_id", "gene_id", "transcript_id_control", "transcript_id_case",
                "chr", "strand", "event_type", cn_adjust)
  
  data.table::setcolorder(out, cols_old)
  out <- out[, ..cols_old]
  data.table::setnames(out, old = cols_old, new = cols_new)
  .return_splice_output(out[], obj = .spi_obj, what = "paired_hits", return_class = return_class)
}
