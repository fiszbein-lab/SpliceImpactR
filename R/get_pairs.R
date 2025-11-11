#' Pair inclusion and exclusion forms of splicing events
#'
#' Builds paired tables of inclusion/exclusion forms for splicing events from
#' rMATS-like or HITindex-like inputs. In rMATS mode, events are paired when
#' both INC and EXC forms exist for a given event ID. In HITindex mode, all
#' positive and negative deltaPSI rows within each event are cross-joined.
#'
#' @param x A data.frame or data.table containing splicing event information.
#' @param source Character string specifying input structure:
#'   \describe{
#'     \item{\code{"paired"}}{(rMATS-like) requires exactly one INC and one EXC
#'       per event ID.}
#'     \item{\code{"multi"}}{(HITindex-like) pairs all positive and negative
#'       \code{delta_psi} values within each event.}
#'   }
#'
#' @return A \link[data.table]{data.table} where each row represents an
#'   inclusion-exclusion pair of the same event.
#' @details
#' In \code{source="paired"} mode, only events with exactly one INC and one EXC
#' row are retained. In \code{source="multi"} mode, all positive deltaPSI rows are
#' joined with all negative deltaPSI rows (cartesian join) within each event.
#'
#' @importFrom data.table as.data.table setkeyv setnames setorderv
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
#' pairs <- get_pairs(x_seq, source="multi")
get_pairs <- function(x,
                      source = c("paired","multi")) {

  source  <- match.arg(source)
  DT      <- as.data.table(x)

  # --- guards ---
  miss <- setdiff("event_id", names(DT))
  if (length(miss)) stop("Missing key columns: ", paste(miss, collapse=", "))
  if (source == "rmats"    && !("form" %in% names(DT)))
    stop("rmats mode requires column '", "form", "'.")
  if (source == "hitindex" && !("delta_psi" %in% names(DT)))
    stop("hitindex mode requires column '", "delta_psi", "'.")

  setkeyv(DT, "event_id")

  if (source == "rmats") {
    # keep keys that appear **exactly twice** (one INC + one EXC)
    cnt <- DT[, .N, by = "event_id"]
    keep_keys <- cnt[N == 2L, "event_id"]
    if (!nrow(keep_keys)) return(DT[0])

    DT2 <- DT[keep_keys, on = event_id]

    INC <- DT2[form == "INC"]
    EXC <- DT2[form == "EXC"]

    # restrict to keys present in BOTH forms
    common_keys <- intersect(
      do.call(paste, c(INC[, "event_id"], sep = "\r")),
      do.call(paste, c(EXC[, "event_id"], sep = "\r"))
    )
    if (!length(common_keys)) return(DT[0])

    INC <- INC[do.call(paste, c(.SD, sep="\r")) %chin% common_keys, .SDcols = event_id]
    EXC <- EXC[do.call(paste, c(.SD, sep="\r")) %chin% common_keys, .SDcols = event_id]

    # drop form col (we’re going to suffix all non-keys)
    INC[, (form) := NULL]
    EXC[, (form) := NULL]

    # suffix and join
    inc_cols <- setdiff(names(INC), "event_id")
    exc_cols <- setdiff(names(EXC), "event_id")
    setnames(INC, inc_cols, paste0(inc_cols, "_inc"))
    setnames(EXC, exc_cols, paste0(exc_cols, "_exc"))

    setkeyv(INC, "event_id")
    setkeyv(EXC, "event_id")
    out <- EXC[INC, on = event_id, nomatch = 0L]
    data.table::setorderv(out, event_id)
    return(out[])
  }

  # ---------- HITindex mode: pair ALL positive with ALL negative within each key ----------
  # keep only keys that have at least one + and one – deltaPSI
  sign_tbl   <- DT[, .(has_pos = any(delta_psi > 0, na.rm=TRUE),
                       has_neg = any(delta_psi < 0, na.rm=TRUE)), by = "event_id"]
  valid_keys <- sign_tbl[has_pos & has_neg, "event_id"]
  if (!nrow(valid_keys)) return(DT[0])

  DT2 <- DT[valid_keys, on = 'event_id']

  POS <- DT2[delta_psi > 0]
  NEG <- DT2[delta_psi < 0]

  # suffix and pair (cartesian join)
  left_cols  <- setdiff(names(POS), "event_id")
  right_cols <- setdiff(names(NEG), "event_id")
  setnames(POS, left_cols,  paste0(left_cols,  "_inc"))
  setnames(NEG, right_cols, paste0(right_cols, "_exc"))

  setkeyv(POS, "event_id")
  setkeyv(NEG, "event_id")
  out <- merge(POS, NEG, by = 'event_id', allow.cartesian = TRUE)

  # order like combine_inc_exc(): by key (then strongest |deltaPSI| on each side if available)
  dA <- paste0('delta_psi', "_inc")
  dB <- paste0('delta_psi', "_exc")
  out[, ordA := -abs(get(dA))]
  out[, ordB := -abs(get(dB))]
  data.table::setorderv(out, c("event_id", "ordA", "ordB"))
  out[, c("ordA","ordB") := NULL]


  out[]
}
