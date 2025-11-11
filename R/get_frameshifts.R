#' Build a per-transcript coding exon index
#'
#' Internal helper that constructs a coding exon index from a full
#' GTF- or annotation-derived data frame. The output table provides
#' exon-level information such as CDS start/end coordinates, reading
#' frames, coding order, and strand direction per transcript.
#'
#' This function filters for exons that have coding sequence coordinates
#' (`cds_gen_start` / `cds_gen_stop`), infers missing frame values, and
#' orders exons by genomic position according to strand. It also
#' computes a `next_exon_id` pointer within each transcript to aid in
#' downstream frame-consistency or “downstream rescue” analyses.
#'
#' @param ann annotations from get_annotations (annotations)
#'
#' @return A `data.table` of coding exons with one row per exon, including:
#'   \describe{
#'     \item{transcript_id}{Transcript identifier.}
#'     \item{exon_id}{Exon identifier.}
#'     \item{strand}{Strand ("+" or "-").}
#'     \item{cds_start, cds_end}{CDS coordinates.}
#'     \item{start_frame, stop_frame}{Reading frames (0–2).}
#'     \item{ord, code_ord, code_class}{Ordering and coding classification.}
#'     \item{next_exon_id}{ID of the next exon in transcript order.}
#'     \item{strand_num}{Numeric strand (+1 or −1).}
#'   }
#'
#' @seealso
#' Used internally by splicing- and frame-preservation analyses.
#'
#' @keywords internal
#' @noRd
build_coding_index <- function(ann) {
  A <- as.data.table(ann)

  # exon id column can be exon_id or exonID
  exon_key <- if ("exon_id" %in% names(A)) "exon_id" else "exonID"
  if (is.null(exon_key)) stop("annotations must have exon_id or exonID")

  E <- A[type %chin% "exon", .(
    transcript_id = as.character(transcript_id),
    exon_id       = as.character(get(exon_key)),
    strand        = as.character(strand),
    cds_start     = as.integer(cds_gen_start),
    cds_end       = as.integer(cds_gen_stop),
    start_frame   = as.integer(start_frame),
    stop_frame    = as.integer(stop_frame),
    ord           = as.integer(absolute_exon_position),
    code_ord      = as.integer(coding_exon_position),
    code_class      = as.character(coding_exon_class)
  )]

  # keep coding exons only
  E <- E[!is.na(cds_start) & !is.na(cds_end)]

  # prefer provided frames; fallback to modulo (0..2) when missing/invalid
  bad_sf <- is.na(E$start_frame) | !(E$start_frame %in% 0:2)
  bad_ef <- is.na(E$stop_frame)  | !(E$stop_frame  %in% 0:2)
  E[bad_sf, start_frame := (cds_start - 1L) %% 3L]
  E[bad_ef, stop_frame  := (cds_end   - 1L) %% 3L]

  # order within transcript; + uses start asc, - uses end desc via sign trick
  E[, t_ord := fifelse(strand == "+", cds_start, -cds_end)]
  data.table::setorder(E, transcript_id, t_ord)

  # next exon pointer (helpful for “downstream rescue” scans)
  E[, next_exon_id := data.table::shift(exon_id, type = "lead"), by = transcript_id]

  # numeric strand helper
  E[, strand_num := fifelse(strand == "+", 1L, -1L)]

  data.table::setkey(E, transcript_id, exon_id)
  E[]
}


#' Identify the terminal coding exon overlap between two transcripts
#'
#' Internal helper that determines the pair of exons (one per transcript)
#' that overlap at the 5′ or 3′ end of coding sequence, depending on
#' event type (AFE vs ALE).
#'
#' @param E `data.table` from [build_coding_index()], containing exon-level
#'   coding coordinates (`cds_start`, `cds_end`, `strand`, etc.).
#' @param tx1,tx2 Character. Transcript IDs to compare.
#' @param mode Character scalar, `"AFE"` or `"ALE"`, indicating whether to
#'   search for the earliest (5′) or latest (3′) coding overlap.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{exon1}{Exon ID from `tx1` participating in terminal overlap.}
#'     \item{exon2}{Exon ID from `tx2` participating in terminal overlap.}
#'   }
#' Returns `NULL` if no valid coding overlap is detected.
#'
#' @seealso [.compare_overlap_frames()], [.find_rescue()]
#' @keywords internal
#' @noRd
.pick_terminal_overlap <- function(E, tx1, tx2, mode = c("AFE","ALE")) {
  mode <- match.arg(mode)

  x <- E[transcript_id == tx1 &
           is.finite(cds_start) & is.finite(cds_end),
         .(exon = exon_id,
           start = as.numeric(cds_start),
           end   = as.numeric(cds_end),
           strand = unique(strand))]

  y <- E[transcript_id == tx2 &
           is.finite(cds_start) & is.finite(cds_end),
         .(exon = exon_id,
           start = as.numeric(cds_start),
           end   = as.numeric(cds_end),
           strand = unique(strand))]

  if (!nrow(x) || !nrow(y) || anyNA(x$start) || anyNA(x$end) || anyNA(y$start) || anyNA(y$end))
    return(NULL)

  # order exons along transcript direction (5'→3')
  if (unique(x$strand) == "+") x[, ord := frank(start, ties.method="first")] else x[, ord := frank(-end, ties.method="first")]
  if (unique(y$strand) == "+") y[, ord := frank(start, ties.method="first")] else y[, ord := frank(-end, ties.method="first")]

  # coding overlap (genomic)
  ov <- outer(x$start, y$end, "<=") & outer(x$end, y$start, ">=")
  if (!any(ov)) return(NULL)

  idx <- which(ov, arr.ind = TRUE)

  # choose terminal pair
  if (mode == "AFE") {
    # earliest overlap toward 5′ end: minimize the later of the two orders
    sel <- idx[ which.min( pmax(x$ord[idx[,1]], y$ord[idx[,2]]) ), ]
  } else {
    # latest overlap toward 3′ end: maximize the earlier of the two orders
    sel <- idx[ which.max( pmin(x$ord[idx[,1]], y$ord[idx[,2]]) ), ]
  }

  list(exon1 = x$exon[sel[1]], exon2 = y$exon[sel[2]])
}

#' Check for one-dimensional genomic overlap
#'
#' @param a1,b1 Integers. Start and end coordinates of the first interval
#'   (must satisfy `a1 <= b1`).
#' @param a2,b2 Integers. Start and end coordinates of the second interval
#'   (must satisfy `a2 <= b2`).
#'
#' @return Logical vector indicating whether the intervals overlap.
#'
#' @keywords internal
#' @noRd
.overlap1d <- function(a1, b1, a2, b2) {
  # assumes a1<=b1 and a2<=b2
  pmax(a1, a2) <= pmin(b1, b2)
}

#' Compare coding reading frames across an overlapping exon pair
#'
#' @param E `data.table` from [build_coding_index()], containing CDS coordinates
#'   (`cds_start`, `cds_end`), strand information, and inferred or annotated
#'   reading frames (`start_frame`, `stop_frame`).
#' @param tx1,tx2 Character. Transcript IDs of the inclusion/exclusion pair.
#' @param e1,e2 Character. Exon IDs defining the overlap to be compared.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{frame_call}{`"PartialMatch"` if frames align at both ends,
#'       or `"FrameShift"` if they do not.}
#'     \item{rescue}{Always `"noRescue"` here; downstream steps may replace this
#'       with a rescue annotation.}
#'   }
#'
#' @details
#' For the `+` strand, frame offsets are computed relative to `cds_start`;
#' for the `−` strand, they are computed relative to `cds_end`, counting
#' backwards. Overlapping coordinates are compared at both the start and end
#' of the shared coding region.
#'
#' This function is strand-aware and tolerant of partial overlaps, returning
#' `"PartialMatch"` if either transcript lacks valid frame data or if the
#' overlap is zero-length.
#'
#' @seealso [.pick_terminal_overlap()], [.find_rescue()]
#' @keywords internal
#' @noRd
.compare_overlap_frames <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]
  y <- E[list(tx2, e2)]
  if (nrow(x) == 0L || nrow(y) == 0L) return(list(frame_call="PartialMatch", rescue="noRescue"))

  # overlap on absolute genomic coordinates
  os <- max(x$cds_start, y$cds_start)
  oe <- min(x$cds_end,   y$cds_end)
  if (!(os <= oe)) return(list(frame_call="PartialMatch", rescue="noRescue"))

  # For + strand: frame at overlap start = (start_frame + (os - cds_start)) %% 3
  # For - strand: frame at overlap start counts from cds_end backwards.
  if (x$strand[1] == "+") {
    f1_start <- (x$start_frame + (os - x$cds_start)) %% 3L
    f1_end   <- (f1_start + (oe - os)) %% 3L
    f2_start <- (y$start_frame + (os - y$cds_start)) %% 3L
    f2_end   <- (f2_start + (oe - os)) %% 3L
  } else {
    f1_start <- (x$start_frame + abs(oe - x$cds_end)) %% 3L
    f1_end   <- (f1_start + abs(os - oe)) %% 3L
    f2_start <- (y$start_frame + abs(oe - y$cds_end)) %% 3L
    f2_end   <- (f2_start + abs(os - oe)) %% 3L
  }

  same_start <- (f1_start == f2_start)
  same_end   <- (f1_end   == f2_end)
  if (same_start && same_end) {
    list(frame_call="PartialMatch", rescue="noRescue")
  } else {
    list(frame_call="FrameShift", rescue="noRescue")  # rescue filled by separate pass if desired
  }
}

#' Identify downstream frame “rescue” between two transcripts
#'
#' @param E `data.table` from [build_coding_index()], containing CDS
#'   coordinates, reading frames (`start_frame`), and coding order
#'   information (`code_ord`, `code_class`).
#' @param tx1,tx2 Character. Transcript identifiers for the pair being compared.
#' @param e1,e2 Character. Exon IDs defining the initial overlapping exons
#'   that were found to be frame‐shifted.
#'
#' @return
#' A character string:
#' \describe{
#'   \item{`"noRescue"`}{if no downstream exons restore frame alignment.}
#'   \item{`"tx1|exon1|tx2|exon2"`}{if a downstream overlapping exon pair
#'     re-enters frame consistency.}
#' }
#'
#' @details
#' The function walks downstream from the input exons (based on coding order)
#' and checks overlapping exons for identical start and end reading frames.
#' If a matching overlap is found before either transcript terminates its
#' CDS (`code_class == 'last'`), it is reported as a potential *rescue* event.
#'
#' Frame offsets are calculated modulo 3 along the CDS, accounting for strand
#' orientation. Only overlapping exons are tested for frame recovery.
#'
#' @seealso [.compare_overlap_frames()], [.incoming_upstream_shift()]
#' @keywords internal
#' @noRd
.find_rescue <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]
  y <- E[list(tx2, e2)]
  if (nrow(x) == 0L || nrow(y) == 0L) return("noRescue")

  # build ordered lists of overlapping coding exons downstream (inclusive)
  X <- E[transcript_id == tx1]
  Y <- E[transcript_id == tx2]
  data.table::setorder(X, code_ord)
  data.table::setorder(Y, code_ord)


  xi <- which(X$exon_id == e1)[1]
  yi <- which(Y$exon_id == e2)[1]
  if (is.na(xi) || is.na(yi)) return("noRescue")


  if (X$code_class[xi] == 'last' | Y$code_class[yi] == 'last') {
    return("noRescue")
  } else {
    xi <- xi+1
    yi <- yi+1
  }

  for (i in xi:nrow(X)) {
    for (j in yi:nrow(Y)) {
      if (anyNA(c(X$cds_start[i], X$cds_end[i], Y$cds_start[j], Y$cds_end[j]))) {
        return("noRescue")
      }
      if (.overlap1d(X$cds_start[i], X$cds_end[i], Y$cds_start[j], Y$cds_end[j])) {
        os <- max(X$cds_start[i], Y$cds_start[j])
        oe <- min(X$cds_end[i],   Y$cds_end[j])

        if (X$strand[i] == "+") {
          f1s <- (X$start_frame[i] + (os - X$cds_start[i])) %% 3L
          f1e <- (f1s + (oe - os)) %% 3L
          f2s <- (Y$start_frame[j] + (os - Y$cds_start[j])) %% 3L
          f2e <- (f2s + (oe - os)) %% 3L
        } else {
          f1s <- (X$start_frame[i] + abs(oe - X$cds_end[i])) %% 3L
          f1e <- (f1s + abs(os - oe)) %% 3L
          f2s <- (Y$start_frame[j] + abs(oe - Y$cds_end[j])) %% 3L
          f2e <- (f2s + abs(os - oe)) %% 3L
        }

        if (f1s == f2s && f1e == f2e) {
          if (X$code_class[i] == 'last' | Y$code_class[j] == 'last') {
            return("noRescue")
          } else {
            return(paste(tx1, X$exon_id[i], tx2, Y$exon_id[j], sep="|"))
          }

        }
      }
    }
  }
  "noRescue"
}


#' Parse a delimited list of exon identifiers
#'
#' @param x Character. String containing one or more exon IDs, separated by
#'   delimiters such as `","`, `";"`, `"|"`, or spaces.
#'
#' @return Character vector of unique, trimmed exon IDs.
#'
#' @keywords internal
#' @noRd
.parse_exon_list <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  unique(trimws(unlist(strsplit(as.character(x), "[,;| ]+"))))
}

#' Select the terminal exon from a set within a transcript
#'
#' @param annotations get_annotations annotations sub data frae
#' @param tx Character. Transcript ID.
#' @param exons_vec Character vector of candidate exon IDs.
#'
#' @return
#' Character scalar of the selected exon ID (the terminal one),
#' or `NA_character_` if no valid match is found.
#'
#' @details
#' This is typically used to identify the last coding exon involved in
#' skipped-exon or mutually exclusive exon events, ensuring frame tests
#' are anchored at the transcript’s downstream boundary.
#'
#' @seealso [.find_rescue()], [.compare_boundary_end()]
#' @keywords internal
#' @noRd
.pick_last_exon <- function(annotations, tx, exons_vec) {
  outE <- exons_vec
  if (length(exons_vec) > 1) {
    if (!length(exons_vec)) return(NA_character_)
    cand <- annotations[transcript_id == tx & exon_id %chin% exons_vec]
    if (!nrow(cand)) return(NA_character_)
    outE <- cand[which.max(absolute_exon_position), exon_id]
  }
  outE
}

#' Compare reading frames at the exon start boundary
#' @param E `data.table` from [build_coding_index()], containing exon-level
#'   frame information (`start_frame`, `stop_frame`, `strand`).
#' @param tx1,tx2 Character. Transcript IDs.
#' @param e1,e2 Character. Exon IDs to compare.
#'
#' @return
#' A list with elements:
#' \describe{
#'   \item{frame_call}{`"PartialMatch"` if frames are consistent, or `"FrameShift"` if not.}
#'   \item{rescue}{Always `"noRescue"` here.}
#' }
#' @seealso [.compare_boundary_end()], [.compare_overlap_frames()]
#' @keywords internal
#' @noRd
.compare_boundary_start <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]; y <- E[list(tx2, e2)]
  if (!nrow(x) || !nrow(y)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  f1 <- if (x$strand[1] == "+") x$start_frame[1] else x$stop_frame[1]
  f2 <- if (y$strand[1] == "+") y$start_frame[1] else y$stop_frame[1]
  if (is.na(f1) || is.na(f2))  return(list(frame_call="PartialMatch", rescue="noRescue"))
  if ((f1 %% 3L) == (f2 %% 3L)) list(frame_call="PartialMatch", rescue="noRescue") else list(frame_call="FrameShift", rescue="noRescue")
}

#' Compare reading frames at the exon end boundary
#' @inheritParams .compare_boundary_start
#'
#' @return
#' A list with `frame_call` (`"PartialMatch"` or `"FrameShift"`) and
#' `rescue = "noRescue"`.
#' @seealso [.compare_boundary_start()], [.compare_at_overlap_end()]
#' @keywords internal
#' @noRd
.compare_boundary_end <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]; y <- E[list(tx2, e2)]
  if (!nrow(x) || !nrow(y)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  f1 <- if (x$strand[1] == "+") x$stop_frame[1] else x$start_frame[1]
  f2 <- if (y$strand[1] == "+") y$stop_frame[1] else y$start_frame[1]
  if (is.na(f1) || is.na(f2))  return(list(frame_call="PartialMatch", rescue="noRescue"))
  if ((f1 %% 3L) == (f2 %% 3L)) list(frame_call="PartialMatch", rescue="noRescue") else list(frame_call="FrameShift", rescue="noRescue")
}

#' Compare reading frames at the start of the coding overlap
#' @inheritParams .compare_boundary_start
#'
#' @return
#' A list with `frame_call` (`"PartialMatch"` or `"FrameShift"`)
#' and `rescue = "noRescue"`.
#' @seealso [.compare_at_overlap_end()], [.compare_overlap_frames()]
#' @keywords internal
#' @noRd
.compare_at_overlap_start <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]; y <- E[list(tx2, e2)]
  if (!nrow(x) || !nrow(y)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  os <- max(x$cds_start, y$cds_start); oe <- min(x$cds_end, y$cds_end)
  if (!(os <= oe)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  if (x$strand[1] == "+") {
    f1 <- (x$start_frame + (os - x$cds_start)) %% 3L
    f2 <- (y$start_frame + (os - y$cds_start)) %% 3L
  } else {
    # “first” base of overlap on −strand is near cds_end
    f1 <- (x$start_frame + abs(oe - x$cds_end)) %% 3L
    f2 <- (y$start_frame + abs(oe - y$cds_end)) %% 3L
  }
  if (f1 == f2) list(frame_call="PartialMatch", rescue="noRescue") else list(frame_call="FrameShift", rescue="noRescue")
}

#' Compare reading frames at the end of the coding overlap
#'
#' @inheritParams .compare_boundary_start
#'
#' @return
#' A list with `frame_call` (`"PartialMatch"` or `"FrameShift"`)
#' and `rescue = "noRescue"`.
#' @seealso [.compare_at_overlap_start()], [.compare_overlap_frames()]
#' @keywords internal
#' @noRd
.compare_at_overlap_end <- function(E, tx1, tx2, e1, e2) {
  x <- E[list(tx1, e1)]; y <- E[list(tx2, e2)]
  if (!nrow(x) || !nrow(y)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  os <- max(x$cds_start, y$cds_start); oe <- min(x$cds_end, y$cds_end)
  if (!(os <= oe)) return(list(frame_call="PartialMatch", rescue="noRescue"))
  if (x$strand[1] == "+") {
    # last base of overlap is oe on +
    f1 <- (x$start_frame + (oe - x$cds_start)) %% 3L
    f2 <- (y$start_frame + (oe - y$cds_start)) %% 3L
  } else {
    # last base of overlap is os on −
    f1 <- (x$start_frame + abs(os - x$cds_end)) %% 3L
    f2 <- (y$start_frame + abs(os - y$cds_end)) %% 3L
  }
  if (f1 == f2) list(frame_call="PartialMatch", rescue="noRescue") else list(frame_call="FrameShift", rescue="noRescue")
}

#' @description
#' `.incoming_upstream_shift()` scans the coding index (`E`) of two transcripts
#' for the first *upstream overlapping exon pair* relative to the event exons (`e1`, `e2`),
#' and checks whether that upstream pair is already out of frame.
#'
#' @details
#' The function walks upstream in transcript order (`t_ord` from [build_coding_index()])
#' using a two-pointer approach until it finds the first overlapping pair of coding exons.
#' Once found, it reuses `.compare_overlap_frames()` to test if the pair is a frame shift.
#'
#' This is primarily used to determine if a local boundary mismatch is a *new* frame shift
#' or inherited from an existing upstream shift.
#'
#' @param E `data.table` produced by [build_coding_index()], containing coding exon metadata
#'   (columns such as `transcript_id`, `exon_id`, `cds_start`, `cds_end`, `t_ord`, etc.).
#' @param tx1,tx2 Character. Transcript IDs to compare.
#' @param e1,e2 Character. Exon IDs corresponding to the event exons in each transcript.
#'
#' @return Logical scalar. `TRUE` if an upstream overlapping exon pair exists *and*
#' that pair is a frame shift according to `.compare_overlap_frames()`.
#' `FALSE` otherwise (including when exons do not exist or no upstream overlap found).
#'
#' @seealso [build_coding_index()], [.compare_overlap_frames()]
#' @keywords internal
#' @noRd
.incoming_upstream_shift <- function(E, tx1, tx2, e1, e2) {
  X <- E[transcript_id == tx1]
  Y <- E[transcript_id == tx2]
  if (!nrow(X) || !nrow(Y)) return(FALSE)

  i <- which(X$exon_id == e1)[1]
  j <- which(Y$exon_id == e2)[1]
  if (is.na(i) || is.na(j)) return(FALSE)

  # Walk upstream in transcript order (t_ord already encodes 5'->3' for both strands)
  ip <- i - 1L
  jp <- j - 1L
  if (ip < 1L || jp < 1L) return(FALSE)

  # Search the nearest upstream overlapping pair
  found <- FALSE
  up_e1 <- up_e2 <- NA_character_

  # two-pointer style walk toward upstream region
  ii <- ip; jj <- jp
  while (ii >= 1L && jj >= 1L && !found) {
    if ( (max(X$cds_start[ii], Y$cds_start[jj]) <= min(X$cds_end[ii], Y$cds_end[jj])) ) {
      found <- TRUE
      up_e1 <- X$exon_id[ii]
      up_e2 <- Y$exon_id[jj]
    } else {
      # advance the exon that is downstream-most (greater t_ord) back toward upstream
      if (X$t_ord[ii] >= Y$t_ord[jj]) {
        ii <- ii - 1L
      } else {
        jj <- jj - 1L
      }
    }
  }

  if (!found) return(FALSE)

  # Reuse your existing overlap-frame comparer on that upstream pair
  cmp <- .compare_overlap_frames(E, tx1, tx2, up_e1, up_e2)
  identical(cmp$frame_call, "FrameShift")
}

#' @title Check for upstream frame shift
#' @keywords internal
#' @description
#' `.needs_upstream_check()` indicates whether a given alternative splicing
#' event type should undergo upstream frame-shift evaluation via
#' `.incoming_upstream_shift()`.
#'
#' @param event_type Character scalar or vector. Event class label(s)
#'   such as `"A5SS"`, `"A3SS"`, `"RI"`, `"SE"`, `"MXE"`, etc.
#'
#' @return Logical vector of same length as input.
#' `TRUE` for event types where an upstream frame-shift scan is relevant.
.needs_upstream_check <- function(event_type) {
  event_type %chin% c("A5SS","A3SS","RI","SE","MXE")
}


#' @title Identify frameshifts and rescue events between transcript pairs
#' @description
#' Function that evaluates whether inclusion and exclusion transcript
#' isoforms for each alternative splicing event maintain the reading frame or
#' induce a frameshift. Optionally identifies "rescue" cases where downstream
#' coding structure re-aligns the frame.
#'
#' @param hits `data.frame` or `data.table` containing splicing event metadata,
#'   typically output from [compare_sequences_alignment()]. Must include
#'   columns such as `event_type_inc`, `transcript_id_inc`, `transcript_id_exc`,
#'   `exons_inc`, `exons_exc`, and `pc_class`.
#' @param annotations from get_annotations (annotations)
#' @param allow_ale_fs Logical (default `FALSE`).
#'   Whether to allow ALE/HLE events to be considered frameshifting.
#'
#' @return A `data.table` identical to `hits` with four appended columns:
#' \describe{
#'   \item{frame_call}{`"FrameShift"` or `"PartialMatch"`.}
#'   \item{rescue}{Rescue classification (e.g. `"noRescue"` or type string).}
#'   \item{frame_check_exon1}{Exon ID used for inclusion isoform comparison.}
#'   \item{frame_check_exon2}{Exon ID used for exclusion isoform comparison.}
#' }
#'
#' @seealso
#' [build_coding_index()],
#' [.compare_boundary_start()],
#' [.compare_boundary_end()],
#' [.compare_at_overlap_start()],
#' [.compare_at_overlap_end()],
#' [.compare_overlap_frames()],
#' [.find_rescue()],
#' [.incoming_upstream_shift()],
#' [.needs_upstream_check()]
#'
#' @keywords internal
compare_frames <- function(hits,
                           annotations,
                           allow_ale_fs = FALSE) {
  print(paste0("[Processing] Identifying frame shifts and rescues"))
  H <- as.data.table(hits)

  E <- build_coding_index(annotations)

  # compute per-row
  res <- H[, {
    if (pc_class == "protein_coding") {

      et <- as.character(event_type_inc)
      tx1 <- as.character(transcript_id_inc)
      tx2 <- as.character(transcript_id_exc)
      e1  <- .parse_exon_list(as.character(exons_inc))
      e2  <- .parse_exon_list(as.character(exons_exc))
      if (et %chin% c("AFE","HFE","ALE","HLE")) {
        mode <- if (et %chin% c("AFE","HFE")) "AFE" else "ALE"
        pick <- .pick_terminal_overlap(E, tx1, tx2, mode)
        e1 <- pick$exon1
        e2 <- pick$exon2
      } else {
        e1 <- .pick_last_exon(annotations, tx1, e1)
        e2 <- .pick_last_exon(annotations, tx2, e2)
      }


      frame_call <- "PartialMatch"
      rescue <- "noRescue"
      if (length(e1) == 1L && length(e2) == 1L &&
          !is.na(e1) && !is.na(e2) &&
          e1 %in% E$exon_id && e2 %in% E$exon_id) {
      # if (e1 %in% E$exon_id & e2 %in% E$exon_id) {
        cmp <- switch(et,
                      "A3SS" = .compare_at_overlap_start(E, tx1, tx2, e1, e2),          # first base of overlap
                      "A5SS" = .compare_boundary_end(E, tx1, tx2, e1, e2),               # exon end (3′ boundary)
                      "RI"   = .compare_at_overlap_end(E, tx1, tx2, e1, e2),             # last base of overlap
                      "SE"   = .compare_boundary_start(E, tx1, tx2, e1, e2),             # start of last exons
                      "MXE"  = .compare_boundary_start(E, tx1, tx2, e1, e2),             # start of last exons
                      "AFE"  = .compare_overlap_frames(E, tx1, tx2, e1, e2),             # terminal overlap (earliest)
                      "HFE"  = .compare_overlap_frames(E, tx1, tx2, e1, e2),
                      "ALE"  = .compare_overlap_frames(E, tx1, tx2, e1, e2),             # terminal overlap (latest)
                      "HLE"  = .compare_overlap_frames(E, tx1, tx2, e1, e2),
                      .compare_boundary_end(E, tx1, tx2, e1, e2)                          # safe default
        )
        frame_call <- cmp$frame_call
        if (frame_call == "FrameShift") {
          rescue <- .find_rescue(E, tx1, tx2, e1, e2)
          upstream_shift <- FALSE
          if (.needs_upstream_check(et)) {
            upstream_shift <- .incoming_upstream_shift(E, tx1, tx2, e1, e2)
            if (upstream_shift == TRUE) {
              frame_call <- "PartialMatch"
              rescue     <- "noRescue"
              # print(upstream_shift)
            }
          }

        }
      }

    } else {
      frame_call <- NA_character_
      rescue <- NA_character_
      e1 <- NA_character_
      e2 <- NA_character_
    }



    .(frame_call = frame_call,
      rescue = rescue,
      frame_check_exon1 = e1,
      frame_check_exon2 = e2)
  }, by = .I]

  # stitch back to input
  out <- cbind(H, res[, .(frame_call, rescue, frame_check_exon1, frame_check_exon2)])
  if (allow_ale_fs == FALSE) {
    out[event_type_exc %chin% c("ALE", "HLE") | event_type_inc %chin% c("ALE", "HLE"), `:=` (frame_call = "PartialMatch",
                                                                                             rescue = "noRescue")]
  }
  print(paste0("[INFO] ", sum(out$frame_call[!is.na(out$frame_call)] == "FrameShift") ," frameshifts (",
               sum(out$rescue[!is.na(out$rescue)] != "noRescue")," rescues) and ",
               sum(out$frame_call[!is.na(out$frame_call)] == "PartialMatch"), " non-frameshifts were identified, "))
  out[]
}


#' @title Compare frame states after sequence alignment
#' @description
#' Wrapper function that performs sequence alignment (via
#' [compare_sequences_alignment()]) and frame-shift analysis (via
#' [compare_frames()]) for a complete set of inclusion/exclusion transcript
#' pairs. It then summarizes each event by a high-level classification label.
#'
#' @details
#' `compare_sequence_frame()` is a convenience function that integrates
#' sequence and frame comparison stages in one call, producing an annotated
#' table suitable for downstream summarization or visualization.
#'
#' The summary label `summary_classification` follows this precedence:
#' 1. `"Match"` - identical protein sequences.
#' 2. `"FrameShift"` - frame disrupted.
#' 3. `"Rescue"` - frame restored downstream.
#' 4. Otherwise, inherited from `pc_class`.
#'
#' @param complete_hits `data.frame` or `data.table` containing complete event
#'   information for inclusion/exclusion transcript pairs, typically from
#'   [get_pairs()] or similar.
#' @param ann Annotation object (output of [get_annotation()]) used for both
#'   sequence alignment and coding index construction.
#'
#' @return A `data.table` containing all columns from `complete_hits`, plus:
#' \describe{
#'   \item{frame_call}{Result from [compare_frames()].}
#'   \item{rescue}{Rescue classification.}
#'   \item{summary_classification}{One of `"FrameShift"`, `"Rescue"`,
#'   `"Match"`, or the original `pc_class`.}
#' }
#'
#' @seealso [compare_frames()], [compare_sequences_alignment()]
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
#' seq_compare <-compare_sequence_frame(pairs, annots$annotations)
#' @export
compare_sequence_frame <- function(complete_hits, ann) {
  hits_compare_sequence <- compare_sequences_alignment(hits = complete_hits, annotations = ann, include_sequences = TRUE, verbose = TRUE)
  hits_compare_frame <- compare_frames(hits = hits_compare_sequence, annotations = ann,allow_ale_fs = FALSE)

  # summarize classifications for plotting
  hits_compare_frame[, summary_classification := pc_class]
  hits_compare_frame[frame_call == 'FrameShift', summary_classification := 'FrameShift']
  hits_compare_frame[rescue != 'noRescue' & !is.na(rescue), summary_classification := 'Rescue']
  hits_compare_frame[protein_seq_exc == protein_seq_inc & !is.na(protein_seq_inc), summary_classification := "Match"]

  return(hits_compare_frame)
}
