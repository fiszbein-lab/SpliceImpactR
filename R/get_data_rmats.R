#' @keywords internal
.clean_pair <- function(s, e) {
  bad <- is.na(s) | is.na(e) | e <= s
  s[bad] <- NA_integer_
  e[bad] <- NA_integer_
  return(list(s = s, e = e))
}

#' @keywords internal
.fmt_pair <- function(s, e) {
  ok <- !(is.na(s) | is.na(e) | e <= s)
  out <- character(length(s))
  out[ok] <- paste0(as.integer(s[ok]), "-", as.integer(e[ok]))
  return(out)
}

#' @keywords internal
.collapse3 <- function(a, b, c) {
  out <- a
  addb <- nzchar(b)
  out[addb & nzchar(out)] <- paste0(out[addb & nzchar(out)], ";", b[addb & nzchar(out)])
  out[addb & !nzchar(out)] <- b[addb & !nzchar(out)]

  addc <- nzchar(c)
  out[addc & nzchar(out)] <- paste0(out[addc & nzchar(out)], ";", c[addc & nzchar(out)])
  out[addc & !nzchar(out)] <- c[addc & !nzchar(out)]

  return(out)
}

#' Compute 1-based coordinates for the tail segment
#' when one interval fully overlaps another but differs
#' at exactly one boundary (start or end).
#'
#' @param longS Integer vector of start coordinates for the longer interval.
#' @param longE Integer vector of end coordinates for the longer interval.
#' @param shortS Integer vector of start coordinates for the shorter interval.
#' @param shortE Integer vector of end coordinates for the shorter interval.
#'
#' @return A list with integer vectors:
#'   \item{start}{1-based start positions of the tail region.}
#'   \item{end}{1-based end positions of the tail region.}
#'
#' @details
#' For intervals sharing either their start or end but not both,
#' returns the coordinates of the extra tail portion on the longer interval.
#' Invalid or negative-length intervals are returned as `NA`.
#'
#' @keywords internal
.tail_coords_1based <- function(longS, longE, shortS, shortE) {
  same_start <- !is.na(longS) & !is.na(shortS) & (longS == shortS)
  same_end   <- !is.na(longE) & !is.na(shortE) & (longE == shortE)

  tS <- ifelse(same_start & !same_end, shortE + 1L,
               ifelse(same_end   & !same_start, longS, NA_integer_))
  tE <- ifelse(same_start & !same_end, longE,
               ifelse(same_end   & !same_start, shortS - 1L, NA_integer_))

  # invalidate non-positive length / nonsense
  bad <- !is.na(tS) & !is.na(tE) & (tE < tS)
  tS[bad] <- NA_integer_; tE[bad] <- NA_integer_

  list(start = tS, end = tE)
}

#' Load rMATS event files into standardized data.tables
#'
#' Parses rMATS output (.MATS.JC.txt or .MATS.JCEC.txt) for multiple event types
#' and returns unified event tables ready for downstream inclusion/exclusion processing.
#'
#' @param paths A data.frame with columns \code{path}, \code{sample_name}, and \code{condition}.

#' @param use Character scalar, one of \code{"JC"} or \code{"JCEC"}.
#' @param event_types Event types to include: one or more of
#'   \code{c("SE", "RI", "A5SS", "A3SS", "MXE")}.
#'
#' @return A `data.table` with unified rMATS event annotations, including columns:
#' \itemize{
#'   \item event_type, event_id, gene_id, chr, strand
#'   \item sample, condition (if applicable)
#'   \item delta_psi, pvalue, fdr (if present)
#'   \item inclusion/exclusion read counts
#' }
#' @importFrom data.table fread rbindlist as.data.table setkey setcolorder
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
#' rmats <- load_rmats(sample_frame, event_types = c("MXE", "SE", "A3SS", "A5SS", "RI"))
#'
#' @export
load_rmats <- function(paths,
                       use = c("JC","JCEC")[2],
                       event_types = c("SE","RI","A5SS","A3SS","MXE")) {

  use <- match.arg(use)

  resolve_files <- function(p, event_types, use) {
    if (file.exists(p) && !dir.exists(p)) return(p)
    if (dir.exists(p)) {
      patt <- sprintf("^(%s)\\.MATS\\.%s\\.txt$", paste(event_types, collapse="|"), use)
      return(list.files(p, pattern = patt, full.names = TRUE))
    }
    character(0)
  }

  read_one <- function(f) {
    ev <- sub("\\.MATS\\..*$", "", basename(f))
    dt <- suppressWarnings(fread(f, na.strings = c("NA","NaN","")))
    dt[, `:=`(event_type = ev, source_file = f)]
    dt <- dt[, .SD, .SDcols = unique(names(dt))]
    setcolorder(dt, c("event_type", setdiff(names(dt), "event_type")))
    dt
  }

  # -------- mode B: data.frame with path/sample_name/condition --------
  if (is.data.frame(paths)) {
    req <- c("path","sample_name","condition")
    miss <- setdiff(req, names(paths))
    if (length(miss)) stop("When passing a data.frame, include columns: ", paste(miss, collapse = ", "))

    # In sample-mode, leave compute_summary as-is (user decides); we still parse lists correctly.
    parts <- lapply(seq_len(nrow(paths)), function(i) {
      pth <- paths$path[i]
      samp <- as.character(paths$sample_name[i])
      cond <- as.character(paths$condition[i])
      files <- resolve_files(pth, event_types, use)
      if (!length(files)) stop("No rMATS files found under: ", pth)
      dt <- rbindlist(lapply(files, function(f1) {
        out <- read_one(f1)
        out <- out[, GeneID := sub("\\.\\d+$", "", GeneID)]
        return(out[IJC_SAMPLE_1 != 0 | SJC_SAMPLE_1 != 0])
      }), use.names = TRUE, fill = TRUE)
      dt[, `:=`(sample = samp, condition = cond)]

      # dt[, event_id := sprintf("%s:%s", event_type, as.character(ID))]
      setcolorder(dt, c("sample","condition","event_type", setdiff(names(dt), c("sample","condition","event_type"))))
      dt[]
    })

    DT <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    need <- c("ID","chr","strand","event_type","sample","condition")
    miss <- setdiff(need, names(DT)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))
    setkey(DT, sample, condition, event_type, ID)
    DT[, source := "rmats"]
    return(DT[])
  }

  stop("`paths` must be a character vector of paths OR a data.frame with columns: path, sample_name, condition.")
}

#' Expand rMATS event tables into scalar exon inclusion/exclusion coordinates.
#'
#' Converts rMATS "event" tables (SE, MXE, A3SS, A5SS, RI)
#' into standardized scalar representations with explicit inclusion/exclusion
#' segments for downstream genomic mapping.
#'
#' @param DT A `data.table` or `data.frame` of rMATS output (merged or per-sample).
#'
#' @return A standardized `data.table` containing:
#' \itemize{
#'   \item event_id unique event identifier.
#'   \item event_type rMATS event type (SE, MXE, etc.).
#'   \item form inclusion/exclusion form.
#'   \item gene_id, chr, strand.
#'   \item inc, exc scalar genomic segments (string: e.g. `"100-200;300-400"`).
#'   \item inclusion_reads, exclusion_reads, psi - numeric metrics.
#'   \item condition, sample, source_file - carried forward if present.
#' }
#'
#' @details
#' Handles all five canonical rMATS event types (SE, MXE, A3SS, A5SS, RI),
#' applying strand-aware logic for MXE and coordinate adjustments for A3/A5.
#' Non-standard columns (e.g. IJC_SAMPLE_1) are checked for presence.
#'
#' @importFrom data.table as.data.table copy fifelse setorder setcolorder rbindlist %chin% :=
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
#' rmats <- get_rmats(load_rmats(sample_frame, event_types = c("MXE", "SE", "A3SS", "A5SS", "RI")))
#' @export
get_rmats <- function(DT) {
  x <- data.table::as.data.table(DT)
  x[, strand := fifelse(strand %chin% c("+","-"), as.character(strand), "+")]

  g <- function(nm) if (nm %in% names(x)) as.integer(x[[nm]]) else rep(NA_integer_, nrow(x))

  # Common coords (0-based, half-open)
  upES   <- g("upstreamES");    upEE   <- g("upstreamEE")
  dnES   <- g("downstreamES");  dnEE   <- g("downstreamEE")

  # SE
  seS <- g("exonStart_0base");  seE <- g("exonEnd")

  # MXE
  m1S <- g("1stExonStart_0base"); m1E <- g("1stExonEnd")
  m2S <- g("2ndExonStart_0base"); m2E <- g("2ndExonEnd")

  # A3/A5
  longS <- g("longExonStart_0base"); longE <- g("longExonEnd")
  shS   <- g("shortES");             shE   <- g("shortEE")
  flS   <- g("flankingES");          flE   <- g("flankingEE")
  t <- .tail_coords_1based(longS, longE, shS, shE)

  # # Pre-allocate result skeleton (two copies of base rows)
  base <- x
  INC  <- copy(base); INC[, form := "INC"]
  EXC  <- copy(base); EXC[, form := "EXC"]


  # ---------- SE ----------
  idx <- x$event_type %chin% "SE"
  if (any(idx)) {
    # INC inc: [upES, upEE), [seS, seE), [dnES, dnEE); INC exc: empty
    p1 <- .clean_pair(upES[idx], upEE[idx])
    p2 <- .clean_pair(seS[idx],  seE[idx])
    p3 <- .clean_pair(dnES[idx], dnEE[idx])
    INC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(seS[idx], seE[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx])),
      exc = ""
    )]

    # EXC: inc = [upES,upEE); [dnES,dnEE] ; exc = [seS,seE]
    EXC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx]),
                      rep("", sum(idx))),
      exc = .fmt_pair(seS[idx], seE[idx])
    )]
  }

  # ---------- MXE (strand-aware) ----------
  idx <- x$event_type %chin% "MXE"
  if (any(idx)) {
    plus  <- idx & x$strand == "+"
    minus <- idx & x$strand == "-"

    # + strand: INC includes exon1; EXC includes exon2; EXC.exc = exon1
    if (any(plus)) {
      INC[plus, `:=`(
        inc = .collapse3(.fmt_pair(upES[plus],upEE[plus]),
                         .fmt_pair(m1S[plus],m1E[plus]),
                         .fmt_pair(dnES[plus],dnEE[plus])),
        exc = .fmt_pair(m2S[plus], m2E[plus])
      )]
      EXC[plus, `:=`(
        inc = .collapse3(.fmt_pair(upES[plus],upEE[plus]),
                         .fmt_pair(m2S[plus],m2E[plus]),
                         .fmt_pair(dnES[plus],dnEE[plus])),
        exc = .fmt_pair(m1S[plus], m1E[plus])
      )]
    }

    # − strand: swap 1st/2nd
    if (any(minus)) {
      INC[minus, `:=`(
        inc = .collapse3(.fmt_pair(upES[minus],upEE[minus]),
                         .fmt_pair(m2S[minus],m2E[minus]),
                         .fmt_pair(dnES[minus],dnEE[minus])),
        exc = .fmt_pair(m1S[minus],m1E[minus])
      )]
      EXC[minus, `:=`(
        inc = .collapse3(.fmt_pair(upES[minus],upEE[minus]),
                         .fmt_pair(m1S[minus],m1E[minus]),
                         .fmt_pair(dnES[minus],dnEE[minus])),
        exc = .fmt_pair(m2S[minus], m2E[minus])
      )]
    }
  }

  # ---------- A3SS ----------
  # A3SS rows
  idx_A3 <- x$event_type %chin% "A3SS"
  if (any(idx_A3)) {
    INC[idx_A3, inc := .collapse3(.fmt_pair(longS[idx_A3],longE[idx_A3]),
                                 rep("", sum(idx_A3)),
                                 rep("", sum(idx_A3)))]
    INC[idx_A3, exc := ""]
    EXC[idx_A3, inc := .collapse3(.fmt_pair(shS[idx_A3],shE[idx_A3]),
                                 rep("", sum(idx_A3)),
                                 rep("", sum(idx_A3)))]
    EXC[idx_A3, exc := .fmt_pair(t$start[idx_A3], t$end[idx_A3])]
  }

  # A5SS rows
  idx_A5 <- x$event_type %chin% "A5SS"
  if (any(idx_A5)) {
    INC[idx_A5, inc := .collapse3(.fmt_pair(longS[idx_A5],longE[idx_A5]),
                                 rep("", sum(idx_A5)),
                                 rep("", sum(idx_A5)))]
    INC[idx_A5, exc := ""]
    EXC[idx_A5, inc := .collapse3(.fmt_pair(shS[idx_A5],shE[idx_A5]),
                                 rep("", sum(idx_A5)),
                                 rep("", sum(idx_A5)))]
    EXC[idx_A5, exc := .fmt_pair(t$start[idx_A5], t$end[idx_A5])]
  }

  # ---------- RI ----------
  idx <- x$event_type %chin% "RI"
  if (any(idx)) {
    # INC: inc = upstream piece + intron + downstream piece ; exc = empty
    pU <- .clean_pair(upES[idx], upEE[idx])
    pI <- .clean_pair(upEE[idx], dnES[idx])  # intron
    pD <- .clean_pair(dnES[idx], dnEE[idx])
    INC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                      .fmt_pair(upEE[idx],dnES[idx]),
                      .fmt_pair(dnES[idx],dnEE[idx])),
      exc = ""
    )]

    # EXC: inc = flanks ; exc = intron
    EXC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                      .fmt_pair(dnES[idx],dnEE[idx]),
                      rep("", sum(idx))),
      exc = .fmt_pair(upEE[idx], dnES[idx])
    )]
  }
  dup_remover <- cbind(INC, EXC)
  data.table::setnames(dup_remover,
           (ncol(dup_remover)/2 + 1):ncol(dup_remover),
           paste0("EXC_", names(dup_remover)[(ncol(dup_remover)/2 + 1):ncol(dup_remover)]))
  dup_remover[is.na(inc),     inc := ""]
  dup_remover[is.na(exc),     exc := ""]
  dup_remover[is.na(EXC_inc), EXC_inc := ""]
  dup_remover[is.na(EXC_exc), EXC_exc := ""]
  dup_remover[, event_id := sprintf("%s:%d", event_type,
                                    as.integer(factor(paste(GeneID, chr, inc, exc, EXC_inc, EXC_exc), levels = unique(paste(GeneID, chr, inc, exc, EXC_inc, EXC_exc))))),
              by = event_type]

  dup_remover[, `:=`(inclusion_reads = as.integer(IJC_SAMPLE_1),
                     exclusion_reads = as.integer(SJC_SAMPLE_1))]
  dup_remover$depth <- dup_remover$inclusion_reads + dup_remover$exclusion_reads
  dup_remover$depth[is.na(dup_remover$depth)] <- -Inf
  data.table::setorder(dup_remover, event_id, sample, -depth)
  removers <- !duplicated(dup_remover[, .(event_id, sample)])
  dup_remover <- dup_remover[removers]
  key_i <- c(which(colnames(dup_remover) == 'sample'), which(colnames(dup_remover) == 'exc'),
             which(colnames(dup_remover) == 'EXC_sample'), which(colnames(dup_remover) == 'EXC_exc'),
             which(colnames(dup_remover) == 'event_id'),
             which(colnames(dup_remover) == 'inclusion_reads'), which(colnames(dup_remover) == 'exclusion_reads'))

  INC <- dup_remover[,.SD, .SDcols = c(key_i[1]:key_i[2], key_i[5], key_i[6], key_i[7])]
  EXC <- dup_remover[, .SD, .SDcols = c(key_i[3]:key_i[4], key_i[5])]
  colnames(EXC) <- gsub("EXC_", "", colnames(EXC))

  INC[,  psi := suppressWarnings(as.numeric(IncLevel1))]
  EXC[, `:=`(inclusion_reads = INC$exclusion_reads,
             exclusion_reads = INC$inclusion_reads,
             psi = 1-INC$psi)]

  # Bind INC/EXC rows
  out <- data.table::rbindlist(list(INC, EXC), use.names = TRUE, fill = TRUE)

  data.table::setnames(out, c("GeneID"), c("gene_id"))
  data.table::setcolorder(out, c("event_id", "event_type", "form", "gene_id", "chr", "strand", "inc", "exc", "inclusion_reads", "exclusion_reads", "psi", "sample", "condition", "source_file"))

  return(out[, .SD, .SDcols = c("event_id","event_type","form","gene_id","chr","strand",
                         "inc","exc","inclusion_reads","exclusion_reads","psi",
                         "sample","condition","source_file")])
}


#' @title Import post-differential-inclusion rMATS results
#'
#' @description
#' This function reads post-DI rMATS results and converts them into
#' the standardized SpliceImpactR long format with one row per
#' event x (INC/EXC) form.
#'
#' Input can be:
#' * a data.frame with columns `path`, `grp1`, `grp2`, `event_type`,
#'   in which case each file is read and processed; or
#' * a single rMATS results data.frame, in which case `event_type`
#'   must be supplied.
#'
#' For each event, the function constructs paired INC and EXC entries:
#' * `inc` contains genomic segments included in the form
#' * `exc` contains the excluded segment(s)
#' * `delta_psi`, `p.value`, and `padj` are assigned using
#'   the rMATS-reported values
#'
#' Event IDs are automatically generated (event_type:N) if not supplied
#'
#' @param input Either:
#'   * a data.frame with columns `path`, `grp1`, `grp2`, `event_type`, or
#'   * a data.frame of rMATS post-DI results.
#' @param event_type Optional event type when `input` contains a single
#'   rMATS data.frame. Ignored when file metadata table is supplied.
#'
#' @return A `data.table` with columns:
#' \describe{
#'   \item{event_id}{unique event identifier}
#'   \item{event_type}{splicing event type}
#'   \item{form}{\code{"INC"} or \code{"EXC"}}
#'   \item{gene_id}{gene ID}
#'   \item{chr}{chromosome}
#'   \item{strand}{strand}
#'   \item{inc}{genomic coordinates of included segment(s)}
#'   \item{exc}{genomic coordinates of excluded segment(s)}
#'   \item{p.value}{rMATS p-value}
#'   \item{padj}{FDR}
#'   \item{delta_psi}{signed PSI change (+INC, -EXC)}
#' }
#'
#' @examples
#' # # Multiple files
#' # input <- data.frame(
#' #   path = c('/path/A3SS.MATS.JC.txt', '/path2/A5SS.MATS.JC.txt'),
#' #   grp1 = c("WT","WT"),
#' #   grp2 = c("KO","KO"),
#' #   event_type = c("A3SS", "A5SS")
#' # )
#' # res <- get_rmats_post_di(meta)
#'
#' # Single rMATS table already loaded as df
#' df <- data.frame(
#'   ID = 1L,
#'   GeneID = "ENSG00000182871",
#'   geneSymbol = "COL18A1",
#'   chr = "chr21",
#'   strand = "+",
#'   longExonStart_0base = 45505834L,
#'   longExonEnd = 45505966L,
#'   shortES = 45505837L,
#'   shortEE = 45505966L,
#'   flankingES = 45505357L,
#'   flankingEE = 45505431L,
#'   ID.2 = 2L,
#'   IJC_SAMPLE_1 = "4,1,0",
#'   SJC_SAMPLE_1 = "9,12,3",
#'   IJC_SAMPLE_2 = "0,4,5",
#'   SJC_SAMPLE_2 = "11,15,15",
#'   IncFormLen = 52L,
#'   SkipFormLen = 49L,
#'   PValue = 0.6967562,
#'   FDR = 1,
#'   IncLevel1 = "0.295,0.073,0.0",
#'   IncLevel2 = "0.0,0.201,0.239",
#'   IncLevelDifference = -0.024,
#'   stringsAsFactors = FALSE
#' )
#' res2 <- get_rmats_post_di(df, event_type = "A3SS")
#' @export
get_rmats_post_di <- function(input,
                               event_type=NULL) {

  if (sum(c("path", "event_type") %in% colnames(input)) == 2) {
    out_list <- lapply(seq_len(nrow(input)), function(i) {
      dt <- data.table::fread(input$path[i])
      dt[, event_type := input$event_type[i]]
      dt[, GeneID := tstrsplit(GeneID, "[.]")[[1]]]
      .get_rmats_di_helper(dt)
    })
    return(unique(data.table::rbindlist(out_list, fill=TRUE)))
  } else {
    dt <- data.table::as.data.table(input)
    dt[, event_type := event_type]
    out <- .get_rmats_di_helper(dt)
    return(data.table::data.table(out))
  }
}


#' Helper to adjust locations and setup INC and EXC from di output from rmats
#' @description
#' Internal helper for `get_rmats_post_di()`.
#' Converts one rMATS post-DI result table into long
#' INC/EXC representation and assigns genomic coordinates,
#' `delta_psi`, `p.value`, and `padj`.
#'
#' @details
#' Uses rMATS genomic coordinate columns to derive:
#' * included coordinates (`inc`)
#' * excluded coordinates (`exc`)
#' for each event type (SE, MXE, RI, A3SS, A5SS).
#'
#' Produces exactly two rows per event: INC and EXC,
#' flipping the sign of `delta_psi` for EXC.
#'
#' @param DT rMATS post-DI data.frame from get_rmats_post_di
#'
#' @return A `data.table` with columns:
#' \code{event_id, event_type, form, gene_id, chr, strand,
#' inc, exc, p.value, delta_psi, padj}
#' @keywords internal
#' @noRd
.get_rmats_di_helper <- function(DT) {
  x <- data.table::as.data.table(DT[,.SD, .SDcols = which(!((colnames(DT) %in% c("IJC_SAMPLE_1",
                                                        "SJC_SAMPLE_1",
                                                        "IJC_SAMPLE_2",
                                                        "SJC_SAMPLE_2",
                                                        "IncLevel1",
                                                        "IncLevel2"))))])
  x[, strand := fifelse(strand %chin% c("+","-"), as.character(strand), "+")]

  g <- function(nm) if (nm %in% names(x)) as.integer(x[[nm]]) else rep(NA_integer_, nrow(x))

  # Common coords (0-based, half-open)
  upES   <- g("upstreamES");    upEE   <- g("upstreamEE")
  dnES   <- g("downstreamES");  dnEE   <- g("downstreamEE")

  # SE
  seS <- g("exonStart_0base");  seE <- g("exonEnd")

  # MXE
  m1S <- g("1stExonStart_0base"); m1E <- g("1stExonEnd")
  m2S <- g("2ndExonStart_0base"); m2E <- g("2ndExonEnd")

  # A3/A5
  longS <- g("longExonStart_0base"); longE <- g("longExonEnd")
  shS   <- g("shortES");             shE   <- g("shortEE")
  flS   <- g("flankingES");          flE   <- g("flankingEE")
  t <- .tail_coords_1based(longS, longE, shS, shE)

  # # Pre-allocate result skeleton (two copies of base rows)
  base <- x
  INC  <- copy(base); INC[, form := "INC"]
  EXC  <- copy(base); EXC[, form := "EXC"]


  # ---------- SE ----------
  idx <- x$event_type %chin% "SE"
  if (any(idx)) {
    # INC inc: [upES, upEE), [seS, seE), [dnES, dnEE); INC exc: empty
    p1 <- .clean_pair(upES[idx], upEE[idx])
    p2 <- .clean_pair(seS[idx],  seE[idx])
    p3 <- .clean_pair(dnES[idx], dnEE[idx])
    INC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(seS[idx], seE[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx])),
      exc = ""
    )]

    # EXC: inc = [upES,upEE); [dnES,dnEE] ; exc = [seS,seE]
    EXC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx]),
                       rep("", sum(idx))),
      exc = .fmt_pair(seS[idx], seE[idx])
    )]
  }

  # ---------- MXE (strand-aware) ----------
  idx <- x$event_type %chin% "MXE"
  if (any(idx)) {
    plus  <- idx & x$strand == "+"
    minus <- idx & x$strand == "-"

    # + strand: INC includes exon1; EXC includes exon2; EXC.exc = exon1
    if (any(plus)) {
      INC[plus, `:=`(
        inc = .collapse3(.fmt_pair(upES[plus],upEE[plus]),
                         .fmt_pair(m1S[plus],m1E[plus]),
                         .fmt_pair(dnES[plus],dnEE[plus])),
        exc = .fmt_pair(m2S[plus], m2E[plus])
      )]
      EXC[plus, `:=`(
        inc = .collapse3(.fmt_pair(upES[plus],upEE[plus]),
                         .fmt_pair(m2S[plus],m2E[plus]),
                         .fmt_pair(dnES[plus],dnEE[plus])),
        exc = .fmt_pair(m1S[plus], m1E[plus])
      )]
    }

    # − strand: swap 1st/2nd
    if (any(minus)) {
      INC[minus, `:=`(
        inc = .collapse3(.fmt_pair(upES[minus],upEE[minus]),
                         .fmt_pair(m2S[minus],m2E[minus]),
                         .fmt_pair(dnES[minus],dnEE[minus])),
        exc = .fmt_pair(m1S[minus],m1E[minus])
      )]
      EXC[minus, `:=`(
        inc = .collapse3(.fmt_pair(upES[minus],upEE[minus]),
                         .fmt_pair(m1S[minus],m1E[minus]),
                         .fmt_pair(dnES[minus],dnEE[minus])),
        exc = .fmt_pair(m2S[minus], m2E[minus])
      )]
    }
  }

  # ---------- A3SS ----------
  # A3SS rows
  idx_A3 <- x$event_type %chin% "A3SS"
  if (any(idx_A3)) {
    INC[idx_A3, inc := .collapse3(.fmt_pair(longS[idx_A3],longE[idx_A3]),
                                  rep("", sum(idx_A3)),
                                  rep("", sum(idx_A3)))]
    INC[idx_A3, exc := ""]
    EXC[idx_A3, inc := .collapse3(.fmt_pair(shS[idx_A3],shE[idx_A3]),
                                  rep("", sum(idx_A3)),
                                  rep("", sum(idx_A3)))]
    EXC[idx_A3, exc := .fmt_pair(t$start[idx_A3], t$end[idx_A3])]
  }

  # A5SS rows
  idx_A5 <- x$event_type %chin% "A5SS"
  if (any(idx_A5)) {
    INC[idx_A5, inc := .collapse3(.fmt_pair(longS[idx_A5],longE[idx_A5]),
                                  rep("", sum(idx_A5)),
                                  rep("", sum(idx_A5)))]
    INC[idx_A5, exc := ""]
    EXC[idx_A5, inc := .collapse3(.fmt_pair(shS[idx_A5],shE[idx_A5]),
                                  rep("", sum(idx_A5)),
                                  rep("", sum(idx_A5)))]
    EXC[idx_A5, exc := .fmt_pair(t$start[idx_A5], t$end[idx_A5])]
  }

  # ---------- RI ----------
  idx <- x$event_type %chin% "RI"
  if (any(idx)) {
    # INC: inc = upstream piece + intron + downstream piece ; exc = empty
    pU <- .clean_pair(upES[idx], upEE[idx])
    pI <- .clean_pair(upEE[idx], dnES[idx])  # intron
    pD <- .clean_pair(dnES[idx], dnEE[idx])
    INC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(upEE[idx],dnES[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx])),
      exc = ""
    )]

    # EXC: inc = flanks ; exc = intron
    EXC[idx, `:=`(
      inc = .collapse3(.fmt_pair(upES[idx],upEE[idx]),
                       .fmt_pair(dnES[idx],dnEE[idx]),
                       rep("", sum(idx))),
      exc = .fmt_pair(upEE[idx], dnES[idx])
    )]
  }
  dup_remover <- cbind(INC, EXC)
  data.table::setnames(dup_remover,
                       (ncol(dup_remover)/2 + 1):ncol(dup_remover),
                       paste0("EXC_", names(dup_remover)[(ncol(dup_remover)/2 + 1):ncol(dup_remover)]))
  dup_remover[is.na(inc),     inc := ""]
  dup_remover[is.na(exc),     exc := ""]
  dup_remover[is.na(EXC_inc), EXC_inc := ""]
  dup_remover[is.na(EXC_exc), EXC_exc := ""]

  IncLevelDifference <- dup_remover$IncLevelDifference
  PValue <- dup_remover$PValue
  FDR <- dup_remover$FDR

  
  dup_remover[, event_id := sprintf("%s:%d", event_type,
                                    as.integer(factor(paste(GeneID, chr, inc, exc, EXC_inc, EXC_exc, IncLevelDifference),
                                                      levels = unique(paste(GeneID, chr, inc, exc, EXC_inc, EXC_exc, IncLevelDifference))))),
              by = event_type]


  removers <- !duplicated(dup_remover[, .(event_id)])
  dup_remover <- dup_remover[removers]
  key_i <- c(which(colnames(dup_remover) == 'exc'),
             which(colnames(dup_remover) == 'EXC_exc'),
             which(colnames(dup_remover) == 'event_id'))

  INC <- dup_remover[,.SD, .SDcols = c(2:key_i[1], key_i[3])]
  EXC <- dup_remover[, .SD, .SDcols = c((key_i[1]+2):key_i[2], key_i[3])]
  colnames(EXC) <- gsub("EXC_", "", colnames(EXC))

  INC[,  `:=` (
    delta_psi = suppressWarnings(as.numeric(IncLevelDifference)),
    p.value = suppressWarnings(as.numeric(PValue)),
    padj = suppressWarnings(as.numeric(FDR))
  )]
  EXC[,  `:=` (
    delta_psi = suppressWarnings(-1*as.numeric(IncLevelDifference)),
    p.value = suppressWarnings(as.numeric(PValue)),
    padj = suppressWarnings(as.numeric(FDR))
  )]

  # Bind INC/EXC rows
  out <- data.table::rbindlist(list(INC, EXC), use.names = TRUE, fill = TRUE)

  data.table::setnames(out, c("GeneID"), c("gene_id"))
  data.table::setcolorder(out, c("event_id", "event_type",
                                 "form", "gene_id",
                                 "chr", "strand",
                                 "inc", "exc",
                                 "delta_psi", "p.value", "padj"))

  return(out[, .SD, .SDcols = c("event_id","event_type",
                                "form","gene_id",
                                "chr","strand",
                                "inc","exc", "p.value",
                                "delta_psi", "padj")])
}

















































