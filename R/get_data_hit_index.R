#' Parse exon coordinate strings of form "chr:start-end" into components.
#'
#' @param exon_col Character vector of exon coordinate strings.
#' @return A list with elements \code{chr}, \code{start}, and \code{end}.
#' @keywords internal
.parse_exon_coords <- function(exon_col) {
  stopifnot(is.character(exon_col))
  pattern <- ".*\\b(chr?[A-Za-z0-9_.]+):([0-9]+)-([0-9]+).*"
  chr <- sub(pattern, "\\1", exon_col, perl = TRUE)
  s   <- as.integer(sub(pattern, "\\2", exon_col, perl = TRUE))
  e   <- as.integer(sub(pattern, "\\3", exon_col, perl = TRUE))
  bad <- !grepl("\\b[0-9]+-[0-9]+\\b", exon_col)
  chr[bad] <- NA_character_; s[bad] <- NA_integer_; e[bad] <- NA_integer_
  list(chr = chr, start = s, end = e)
}

#' Combine start/end coordinates into "start-end" strings (1-based).
#'
#' @param start,end Integer vectors of equal length.
#' @return Character vector of "start-end" strings (NA for invalid pairs).
#' @keywords internal
.coords_to_string_1b <- function(start, end) {
  stopifnot(length(start) == length(end))
  ok <- !(is.na(start) | is.na(end) | end < start)
  out <- rep(NA_character_, length(start))
  out[ok] <- paste0(as.integer(start[ok]), "-", as.integer(end[ok]))
  out
}


#' Locate HIT index PSI files within a directory or return the file if directly provided.
#'
#' @param p Directory path or file path.
#' @return Character vector of full paths to matching files.
#' @keywords internal
.find_hitindex_files <- function(p) {
  if (dir.exists(p)) {
    list.files(
      p,
      pattern = "(^|\\.)((AFE|ALE|HFE|HLE)PSI)(\\.|$)",
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else if (file.exists(p)) {
    p
  } else {
    character(0)
  }
}

#' Read exon-level counts or annotations used in HIT index files.
#'
#' @param path Base path to HIT index outputs (excluding suffix like ".exon").
#' @param columns Columns to select from the exon file.
#' @return A data.table with the specified columns.
#' @keywords internal
.read_exon_files <- function(path, columns = c("gene", "exon", "ID")) {
  exon_path <- paste0(path, "exon")
  if (!file.exists(exon_path))
    stop("Expected exon file not found at: ", exon_path)
  data.table::fread(exon_path, select = columns, showProgress = FALSE)
}

#' Read and merge a single HIT index PSI file with its exon annotations.
#'
#' @param f Path to a .AFEPSI/.ALEPSI/etc. file.
#' @param sample Sample name.
#' @param condition Experimental condition.
#' @return A standardized data.table with inclusion/exclusion metrics and metadata.
#' @keywords internal
.read_one_hit <- function(f, sample, condition) {
  # --- Determine event type from filename ---
  ev_guess <- toupper(basename(f))
  ev <- if (grepl("AFE", ev_guess)) {
    "AFE"
  } else if (grepl("ALE", ev_guess)) {
      "ALE"
    } else if (grepl("HFE", ev_guess)) {
    "HFE"
    } else if (grepl("HLE", ev_guess)) {
    "HLE"
    } else {"AFE"}

  dt_init <- suppressWarnings(data.table::fread(f, na.strings = c("NA", "NaN"), quote = ""))

  exon_limited <- .read_exon_files(gsub(paste0(ev, "PSI"), "", f, ignore.case = TRUE))
  dt <- exon_limited[dt_init, on = .(gene, exon)]

  # rename PSI column
  psi_col <- grep("PSI$", names(dt), value = TRUE, ignore.case = TRUE)
  if (length(psi_col)) data.table::setnames(dt, psi_col[1], "PSI_raw")

  coord <- .parse_exon_coords(dt$exon)
  included <- .coords_to_string_1b(coord$start, coord$end)

  nUP <- suppressWarnings(as.numeric(dt$nUP))
  nDOWN <- suppressWarnings(as.numeric(dt$nDOWN))
  inc_r <- abs(nUP - nDOWN)

  total_col <- if ("sumR-L" %in% names(dt)) {
    "sumR-L"
  } else if ("sumL-R" %in% names(dt)) {
      "sumL-R"
    } else {stop("Neither 'sumR-L' nor 'sumL-R' present in file: ", f)}
  total_r <- dt[[total_col]]
  exclusion_r <- total_r - inc_r

  class_map <- c(
    "first" = "first", "FirstInternal_high" = "hybrid",
    "FirstInternal_medium" = "hybrid", "last" = "last",
    "LastInternal_high" = "hybrid", "LastInternal_medium" = "hybrid"
  )

  gene_stripped <- sub("\\.\\d+$", "", dt$gene)
  class_i <- class_map[dt$ID]

  out <- data.table::data.table(
    event_type = ev,
    form = "SITE",
    chr = coord$chr,
    strand = as.character(dt$strand %||% "+"),
    inc = included,
    exc = "",
    gene_id = gene_stripped,
    sample = sample,
    condition = condition,
    inclusion_reads = inc_r,
    exclusion_reads = exclusion_r,
    total_reads = total_r,
    psi_original = suppressWarnings(as.numeric(dt$PSI_raw)),
    nFE = suppressWarnings(as.numeric(dt$nFE %||% 0)),
    nLE = suppressWarnings(as.numeric(dt$nLE %||% 0)),
    nUP = nUP,
    nDOWN = nDOWN,
    nTXPT = suppressWarnings(as.numeric(dt$nTXPT)),
    HITindex = suppressWarnings(as.numeric(dt$HITindex)),
    source = "hitindex",
    source_file = f,
    class = class_i
  )
  out[]
}



#' Load HIT index PSI files for one or more samples/conditions.
#'
#' @param paths_df Data.frame with columns: \code{path}, \code{condition}, and optionally \code{sample_name}.
#' @param keep_annotated_first_last Logical; if TRUE, retain only annotated first/last exons and normalize PSI.
#' @return A standardized `data.table` of HIT index PSI values with inclusion/exclusion and metadata.
#'
#' @importFrom data.table rbindlist setcolorder
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
#' hit_index <- get_hitindex(sample_frame, keep_annotated_first_last = TRUE)
#'
#' @export
get_hitindex <- function(paths_df, keep_annotated_first_last = FALSE) {
  stopifnot(is.data.frame(paths_df))
  req <- c("path", "condition")
  miss <- setdiff(req, names(paths_df))
  if (length(miss)) stop("Data frame must include: ", paste(req, collapse = ", "))

  if (!"sample_name" %in% names(paths_df)) {
    paths_df$sample_name <- basename(normalizePath(paths_df$path, mustWork = FALSE))
  }

  parts <- lapply(seq_len(nrow(paths_df)), function(i) {
    files <- .find_hitindex_files(paths_df$path[i])
    if (!length(files))
      stop("No .AFEPSI/.ALEPSI files found under: ", paths_df$path[i])

    data.table::rbindlist(
      lapply(files, function(f) {
        .read_one_hit(f, paths_df$sample_name[i], paths_df$condition[i])
      }),
      use.names = TRUE, fill = TRUE
    )
  })

  DT <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

  if (isTRUE(keep_annotated_first_last)) {
    event_count <- nrow(DT)
    DT <- DT[
      (event_type %chin% c("ALE", "HLE") & nLE != 0) |
        (event_type %chin% c("AFE", "HFE") & nFE != 0)
    ]
    DT[, total := sum(psi_original, na.rm = TRUE), by = .(gene_id, sample, event_type)]
    DT[, psi := psi_original / total]
    message(sprintf("[INFO] Filtered %d unannotated ALE/AFE events (%d -> %d)",
                    event_count - nrow(DT), event_count, nrow(DT)))
  } else {
    DT[, psi := psi_original]
  }

  DT[, event_id := paste0(gene_id, ":", event_type)]
  data.table::setcolorder(DT, c(
    "event_id", "event_type", "form", "gene_id", "chr", "strand",
    "inc", "exc", "inclusion_reads", "exclusion_reads", "psi",
    "sample", "condition", "source_file", "HITindex", "class"
  ))

  DT[]
}


#' Wrapper function to get both rmats and hit index cleanly
#' @param sample_frame Data.frame with columns: \code{path}, \code{condition}, and \code{sample_name}.
#' @param event_types event types to load from rMATS
#' @param use Character scalar, one of \code{"JC"} or \code{"JCEC"}.
#' @param keep_annotated_first_last Logical; if TRUE, retain only annotated first/last exons and normalize PSI.
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
#' data <- get_rmats_hit(sample_frame, event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI"))
#' @return a `data.table` for all the event types desired from the paths supplied
#' contains: event_id (unique id for event), event_type (AS event type), form
#' (INC/EXC), gene_id (ensembl id), strand, inc, exc (inclusion and exclusion coords)
#' inclusion reads, exclusion reads, psi, sample, condition, source file
#' @export
get_rmats_hit <- function(sample_frame,
                          event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI"),
                          use = 'JCEC',
                          keep_annotated_first_last = TRUE) {
  if ("ALE" %in% event_types | "AFE" %in% event_types) {
    hit_index <- get_hitindex(sample_frame, keep_annotated_first_last)
    if (sum(c("MXE", "SE", "A3SS", "A5SS", "RI") %in% event_types) > 0) {
      rmats <- get_rmats(load_rmats(sample_frame, use, event_types))
      data <- rbind(rmats, hit_index[, .SD, .SDcols = seq(1, ncol(rmats))])
    } else {
      data <- hit_index
    }

  } else if (sum(c("MXE", "SE", "A3SS", "A5SS", "RI") %in% event_types) > 0) {
    data <- get_rmats(load_rmats(sample_frame, use, event_types))
  }
  return(data)
}








































