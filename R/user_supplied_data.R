#' @title Get User-Supplied Splicing Event Data
#'
#' @description
#'
#' SpliceImpactR accepts user-supplied splicing data in the same structure
#' produced by \code{get_rmats_hit()}. Each event must be represented at the
#' event-form and sample level.
#'
#' **Event representation**
#' - Each splicing event must be split into two forms:
#'   - \code{INC}: the inclusion isoform
#'   - \code{EXC}: the exclusion isoform
#' - Alternative first/last exon events (AFE/ALE) or more abstract events may
#'  instead provide a single \code{SITE} form.
#'
#' **Coordinates**
#' - \code{inc} column: genomic coordinates included in the given form
#' - \code{exc} column: genomic coordinates excluded in the given form
#' - Coordinates may be one or multiple ranges (e.g., "100-200" or "100-150;300-350")
#' - User must supply at least an inc and if supplying an exc, accompany with an inc coord
#'
#' **Counts and PSI**
#' - \code{inclusion_reads} and \code{exclusion_reads} must be provided per form
#' - \code{psi} must be provided per sample (range 0-1)
#' If psi isn't given, it will be extracted through
#' \code{inclusion_reads}/\code{exclusion_reads}
#'
#' **Sample structure**
#' - Each event must have \code{INC} and \code{EXC} rows (or just \code{SITE})
#' - Each event must have >1 sample per condition (e.g., case vs control)
#' - Required sample annotations: \code{sample}, \code{condition}
#'
#' **Required columns**
#' \code{event_id, event_type, form, gene_id, chr, strand, inc, exc,
#'  inclusion_reads, exclusion_reads, psi, sample, condition}
#'
#' **Defaults**
#' - If \code{event_type} not supplied: filled as \code{"unknown"}
#' - If \code{source_file} not supplied: filled with empty string
#'
#' This format enables downstream functionality
#' including PSI modeling, annotation integration, and protein consequence
#' prediction.
#' @param df Data frame with splicing events. Detailed in description
#'
#' @return data.table with cols, detailed above: "event_id","event_type","form",
#' "gene_id","chr", "strand", "inc","exc","inclusion_reads","exclusion_reads",
#' "psi", "sample", "condition","source_file" -- designed to match get_rmats_hit
#' output
#'
#' @examples
#' example_df <- data.frame(
#'   event_id = rep("A3SS:1", 8),
#'   event_type = "A3SS",
#'   form = rep(c("INC","EXC"), each = 4),
#'   gene_id = "ENSG00000158286",
#'   chr = "chrX",
#'   strand = "-",
#'   inc = c(rep("149608626-149608834",4), rep("149608626-149608829",4)),
#'   exc = c(rep("",4), rep("149608830-149608834",4)),
#'   inclusion_reads = c(30,32,29,31, 2,3,4,3),
#'   exclusion_reads = c(1,1,2,1, 28,27,26,30),
#'   sample = c("S1","S2","S3","S4","S1","S2","S3","S4"),
#'   condition = rep(c("case","case","control","control"), 2),
#'   stringsAsFactors = FALSE
#' )
#' example_df$psi < example_df$inclusion_reads / example_df$exclusion_reads
#' user_data <- get_user_data(example_df)
#' print(user_data)
#'
#' @export
get_user_data <- function(df) {
  dt <- data.table::as.data.table(df)

  required_cols <- c(
    "event_id","form","gene_id","chr","strand",
    "inc","exc","inclusion_reads","exclusion_reads","sample","condition"
  )

  missing <- setdiff(required_cols, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse=", "))
  }

  # Defaults
  if (!"source_file" %in% names(dt)) dt[, source_file := ""]
  if (!"event_type" %in% names(dt)) dt[, event_type := "unknown"]

  # Check INC/EXC or SITE
  ok_events <- dt[, {
    forms <- unique(form)
    valid <- ("INC" %in% forms & "EXC" %in% forms) || ("SITE" %in% forms)
    .(ok = valid)
  }, by = event_id]$ok

  if (!all(ok_events)) stop("Some events do not have INC+EXC or SITE")

  # Disallow events that mix INC/EXC with SITE
  mixed <- dt[, {
    f <- unique(form)
    .(bad = ("SITE" %in% f) & (any(f %in% c("INC","EXC"))))
  }, by = event_id]$bad

  if (any(mixed)) stop("Some events contain both SITE and INC/EXC forms")

  # Check >1 sample per condition
  cond_check <- dt[, .N, by = .(event_id, condition)][, all(N > 1), by = event_id]$V1
  if (!all(cond_check)) stop("Some events do not have >1 sample per condition")

  if (!"psi" %in% names(dt)) {
    dt[, psi := inclusion_reads / (inclusion_reads + exclusion_reads)]
  } else {
    # If psi present but has NAs, fill those too
    dt[is.na(psi), psi := inclusion_reads / (inclusion_reads + exclusion_reads)]
  }
  final_order <- c(
    "event_id","event_type","form","gene_id","chr","strand",
    "inc","exc","inclusion_reads","exclusion_reads","psi",
    "sample","condition","source_file"
  )

  dt <- dt[, c(final_order, setdiff(names(dt), final_order)), with = FALSE]
  return(dt)
}


#' @title Format user-supplied post-DI (post-differential-inclusion) splicing results
#'
#' @description
#' This function converts user-supplied event results into the internal
#' SpliceImpactR DI format. It accepts per-event statistics and ensures
#' each splicing event contains valid inclusion/exclusion structure.
#'
#' ## Requirements
#' **Event representation**
#' - Each splicing event must be split into two forms:
#'   - \code{INC}: the inclusion isoform
#'   - \code{EXC}: the exclusion isoform
#' - Alternative first/last exon events (AFE/ALE) or more abstract events may
#'  instead provide a single \code{SITE} form.
#'
#' **Coordinates**
#' - \code{inc} column: genomic coordinates included in the given form
#' - \code{exc} column: genomic coordinates excluded in the given form
#' - Coordinates may be one or multiple ranges (e.g., "100-200" or "100-150;300-350")
#'
#' - User **must supply `event_id`** (unique ID per splicing event, event_id =
#' event_type:x for form != SITE. event_id = gene:event_type for form = SITE.
#' Look at example output from get_differential_inclusion using test data for
#' more examples)
#' - Each event must be either:
#'   - **INC + EXC** forms (paired isoforms), or
#'   - **SITE** (single isoform)
#' - Required columns:
#'   `gene_id, chr, strand, inc, exc, form, event_id`
#'
#' ## Behavior
#' - Does **not** generate event IDs: user must provide them
#' - Constructs `site_id = event_type|gene_id|chr|inc|exc|form`
#' - `delta_psi`:
#'   - If missing: INC = +1, EXC = -1, SITE expands to +1 and -1 each to get
#'   all relevant comparisons
#' - `p.value`, `padj`:
#'   - If missing: set to 0
#' - All diagnostic fields
#'   (`cooks_max, n, n_used, n_samples, n_case,
#'     mean_psi_ctrl, mean_psi_case, n_control`)
#'    set to -1 if missing
#'
#' ## Validation
#' - Throws error if:
#'   - Any event lacks **INC+EXC** or **SITE**
#'   - An event mixes SITE with INC/EXC
#'
#' ## Output
#' Returns a `data.table` formatted like SpliceImpactR DI output.
#' Ready for annotation, pairing, enrichment, and PPI analysis.
#'
#' @param df Data frame of post-DI results, detailed in description
#' @examples
#' example_user_data <- data.frame(
#'   event_id = rep("A3SS:1", 8),
#'   event_type = "A3SS",
#'   gene_id = "ENSG00000158286",
#'   chr = "chrX",
#'   strand = "-",
#'   form = rep(c("INC","EXC"), each = 4),
#'   inc = c(
#'     rep("149608626-149608834", 4),
#'     rep("149608626-149608829", 4)
#'   ),
#'   exc = c(
#'     rep("", 4),
#'     rep("149608830-149608834", 4)
#'   ),
#'   inclusion_reads = c(30, 28, 25, 32,  2, 3, 4, 3),
#'   exclusion_reads = c(1, 2, 1, 1, 28, 27, 26, 30),
#'   sample = c("S1","S2","S3","S4","S1","S2","S3","S4"),
#'   condition = rep(c("case","case","control","control"), 2),
#'   stringsAsFactors = FALSE
#' )
#'
#' # compute psi if missing, just for demo
#'
#' post_di_user_data <- get_user_data_post_di(example_user_data)
#' print(post_di_user_data)
#'
#' @return A data.table formatted like SpliceImpactR DI output
#' (get_differential_inclusion)
#' @importFrom data.table copy as.data.table
#' @export
get_user_data_post_di <- function(df) {
  .spi_in <- .resolve_splice_input(df, what = "raw_events")
  dt <- data.table::as.data.table(.spi_in$dt)

  required_cols <- c("gene_id","chr","strand","inc","exc","form")
  missing <- setdiff(required_cols, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse=", "))
  }

  # Fill event_type if missing
  if (!"event_type" %in% names(dt))
    dt[, event_type := "unknown"]

  # Assign event_id if missing
  if (!"event_id" %in% names(dt)) {
    dt[, event_id := paste0(gene_id, ":", event_type, ":", .I)]
  }

  # Site ID
  dt[, site_id := paste(event_type, gene_id, chr, inc, exc, form, sep="|")]

  # Default DI stats
  if (!"delta_psi" %in% names(dt)) {
    dt[, delta_psi := ifelse(form == "INC", 1,
                             ifelse(form == "EXC", -1, 1))]
    # Expand SITE events into +/- delta_psi forms
    if (any(dt$form == "SITE")) {

      site_events <- dt[form == "SITE"]

      # original copy (+1)
      site_plus  <- copy(site_events)
      site_plus[, delta_psi := 1]

      # mirrored copy (-1)
      site_minus <- copy(site_events)
      site_minus[, delta_psi := -1]

      # keep form as SITE for both, but mark direction
      site_plus[, form := "SITE"]
      site_minus[, form := "SITE"]

      # combine back
      dt <- rbind(dt[form != "SITE"], site_plus, site_minus, fill = TRUE)
    }
  }
  if (!"p.value" %in% names(dt)) dt[, p.value := 0]
  if (!"padj" %in% names(dt))   dt[, padj := 0]

  # Default diagnostic stats
  defaults <- c(
    "cooks_max","n","n_used","n_samples","n_case",
    "mean_psi_ctrl","mean_psi_case"
  )
  for (col in defaults) {
    if (!col %in% names(dt)) dt[, (col) := -1]
  }

  # Final column order (same as DI output)
  final_cols <- c(
    "site_id","event_type","event_id","gene_id","chr","strand",
    "inc","exc","n_samples","n_control","n_case",
    "mean_psi_ctrl","mean_psi_case","delta_psi","p.value",
    "padj","cooks_max","form","n","n_used"
  )

  # Fill missing n_control if not present
  if (!"n_control" %in% names(dt)) dt[, n_control := -1]

  # Keep any extras trailing
  dt <- dt[, c(final_cols, setdiff(names(dt), final_cols)), with = FALSE]
  return(dt)
}

#' @title Compare user-selected transcript pairs
#'
#' @description
#' Builds a matched-like table for pairs of transcripts by extracting
#' all coding exon coordinates from annotations.
#'
#' @param transcript_pairs data.frame with columns `transcript1`, `transcript2`
#' @param annotations flattened GTF-style data.frame or data.table (from get_annotation)
#'
#' @examples
#' annotation_df <- load_example_data("annotation_df")$annotation_df
#' pairs <- data.frame(
#'     transcript1 = c("ENST00000337907", "ENST00000426559"),
#'     transcript2 = c("ENST00000400908", "ENST00000399728")
#' )
#' matched <- compare_transcript_pairs(pairs, annotation_df$annotations)
#' print(matched)
#'
#' @return data.table mimicking `matched` structure, ready for downstream comparison.
#' @export
compare_transcript_pairs <- function(transcript_pairs, annotations) {
  ann <- build_from_annotations(annotations)
  EX  <- ann$exons
  TX  <- ann$transcripts
  matched_cols <- c(
    "event_id", "event_type", "form", "gene_id", "chr", "strand",
    "inc", "exc", "delta_psi", "p.value", "padj",
    "n_samples", "n_control", "n_case", "transcript_id", "exons"
  )
  empty_matched <- function() {
    data.table::data.table(
      event_id = character(),
      event_type = character(),
      form = character(),
      gene_id = character(),
      chr = character(),
      strand = character(),
      inc = character(),
      exc = character(),
      delta_psi = numeric(),
      p.value = numeric(),
      padj = numeric(),
      n_samples = integer(),
      n_control = integer(),
      n_case = integer(),
      transcript_id = character(),
      exons = character()
    )
  }

  if (methods::is(transcript_pairs, "SpliceImpactResult")) {
    .spi_in <- .resolve_splice_input(transcript_pairs, what = "paired_hits")
    ph <- data.table::as.data.table(.spi_in$dt)
    need_tx <- c("transcript_id_control", "transcript_id_case")
    miss_tx <- setdiff(need_tx, names(ph))
    if (length(miss_tx)) {
      stop("compare_transcript_pairs: SpliceImpactResult paired_hits missing required columns: ", paste(miss_tx, collapse = ", "))
    }
    pairs <- unique(ph[, .(transcript1 = transcript_id_control, transcript2 = transcript_id_case)])
  } else {
    pairs <- as.data.table(transcript_pairs)
  }
  stopifnot(all(c("transcript1", "transcript2") %in% names(pairs)))
  pairs[, transcript1 := trimws(as.character(transcript1))]
  pairs[, transcript2 := trimws(as.character(transcript2))]
  pairs <- unique(
    pairs[
      !is.na(transcript1) & !is.na(transcript2) &
        nzchar(transcript1) & nzchar(transcript2)
    ]
  )
  if (!nrow(pairs)) return(empty_matched())

  keepers <- pairs$transcript1 %in% TX$transcript_id &
    pairs$transcript2 %in% TX$transcript_id
  outMessage <- paste0(sum(keepers), " out of ", nrow(pairs), " transcript pairs",
         " contained within annotations")
  message(outMessage)
  if (sum(keepers) == 0) return(empty_matched())
  pairs <- pairs[keepers]

  # subset exon + tx info
  tx_info <- TX[, .(
    transcript_id = as.character(transcript_id),
    gene_id = as.character(gene_id),
    chr = as.character(chr),
    strand = as.character(strand)
  )]
  setkey(EX, transcript_id)
  setkey(tx_info, transcript_id)

  get_exon_info <- function(txid) {
    ex_sub <- EX[transcript_id == txid]
    if (!nrow(ex_sub))
      return(list(coords = "", exons = ""))
    ex_sub <- ex_sub[order(start)]
    coords_str <- paste(sprintf("%d-%d", ex_sub$start, ex_sub$stop), collapse = ";")
    exons_str  <- paste(ex_sub$exon_id, collapse = ";")
    list(coords = coords_str, exons = exons_str)
  }

  # ---- build coordinate + exon strings for both transcripts ----
  pairs[, c("inc", "inc_exons") := {
    tmp <- lapply(transcript1, get_exon_info)
    list(vapply(tmp, `[[`, character(1), "coords"),
         vapply(tmp, `[[`, character(1), "exons"))
  }]

  pairs[, c("exc", "exc_exons") := {
    tmp <- lapply(transcript2, get_exon_info)
    list(vapply(tmp, `[[`, character(1), "coords"),
         vapply(tmp, `[[`, character(1), "exons"))
  }]

  # add gene/chr/strand from both sides; prefer transcript1 metadata
  tx1 <- copy(tx_info)
  setnames(tx1, c("transcript_id", "gene_id", "chr", "strand"),
           c("transcript1", "gene_id_1", "chr_1", "strand_1"))
  tx2 <- copy(tx_info)
  setnames(tx2, c("transcript_id", "gene_id", "chr", "strand"),
           c("transcript2", "gene_id_2", "chr_2", "strand_2"))
  pairs <- merge(pairs, tx1, by = "transcript1", all.x = TRUE)
  pairs <- merge(pairs, tx2, by = "transcript2", all.x = TRUE)
  pairs[, gene_id := ifelse(!is.na(gene_id_1) & nzchar(gene_id_1), gene_id_1, gene_id_2)]
  pairs[, chr := ifelse(!is.na(chr_1) & nzchar(chr_1), chr_1, chr_2)]
  pairs[, strand := ifelse(!is.na(strand_1) & nzchar(strand_1), strand_1, strand_2)]
  pairs[, c("gene_id_1", "chr_1", "strand_1", "gene_id_2", "chr_2", "strand_2") := NULL]

  # assign synthetic IDs
  pairs[, event_id := paste0("TXCMP:", seq_len(.N))]
  pairs[, event_type := "TXCMP"]
  inc_dt <- pairs[, .(
    event_id, event_type,
    form = "INC",
    gene_id, chr, strand,
    inc, exc,
    delta_psi = 1, p.value = 0, padj = 0,
    n_samples = 0, n_control = 0, n_case = 0,
    transcript_id = transcript1,
    exons = inc_exons
  )]

  exc_dt <- pairs[, .(
    event_id, event_type,
    form = "EXC",
    gene_id, chr, strand,
    inc, exc,
    delta_psi = -1,  # mirror sign if you ever compute it
    p.value = 0, padj = 0,
    n_samples = 0, n_control = 0, n_case = 0,
    transcript_id = transcript2,
    exons = exc_exons
  )]

  matched_like <- rbindlist(list(inc_dt, exc_dt), use.names = TRUE, fill = TRUE)
  matched_like <- matched_like[, ..matched_cols]
  setorder(matched_like, event_id, form)

  # Keep only complete, non-empty event pairs.
  bad <- matched_like[
    is.na(gene_id) | !nzchar(gene_id) |
      is.na(transcript_id) | !nzchar(transcript_id),
    unique(event_id)
  ]
  if (length(bad)) {
    matched_like <- matched_like[!event_id %chin% bad]
  }
  if (!nrow(matched_like)) return(empty_matched())

  return(matched_like)
}



