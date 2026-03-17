#' Internal helper that constructs per-exon (`Dexon`) and per-transcript (`Dtx`)
#' mappings of domain identifiers from an exon-feature annotation table.
#'
#' @param exon_features `data.frame` containing exon-feature relationships with
#'   amino acid start/end positions.
#'
#' @return A named `list` with:
#' \describe{
#'   \item{`Dexon`}{`data.table` keyed by `transcript_id, exon_id`}
#'   \item{`Dtx`}{`data.table` keyed by `transcript_id`}
#' }
#'
#' @keywords internal
build_domain_lookup <- function(exon_features) {
  D <- as.data.table(exon_features)[
    , .(transcript_id = as.character(ensembl_transcript_id),
        peptide_id    = as.character(ensembl_peptide_id),
        exon_id       = as.character(exon_id),
        database      = as.character(database),
        feature_id    = as.character(feature_id),
        name          = gsub("[|]", " ", gsub(';', " ", as.character(name))),
        aa_start      = as.integer(overlap_aa_start),
        aa_end        = as.integer(overlap_aa_end))
  ]
  D <- D[is.finite(aa_start) & is.finite(aa_end) & aa_start <= aa_end]

  # Instance key keeps repeats distinct
  D[, dom_key := sprintf("%s;%s", database, name)]

  # per-exon lookup
  Dexon <- D[, .(doms = list(dom_key)), by = .(transcript_id, exon_id)]
  setkey(Dexon, transcript_id, exon_id)

  # per-transcript (protein) lookup: keep all instances
  Dtx <- D[, .(doms = list(dom_key)), by = .(transcript_id)]
  setkey(Dtx, transcript_id)

  list(Dexon = Dexon, Dtx = Dtx)
}


#' Internal helper that returns all domain identifiers overlapping a given
#' set of exons within a transcript.
#'
#' @param Dexon Output of [build_domain_lookup()] (`Dexon` element).
#' @param tx Character scalar; transcript ID.
#' @param exons_vec Character vector of exon IDs.
#'
#' @return Character vector of unique domain identifiers for those exons.
#' @keywords internal
domains_on_exons <- function(Dexon, tx, exons_vec) {
  if (!length(exons_vec)) return(character(0))
  got <- Dexon[list(tx)][exon_id %chin% exons_vec]
  if (!nrow(got)) return(character(0))
  unique(unlist(got$doms, use.names = FALSE))
}

#' Internal helper returning all domains mapped to a given transcript's protein.
#'
#' @param Dtx Output of [build_domain_lookup()] (`Dtx` element).
#' @param tx Character scalar; transcript ID.
#'
#' @return Character vector of domain identifiers.
#' @keywords internal
domains_on_protein <- function(Dtx, tx) {
  row <- Dtx[list(tx)]
  if (!nrow(row)) return(character(0))
  unique(unlist(row$doms, use.names = FALSE))
}

#' Internal helper that collapses a character vector of domain identifiers into
#' a single, stable | delimited string.
#'
#' @param v Character vector of domain identifiers.
#' @return Single string joining unique sorted identifiers.
#' @keywords internal
collapse_domains <- function(v) {
  if (!length(v)) return("")
  paste(sort(unique(v)), collapse = "|")
}


#' @description helper fxn to prep domain coords for domain name + overlap check
#' @return data.table with domains, chr, start, end
#' @keywords internal
#' @noRd
.parse_domain_coords <- function(domains) {
  if (length(domains) == 0L)
    return(data.table(domain = character(), chr = character(),
                      start = integer(), end = integer()))

  # parts <- sub("^(.*) chr([^ ]+:[0-9]+-[0-9]+)$", "\\1\tchr\\2", domains)
  parts <- sub("^(.*) chr([^ ]+)$", "\\1\tchr\\2", domains)
  parts <- tstrsplit(parts, "\t", fixed = TRUE)
  parts[[2]] <- paste0("chr", gsub("chr", "", parts[[2]]))
  coords <- tstrsplit(parts[[2]], "[:-]")
  chr    <- coords[[1]]
  start  <- suppressWarnings(as.integer(coords[[2]]))
  end    <- suppressWarnings(as.integer(coords[[3]]))

  start <- ifelse(is.na(start), end, start)
  end   <- ifelse(is.na(end),   start, end)
  both_na <- is.na(start) & is.na(end)
  start[both_na] <- 0L
  end[both_na]   <- 0L
  data.table(
    domain = trimws(parts[[1]]),
    chr    = trimws(chr),
    start  = start,
    end    = end
  )[!is.na(chr) & !is.na(start) & !is.na(end)]
}


#' @description helper fxn to use foverlaps to get domain differences
#' @return data.table with domains, chr, start, end
#' @keywords internal
#' @noRd
.diff_domains <- function(a, b) {
  if (nrow(a) == 0L || nrow(b) == 0L) return(a)

  a <- a[!is.na(start) & !is.na(end)]
  b <- b[!is.na(start) & !is.na(end)]
  if (nrow(a) == 0L || nrow(b) == 0L) return(a)

  a[, c("start", "end") := .(pmin(start, end), pmax(start, end))]
  b[, c("start", "end") := .(pmin(start, end), pmax(start, end))]

  setkey(a, domain, chr, start, end)
  setkey(b, domain, chr, start, end)

  ol <- foverlaps(
    a, b,
    by.x = c("domain","chr","start","end"),
    by.y = c("domain","chr","start","end"),
    nomatch = 0
  )

  # Get non-overlapping entries from a
  out <- a[!unique(ol, by = c("domain","chr","i.start","i.end"))
           [, .(domain, chr, start = i.start, end = i.end)]]

  # Return as data.table (even if empty)
  out[]
}



#' @title Add protein domain annotations to splicing events
#' @description
#' Annotates each splicing event with protein domains that are gained,
#' lost, or uniquely present in inclusion or exclusion isoforms.
#'
#' @details
#' Internally, this function builds a domain lookup table from an exon feature
#' annotation (e.g. InterPro, Pfam) and extracts per-exon and per-transcript
#' domain lists for each isoform in `hits`. Differences between the inclusion
#' (`*_case`) and exclusion (`*_control`) isoforms are then summarized as:
#'
#' * `case_only_domains`: domains unique to the inclusion isoform
#' * `control_only_domains`: domains unique to the exclusion isoform
#' * `diff_n`: total number of non-shared domains
#'
#' If `show_protein_domains = TRUE`, additional columns report full domain
#' sets across the entire inclusion/exclusion proteins.
#'
#' @param hits `data.frame`, `data.table`, or `SpliceImpactResult` containing
#'   transcript pairs with at least `transcript_id_case`,
#'   `transcript_id_control`, `exons_case`, `exons_control`, and `event_type`.
#' @param exon_features `data.frame` of exon-domain annotations with columns
#'   `ensembl_transcript_id`, `ensembl_peptide_id`, `exon_id`, `database`,
#'   `feature_id`, `name`, `overlap_aa_start`, `overlap_aa_end`.
#' @param show_protein_domains Logical; if `TRUE`, include full protein-level
#'   domain sets (`domains_protein_case` / `domains_protein_control`).
#' @param return_class Character. Output mode: `"data.table"`, `"S4"`, or
#'   `"auto"` (default). In `auto`, S4 input returns updated S4 output.
#'
#' @return
#' The input `hits` table with added columns (or updated
#' `SpliceImpactResult` when `return_class` resolves to S4):
#' \itemize{
#'   \item `domains_exons_case`, `domains_exons_control` domains found on event exons
#'   \item `case_only_domains`, `control_only_domains` domains unique to each isoform
#'   \item `case_only_domains_list`, `control_only_domains_list`, `either_domains_list` list-columns
#'   \item `case_only_n`, `control_only_n`, `diff_n` domain counts
#'   \item optionally, `domains_protein_case` / `domains_protein_control`
#' }
#'
#' @examples
#' ex <- load_example_data("sample_frame")
#' sample_frame <- ex$sample_frame
#' hit_index <- get_hitindex(sample_frame)
#' res <- get_differential_inclusion(hit_index)
#' annotation_df <- load_example_data("annotation_df")$annotation_df
#' matched <- get_matched_events_chunked(res, annotation_df$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annotation_df$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
#' interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#' print(hits_domain)
#'
#' @export
get_domains <- function(hits, exon_features, show_protein_domains = FALSE, return_class = c("auto", "data.table", "S4")) {
  return_class <- match.arg(return_class)
  .spi_in <- .resolve_splice_input(hits, what = "paired_hits")
  .spi_obj <- .spi_in$obj
  LU <- build_domain_lookup(exon_features)
  Dexon <- LU$Dexon; Dtx <- LU$Dtx

  H <- as.data.table(.spi_in$dt)

  res <- H[, {

    txi <- as.character(transcript_id_case)
    txe <- as.character(transcript_id_control)
    exi <- .parse_exon_ids(exons_case)
    exe <- .parse_exon_ids(exons_control)

    dpi <- domains_on_protein(Dtx,  txi)
    dpe <- domains_on_protein(Dtx,  txe)
    dei <- domains_on_exons(  Dexon, txi, exi)
    dee <- domains_on_exons(  Dexon, txe, exe)

    # event class
    et <- as.character(event_type %||% event_type %||% "")

    term <- et %chin% c("AFE","ALE","HFE","HLE")
    if (term) {
      dei_dt <- .parse_domain_coords(dei)
      dpe_dt <- .parse_domain_coords(dpe)
      dee_dt <- .parse_domain_coords(dee)
      dpi_dt <- .parse_domain_coords(dpi)

      case_only_dt <- .diff_domains(dei_dt, dpe_dt)
      control_only_dt <- .diff_domains(dee_dt, dpi_dt)

      case_only <- if (nrow(case_only_dt) == 0) character(0) else case_only_dt$domain
      control_only <- if (nrow(control_only_dt) == 0) character(0) else control_only_dt$domain

    } else {
      dei_dt <- .parse_domain_coords(dei)  # inc event exons
      dee_dt <- .parse_domain_coords(dee)  # exc event exons
      dpi_dt <- .parse_domain_coords(dpi)  # inc full protein
      dpe_dt <- .parse_domain_coords(dpe)  # exc full protein
      
      case_only_dt <- .diff_domains(dei_dt, dpe_dt)  # inc event vs exc protein
      control_only_dt <- .diff_domains(dee_dt, dpi_dt)
      
      # dei_dt <- .parse_domain_coords(dei)
      # dee_dt <- .parse_domain_coords(dee)
      # 
      # inc_only_dt <- .diff_domains(dei_dt, dee_dt)
      # exc_only_dt <- .diff_domains(dee_dt, dei_dt)

      case_only <- if (nrow(case_only_dt) == 0) character(0) else case_only_dt$domain
      control_only <- if (nrow(control_only_dt) == 0) character(0) else control_only_dt$domain
    }

    dpi <- sub(" chr[0-9]+.*$", "", dpi)
    dpe <- sub(" chr[0-9]+.*$", "", dpe)
    dee <- sub(" chr[0-9]+.*$", "", dee)
    dei <- sub(" chr[0-9]+.*$", "", dei)

    c(
      if (isTRUE(show_protein_domains)) list(
        domains_protein_case = collapse_domains(dpi),
        domains_protein_control = collapse_domains(dpe)
      ),
      list(
        domains_exons_case = collapse_domains(dei),
        domains_exons_control = collapse_domains(dee),

        case_only_domains  = collapse_domains(case_only),
        control_only_domains  = collapse_domains(control_only),
        case_only_domains_list  = list(case_only),
        control_only_domains_list  = list(control_only),
        either_domains_list = list(unlist(c(case_only, control_only))),
        case_only_n        = uniqueN(case_only[case_only != ""]),
        control_only_n        = uniqueN(control_only[control_only != ""]),
        diff_n            = uniqueN(case_only[case_only != ""]) +
          uniqueN(control_only[control_only != ""])
      )
    )
  }, by = .I]

  out <- cbind(H, res[])
  if ("I" %in% names(out)) out[, I := NULL]
  .return_splice_output(out, obj = .spi_obj, what = "paired_hits", return_class = return_class)
}
