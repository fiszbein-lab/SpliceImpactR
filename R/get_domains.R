#' Internal helper that constructs per-exon (`Dexon`) and per-transcript (`Dtx`)
#' mappings of domain identifiers from an exon-feature annotation table.
#'
#' @param exon_features `data.frame` containing exon–feature relationships with
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

#' Internal helper returning all domains mapped to a given transcript’s protein.
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
#' a single, stable “|”-delimited string.
#'
#' @param v Character vector of domain identifiers.
#' @return Single string joining unique sorted identifiers.
#' @keywords internal
collapse_domains <- function(v) {
  if (!length(v)) return("")
  paste(sort(unique(v)), collapse = "|")
}

#' @title Add protein domain annotations to splicing events
#' @description
#' Annotates each splicing event with protein domains that are gained,
#' lost, or uniquely present in inclusion or exclusion isoforms.
#'
#' @details
#' Internally, this function builds a domain lookup table from an exon–feature
#' annotation (e.g. InterPro, Pfam) and extracts per-exon and per-transcript
#' domain lists for each isoform in `hits`. Differences between the inclusion
#' (`*_inc`) and exclusion (`*_exc`) isoforms are then summarized as:
#'
#' * `inc_only_domains`: domains unique to the inclusion isoform
#' * `exc_only_domains`: domains unique to the exclusion isoform
#' * `diff_n`: total number of non-shared domains
#'
#' If `show_protein_domains = TRUE`, additional columns report full domain
#' sets across the entire inclusion/exclusion proteins.
#'
#' @param hits `data.frame` or `data.table` containing transcript pairs with
#'   at least `transcript_id_inc`, `transcript_id_exc`, `exons_inc`,
#'   `exons_exc`, and `event_type_inc` (or `event_type_exc`).
#' @param exon_features `data.frame` of exon–domain annotations with columns
#'   `ensembl_transcript_id`, `ensembl_peptide_id`, `exon_id`, `database`,
#'   `feature_id`, `name`, `overlap_aa_start`, `overlap_aa_end`.
#' @param show_protein_domains Logical; if `TRUE`, include full protein-level
#'   domain sets (`domains_protein_inc` / `domains_protein_exc`).
#'
#' @return
#' The input `hits` table with added columns:
#' \itemize{
#'   \item `domains_exons_inc`, `domains_exons_exc` — domains found on event exons
#'   \item `inc_only_domains`, `exc_only_domains` — domains unique to each isoform
#'   \item `inc_only_domains_list`, `exc_only_domains_list`, `either_domains_list` — list-columns
#'   \item `inc_only_n`, `exc_only_n`, `diff_n` — domain counts
#'   \item optionally, `domains_protein_inc` / `domains_protein_exc`
#' }
#'
#' @seealso
#' Internal helpers:
#' * [build_domain_lookup()] — builds per-exon / per-transcript domain tables
#' * [domains_on_exons()] — extract domains on a given exon set
#' * [domains_on_protein()] — extract domains on a full transcript
#' * [collapse_domains()] — canonicalize domain vectors
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
#' annotation_df <- get_annotation(load = "test")
#' matched <- get_matched_events_chunked(res, annotation_df$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annotation_df$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' @export
get_domains <- function(hits, exon_features, show_protein_domains = FALSE) {
  LU <- build_domain_lookup(exon_features)
  Dexon <- LU$Dexon; Dtx <- LU$Dtx

  H <- as.data.table(hits)

  res <- H[, {

    txi <- as.character(transcript_id_inc)
    txe <- as.character(transcript_id_exc)
    exi <- .parse_exon_ids(exons_inc)
    exe <- .parse_exon_ids(exons_exc)

    dpi <- domains_on_protein(Dtx,  txi)
    dpe <- domains_on_protein(Dtx,  txe)
    dei <- domains_on_exons(  Dexon, txi, exi)
    dee <- domains_on_exons(  Dexon, txe, exe)

    # event class
    et <- as.character(event_type_inc %||% event_type_exc %||% "")

    term <- et %chin% c("AFE","ALE","HFE","HLE")
    if (term) {
      inc_only <- setdiff(dei, dpe)   # event domains on INC not anywhere in EXC protein
      exc_only <- setdiff(dee, dpi)
    } else {
      inc_only <- setdiff(dei, dee)   # event domains differ between event regions
      exc_only <- setdiff(dee, dei)
    }

    c(
      if (isTRUE(show_protein_domains)) list(
        domains_protein_inc = collapse_domains(dpi),
        domains_protein_exc = collapse_domains(dpe)
      ),
      list(
        domains_exons_inc = collapse_domains(dei),
        domains_exons_exc = collapse_domains(dee),

        inc_only_domains  = collapse_domains(inc_only),
        exc_only_domains  = collapse_domains(exc_only),
        inc_only_domains_list  = list(inc_only),
        exc_only_domains_list  = list(exc_only),
        either_domains_list = list(unlist(c(inc_only, exc_only))),
        inc_only_n        = uniqueN(inc_only),
        exc_only_n        = uniqueN(exc_only),
        diff_n            = uniqueN(inc_only) + uniqueN(exc_only)
      )
    )
  }, by = .I]

  cbind(H, res[])
}


