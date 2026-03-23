#' End-to-end SpliceImpactR wrapper with selectable output class
#'
#' Runs the core SpliceImpactR pipeline from raw event table (or sample paths)
#' through paired domain/PPI calls, then returns either a compact `data.table`
#' bundle (`data`, `res`, `hits_final`) or an S4 [SpliceImpactResult].
#'
#' @param sample_frame Optional sample manifest for [get_hitindex()] / [get_rmats_hit()].
#'   Must include `path`, `condition`, and optional `sample_name`.
#' @param data Optional precomputed raw event-level table.
#' @param res Optional precomputed differential inclusion table.
#' @param annotation_df Optional annotation list from [get_annotation()] with
#'   `annotations` and `sequences`.
#' @param protein_feature_total Optional protein feature table from
#'   [get_comprehensive_annotations()]. Required if `exon_features` is not supplied
#'   and domain/PPI steps are run.
#' @param exon_features Optional precomputed exon-feature overlap table from
#'   [get_exon_features()].
#' @param ppi Optional preloaded PPI table. If `NULL`, [get_ppi_interactions()]
#'   is used when needed.
#' @param source_data Which ingestion route to use when `data` is `NULL`.
#'   One of `"hitindex"`, `"rmats"`, `"both"`, or legacy alias `"rmats_hit"`.
#' @param event_types Event types for [get_rmats_hit()] / [load_rmats()].
#'   Use both terminal and non-terminal types for `source_data = "both"`.
#' @param use Junction count mode for rMATS ingestion (`"JC"` or `"JCEC"`).
#' @param keep_annotated_first_last Passed to [get_hitindex()] for terminal events.
#' @param min_total_reads Passed to [get_differential_inclusion()].
#' @param minimum_proportion_containing_event Passed to [get_differential_inclusion()].
#' @param terminal_fill Passed to [get_differential_inclusion()].
#' @param cooks_cutoff Passed to [get_differential_inclusion()].
#' @param adjust_method Passed to [get_differential_inclusion()].
#' @param fdr_threshold Passed to [keep_sig_pairs()] as the adjusted p-value cutoff.
#' @param delta_psi_threshold Passed to [keep_sig_pairs()] as the absolute
#'   delta-psi cutoff.
#' @param parallel_glm Passed to [get_differential_inclusion()].
#' @param chunk_size_glm Passed to [get_differential_inclusion()].
#' @param BPPARAM Passed to [get_differential_inclusion()]. Use a
#'   [BiocParallel::BiocParallelParam-class] object (for example
#'   [BiocParallel::SerialParam()], [BiocParallel::SnowParam()], or
#'   [BiocParallel::MulticoreParam()]).
#' @param chunk_size_match Chunk size for [get_matched_events_chunked()].
#' @param source_pairs Pairing mode for [get_pairs()] (`"multi"` or `"paired"`).
#' @param show_protein_domains Passed to [get_domains()].
#' @param return_class One of `"data.table"` or `"S4"`.
#' @param debug_steps Logical; if `TRUE`, includes intermediates (`matched`,
#'   `hits_sequences`, `pairs`, `seq_compare`, `hits_domain`) in data.table mode.
#' @param metadata Optional list attached to `SpliceImpactResult@metadata`.
#' @param verbose Logical; emit progress messages.
#'
#' @return
#' If `return_class = "data.table"`, returns a named list with
#' `data`, `res`, and `hits_final`.
#'
#' If `return_class = "S4"`, returns a [SpliceImpactResult] containing
#' `raw_events`, `di_events`, and `paired_hits` slots.
#'
#' @examples
#' ex <- load_example_data(
#'   c("sample_frame", "annotation_df", "protein_feature_total", "ppi")
#' )
#' out <- get_splicing_impact(
#'   sample_frame = ex$sample_frame,
#'   annotation_df = ex$annotation_df,
#'   protein_feature_total = ex$protein_feature_total,
#'   ppi = ex$ppi,
#'   source_data = "rmats",
#'   event_types = c("SE"),
#'   use = "JCEC",
#'   parallel_glm = FALSE,
#'   BPPARAM = BiocParallel::SerialParam(),
#'   verbose = FALSE
#' )
#' print(names(out))
#'
#' @export
get_splicing_impact <- function(
    sample_frame = NULL,
    data = NULL,
    res = NULL,
    annotation_df = NULL,
    protein_feature_total = NULL,
    exon_features = NULL,
    ppi = NULL,
    source_data = c("rmats_hit", "hitindex", "rmats", "both"),
    event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI", "HFE", "HLE"),
    use = "JCEC",
    keep_annotated_first_last = FALSE,
    min_total_reads = 10L,
    minimum_proportion_containing_event = 0.5,
    terminal_fill = "event_max",
    cooks_cutoff = "Inf",
    adjust_method = "fdr",
    fdr_threshold = 0.05, 
    delta_psi_threshold = 0.10,
    parallel_glm = TRUE,
    chunk_size_glm = 1000L,
    BPPARAM = BiocParallel::SerialParam(),
    chunk_size_match = 2000L,
    source_pairs = c("multi", "paired"),
    show_protein_domains = FALSE,
    return_class = c("data.table", "S4"),
    debug_steps = FALSE,
    metadata = list(),
    verbose = TRUE
) {
  source_data <- match.arg(source_data)
  source_pairs <- match.arg(source_pairs)
  return_class <- match.arg(return_class)
  si <- NULL

  if (identical(source_data, "rmats_hit")) source_data <- "both"

  if (methods::is(data, "SpliceImpactResult")) {
    si <- data
    data <- as_dt_from_s4(si, "raw_events")
    if (is.null(res)) {
      res <- as_dt_from_s4(si, "di_events")
      if (!nrow(res)) res <- NULL
    }
    if (is.null(sample_frame)) {
      sf0 <- as_dt_from_s4(si, "sample_frame")
      if (nrow(sf0)) sample_frame <- sf0
    }
  }

  have_data <- !is.null(data) && nrow(data.table::as.data.table(data)) > 0L
  if (!have_data) {
    if (is.null(sample_frame)) {
      stop("get_splicing_impact: provide either `data` or `sample_frame`.")
    }
    if (!is.null(si) && !is.null(sample_frame)) {
      si <- add_splice_part(si, sample_frame = sample_frame)
    }

    if (verbose) message("[STEP] Loading raw event table")

    if (identical(source_data, "hitindex")) {
      if (!is.null(si)) {
        si <- get_hitindex(si, keep_annotated_first_last = keep_annotated_first_last)
        data <- as_dt_from_s4(si, "raw_events")
      } else {
        data <- get_hitindex(
          paths_df = sample_frame,
          keep_annotated_first_last = keep_annotated_first_last
        )
      }

    } else if (identical(source_data, "rmats")) {
      rmats_event_types <- intersect(event_types, c("MXE", "SE", "A3SS", "A5SS", "RI"))
      if (!length(rmats_event_types)) {
        stop("get_splicing_impact: `source_data='rmats'` requires at least one rMATS event type in `event_types` (MXE/SE/A3SS/A5SS/RI).")
      }
      if (!is.null(si)) {
        si <- load_rmats(si, use = use, event_types = rmats_event_types)
        si <- get_rmats(si)
        data <- as_dt_from_s4(si, "raw_events")
      } else {
        data <- get_rmats(load_rmats(sample_frame, use = use, event_types = rmats_event_types))
      }

    } else {
      if (!is.null(si)) {
        si <- get_rmats_hit(
          sample_frame = si,
          event_types = event_types,
          use = use,
          keep_annotated_first_last = keep_annotated_first_last
        )
        data <- as_dt_from_s4(si, "raw_events")
      } else {
        data <- get_rmats_hit(
          sample_frame = sample_frame,
          event_types = event_types,
          use = use,
          keep_annotated_first_last = keep_annotated_first_last
        )
      }
    }

  } else {
    data <- data.table::as.data.table(data)
  }

  have_res <- !is.null(res) && nrow(data.table::as.data.table(res)) > 0L
  if (!have_res) {
    if (verbose) message("[STEP] Differential inclusion")
    res <- get_differential_inclusion(
      DT = data,
      min_total_reads = min_total_reads,
      minimum_proportion_containing_event = minimum_proportion_containing_event,
      terminal_fill = terminal_fill,
      cooks_cutoff = cooks_cutoff,
      adjust_method = adjust_method,
      parallel_glm = parallel_glm,
      chunk_size_glm = chunk_size_glm,
      BPPARAM = BPPARAM,
      verbose = verbose
    )
  } else {
    res <- data.table::as.data.table(res)
  }
    res_di <- keep_sig_pairs(res, fdr_threshold, delta_psi_threshold)
  
  if (is.null(annotation_df) || !all(c("annotations", "sequences") %in% names(annotation_df))) {
    stop("get_splicing_impact: `annotation_df` must be provided and contain `annotations` and `sequences`.")
  }

  if (verbose) message("[STEP] Match -> sequence attach -> pairing")
  matched <- get_matched_events_chunked(
    events = res_di,
    annotations = annotation_df$annotations,
    chunk_size = chunk_size_match
  )
  hits_sequences <- attach_sequences(matched, annotation_df$sequences)
  pairs <- get_pairs(hits_sequences, source = source_pairs)

  if (verbose) message("[STEP] Sequence/frame comparison")
  seq_compare <- compare_sequence_frame(pairs, annotation_df$annotations)

  if (is.null(exon_features)) {
    if (is.null(protein_feature_total)) {
      stop("get_splicing_impact: provide `exon_features` or `protein_feature_total`.")
    }
    if (verbose) message("[STEP] Exon-feature mapping")
    exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
  }

  if (verbose) message("[STEP] Domain calls")
  hits_domain <- get_domains(
    hits = seq_compare,
    exon_features = exon_features,
    show_protein_domains = show_protein_domains
  )

  if (is.null(ppi)) {
    if (verbose) message("[STEP] Loading PPI interactions")
    ppi <- get_ppi_interactions()
  }
  if (is.null(protein_feature_total)) {
    stop("get_splicing_impact: `protein_feature_total` is required for get_ppi_switches().")
  }

  if (verbose) message("[STEP] PPI switch calls")
  hits_final <- get_ppi_switches(
    hits_domain = hits_domain,
    ppi = ppi,
    protein_feature_total = protein_feature_total
  )

  if (identical(return_class, "S4")) {
    md <- c(
      metadata,
      list(
        source_data = source_data,
        source_pairs = source_pairs,
        has_raw = TRUE,
        has_di = TRUE,
        has_res_di = TRUE,
        has_matched = TRUE,
        has_sample_frame = !is.null(sample_frame),
        has_hits = TRUE
      )
    )
    if (!is.null(si)) {
      si <- add_splice_part(si, data = data)
      si <- add_splice_part(si, res = res)
      si <- add_splice_part(si, res_di = res_di)
      si <- add_splice_part(si, matched = matched)
      if (!is.null(sample_frame)) si <- add_splice_part(si, sample_frame = sample_frame)
      si <- add_splice_part(si, hits_final = hits_final)
      si@metadata <- c(si@metadata, md)
      return(si)
    } else {
      return(as_splice_impact_result(
        data = data,
        res = res,
        res_di = res_di,
        matched = matched,
        sample_frame = sample_frame,
        hits_final = hits_final,
        metadata = md
      ))
    }
  }

  out <- list(
    data = data,
    res = res_di,
    hits_final = hits_final
  )

  if (isTRUE(debug_steps)) {
    out$matched <- matched
    out$hits_sequences <- hits_sequences
    out$pairs <- pairs
    out$seq_compare <- seq_compare
    out$hits_domain <- hits_domain
  }

  out
}
