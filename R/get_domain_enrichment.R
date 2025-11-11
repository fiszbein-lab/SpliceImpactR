#' Internal helper that normalizes list columns or delimited strings
#' of domain identifiers into unique trimmed character vectors.
#'
#' @param x List or character vector of domain entries.
#' @param delim Regular expression used to split delimited strings.
#'
#' @return A list of character vectors, one per input element.
#' @keywords internal
.parse_listcol <- function(x, delim = "[,;|[:space:]]+") {
  # Accept list-col of character vectors OR a single delimited string
  if (is.list(x)) {
    lapply(x, function(v) unique( trimws( as.character(v[!is.na(v) & nzchar(v)]) ) ))
  } else {
    lapply(as.character(x), function(s) {
      if (is.na(s) || !nzchar(s)) character(0)
      else unique(trimws(unlist(strsplit(s, delim))))
    })
  }
}

#' Removes the part of a domain identifier after the first semicolon.
#'
#' @param id Character vector of domain identifiers (e.g. `"Pfam;Kinase"`).
#' @return Character vector containing only the database name (e.g. `"Pfam"`).
#' @keywords internal
.db_remove_domain <- function(id) {
  sub(";.*$", "", id)
}

#' Removes the database portion from a semicolon-delimited domain ID.
#'
#' @param id Character vector of domain identifiers (e.g. `"Pfam;Kinase"`).
#' @return Character vector of domain names (e.g. `"Kinase"`).
#' @keywords internal
.domain_remove_db <- function(id) {
  sub("^[^;]*;", "", id)
}

#' @title Domain-level enrichment via hypergeometric test
#' @description
#' Tests whether particular protein domains are overrepresented among
#' inclusion/exclusion transcript pairs (foreground) relative to a
#' matched background set, using the hypergeometric test.
#'
#' @details
#' Each domain identifier is counted once per transcript pair based on
#' its presence in a list column (e.g. `either_domains_list`). The
#' probability of observing at least `k` such pairs is computed under
#' the hypergeometric distribution
#' \deqn{P(X \ge k), \quad X \sim \mathrm{Hypergeom}(M, B-M, K)}
#' where:
#' \itemize{
#'   \item \code{K} = number of foreground pairs,
#'   \item \code{B} = number of background pairs,
#'   \item \code{M} = background count of pairs containing the domain,
#'   \item \code{k} = foreground count of pairs containing the domain.
#' }
#' P-values are Benjamini-Hochberg adjusted (`padj`).
#'
#' Optionally, analyses can be restricted by event type or database
#' prefix (e.g. `"Pfam"`, `"SMART"`) and domains with fewer than
#' `min_fg_count` foreground occurrences are skipped.
#'
#' @param hits `data.frame` or `data.table` containing the foreground
#'   transcript pairs (typically the significant inclusion/exclusion
#'   events). Must include a list column of domain IDs (default:
#'   `"either_domains_list"`).
#' @param background `data.frame` or `data.table` representing the
#'   matched background pairs. Must include a list column of domain
#'   IDs (default: `"total_sd_domains"`).
#' @param domain_col_fg Name of the domain list column in `hits`.
#' @param domain_col_bg Name of the domain list column in `background`.
#' @param event_col Name of the column giving event type (default:
#'   `"event_type_inc"`). Set `NULL` to skip event filtering.
#' @param event_filter Character vector of event types to include
#'   (e.g. `"A5SS"`, `"A3SS"`). If `NULL`, all events are used.
#' @param db_filter Character vector of database prefixes to retain
#'   (e.g. `"Pfam"`, `"SMART"`). If `NULL`, all domains are used.
#' @param min_fg_count Minimum number of foreground hits required to
#'   test a domain (default `2`).
#' @param delim Regular expression describing the delimiters in string
#'   list columns (default `[,;|[:space:]]+`).
#'
#' @return A `data.table` with one row per domain, including:
#' \describe{
#'   \item{`domain_id`}{Domain identifier (database prefix removed).}
#'   \item{`db`}{Database prefix (e.g. `"Pfam"`).}
#'   \item{`k`, `K`}{Foreground domain count and total foreground pairs.}
#'   \item{`M`, `B`}{Background domain count and total background pairs.}
#'   \item{`fg_prop`, `bg_prop`}{Proportion of pairs with the domain.}
#'   \item{`OR`}{Odds ratio (Haldane-Anscombe corrected).}
#'   \item{`pval`, `padj`}{Raw and BH-adjusted p-values.}
#'   \item{`events`}{Event IDs contributing to the domain count.}
#' }
#'
#' @seealso
#' * [enrich_by_event()] run per event type
#' * [enrich_by_db()] run per database
#' * [add_domain_columns()] attach domain lists to hits
#'
#' @importFrom stats phyper start end setNames
#' @importFrom methods is
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
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#'
#' enriched_domains <- enrich_domains_hypergeo(hits_domain, bg, db_filter = 'interpro')
#'
#' @export
enrich_domains_hypergeo <- function(
    hits,
    background,
    domain_col_fg = "either_domains_list",
    domain_col_bg = "total_sd_domains",
    event_col     = "event_type_inc",      # set NULL if you don't want event filtering
    event_filter  = NULL,                  # e.g. c("A5SS","A3SS") or "AFE"
    db_filter     = NULL,                  # e.g. "Pfam" | c("Pfam","SMART")
    min_fg_count  = 2,                     # minimum foreground hits to test a domain
    delim         = "[,;|[:space:]]+"
) {
  FG <- as.data.table(hits)
  BG <- as.data.table(background)

  # Optional event filtering
  if (!is.null(event_col) && !is.null(event_filter) && event_col %in% names(FG)) {
    FG <- FG[get(event_col) %chin% as.character(event_filter)]
  }

  # Parse list columns -> list of unique domain ids per row (presence/absence)
  FG[, .fg_domains := .parse_listcol(either_domains_list, delim = delim)]
  BG[, .bg_domains := .parse_listcol(total_sd_domains, delim = delim)]

  # Optional DB filter (keep only domains whose prefix matches db_filter)
  if (!is.null(db_filter)) {
    keep_db <- function(v) {
      if (!length(v)) return(v)
      db <- .db_remove_domain(v)
      v[!is.na(db) & db %chin% as.character(db_filter)]
    }
    FG[, .fg_domains := lapply(.fg_domains, keep_db)]
    BG[, .bg_domains := lapply(.bg_domains, keep_db)]
  }

  # Foreground totals
  K <- nrow(FG)  # number of FG pairs
  # Background totals
  B <- nrow(BG)  # number of BG pairs

  # Explode to long (presence), counting once per pair
  fg_long <- FG[, .(domain_id = unlist(.fg_domains),
                    event_id = event_id), by = .I][!is.na(domain_id)]
  bg_long <- BG[, .(domain_id = unlist(.bg_domains)), by = .I]

  # Tally presence counts
  fg_counts <- fg_long[ , .(.N, events = paste(event_id, collapse="|")), by = .(domain_id)][, .(domain_id, events, k = N)]
  bg_counts <- bg_long[ , .N, by = domain_id][, .(domain_id, M = N)]

  # Merge; keep only domains that occur somewhere
  tallies <- merge(bg_counts, fg_counts, by = "domain_id", all = TRUE)
  tallies[is.na(M), M := 0L]
  tallies[is.na(k), k := 0L]

  # Drop ultra-rare in FG if requested
  tallies <- tallies[k >= as.integer(min_fg_count)]

  if (!nrow(tallies)) {
    return(data.table(
      domain_id = character(), db = character(), K = integer(), B = integer(),
      k = integer(), M = integer(), fg_prop = numeric(), bg_prop = numeric(),
      OR = numeric(), pval = numeric(), padj = numeric()
    ))
  }

  # Hypergeometric: P(X >= k) with X~Hyper(M, B-M, K)
  # phyper(q, m, n, k, lower.tail=FALSE) uses q = k-1
  tallies[, `:=`(
    K = K,
    B = B,
    fg_prop = ifelse(K > 0, k / K, NA_real_),
    bg_prop = ifelse(B > 0, M / B, NA_real_)
  )]

  # Haldane-Anscombe for OR to avoid div-by-zero
  tallies[, `:=`(
    a = k + 0.5,
    b = (K - k) + 0.5,
    c = M + 0.5,
    d = (B - M) + 0.5
  )]

  tallies[, `:=`(
    OR = (a * d) / (b * c),
    pval = stats::phyper(pmax(k - 1L, 0L), M, pmax(B - M, 0L), K, lower.tail = FALSE)
  )]
  tallies[, padj := p.adjust(pval, method = "BH")]
  # Adjust p-values (BH)


  # Annotate db prefix for convenience
  tallies[, `:=` (
    db = .db_remove_domain(domain_id),
    domain_id = .domain_remove_db(domain_id)
  )
  ]

  # Nice ordering
  data.table::setorder(tallies, padj, -OR)

  # Final columns
  tallies[, .(domain_id, db, k, K, M, B, fg_prop, bg_prop, OR, pval, padj, events)]
}

#' @title Run domain enrichment by event type
#' @description
#' Convenience wrapper for [enrich_domains_hypergeo()] that runs the
#' enrichment test separately for each event type and combines results.
#'
#' @param hits,background See [enrich_domains_hypergeo()].
#' @param events Character vector of event types to analyze.
#' @param ... Additional arguments passed to [enrich_domains_hypergeo()].
#'
#' @return Combined `data.table` with an added `event_type` column.
#' @seealso [enrich_by_db()]
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
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#'
#' enriched_domains <- enrich_by_event(hits_domain, bg, events = 'AFE', db_filter = 'interpro')
#'
#'
#' @export
enrich_by_event <- function(hits, background, events, ...) {
  rbindlist(lapply(events, function(ev) {
    out <- enrich_domains_hypergeo(hits, background, event_filter = ev, ...)
    out[, event_type := ev][]
  }), use.names = TRUE, fill = TRUE)
}

#' @title Run domain enrichment by database
#' @description
#' Convenience wrapper for [enrich_domains_hypergeo()] that runs the
#' enrichment test separately for each database (e.g. `"Pfam"`,
#' `"SMART"`) and combines the results.
#'
#' @param hits,background See [enrich_domains_hypergeo()].
#' @param dbs Character vector of database prefixes to test.
#' @param ... Additional arguments passed to [enrich_domains_hypergeo()].
#'
#' @return Combined `data.table` with an added `database` column.
#' @seealso [enrich_by_event()]
#' @examples
#'
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
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#'
#' enriched_domains <- enrich_by_db(hits_domain, bg, dbs = 'interpro')
#'
#'
#' @export
enrich_by_db <- function(hits, background, dbs, ...) {
  rbindlist(lapply(dbs, function(db) {
    out <- enrich_domains_hypergeo(hits, background, db_filter = db, ...)
    out[, database := db][]
  }), use.names = TRUE, fill = TRUE)
}


#' @title Plot enriched domains by associated event count
#'
#' @description
#' Visualizes the top enriched protein domains based on the number of
#' events contributing to each domain's enrichment, optionally coloring
#' by -log10 adjusted p-value.
#'
#' @details
#' The function expects the output of [enrich_domains_hypergeo()],
#' typically a `data.table` or `data.frame` containing `domain_id`,
#' `events`, and optionally `padj` and `OR`.
#'
#' Each bar corresponds to one domain, with height proportional to the
#' number of unique `event_id`s contributing to that domain.
#' Bars are ordered by ascending adjusted p-value (`padj`), and colored
#' by -log10(padj) if available. When no `padj` column is present,
#' the bars are shown in a uniform fill color.
#'
#' @param enriched_domains `data.frame` or `data.table`
#'   Output table from [enrich_domains_hypergeo()], including at least
#'   columns `domain_id` and `events`. Optional columns `padj` and `OR`
#'   are used for coloring and labeling.
#' @param top_n Integer (default `25`)
#'   Number of top domains to display, ranked by increasing `padj`.
#'
#' @return A `ggplot` object showing bars of domain counts colored by
#'   enrichment significance.
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
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#' enriched_domains <- enrich_domains_hypergeo(hits_domain, bg, db_filter = 'interpro')
#' plot_enriched_domains_counts(enriched_domains, top_n = 20)
#'
#' @seealso [enrich_domains_hypergeo()], [enrich_by_event()],
#'   [enrich_by_db()]
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_gradient labs
#'   theme_minimal theme element_blank element_text geom_text expand_limits
#' @export
plot_enriched_domains_counts <- function(enriched_domains,
                                         top_n      = 25) {
  DT <- as.data.table(enriched_domains)

  # count event_ids from the pipe-separated 'events' column
  DT[, n_events := vapply(strsplit(as.character(events), "\\|"),
                          function(v) sum(nzchar(v)), integer(1))]

  # convenience metrics for color/tooltip
  if ("padj" %in% names(DT)) DT[, ml10 := -log10(padj)]
  if (!("OR" %in% names(DT))) DT[, OR := NA_real_]

  # choose order and take top_n
  ORD <- order(DT$padj, decreasing = FALSE, na.last = NA)

  keep <- DT[ORD][seq_len(min(top_n, .N))]

  # factor for display order
  keep[, label := as.character(domain_id)]
  keep[, label := ifelse(nchar(label) > 50, paste0(substr(label, 1, 47), "..."), label)]
  keep[, label := factor(label, levels = rev(label))]

  # plot: bars by event count, color by -log10(padj) if available
  p <- ggplot(keep, aes(x = label, y = n_events, fill = ml10)) +
    geom_col(width = 0.7, color = "white") +
    coord_flip() +
    scale_fill_gradient(name = expression(-log[10](padj)),
                        low = "#f7f7f7", high = "deeppink4",
                        na.value = "grey80") +
    labs(x = NULL, y = "Number of event_id",
         title = "Enriched domains by associated events") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold"))

  p <- p + geom_text(aes(label = n_events),
                     hjust = -0.1, size = 3,
                     color = "black")

  # keep bars fully visible when labels are drawn
  p <- p + expand_limits(y = max(keep$n_events) * 1.1)

  return(p)
}
