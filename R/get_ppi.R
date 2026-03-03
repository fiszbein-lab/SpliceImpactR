#' Pull PPI from SpliceImpactR's data
#' 
#' Generation details in inst/scripts
#' @examples
#' ppi_int <- get_ppi_interactions()
#' 
#' @export
get_ppi_interactions <- function() {
  path <- system.file(
    "extdata",
    "ppi.RDS",
    package = "SpliceImpactR"
  )
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rb") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  as.data.table(readRDS(con))
}



#' Convert InterPro back to PFAM
#'
#' @param hits_domain data.table with gene_id and list-cols case_only_domains_list / control_only_domains_list
#' @param ppi wide interaction table from saved data (get_ppi)
#' @param protein_feature_total table with database/clean_name/feature_id for interpro mapping
#' @return hits_domain with added columns
#' 
#' @keywords internal
ipr_to_pfam <- function(ipr_ids) {
  ipr_ids <- unique(ipr_ids[!is.na(ipr_ids) & nzchar(ipr_ids)])
  if (!length(ipr_ids)) return(character())
  
  x  <- PFAM.db::PFAMINTERPRO2AC
  mk <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mk])
  
  pf <- unique(unlist(xx[ipr_ids], use.names = FALSE))
  pf[!is.na(pf) & nzchar(pf)]
}
#' Helper to annotate PPI changes from DDI and DMI
#' @keywords internal
mark_changing_partners_split <- function(ppi,
                                         gene_id,
                                         changed_pfam_case,
                                         changed_pfam_control,
                                         changed_motif_case = character(),
                                         changed_motif_control = character()) {
  ppi  <- as.data.table(ppi)
  sub <- ppi[geneA == gene_id | geneB == gene_id]
  
  any_in <- function(x, set) {
    if (!length(set)) return(FALSE)
    if (is.list(x)) any(unlist(x, use.names = FALSE) %chin% set)
    else any(as.character(x) %chin% set)
  }
  
  # ensure expected output cols exist even when sub is empty
  sub[, `:=`(
    partner_gene = if (nrow(sub)) fifelse(geneA == gene_id, geneB, geneA) else character(),
    DDI_changed_case = FALSE, DDI_changed_control = FALSE,
    DMI_changed_case = FALSE, DMI_changed_control = FALSE,
    interaction_changed_case = FALSE,
    interaction_changed_control = FALSE
  )]
  
  if (!nrow(sub)) return(sub)
  
  if (all(c("DDI","DDI_A","DDI_B") %in% names(sub))) {
    sub[, DDI_changed_case := DDI & (any_in(DDI_A, changed_pfam_case) | any_in(DDI_B, changed_pfam_case)), by = .I]
    sub[, DDI_changed_control := DDI & (any_in(DDI_A, changed_pfam_control) | any_in(DDI_B, changed_pfam_control)), by = .I]
  }
  
  if (all(c("DMI","DMI_A","DMI_B") %in% names(sub))) {
    # convention: DMI_A is PFAM domain; DMI_B is motif/feature id (e.g., ELM)
    sub[, DMI_changed_case := DMI & (
      any_in(DMI_A, changed_pfam_case) |
        (length(changed_motif_case) > 0L && any_in(DMI_B, changed_motif_case))
    ), by = .I]
    
    sub[, DMI_changed_control := DMI & (
      any_in(DMI_A, changed_pfam_control) |
        (length(changed_motif_control) > 0L && any_in(DMI_B, changed_motif_control))
    ), by = .I]
  }
  
  sub[, interaction_changed_case := DDI_changed_case | DMI_changed_case]
  sub[, interaction_changed_control := DDI_changed_control | DMI_changed_control]
  sub[]
}

#' Annotate hits_domain with PPI changes for inclusion vs exclusion forms
#'
#' Adds list-cols case_ppi/control_ppi (partner genes) plus counts.
#' Also returns (optionally useful) per-event token sets in PFAM + ELM forms.
#'
#' @param hits_domain data.table with gene_id and list-cols case_only_domains_list / control_only_domains_list
#' @param ppi wide interaction table from saved data (get_ppi)
#' @param protein_feature_total table with database/clean_name/feature_id for interpro mapping
#' @return hits_domain with added columns
#' @return A `data.table` identical to `hits_domain` with added columns:
#' \describe{
#'   \item{`case_ppi`, `control_ppi`}{Lists of partner transcripts unique to
#'     inclusion or exclusion isoforms.}
#'   \item{`n_control_ppi`, `n_control_ppi`}{Counts of gained/lost interactions.}
#'   \item{`n_ppi`}{Total PPI changes (sum of both directions).}
#' }
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
#'
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
#' ppi <- get_ppi_interactions()             
#' hits_ppi <- get_ppi_switches(hits_domain, ppi, protein_feature_total)
#' @examples
#' 
#' hits_domain[n_ppi > 0, .(event_id, gene_id, n_case_ppi, n_control_ppi, n_ppi, case_ppi, control_ppi)]
#' 
#' @export
get_ppi_switches <- function(hits_domain, ppi, protein_feature_total) {
  hd <- as.data.table(hits_domain)
  
  ipr_map <- unique(as.data.table(protein_feature_total)[
    database == "interpro",
    .(clean_name, ipr = feature_id)
  ])
  ipr_map <- ipr_map[!is.na(clean_name) & nzchar(clean_name) &
                       !is.na(ipr) & nzchar(ipr)]
  
  # Core parser you provided, generalized to accept a *single list-cell* input.
  # Returns a list with:
  # - pfam: PFAM IDs (including those derived from InterPro names)
  # - elm : ELM IDs (raw ELM tokens)
  parse_tokens <- function(names_vec) {
    names_vec <- unlist(names_vec, use.names = FALSE)
    names_vec <- names_vec[!is.na(names_vec) & nzchar(names_vec)]
    if (!length(names_vec)) return(list(pfam = character(), elm = character()))
    
    sp <- tstrsplit(names_vec, ";", fixed = TRUE)
    src <- sp[[1]]
    val <- sp[[2]]
    
    preconvert_pfam <- unique(val[src == "pfam"])
    preconvert_elm  <- unique(val[src == "elm"])
    
    ip_names <- unique(val[src == "interpro"])
    if (length(ip_names)) {
      ipr <- unique(ipr_map[list(ip_names), on = .(clean_name)][, ipr])
      ipr <- ipr[!is.na(ipr) & nzchar(ipr)]
      pf_from_ip <- ipr_to_pfam(ipr)
      pfam <- unique(c(preconvert_pfam, pf_from_ip))
    } else {
      pfam <- preconvert_pfam
    }
    
    list(
      pfam = pfam[!is.na(pfam) & nzchar(pfam)],
      elm  = preconvert_elm[!is.na(preconvert_elm) & nzchar(preconvert_elm)]
    )
  }
  
  case_ppi <- vector("list", nrow(hd))
  control_ppi <- vector("list", nrow(hd))
  n_case   <- integer(nrow(hd))
  n_control   <- integer(nrow(hd))
  n_all   <- integer(nrow(hd))
  
  # optional: keep token sets (useful for debugging / downstream)
  case_pfam <- vector("list", nrow(hd))
  control_pfam <- vector("list", nrow(hd))
  case_elm  <- vector("list", nrow(hd))
  control_elm  <- vector("list", nrow(hd))
  
  for (i in seq_len(nrow(hd))) {
    gene_id <- as.character(hd$gene_id[i])
    if (is.na(gene_id) || !nzchar(gene_id)) next
    
    tok_case <- parse_tokens(hd$case_only_domains_list[i])
    tok_control <- parse_tokens(hd$control_only_domains_list[i])
    
    case_pfam[[i]] <- tok_case$pfam
    control_pfam[[i]] <- tok_control$pfam
    case_elm[[i]]  <- tok_case$elm
    control_elm[[i]]  <- tok_control$elm
    
    edges <- mark_changing_partners_split(
      ppi = ppi,
      gene_id = gene_id,
      changed_pfam_case = tok_case$pfam,
      changed_pfam_control = tok_control$pfam,
      changed_motif_case = tok_case$elm,
      changed_motif_control = tok_control$elm
    )
    
    case_genes <- unique(edges[interaction_changed_case == TRUE, partner_gene])
    control_genes <- unique(edges[interaction_changed_control == TRUE, partner_gene])
    
    case_ppi[[i]] <- case_genes
    control_ppi[[i]] <- control_genes
    n_case[i]     <- length(case_genes)
    n_control[i]     <- length(control_genes)
    n_all[i]     <- length(unique(c(case_genes, control_genes)))
  }
  
  hd[, `:=`(
    case_ppi   = case_ppi,
    control_ppi   = control_ppi,
    n_case_ppi = n_case,
    n_control_ppi = n_control,
    n_ppi     = n_all,
    case_pfam_changed = case_pfam,
    control_pfam_changed = control_pfam,
    case_elm_changed  = case_elm,
    control_elm_changed  = control_elm
  )]
  
  hd[]
}



#' @title Plot summary of altered PPI interactions
#' @description
#' Visualizes the frequency and magnitude of gained/lost PPIs per event,
#' using a dual-panel layout:
#' - left: proportion of events with any PPI change
#' - right: histograms of CASE and CONTROL partner counts (non-zero only)
#'
#' @param df `data.table` or `data.frame` with PPI counts per event,
#'   as returned by [ppi_switches_for_hits()].
#' @param bins Integer; number of histogram bins (default `30`).
#' @param palette Named character vector of fill colors for the plot
#'   (default includes `"no"`, `"yes"`, `"CASE"`, `"CONTROL"`).
#' @param output_file Optional path to save the figure (`.png` or `.pdf`).
#' @param width,height Numeric dimensions (in inches) for saved plot.
#'
#' @return A `ggplot` object combining two panels (using `patchwork`).
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
#'
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
#' ppi <- get_ppi_interactions()             
#' hits_final <- get_ppi_switches(hits_domain, ppi, protein_feature_total)
#' ppi_plot <- plot_ppi_summary(hits_final)
#'
#' @seealso [ppi_switches_for_hits()]
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_histogram facet_wrap
#'   scale_fill_manual scale_x_discrete labs theme_classic theme_bw theme element_blank
#'   element_text expand_limits
#' @importFrom patchwork plot_layout
#' @importFrom scales percent
#' @importFrom ggplot2 margin
#' @export
plot_ppi_summary <- function(df,
                             bins = 30,
                             palette = c("no" = "grey80", "yes" = "deeppink4",
                                         "CASE" = "#2b8cbe", "CONTROL" = "#e34a33"),
                             output_file = NULL, width = 9, height = 4.8) {
  
  
  DT <- as.data.table(df)
  
  # ----- Left panel: binary any-ppi -----
  left_dt <- DT[, .(has_ppi = ifelse(n_ppi > 0, "yes", "no"))][, .N, by = has_ppi]
  left_dt[, frac := N / sum(N)]
  p_left <- ggplot(left_dt, aes(x = has_ppi, y = N, fill = has_ppi)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(N, " (", scales::percent(frac, accuracy = 1), ")")),
              vjust = -0.25, size = 3.6) +
    scale_fill_manual(values = palette[c("no","yes")], guide = "none") +
    scale_x_discrete(labels = c(no = "No changed PPIs", yes = "Changed PPIs")) +
    labs(x = NULL, y = "Events") +
    theme_classic(base_size = 11) +
    theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
  
  # ----- Right panel: histograms of non-zero INC/EXC PPI counts -----
  long_dt <- rbind(
    DT[, .(type = "CASE", value = as.integer(n_case_ppi))],
    DT[, .(type = "CONTROL", value = as.integer(n_control_ppi))]
  )[value > 0]  # drop zeros as requested
  
  p_right <- ggplot(long_dt, aes(x = value, fill = type)) +
    geom_histogram(bins = bins, color = "white", linewidth = 0.2, show.legend = FALSE) +
    facet_wrap(~type, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = palette[c("CASE","CONTROL")]) +
    labs(x = "PPI partners (non-zero)", y = "Count") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 10))
  
  # ----- Assemble -----
  plt <- p_left + p_right + plot_layout(widths = c(1, 2))
  
  if (!is.null(output_file)) {
    ggsave(output_file, plt, width = width, height = height, dpi = 300)
  }
  plt
}
