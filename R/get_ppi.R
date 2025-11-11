#' Map Pfam accessions to InterPro identifiers
#'
#' Internal helper that joins Pfam accessions to InterPro IDs using
#' mappings from the \pkg{PFAM.db} annotation package.
#'
#' @param df A \code{data.table} or \code{data.frame} containing at
#'   least two columns \code{D1} and \code{D2} with Pfam accessions.
#'
#' @return A \code{data.table} with additional columns
#'   \code{I1} and \code{I2} containing corresponding InterPro IDs.
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom PFAM.db PFAMINTERPRO
pfam_to_ip <- function(df) {
  pf2ipr <- data.table::setDT(as.data.frame(PFAM.db::PFAMINTERPRO))[, .(IPR = interpro), by = .(PF = ac)]
  df[, I1 := pf2ipr[.SD, on = .(PF = D1), IPR]]
  df[, I2 := pf2ipr[.SD, on = .(PF = D2), IPR]]
  return(df[!is.na(I1) & !is.na(I2)])
}

#' Download and process the PPIDM dataset (Pfam-Pfam domain interactions)
#'
#' Retrieves the 2021 PPIDM dataset from Zenodo, extracts only
#' \code{OutputData/} contents, merges them into a single table,
#' and maps Pfam accessions to InterPro IDs.
#'
#' @param download Logical; if TRUE, the PPIDM zip archive is
#'   automatically fetched via \pkg{BiocFileCache}. If FALSE, assumes
#'   the processed file already exists in \code{load_dir}.
#' @param load_dir Directory where files will be downloaded or read
#'   from. The directory will be created if missing. If \code{NULL}, creates temp
#'   dir.
#' @param test Logical; if TRUE, loads the minimal test set. If FALSE,
#'   normal behavior
#'
#' @return A \code{data.table} containing Pfam-InterPro domain pairs.
#' @export
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcdownload bfcrpath
#' @importFrom utils unzip
#' @importFrom data.table fread fwrite setDT rbindlist
#' @examples
#' ppidm <- get_ppidm(test=TRUE)
#' # ppidm <- get_ppidm(download = TRUE)
#'
get_ppidm <- function(
    download = TRUE,
    load_dir = NULL,
    test = FALSE
) {
  if (test == TRUE) {
    return(fread(get_example_data("test_ppidm.csv")))
  }
  if (download == TRUE) {
    if (is.null(load_dir)) {
      load_dir <- file.path(tempdir(), "SpliceImpactR_PPIDM")
      message("[INFO] Using temporary directory: ", load_dir)
    }
    dir.create(load_dir, recursive = TRUE, showWarnings = FALSE)

    bfc  <- BiocFileCache::BiocFileCache(cache = load_dir, ask = FALSE)
    url  <- "https://zenodo.org/records/4880347/files/PPIDMdata_2021.zip?download=1"
    name <- "PPIDMdata_2021.zip"

    bfc <- BiocFileCache::BiocFileCache(cache = load_dir, ask = FALSE)

    q <- BiocFileCache::bfcquery(bfc, name, field="rname")

    if (nrow(q)) {
      rid <- q$rid[1]   # <-- THIS, not q[[1]]
    } else {
      rid <- BiocFileCache::bfcadd(bfc, rname=name, fpath=url)
    }

    BiocFileCache::bfcdownload(bfc, rid, ask=FALSE)
    path <- BiocFileCache::bfcrpath(bfc, rids=rid)

    zip_path <- BiocFileCache::bfcrpath(bfc, rids = rid)

    z <- utils::unzip(zip_path, list = TRUE)
    in_outdir <- grepl("(^|/)(OutputData|outputFiles)(/|$)", z$Name, perl = TRUE, ignore.case = TRUE)
    z_keep <- z$Name[in_outdir]

    out_dir <- file.path(load_dir, "PPIDM_OutputData")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(zip_path, files = z_keep, exdir = out_dir)

    ## 4) Build absolute paths & keep only files (not directories)
    paths <- file.path(out_dir, z_keep)
    paths <- paths[file.exists(paths) & !dir.exists(paths)]
    ppidm <- do.call(rbind, lapply(paths, fread))
    fwrite(ppidm, paste0(load_dir, "ppidm.csv"))
    return(pfam_to_ip(data.table::setDT(do.call(rbind, lapply(paths, fread)))))
  } else {
    csv_path <- file.path(load_dir, "ppidm.csv")
    if (!file.exists(csv_path)) {
      stop(sprintf("File '%s' does not exist or is non-readable", csv_path))
    }
    return(pfam_to_ip(fread(csv_path)))
  }
}



#' Identify Interacting Isoform Pairs Using InterPro-Derived Domain Interactions
#'
#' Maps pairs of transcript isoforms (and their corresponding proteins)
#' whose InterPro domains are known to interact according to PPIDM.
#'
#' @param protein_features A data.frame or data.table containing
#'   InterPro domain annotations from \code{get_protein_features()}.
#' @param interpro_pairs A data.table containing interacting InterPro pairs
#'   (e.g., from \code{get_ppidm()}).
#' @param init Logical; if TRUE, build the interaction table fresh. If FALSE,
#'   load a previously saved version from \code{load_dir}.
#' @param load_dir Directory for saving/loading intermediate results.
#'   If \code{NULL}, a temporary directory will be used.
#' @param save Logical; if TRUE, save results to
#'   \code{file.path(load_dir, "protein_protein_interactions_ppidm.csv")}.
#'
#' @return A data.table with columns:
#'   \code{txA, txB, IA, IB, pepA, pepB}.
#'
#' @examples
#' ppidm <- get_ppidm(test=TRUE)
#' annotation_df <- get_annotation(load = "test")
#' interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' ppi <- get_isoform_interactions(protein_feature_total[ensembl_transcript_id %in%
#'                                                         annotation_df$annotations[, unique(transcript_id)]],
#'                                 ppidm, save = FALSE, load_dir = '/projectnb2/evolution/zwakefield/proteinImpacts/', init = TRUE)
#'
#'
#' @importFrom data.table as.data.table fwrite fread setDT
#' @importFrom dplyr distinct
#' @export
get_isoform_interactions <- function(protein_features, interpro_pairs, init = TRUE, load_dir = NULL, save = FALSE) {

  if (init) {
    L <- data.table::as.data.table(protein_features)
    IPI <- data.table::as.data.table(interpro_pairs)[,.(I1, I2)] %>% dplyr::distinct()

    I1 <- as.character(IPI[[1]]); I1[is.na(I1)] <- "<NA>"
    I2 <- as.character(IPI[[2]]); I2[is.na(I2)] <- "<NA>"
    lo <- pmin(I1, I2)
    hi <- pmax(I1, I2)
    key <- paste0(lo, "||", hi)
    is_unique_row <- !duplicated(key) & !duplicated(key, fromLast = TRUE)
    IPI <- IPI[is_unique_row]

    L <- L[tolower(database) == "interpro"]

    L <- L[, .(transcript_id = ensembl_transcript_id,
               IPR = feature_id)] %>% dplyr::distinct()


    # ---- join: I1 → tx1, then I2 → tx2 (cartesian allowed) -------------------
    # I1 with transcripts
    j1 <- IPI[L, on = .(I1 = IPR), allow.cartesian = TRUE, nomatch = 0L] %>% distinct()
    j2 <- L[j1,  on = .(IPR = I2), allow.cartesian = TRUE, nomatch = 0L] %>% distinct()

    lo <- pmin(j2$transcript_id, j2$i.transcript_id)
    hi <- pmax(j2$transcript_id, j2$i.transcript_id)
    key <- paste0(lo, "||", hi)
    is_unique_row <- !duplicated(key) & !duplicated(key, fromLast = TRUE)
    j2 <- j2[is_unique_row]
    pairs <- j2[, .(txA = transcript_id,
                    txB = i.transcript_id,
                    IA = I1,
                    IB = IPR)]
    pf <- protein_features[!is.na(ensembl_transcript_id) & !is.na(ensembl_peptide_id), .(ensembl_peptide_id = first(ensembl_peptide_id)), by = ensembl_transcript_id]

    pairs <- pf[pairs, on = .(ensembl_transcript_id = txA)][, `:=`(pepA = ensembl_peptide_id)][, ensembl_peptide_id := NULL]

    pairs <- pf[pairs, on = .(ensembl_transcript_id = txB)][, `:=`(pepB = ensembl_peptide_id)][, ensembl_peptide_id := NULL]

    if (save) {
      data.table::fwrite(pairs, paste0(load_dir, "/protein_protein_interactions_ppidm.csv"))
    }
    return(pairs)
  } else {
    pairs <- data.table::fread(paste0(load_dir, "/protein_protein_interactions_ppidm.csv"))
    return(data.table::setDT(pairs))
  }

}



#' @title Build transcript-level PPI interaction map
#' @description
#' Converts a protein-protein interaction (PPI) table into a symmetric
#' transcript-transcript edge list, including simple evidence summaries.
#'
#'@param ppi `data.frame` or `data.table` containing at least
#'   `ensembl_transcript_id` and `i.ensembl_transcript_id` columns.
#'
#' @return A `data.table` with one row per transcript-partner pair,
#' including:
#' \describe{
#'   \item{`tx`, `partner_tx`}{Transcript pair identifiers.}
#'   \item{`evidence_n`}{Number of evidence sources supporting the edge.}
#'   \item{`evidence_id`}{Concatenated evidence accession(s).}
#' }
#'
#' @keywords internal
build_ppi_maps <- function(ppi) {
  P <- as.data.table(ppi)

  # normalize column names (robust to accidental factors)
  setnames(P,
           old = c("ensembl_transcript_id", "i.ensembl_transcript_id"),
           new = c("txA",                   "txB"),
           skip_absent = TRUE)

  # keep only transcript–transcript rows
  P <- P[ nzchar(txA) & nzchar(txB) & !is.na(txA) & !is.na(txB) ]

  # simple evidence summary (count non-NA fields among IA/IB/pepA/pepB)
  ev_cols <- intersect(c("IA","IB","pepA","pepB"), names(P))
  if (length(ev_cols)) {
    P[, evidence_n := rowSums(!is.na(.SD) & .SD != "", na.rm=TRUE), .SDcols = ev_cols]
    P[, evidence_id := do.call(paste, c(.SD, sep="|")), .SDcols = intersect(c("IA","IB"), names(P))]
  } else {
    P[, `:=`(evidence_n = NA_integer_, evidence_id = NA_character_)]
  }

  # symmetric edges: A->B and B->A
  S1 <- P[, .(tx = txA, partner_tx = txB, evidence_n, evidence_id)]
  S2 <- P[, .(tx = txB, partner_tx = txA, evidence_n, evidence_id)]
  S  <- rbindlist(list(S1, S2), use.names = TRUE)

  # de-dup (keep strongest evidence)
  setorder(S, tx, partner_tx, -evidence_n)
  S <- S[!duplicated(S[, .(tx, partner_tx)])]

  S[]
}


#' @title Attach gene IDs to transcript-partner pairs
#' @description
#' Annotates a transcript-transcript PPI map with corresponding gene IDs.
#'
#' @param x `data.table` or `data.frame` with columns `tx` and `partner_tx`.
#' @param annotations from get_annotations
#'
#' @return The same table with added columns:
#'   `tx_gene`, `partner_gene`, and original columns retained.
#' @keywords internal
attach_gene_to_tx <- function(x, annotations) {
  A <- as.data.table(annotations)[, .(transcript_id = as.character(transcript_id),
                                      gene_id      = as.character(gene_id))]
  A <- A[!is.na(transcript_id) & nzchar(transcript_id)]
  setkey(A, transcript_id)
  x  <- as.data.table(x)
  # left-join for tx and partner
  setnames(x, c("tx","partner_tx"), c(".__tx",".__partner"))
  x <- A[x, on = .(transcript_id = .__tx)][
    A, on = .(transcript_id = .__partner),
    `:=`(partner_tx = i.transcript_id, tx = transcript_id,
         tx_gene = gene_id, partner_gene = i.gene_id)][
           , c("transcript_id","i.transcript_id","gene_id","i.gene_id") := NULL]
  setcolorder(x, c("tx","tx_gene","partner_tx","partner_gene",
                   setdiff(names(x), c("tx","tx_gene","partner_tx","partner_gene"))))
  x[]
}


#' @title Summarize PPI partners per transcript
#' @description
#' Collapses a transcript-partner PPI map into one row per transcript,
#' listing all unique partners and evidence.
#'
#' @param ppi_map_tx_gene `data.table` output of
#'   [attach_gene_to_tx()] or similar.
#'
#' @return A `data.table` with one row per transcript containing:
#' \describe{
#'   \item{`partners_tx`}{List of partner transcript IDs.}
#'   \item{`partners_gene`}{List of partner gene IDs.}
#'   \item{`partners_tbl`}{Compact table of all partners and evidence.}
#' }
#' @keywords internal
summarize_partners <- function(ppi_map_tx_gene) {
  M <- as.data.table(ppi_map_tx_gene)

  partners <- M[, .(
    partners_tx   = list(unique(partner_tx)),
    partners_gene = list(unique(partner_gene)),
    # keep a compact evidence table per tx -> partner
    partners_tbl  = list(.SD[, .(partner_tx, partner_gene, evidence_n, evidence_id)])
  ), by = .(tx, tx_gene)]

  setkey(partners, tx)
  partners[]
}

#' @title Compute PPI gains and losses for inclusion/exclusion transcripts
#' @description
#' Identifies gained and lost protein protein interactions (PPIs)
#' between inclusion (INC) and exclusion (EXC) isoforms for each
#' splicing event.
#'
#' @param hits_all `data.frame` or `data.table` of event-level hits,
#'   including `transcript_id_inc`, `transcript_id_exc`, and `pc_class`.
#' @param ppi PPI map (e.g. IntAct or STRING) containing columns
#'   `ensembl_transcript_id` and `i.ensembl_transcript_id`.
#'
#' @return A `data.table` identical to `hits_all` with added columns:
#' \describe{
#'   \item{`inc_ppi`, `exc_ppi`}{Lists of partner transcripts unique to
#'     inclusion or exclusion isoforms.}
#'   \item{`n_inc_ppi`, `n_exc_ppi`}{Counts of gained/lost interactions.}
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
#' ppidm <- get_ppidm(test=TRUE)
#'
#' ppi <- get_isoform_interactions(protein_feature_total[ensembl_transcript_id %in%
#'                                                         annotation_df$annotations[, unique(transcript_id)]],
#'                                 ppidm, save = FALSE, load_dir = '/projectnb2/evolution/zwakefield/proteinImpacts/', init = TRUE)
#'
#'
#' hits_final <- get_ppi_switches(hits_domain, ppi)
#' @export
get_ppi_switches <- function(hits_all, ppi) {
  H <- as.data.table(hits_all)[, I := NULL]
  res <- H[, {

    if (pc_class == 'protein_coding') {
      inc_ppi <- unique(c(
        ppi[ensembl_transcript_id == transcript_id_inc,
            i.ensembl_transcript_id],
        ppi[i.ensembl_transcript_id == transcript_id_inc,
            ensembl_transcript_id]))

      exc_ppi <- unique(c(
        ppi[ensembl_transcript_id == transcript_id_exc,
            i.ensembl_transcript_id],
        ppi[i.ensembl_transcript_id == transcript_id_exc,
            ensembl_transcript_id]))

      inc_only <- setdiff(inc_ppi[inc_ppi != transcript_id_inc],
                          exc_ppi[exc_ppi != transcript_id_exc])

      exc_only <- setdiff(exc_ppi[exc_ppi != transcript_id_exc],
                          inc_ppi[inc_ppi != transcript_id_inc])

      .(
        inc_ppi = list(inc_only),
        exc_ppi = list(exc_only),
        n_inc_ppi = length(inc_only),
        n_exc_ppi = length(exc_only),
        n_ppi = length(exc_only)+length(inc_only)
      )
    } else {
      .(
        inc_ppi = list(NA_character_),
        exc_ppi = list(NA_character_),
        n_inc_ppi = 0L,
        n_exc_ppi = 0L,
        n_ppi = 0L
      )
    }

  },
  by = .I]

  cbind(H, res)[]
}

#' @title Plot summary of altered PPI interactions
#' @description
#' Visualizes the frequency and magnitude of gained/lost PPIs per event,
#' using a dual-panel layout:
#' - left: proportion of events with any PPI change
#' - right: histograms of INC and EXC partner counts (non-zero only)
#'
#' @param df `data.table` or `data.frame` with PPI counts per event,
#'   as returned by [ppi_switches_for_hits()].
#' @param bins Integer; number of histogram bins (default `30`).
#' @param palette Named character vector of fill colors for the plot
#'   (default includes `"no"`, `"yes"`, `"INC"`, `"EXC"`).
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
#' ppidm <- get_ppidm(test=TRUE)
#'
#' ppi <- get_isoform_interactions(protein_feature_total[ensembl_transcript_id %in%
#'                                                         annotation_df$annotations[, unique(transcript_id)]],
#'                                 ppidm, save = FALSE, load_dir = '/projectnb2/evolution/zwakefield/proteinImpacts/', init = TRUE)
#'
#'
#' hits_final <- get_ppi_switches(hits_domain, ppi)
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
                                         "INC" = "#2b8cbe", "EXC" = "#e34a33"),
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
    DT[, .(type = "INC", value = as.integer(n_inc_ppi))],
    DT[, .(type = "EXC", value = as.integer(n_exc_ppi))]
  )[value > 0]  # drop zeros as requested

  p_right <- ggplot(long_dt, aes(x = value, fill = type)) +
    geom_histogram(bins = bins, color = "white", linewidth = 0.2, show.legend = FALSE) +
    facet_wrap(~type, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = palette[c("INC","EXC")]) +
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


































































