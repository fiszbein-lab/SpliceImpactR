#' Visualize PSI values for a single splicing event
#'
#' @description
#' Creates a per-event PSI summary across samples. For **AFE** and **ALE**
#' events, the plot separates PSI by `inc` entry to distinguish alternative
#' terminal exons; other event types are summarized by `event_id`.
#'
#' @param data `data.frame` or `data.table` containing at least the columns
#'   `event_id`, `event_type`, `psi`, `condition`, and `sample`. For terminal
#'   events, an `inc` column is also required.
#' @param event Character scalar specifying the `event_id` to visualize.
#' @param fill_zeros Logical indicating whether to fill missing PSI values
#'   with zeros for samples that contain any observation of the event
#'   (default: `TRUE`).
#'
#' @return A list with elements:
#' \itemize{
#'   \item `plot`: `ggplot` box/point plot of PSI per sample group.
#'   \item `data`: `data.table` used to build the plot.
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
#' event_probe <- "ENSG00000117632:AFE"
#' probe_individual_event(hit_index, event = event_probe)
#'
#' @importFrom data.table setnames CJ 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point labs theme_minimal
#'   scale_y_continuous position_dodge position_jitterdodge
#' @export
probe_individual_event <- function(data, event, fill_zeros = TRUE) {
  dt <- data.table::as.data.table(data)
  
  required_cols <- c("event_id", "event_type", "psi", "condition", "sample")
  if (!all(required_cols %in% names(dt))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, names(dt)), collapse = ", "))
  }
  
  df <- dt[event_id == event]
  if (!nrow(df)) {
    stop("No rows found for event_id = '", event, "'.")
  }
  
  event_type <- unique(df$event_type)
  if (length(event_type) > 1) {
    warning("Multiple event_type values found; using the first: ", event_type[1])
  }
  
  grouping_col <- if (event_type[1] %chin% c("AFE", "ALE")) "inc" else "event_id"
  
  if (grouping_col == "inc" && !"inc" %in% names(df)) {
    stop("Column 'inc' is required to plot terminal exon events.")
  }
  
  grouping_values <- if (grouping_col == "inc") unique(df[[grouping_col]]) else event
  
  x <- df[, .(psi = suppressWarnings(as.numeric(psi)), condition, sample, grouping = get(grouping_col))]
  setnames(x, "grouping", grouping_col)
  
  if (fill_zeros) {
    all_samples <- unique(dt[, .(sample, condition)])
    samples_to_keep <- x[, .(has_value = any(!is.na(psi))), by = sample][has_value == TRUE, sample]
    
    if (length(samples_to_keep)) {
      all_samples <- all_samples[sample %in% samples_to_keep]
    }
    
    template <- CJ(temp_group = grouping_values, sample = all_samples$sample)
    setnames(template, "temp_group", grouping_col)
    template <- template[all_samples, on = .(sample)]
    
    x <- template[x, on = c("sample", grouping_col)]
    x[!is.finite(psi), psi := 0]
  } else {
    x <- x[is.finite(psi)]
  }
  
  plot_data <- data.table::copy(x)
  setnames(plot_data, grouping_col, "event_group")
  
  p <- ggplot(plot_data, aes(x = event_group, y = psi)) +
    geom_boxplot(aes(fill = condition), position = position_dodge(width = 0.8)) +
    geom_point(aes(color = condition),
               position = position_jitterdodge(dodge.width = 0.8),
               alpha = 0.6) +
    labs(
      x = grouping_col,
      y = "PSI (Percent Spliced In)",
      title = paste("PSI Distribution for Event:", event)
    ) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1))
  
  return(list(plot = p, data = x))
}

#' @title Split genomic coordinate strings into integer ranges
#' @description
#' Internal helper that converts a string like `"12345-12456"`
#' into an integer vector `c(12345, 12456)`.
#'
#' @param coord Character scalar of the form `"start-end"`.
#'
#' @return Integer vector of length two (`start`, `end`).
#' @keywords internal
.split_coord <- function(coord) {
  as.integer(strsplit(coord, "[-]")[[1]])
}

#' @title Classify splicing events as proximal or distal
#'
#' @description
#' Determines whether inclusion/exclusion events correspond to
#' proximal or distal terminal exon usage for **AFE** (Alternative
#' First Exon) and **ALE** (Alternative Last Exon) events, based on
#' genomic coordinates and strand orientation.
#'
#' @details
#' The function compares the genomic coordinates of the inclusion
#' (`inc_inc`) and exclusion (`inc_exc`) segments per event:
#' \itemize{
#'   \item For **AFE** events, proximal = exon with smaller start
#'         coordinate on the `+` strand (or larger end on `-` strand).
#'   \item For **ALE** events, proximal = exon with smaller start
#'         coordinate on the `+` strand (or larger end on `-` strand).
#' }
#' Events outside these types are labeled `"nonTerminal"`.
#'
#' If \code{plot = TRUE}, a summary donut chart is printed showing
#' the proportion of proximal vs distal usage per event type.
#'
#' @param hits `data.frame` or `data.table` containing at least:
#'   \itemize{
#'     \item `event_id`
#'     \item `event_type_inc`
#'     \item `strand_inc`
#'     \item `inc_inc`, `inc_exc`
#'     \item `delta_psi_inc`, `delta_psi_exc`
#'   }
#'
#' @return A `data.table` identical to `hits` with an additional
#' column (named `V1`) specifying `"proximal"`, `"distal"`,
#' `"overlap"`, or `"nonTerminal"`.
#'
#' @seealso [plot_prox_dist()]
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
#' proximal_output <- get_proximal_shift_from_hits(pairs)
#' @export
get_proximal_shift_from_hits <- function(hits) {
  H <- data.table::as.data.table(hits)[event_type_inc %in% c('AFE', 'ALE'), .(event_id, event_type = event_type_inc, strand = strand_inc, pos = inc_inc, delta_psi_inc, neg = inc_exc, delta_psi_exc)]
  res <- H[, {
    pos <- .split_coord(pos)
    neg <- .split_coord(neg)
    if (event_type  == "AFE") {
      if (strand == "+") {
        if (max(pos) > max(neg)) {
          "proximal"
        } else if (max(pos) < max(neg)) {
          "distal"
        } else {
          "overlap"
        }
      } else {
        if (min(pos) < min(neg)) {
          "proximal"
        } else if (min(pos) > min(neg)) {
          "distal"
        } else {
          "overlap"
        }
      }
    } else if (event_type == "ALE") {
      if (strand == "+") {
        if (min(pos) < min(neg)) {
          "proximal"
        } else if (min(pos) > min(neg)) {
          "distal"
        } else {
          "overlap"
        }
      } else {
        if (max(pos) > max(neg)) {
          "proximal"
        } else if (max(pos) < max(neg)) {
          "distal"
        } else {
          "overlap"
        }
      }
    } else {
      "nonTerminal"
    }
  }, by = .I
  ]
  H_prox <- cbind(H, res)
  plot <- plot_prox_dist(H_prox)
  return(list(data = H_prox,
              plot = plot))
}

#' @title Plot proximal vs distal exon usage
#'
#' @description
#' Creates a donut-style summary plot showing the number of events
#' classified as proximal or distal for each event type (AFE, ALE).
#'
#' @details
#' The plot shows each event type as a donut chart:
#' \itemize{
#'   \item Inner label = total number of events.
#'   \item Slices = counts of proximal and distal events.
#' }
#' Only event types with at least one valid classification are shown.
#'
#' @param res `data.table` output from [get_proximal_shift_from_hits()],
#'   containing at least columns `event_type` and `V1` (proximal/distal label).
#'
#' @return A `ggplot` object (donut-style bar chart).
#' @seealso [get_proximal_shift_from_hits()]
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_col coord_polar xlim geom_text
#'   scale_fill_manual facet_wrap theme_void labs theme element_text
#'   position_stack
#' @keywords internal
plot_prox_dist <- function(res) {

  DT <- as.data.table(res)
  DT[, event_type := toupper(as.character(event_type))]
  DT[, group := factor(tolower(V1), levels = c("proximal","distal"))]

  # tidy counts per type
  C <- DT[, .(N = .N), by = .(event_type, group)]
  # ensure both groups exist for each type
  CJT <- CJ(event_type = unique(DT$event_type),
            group = factor(c("proximal","distal"), levels = c("proximal","distal")))
  C <- CJT[C, on = .(event_type, group)]
  C[is.na(N), N := 0L]
  C[, total := sum(N), by = event_type]
  C[, pct := ifelse(total > 0, N / total, 0)]
  C[, N_adj := N*max(total)/total]
  p <- ggplot(C, aes(x = 2, y = N_adj, fill = group)) +
    geom_col(width = 1, color = NA) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +                 # makes the donut hole
    # per-slice labels positioned at stacked midpoints
    geom_text(aes(label = paste0(N)),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 3.8, lineheight = 0.95) +
    # per-facet totals printed in the center
    geom_text(data = unique(C[, .(event_type, total)]),
              aes(x = 0.5, y = 0, label = paste0(total)),
              inherit.aes = FALSE, size = 4) +
    scale_fill_manual(values = c(proximal = "cadetblue4", distal = "deeppink4")) +
    facet_wrap(~ event_type, nrow = 1) +
    theme_void() +
    labs(title = "Proximal vs Distal (AFE vs ALE)") +
    theme(strip.text = element_text(face = "bold"))
  return(p)

}


#' @title Compare inclusion vs. exclusion isoform lengths
#'
#' @description
#' Generates a multi-panel summary comparing isoform lengths between inclusion
#' (INC) and exclusion (EXC) events for either protein or transcript modes.
#'
#' @details
#' Produces three coordinated panels:
#' \enumerate{
#'   \item Paired boxplot showing INC vs EXC lengths per event with Wilcoxon
#'         paired test annotation.
#'   \item Density plot of delta length (INC - EXC).
#'   \item Barplot summarizing protein-coding categories.
#' }
#'
#' @param hits `data.frame` or `data.table` containing event-level results,
#'   including at least `prot_len_inc`, `prot_len_exc`, and `prot_len_diff`
#'   for `mode = "protein"`, or their transcript analogues
#'   `tx_len_inc`, `tx_len_exc`, and `tx_len_diff`.
#' @param phenotypes Named character vector of length 2 giving labels for
#'   control and experimental phenotypes. Must have names
#'   `"control"` and `"experimental"`.
#' @param mode Character, one of `"protein"` or `"transcript"`.
#'   Determines which length columns to use and whether to derive
#'   protein-coding categories.
#' @param output_file Optional path for saving the combined plot
#'   (e.g., `"length_summary.png"`).
#'
#' @return A composite `ggplot` object (assembled with `patchwork`) showing
#' paired boxplots, delta-length density, and protein-coding class distribution.
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
#' proximal_output <- plot_length_comparison(seq_compare)
#'
#' @seealso [plot_alignment_summary()]
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_boxplot geom_density
#'   geom_col geom_text geom_vline scale_fill_manual scale_x_discrete labs theme_bw
#'   theme annotate element_blank element_text position_dodge
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom tools toTitleCase
#' @importFrom stats wilcox.test
#' @export
plot_length_comparison <- function(
    hits,
    phenotypes = c(control = "control", experimental = "case"),
    mode = c("protein", "transcript")[1],
    output_file = NULL
) {
  stopifnot(length(phenotypes) == 2, all(c("control","experimental") %in% names(phenotypes)))
  mode <- match.arg(mode)

  H <- as.data.table(hits)

  # ---- choose columns --------------------------------------------------------
  if (mode == "protein") {
    cols <- c(inc = "prot_len_inc", exc = "prot_len_exc", diff = "prot_len_diff")
  } else {
    cols <- c(inc = "tx_len_inc",   exc = "tx_len_exc",   diff = "tx_len_diff")
  }
  if (!all(cols %in% names(H))) {
    stop("Missing required columns for mode='", mode, "'. Looking for: ",
         paste(cols, collapse=", "))
  }

  # ---- tidy for plots --------------------------------------------------------
  # All rows we can use for per-event diff
  Ddiff <- H[is.finite(get(cols["diff"])), .(diff = get(cols["diff"]))]

  # Rows where BOTH values exist (for paired plot & Wilcoxon)
  Dpair <- H[
    is.finite(get(cols["inc"])) & is.finite(get(cols["exc"])),
    .(value_inc = get(cols["inc"]), value_exc = get(cols["exc"]))
  ][, id := .I]

  # For protein mode, the classic PC categories (using > 0 length = coding)
  if (mode == "protein") {
    pc_cat <- H[, {
      inc_len <- suppressWarnings(as.numeric(get(cols["inc"])))
      exc_len <- suppressWarnings(as.numeric(get(cols["exc"])))
      fcase(
        (inc_len > 0) & (exc_len > 0),               "bothPC",
        (inc_len > 0) & (exc_len <= 0 | is.na(exc_len)), "experimental_only",
        (inc_len <= 0 | is.na(inc_len)) & (exc_len > 0), "control_only",
        default = "noPC"
      )
    }]
  } else {
    # transcript lengths are usually >0; we still offer analogous buckets
    pc_cat <- rep("bothPC", nrow(H))
  }
  PC <- data.table(cat = factor(pc_cat,
                                levels = c("bothPC","experimental_only","control_only","noPC")))
  PC_counts <- PC[, .N, by = cat]

  # ---- stats: Wilcoxon paired (like your old stat_compare_means) ------------
  w_p <- NA_real_
  if (nrow(Dpair)) {
    # If you want "bothPC only" for protein mode (matches your old workflow):
    if (mode == "protein") {
      both_idx <- which(pc_cat == "bothPC")
      if (length(both_idx)) {
        Dpair_use <- Dpair[id %in% both_idx]
      } else {
        Dpair_use <- Dpair[0]
      }
    } else {
      Dpair_use <- Dpair
    }
    if (nrow(Dpair_use)) {
      w_p <- tryCatch(
        suppressWarnings(
          wilcox.test(Dpair_use$value_inc, Dpair_use$value_exc, paired = TRUE)$p.value,
        ),
        error = function(e) NA_real_
      )
    }
  }
  p_label <- if (is.finite(w_p)) sprintf("Wilcoxon, p = %.3g", w_p) else "Wilcoxon, p = NA"

  # ---- plot 1: paired boxplot (INC vs EXC) ----------------------------------
  # Long form for plotting; keep only the rows used in test if protein mode
  if (mode == "protein") {
    keep_ids <- which(pc_cat == "bothPC")
    Dplot <- Dpair[id %in% keep_ids]
  } else {
    Dplot <- Dpair
  }
  Dlong <- rbind(
    Dplot[, .(id, condition = "control",      value = value_exc)],
    Dplot[, .(id, condition = "experimental", value = value_inc)]
  )
  Dlong[, condition := factor(condition, levels = c("control","experimental"))]

  p_box <- ggplot(Dlong, aes(condition, value, fill = condition)) +
    # paired lines
    geom_line(aes(group = id), color = "grey70", linewidth = 0.2, alpha = 0.6,
              position = position_dodge(width = 0)) +
    geom_point(shape = 21, color = "black", size = 1.2, alpha = 0.7) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, color = "black") +
    scale_fill_manual(values = c(control = "#B04A58", experimental = "#2C8C5A")) +
    scale_x_discrete(labels = phenotypes) +
    labs(x = "Phenotype",
         y = sprintf("%s Length", tools::toTitleCase(mode)),
         title = NULL) +
    annotate("text", x = 1.05, y = Inf, label = p_label,
             vjust = 1.5, hjust = 0, size = 4) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  # ---- plot 2: density of delta length (INC - EXC) ------------------------------
  p_density <- ggplot(Ddiff, aes(diff)) +
    geom_density(fill = "#7F5A83", alpha = 0.5, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = sprintf("delta %s Length (INC - EXC)", tools::toTitleCase(mode)),
         y = "Density",
         title = NULL) +
    theme_bw(base_size = 12)

  # ---- plot 3: coding-class counts ------------------------------------------
  nice_lab <- c(
    bothPC = "bothPC",
    experimental_only = phenotypes[["experimental"]],
    control_only = phenotypes[["control"]],
    noPC = "noPC"
  )
  p_pc <- ggplot(PC_counts, aes(cat, N, fill = cat)) +
    geom_col(width = 0.7, color = "black") +
    geom_text(aes(label = N), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c(
      bothPC = "#7E2D5B",
      experimental_only = "#2C8C5A",
      control_only = "#B04A58",
      noPC = "grey70"
    )) +
    scale_x_discrete(labels = nice_lab) +
    labs(x = NULL, y = "Count",
         title = "Protein-coding class (per event)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  # ---- assemble --------------------------------------------------------------
  title_txt <- sprintf("%s vs %s - %s length & coding summary",
                       phenotypes[["experimental"]],
                       phenotypes[["control"]],
                       tools::toTitleCase(mode))
  final_plot <- (p_box | p_density | p_pc) +
    plot_layout(widths = c(1.2, 1, 1.1)) +
    plot_annotation(title = title_txt)

  if (!is.null(output_file)) {
    ggsave(output_file, plot = final_plot, width = 12, height = 5.5)
  }

  final_plot
}


#' @title Plot alignment score distribution and coding summary
#'
#' @description
#' Visualizes alignment scores (`prot_pid` or `dna_pid`) and coding class
#' composition of hits, showing both the categorical composition and
#' histogram of alignment identity values.
#'
#' @param hits `data.frame` or `data.table` containing
#'   `prot_pid` or `dna_pid` and `summary_classification` columns.
#' @param mode Character, either `"protein"` (uses `prot_pid`)
#'   or `"transcript"` (uses `dna_pid`).
#' @param output_file Optional path to save the combined plot.
#'
#' @return A composite `ggplot` object (from `ggpubr::ggarrange`)
#' showing stacked bar counts by coding class and histogram of
#' alignment identity scores.
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
#' alignment_summary <- plot_alignment_summary(seq_compare)
#' @seealso [plot_length_comparison()]
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_bar geom_histogram scale_fill_manual
#'   theme_classic xlab ylab element_blank after_stat
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 facet_wrap
#' @export
plot_alignment_summary <- function(hits,
                                   mode = c("protein", "transcript")[1],
                                   output_file = NULL) {
  if (mode == "protein") pid_col <- 'prot_pid' else pid_col <- 'dna_pid'
  A <- hits[, .(score = get(pid_col), summary_classification, event_type_exc)]
  C <- A[, .N, by = .(summary_classification, event_type_exc)]
  A <- A[!is.na(score)]

  lvl <- rev(c("noPC", "onePC", "protein_coding", "FrameShift", "Rescue", "Match") )
  #                         ^ bottom         ^ middle       ^ top

  A[, summary_classification :=
      factor(summary_classification, levels = lvl)]
  C[, summary_classification :=
      factor(summary_classification, levels = lvl)]

  propCoding <- ggplot2::ggplot(C, ggplot2::aes(fill=summary_classification, x = 1, y = N)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'onePC' = "azure2", "Match" = "#E69F00", 'Rescue' = "#56B4E9", 'FrameShift' = "pink", 'protein_coding' = "deeppink4")) +
    ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("Count") + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
    facet_wrap(~event_type_exc, ncol = 1, scales = 'free_y')

  gdf <- ggplot2::ggplot(A, ggplot2::aes(x = score, fill = summary_classification)) +
    ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count/sum(count))), colour = 1,
                            bins = 20) +
    ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'onePC' = "azure2", "Match" = "#E69F00", 'Rescue' = "#56B4E9", 'FrameShift' = "pink", 'protein_coding' = "deeppink4")) +
    ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction")+
    facet_wrap(~event_type_exc, ncol = 1, scales = 'free_y')

  gdf_comp <- ggpubr::ggarrange(propCoding, gdf, nrow = 1, widths = c(1, 3),
                                common.legend = TRUE)
  if (!is.null(output_file)) {
    ggsave(output_file, plot = final_plot, width = 12, height = 5.5)
  }
  return(gdf_comp)
}

#' Compute library size factors from exon-level read counts
#'
#' @description
#' Estimates per-sample normalization factors (size factors) using
#' the **median-of-ratios** method based on exon-level read coverage
#' (nUP + nDOWN). This approach provides robust depth normalization
#' for splicing event comparisons.
#'
#' @param df `data.frame` or `data.table` with at least columns:
#'   `sample_name`, `exon`, `nUP`, `nDOWN`.
#'
#' @details
#' The function aggregates reads per exon across samples,
#' constructs an exon x sample matrix, and computes size factors:
#' \deqn{SF_j = median_i ((counts_{ij} + 1) / geoMean_i)}
#' where the geometric mean is computed across samples for each exon.
#'
#' @return A `data.frame` identical to the input with an added `sizeFactor`
#'   column, suitable for downstream normalization.
#'
#' @seealso [.get_size_factors_from_exons()]
#' @importFrom dplyr group_by summarize left_join
#' @importFrom tidyr pivot_wider
#' @keywords internal
getSizeFactors <- function(df) {
  df$reads <- df$nUP + df$nDOWN
  counts_wide <- df %>%
    dplyr::group_by(sample_name, exon) %>%
    dplyr::summarize(reads = sum(reads), .groups = "drop") %>%  # aggregate duplicates
    tidyr::pivot_wider(
      names_from  = sample_name,   # columns become sample names
      values_from = reads,
      values_fill = 0              # fill missing combos with 0
    )
  exon_ids <- counts_wide$exon

  # Drop the 'exon' column
  counts_mat <- as.matrix(counts_wide[, -1])

  rownames(counts_mat) <- exon_ids

  geoMean <- exp(rowMeans(log(counts_mat+1)))

  sizeFactors_list <- apply(counts_mat, 2, function(col_j) {
    stats::median((col_j+1)/geoMean)
  })

  df <- dplyr::left_join(df, data.frame(sample_name = names(sizeFactors_list),
                                        sizeFactor = as.numeric(sizeFactors_list)))
  return(df)
}

#' Compute size factors directly from exon count files
#'
#' @description
#' Wrapper around `getSizeFactors()` that can either compute
#' normalization factors from raw exon count files (`exon_files`)
#' or use user-provided values (`user_given`).
#'
#' @param sample_df `data.frame` containing sample metadata with
#'   columns `sample_name`, `path`, and optionally `condition`.
#' @param method Either `'exon_files'` (compute) or `'user_given'` (join).
#'
#' @details
#' - If `method = 'exon_files'`, the function loads exon-level count tables
#'   using `.read_exon_files()` and computes per-sample size factors.
#' - If `method = 'user_given'`, it merges the user-supplied `size_factors`
#'   data frame into `sample_df`.
#'
#' Size factors are normalized to have a geometric mean of 1
#' (DESeq2-style).
#'
#' @return A `data.frame` equal to `sample_df` with an appended
#'   numeric column `sizeFactor`.
#'
#' @seealso [getSizeFactors()]
#' @importFrom data.table dcast rbindlist setDT
#' @importFrom stats median
#' @keywords internal
.get_size_factors_from_exons <- function(sample_df, method = c('exon_files', 'user_given')) {
  if (method == 'user-given') {
    sf <- merge(sample_df, size_factors, by = sample_name)
    return(sf)
  } else {
    counts_long <- data.table::rbindlist(lapply(seq_len(nrow(sample_df)), function(i) {
      x <- sample_df$path[i]
      s <- sample_df$sample_name[i]

      exon_file <- .read_exon_files(
        paste0(x, basename(x), "."),           # keep your pattern
        columns = c("gene", "exon", "nUP", "nDOWN")
      )
      data.table::setDT(exon_file)
      exon_file[, reads := nUP + nDOWN]
      exon_file <- exon_file[reads > 10]
      exon_file[, .(reads = sum(reads)), by = .(exon)][, sample_name := s]
    }), use.names = TRUE, fill = TRUE)

    # 2) Cast to exon x sample matrix
    counts_wide <- data.table::dcast(
      counts_long, exon ~ sample_name, value.var = "reads", fill = 0
    )
    exon_ids    <- counts_wide$exon
    counts_mat  <- as.matrix(counts_wide[, -1])
    rownames(counts_mat) <- exon_ids

    gmean <- function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(NA_real_)
      exp(mean(log(x)))
    }
    counts_mat <- counts_mat + 1
    geoMean <- apply(counts_mat, 1, gmean)

    # Keep informative exons only
    keep <- is.finite(geoMean) & geoMean > 0
    cm   <- counts_mat[keep, , drop = FALSE]
    gm   <- geoMean[keep]

    # 4) Median of ratios per sample (robust to a few DE exons)
    sizeFactors <- apply(cm, 2, function(col_j) stats::median(col_j / gm, na.rm = TRUE))

    # Optional: normalize to geometric mean = 1 (DESeq2-style)
    sizeFactors <- sizeFactors / exp(mean(log(sizeFactors)))

    # 5) Return as a data.frame (ready to join back)
    return(merge(sample_df, data.frame(
      sample_name = names(sizeFactors),
      sizeFactor  = as.numeric(sizeFactors),
      row.names   = NULL,
      check.names = FALSE
    ), by = 'sample_name'))
  }
}


#' Overview of global splicing event distributions between conditions
#'
#' @description
#' Generates a multi-panel comparison summarizing global splicing
#' characteristics (event counts, events-per-gene, PSI distributions)
#' between experimental and control conditions, normalized by sequencing depth.
#'
#' @param events `data.frame` or `data.table` of splicing events with at least:
#'   `sample`, `condition`, `gene_id`, `inclusion_reads`, `exclusion_reads`, `psi`,
#'   `inc`, `exc`.
#' @param sample_df Metadata `data.frame` with sample-level paths and conditions.
#' @param depth_norm Normalization mode: `'exon_files'` (compute from read files)
#'   or `'user-given'` (use provided size factors).
#' @param event_type Character string specifying the event class (e.g., `"AFE"`, `"ALE"`).
#' @param conditions Named character vector mapping `control` and `experimental`
#'   condition labels.
#' @param minReads Minimum reads required for inclusion or exclusion
#'   (default = 10).
#' @param output_file Optional path to save a combined summary figure.
#'
#' @details
#' The output combines:
#' - **Panel A:** Depth-normalized event counts per sample (Wilcoxon test)
#' - **Panel B:** Mean events per gene per sample (Wilcoxon test)
#' - **Panel D:** PSI cumulative distribution comparison (K-S test)
#'
#' Size factors are computed internally using
#' [.get_size_factors_from_exons()], enabling robust normalization.
#'
#' @return Invisibly returns a combined [`patchwork`] plot object.
#'
#' @seealso [.get_size_factors_from_exons()], [getSizeFactors()]
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter scale_fill_manual
#'   labs annotate theme_bw theme stat_ecdf scale_color_manual
#' @importFrom patchwork plot_annotation
#' @importFrom stats wilcox.test ks.test
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
#' ov <- overview_splicing_comparison_fixed(hit_index, sample_frame, 'exon_files')
#' @export
overview_splicing_comparison_fixed <- function(
    events,
    sample_df,
    depth_norm = c('exon_files', 'user-given'),
    event_type = "AFE",                         # or "ALE", etc.
    conditions = c(control = "control", experimental = "case"),
    minReads = 10,
    output_file = NULL
) {
  # ----------------------- validate & normalize -----------------------
  E <- as.data.table(events)

  E <- E[(as.numeric(inclusion_reads) >= minReads | as.numeric(exclusion_reads) >= minReads) & is.finite(psi) & psi >= 0 & psi <= 1]

  # canonical labels
  E[, cond_label := fcase(
    condition == conditions[["control"]],       conditions[["control"]],
    condition == conditions[["experimental"]],  conditions[["experimental"]],
    default = as.character(condition)
  )]

  # unique event key
  E[, event_key := paste(gene_id, inc, exc, sep="|")]

  # ----------------------- Panel A: sample event counts -----------------------
  sample_event_counts <- unique(E[, .(sample, cond_label, event_key)])[, .N, by=.(sample, cond_label)]
  setnames(sample_event_counts, "N", "event_count")

  sample_df <- .get_size_factors_from_exons(sample_df, method = depth_norm)
  s <- sample_event_counts[sample_df, on = .(sample = sample_name, cond_label=condition)]
  s[, event_count_depth_corrected := event_count/sizeFactor]

  wA <- tryCatch(
    wilcox.test(event_count_depth_corrected ~ cond_label, data = as.data.frame(s))$p.value,
    error = function(e) NA_real_
  )
  pA <- ggplot(s, aes(cond_label, event_count_depth_corrected, fill = cond_label)) +
    geom_violin(trim = FALSE, alpha = .75) +
    geom_jitter(width = .12, size = 1, alpha = .35) +
    scale_fill_manual(values = c("#8B2A3C","#2E8B57")) +
    labs(x = "Group", y = paste("Count of", event_type, " (depth normalized)")) +
    annotate("text", x = 1, y = Inf,
             label = sprintf("Wilcoxon p = %s",
                             ifelse(is.finite(wA), formatC(wA, format="g", digits=3), "NA")),
             vjust = 1.4, hjust = 0, size = 4) +
    theme_bw() + theme(legend.position = "none")

  # ----------------------- Panel B: mean events-per-gene per sample -----------
  # events per gene within each sample
  epg_sample_gene <- unique(E[, .(sample, cond_label, gene_id, event_key)])[,
                                                                            .N, by=.(sample, cond_label, gene_id)][
                                                                              sample_df, on = .(sample=sample_name, cond_label=condition)]
  mean_epg_sample <- epg_sample_gene[, .(mean_events_per_gene = mean(N/sizeFactor)), by=.(sample, cond_label)]

  wB <- tryCatch(
    wilcox.test(mean_events_per_gene ~ cond_label, data = as.data.frame(mean_epg_sample))$p.value,
    error = function(e) NA_real_
  )
  pB <- ggplot(mean_epg_sample, aes(cond_label, mean_events_per_gene, fill = cond_label)) +
    geom_violin(trim = FALSE, alpha = .75) +
    geom_jitter(width = .12, size = 1, alpha = .35) +
    scale_fill_manual(values = c("#8B2A3C","#2E8B57")) +
    labs(x = "Group", y = paste("Mean", event_type, "\nper gene (depth normalized)")) +
    annotate("text", x = 1, y = Inf,
             label = sprintf("Wilcoxon p = %s",
                             ifelse(is.finite(wB), formatC(wB, format="g", digits=3), "NA")),
             vjust = 1.4, hjust = 0, size = 4) +
    theme_bw() + theme(legend.position = "none")

  # ----------------------- Panel D: PSI ECDF + KS test ------------------------
  psi_df <- E[, .(psi, cond_label)]
  ks_p <- tryCatch(
    suppressWarnings(ks.test(psi_df$psi[psi_df$cond_label == conditions[["control"]]],
                             psi_df$psi[psi_df$cond_label == conditions[["experimental"]]])$p.value),
    error = function(e) NA_real_
  )
  pD <- ggplot(psi_df, aes(psi, colour = cond_label)) +
    stat_ecdf(geom = "step", linewidth = .8) +
    scale_color_manual(values = c("#8B2A3C","#2E8B57")) +
    labs(x = "PSI", y = paste(event_type, "PSI eCDF")) +
    theme_bw() +
    annotate("text", x = 0.70, y = 0.30,
             label = paste0("K-S p-value = ", ifelse(is.finite(ks_p), formatC(ks_p, format="g", digits=3), "NA")))

  # ----------------------- assemble A/B/C/D -----------------------------------
  # combined <- ((pA + pB) / (pD) +
  #   plot_annotation(tag_levels = "A")) &
  #   theme(plot.tag = element_text(face = "bold", size = 14))
  #
  combined <- ((pA + pB) / (pD)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(face = "bold", size = 14))
    )

  if (!is.null(output_file)) {
    ggsave(output_file, plot = combined, width = 9, height = 8)
  }

  invisible(combined)
}

#' Read HIT index exon tables for all samples
#'
#' @description
#' Internal helper that loads exon-level HIT index tables for each sample
#' and appends corresponding `sample_name` and `condition` metadata.
#'
#' @param sample_df `data.frame` or `data.table` containing at least:
#'   \itemize{
#'     \item `path` directory path for each sample
#'     \item `sample_name` unique sample identifier
#'     \item `condition` experimental group
#'   }
#'
#' @return A `data.table` with columns:
#'   `gene`, `exon`, `HITindex`, `sample`, and `condition`.
#'
#' @keywords internal
#' @importFrom data.table rbindlist
.getHITindex <- function(sample_df) {
  data.table::rbindlist(lapply(seq_len(nrow(sample_df)), function(i) {
    x <- sample_df$path[i]
    s <- sample_df$sample_name[i]

    exon_file <- .read_exon_files(
      paste0(x, basename(x), "."),           # keep your pattern
      columns = c("gene", "exon", "HITindex")
    )
    exon_file[, `:=` (sample = sample_df$sample_name[i],
                      condition = sample_df$condition[i])]
  }), use.names = TRUE, fill = TRUE)
}

#' Compare HIT index values between conditions
#'
#' @description
#' Computes per-event differences in the **HIT index** (Hybrid Intron Tolerance)
#' between control and test conditions, performs t test,
#' adjusts for multiple testing, and produces summary visualizations.
#'
#' @param sample_df `data.frame` or `data.table` containing sample metadata
#'   with columns: `path`, `sample_name`, and `condition`.
#' @param condition_map Named character vector mapping experimental groups,
#'   e.g. `c(control = "control", test = "case")`.
#' @param minimum_proportion Minimum proportion of samples per group
#'   that must have valid HIT values for a test to be performed
#'   (default = 0.5).
#' @param top_n_heatmap Integer; number of events to display in the heatmap
#'   ordered by absolute deltaHIT (default = 5000).
#' @param sig_delta Absolute deltaHIT threshold for volcano highlighting
#'   (default = 0.2).
#' @param fdr_thresh FDR threshold for volcano significance lines
#'   (default = 0.05).
#'
#' @details
#' For each exon-level HIT index event:
#' \itemize{
#'   \item Computes mean HIT per group (`control`, `test`)
#'   \item Performs a t test if both groups
#'         have >50% of samples with valid values
#'   \item Adjusts p-values via Benjamini-Hochberg FDR
#'   \item Computes signed (`diff_HIT`) and absolute (`delta_HIT`) changes
#' }
#'
#' The function returns both the full per-event results table
#' and a combined 2x2 patchwork plot containing:
#' \enumerate{
#'   \item Control vs Test mean HIT scatter
#'   \item Top |deltaHIT| heatmap (control/test columns)
#'   \item Volcano plot of |deltaHIT| vs -log10(FDR)
#'   \item Density of deltaHIT distribution
#' }
#'
#' @return A named list with:
#' \describe{
#'   \item{`results`}{Data.table of per-event statistics (means, p, FDR, deltaHIT)}
#'   \item{`plot`}{Patchwork object with 4-panel summary visualization}
#' }
#'
#' @seealso [.getHITindex()]
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
#' hit_compare <- compare_hit_index(sample_frame, condition_map = c(control = "control", test = "case"))
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient geom_abline
#'   labs theme_minimal element_text geom_tile scale_fill_gradientn
#'   geom_hline geom_vline scale_color_manual geom_density geom_rug
#' @importFrom stats wilcox.test p.adjust t.test
#' @importFrom patchwork plot_annotation
#' @importFrom scales rescale
#' @export
compare_hit_index <- function(
    sample_df,
    condition_map = c(control = "control", test = "case"),
    minimum_proportion = 0.5,
    top_n_heatmap = 5000,
    sig_delta = 0.2,
    fdr_thresh = 0.05
) {

  DT <- .getHITindex(sample_df)

  # event key
  DT[, event_key := do.call(paste, c(.SD, list(sep="|"))), .SDcols = c("gene", "exon")]

  # keep only two groups we care about, rename consistently
  ctrl_lab <- condition_map[["control"]]
  test_lab <- condition_map[["test"]]
  DT <- DT[condition %chin% c(ctrl_lab, test_lab)]
  DT[, group := fifelse(condition == ctrl_lab, "control", "test")]

  # numeric HIT, drop NA
  DT[, hit := suppressWarnings(as.numeric(HITindex))]
  DT <- DT[is.finite(hit)]

  # -------------------- per-event statistics -------------------
  # summary means by event & group
  means <- DT[, .(mean_hit = mean(hit, na.rm = TRUE)),
              by = .(event_key, group)]
  means <- dcast(means, event_key ~ group, value.var = "mean_hit")

  cC <- nrow(DT[group == 'control',unique(sample), by = group])
  cT <- nrow(DT[group == 'test',unique(sample), by = group])

  ev_ok <- DT[, .(
    n_control = sum(group=="control" & !is.na(hit)),
    n_test = sum(group=="test"    & !is.na(hit)),
    uniq_hits = length(unique(hit)),
    maxdiff = abs(max(hit)-min(hit))

  ), by = event_key][
    n_control >= 0.5*cC & n_test >= 0.5*cT & uniq_hits > 1 & maxdiff > .05
  ]
  DTf <- DT[event_key %chin% ev_ok$event_key]

  wcx <- DTf[, {
    ctrl_hits <- .SD[group == "control", hit]
    test_hits <-.SD[group == "test", hit]
    p.val <- 1
    if (length(ctrl_hits) > 2 & length(test_hits) > 2) {
      p.val <- suppressWarnings(t.test(ctrl_hits, test_hits)$p.value)
    }
    .(
      p_value = p.val,
      n_control = length(ctrl_hits),
      n_test = length(test_hits)
    )
  }, by = event_key]
  wcx <- wcx[!is.na(p_value)]
  # combine & compute deltas + FDR
  res <- means[wcx, on = "event_key"]
  # control/test means might be missing if a group has zero data for that event
  res[, `:=`(
    delta_HIT = abs(test - control),
    diff_HIT  = (test - control)  # signed change (test minus control)
  )]
  res[, fdr := p.adjust(p_value, method = "fdr")]

  # tidy, keep useful columns
  setcolorder(res, c("event_key","control","test","diff_HIT","delta_HIT",
                     "n_control","n_test","p_value","fdr"))

  # ----------------------- plotting data ----------------------
  # scatter
  scatter_dt <- res[is.finite(control) & is.finite(test)]

  # volcano
  volcano_dt <- res[!is.na(fdr) & is.finite(delta_HIT)]

  # heatmap (top |deltaHIT| events, two columns: control, test)
  hm_dt <- res[is.finite(control) & is.finite(test)]
  hm_dt <- hm_dt[order(-delta_HIT)]
  # if (nrow(hm_dt) > top_n_heatmap) hm_dt <- hm_dt[1:top_n_heatmap]
  # long format for ggplot heatmap
  hm_long <- melt(hm_dt[, .(event_key, control, test)], id.vars = "event_key",
                  variable.name = "group", value.name = "mean_hit")
  hm_long[, group := factor(group, levels = c("control","test"))]
  # order events by control mean (like your old code)
  ord <- hm_dt[order(-control), event_key]
  hm_long[, event_key := factor(event_key, levels = ord)]

  # delta distribution
  dens_dt <- res[is.finite(diff_HIT)]

  # ---------------------------- plots -------------------------
  # 1) Scatter: control vs test mean HIT, color by |deltaHIT|
  p_scatter <- ggplot(scatter_dt,
                      aes(x = control, y = test, color = delta_HIT)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_gradient(low = "white", high = "red",
                         name = expression(paste("|", Delta, " HIT|"))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Average HIT (Control)", y = "Average HIT (Test)",
         title = "Control vs Test Mean HIT") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # 2) Heatmap: two columns (control/test) for top |deltaHIT| events
  p_heat <- ggplot(hm_long, aes(x = group, y = event_key, fill = mean_hit)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("purple", "orange", "grey80", "turquoise", "red"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      name = "Mean HIT"
    ) +
    labs(x = NULL, y = NULL, title = "Top |deltaHIT| events (means)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  # 3) Volcano: |deltaHIT| vs -log10(FDR)
  p_volcano <- ggplot(volcano_dt,
                      aes(x = delta_HIT, y = -log10(fdr))) +
    geom_point(aes(color = fdr < fdr_thresh & delta_HIT > sig_delta),
               alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey70"),
                       guide = "none") +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = sig_delta, linetype = "dotted") +
    labs(x = expression(paste("|", Delta, " HIT|")),
         y = expression(-log[10]~FDR),
         title = "Volcano") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # 4) deltaHIT distribution (signed)
  p_delta <- ggplot(dens_dt, aes(x = diff_HIT)) +
    geom_density(fill = "deeppink4", alpha = .3) +
    geom_rug(alpha = .15) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = expression(Delta~"HIT (Test - Control)"), y = "Density",
         title = expression(Delta~"HIT distribution")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # combined 2x2
  combined <- (p_scatter + p_heat) / (p_volcano + p_delta)

  list(
    results = res[],
    plot = combined
  )
}


#' @title Integrated summary of event classification, alignment, and domain changes
#'
#' @description
#' Provides a comprehensive visualization and summary of event classification outcomes
#' (frame shifts, rescues, matches, etc.), alignment quality, and domain change prevalence
#' across alternative splicing event types.
#' Integrates results from [`compare_sequence_frame()`] with pre-filter event tables
#' to display pre- vs post-filter usage, event coordination, and gene/domain overlaps.
#'
#' @param hits `data.frame` or `data.table` output from [`compare_sequence_frame()`],
#'   containing per-event alignment and domain information (e.g. `summary_classification`,
#'   `prot_pid`, `event_type_inc`, `inc_only_n`, `exc_only_n`, etc.).
#' @param pre_filter_hits `data.frame` or `data.table` representing the unfiltered
#'   input events (e.g. raw differential inclusion table prior to sequence comparison).
#'
#' @details
#' The function integrates multiple layers of event-level characterization:
#' \itemize{
#'   \item Event-type composition and class proportions
#'   \item Protein alignment quality distribution
#'   \item Domain-change prevalence (INC/EXC/both)
#'   \item Gene and domain cross-event intersections (UpSet plots)
#'   \item Event-type coordination heatmap (Jaccard similarity across genes)
#'   \item Relative event retention pre- vs post-filtering
#' }
#'
#' @return A named list with:
#' \describe{
#'   \item{`summaries`}{List containing per-type tables:
#'     \code{by_type}, \code{class_counts}, \code{score_summary},
#'     \code{domain_prevalence}, and \code{relative_use}.}
#'   \item{`plot`}{A multi-panel [`patchwork`] composite summarizing
#'     classification, alignment, domain changes, and coordination.}
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
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' hits_domain <- get_domains(seq_compare, exon_features)
#' ppidm <- get_ppidm(test=TRUE)
#'
#' restrict_protein_features <- protein_feature_total[ensembl_transcript_id %in%
#'                                                      c(hits_domain$transcript_id_exc, hits_domain$transcript_id_exc)]
#' ppi <- get_isoform_interactions(restrict_protein_features, ppidm, save = FALSE, load_dir = NULL, init = TRUE)
#' hits_final <- get_ppi_switches(hits_domain, ppi)
#' int_summary <- integrated_event_summary(hits_final, res)
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_col geom_violin geom_tile geom_text geom_abline
#'   geom_point scale_fill_manual scale_y_continuous labs theme_minimal theme
#'   element_text ggtitle scale_color_gradient scale_fill_gradient
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom scales percent percent_format rescale
#' @importFrom stats median quantile
#' @importFrom ggplot2 stat_summary ylim xlim
#' @importFrom ComplexUpset upset upset_set_size intersection_size
#' @export
integrated_event_summary <- function(
    hits,
    pre_filter_hits
) {
  DT <- as.data.table(hits)

  DT[, event_type := event_type_inc]
  DT <- DT[!is.na(event_type) & nzchar(event_type)]

  DT[, class_raw := as.character(summary_classification)]
  # map any synonyms to our palette keys
  map_class <- c(
    "Match"="Match", "PartialMatch"="PartialMatch", "FrameShift"="FrameShift",
    "onePC"="onePC", "noPC"="noPC", "Rescue"="Rescue",
    "protein_coding"="Match"  # treat as Match color if this label appears
  )
  DT[, class := fifelse(class_raw %in% names(map_class), map_class[class_raw], class_raw)]

  ## ---------- color palette (same as before) ----------
  class_pal <- c(
    noPC         = "azure4",
    Match        = "#E69F00",
    onePC        = "azure2",
    FrameShift   = "pink",
    protein_coding = "deeppink4",
    Rescue       = "#56B4E9"
  )

  DT[, score := as.numeric(prot_pid)]

  ## ---------- domain-change fields ----------
  DT[, `:=`(
    inc_only_flag = as.integer(inc_only_n > 0),
    exc_only_flag = as.integer(exc_only_n > 0),
    any_dom_flag  = as.integer((inc_only_n > 0) | (exc_only_n > 0))
  )]
  DT[, dom_kind := fifelse(inc_only_flag==1 & exc_only_flag==0, "inc_only",
                           fifelse(inc_only_flag==0 & exc_only_flag==1, "exc_only",
                                   fifelse(inc_only_flag==1 & exc_only_flag==1, "both", "none")))]

  ## =====================  summaries  =====================
  by_type <- DT[, .(n = .N), by = event_type][order(-n)]

  class_counts <- DT[, .N, by = .(event_type, summary_classification)]
  class_props  <- class_counts[, .(prop = N/sum(N)), by = summary_classification]
  lvl <- rev(c("noPC", "onePC", "protein_coding", "FrameShift", "Rescue", "Match") )
  #                         ^ bottom         ^ middle       ^ top
  class_counts[, summary_classification :=
                 factor(summary_classification, levels = lvl)]
  class_props[, summary_classification :=
                factor(summary_classification, levels = lvl)]

  score_summary <- DT[is.finite(prot_score),
                      .(n = .N,
                        median = median(prot_score),
                        q25 = quantile(prot_score, .25),
                        q75 = quantile(prot_score, .75)),
                      by = event_type][order(-n)]

  dom_prev <- DT[, .(
    prop_any  = mean(any_dom_flag==1, na.rm = TRUE),
    prop_inc  = mean(dom_kind=="inc_only", na.rm = TRUE),
    prop_exc  = mean(dom_kind=="exc_only", na.rm = TRUE),
    prop_both = mean(dom_kind=="both",     na.rm = TRUE)
  ), by = event_type]

  ## =====================  plots  =====================


  # P1: stacked class proportions per event type
  p1 <- ggplot(class_counts, aes(event_type, N, fill = summary_classification)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = class_pal, drop = FALSE) +
    scale_y_continuous(labels = percent) +
    labs(x = NULL, y = "Proportion", fill = "Class",
         title = "Event classification by type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # P2: alignment score violin per event type (with median)
  p2 <- if (DT[, sum(is.finite(score))] > 0) {
    ggplot(DT[is.finite(score)], aes(event_type, score)) +
      geom_violin(fill = "grey92", color = "grey50") +
      stat_summary(fun = median, geom = "point", size = 2, color = "deeppink4") +
      labs(x = NULL, y = paste("Protein", "alignment (%)"),
           title = "Alignment score distribution") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else ggplot() + theme_void() + ggtitle("Alignment scores not available")

  # P3: domain change prevalence (any/inc/exc/both) per type
  dom_long <- melt(dom_prev, id.vars = "event_type",
                   variable.name = "kind", value.name = "prop")

  p3 <- ggplot(dom_long, aes(event_type, prop, fill = kind)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c(
      prop_any  = "cadetblue4",
      prop_inc  = "deeppink4",
      prop_exc  = "#56B4E9",
      prop_both = "#E69F00"
    ), labels = c("Any","INC-only","EXC-only","Both")) +
    labs(x = NULL, y = "Proportion", fill = NULL,
         title = "Domain-change prevalence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  have_upset <- requireNamespace("ComplexUpset", quietly = TRUE)

  DT2 <- as.data.table(DT)[pc_class == "protein_coding"]

  DT2[, any_ppi := as.integer(n_ppi > 0)]
  DT2[, event_type := factor(event_type,
                             levels = c("AFE","ALE","A3SS","A5SS","SE","RI","MXE"))]


  ppi1 <- DT2[, .(prop_any_ppi = mean(any_ppi)), by = event_type] %>%
    ggplot(aes(x = event_type, y = prop_any_ppi)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = scales::percent(prop_any_ppi, .1)),
              vjust = -0.3) +
    labs(y = "Proportion with PPI change", x = "",
         title = "Fraction of events with altered PPIs") +
    theme_bw(base_size = 12) +
    ylim(0, 1)

  DT_long <- melt(
    DT2,
    id.vars = "event_type",
    measure.vars = c("n_inc_ppi", "n_exc_ppi"),
    variable.name = "ppi_direction",
    value.name  = "ppi_count"
  )

  DT_long[, ppi_direction := factor(ppi_direction,
                                    labels = c("INC gained","EXC gained"))]

  ppi2 <- ggplot(DT_long[ppi_count > 0],
               aes(x = event_type, y = ppi_count, color = ppi_direction)) +
    # geom_jitter(alpha = 0.4, width = 0.2) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    labs(x = "Event type", y = "# PPI gained",
         title = "PPI gain distributions by event type",
         color = "") +
    theme_bw(base_size = 12)

  ## ---------- 1) Gene x Event-Type UpSet ----------
  # Coalesce gene id (prefers *_inc, falls back to *_exc, then NA removed)
  gene_col <- c("gene_id_inc","gene_id_exc")[c("gene_id_inc","gene_id_exc") %in% names(DT)][1]
  if (is.na(gene_col)) {
    gene_col <- c("gene_id","gene")[c("gene_id","gene") %in% names(DT)][1]
  }
  if (!is.na(gene_col)) {
    G <- unique(DT[!is.na(get(gene_col)) & nzchar(get(gene_col)),
                   .(gene = as.character(get(gene_col)), event_type)])
    gene_wide <- dcast(G[, .(present = TRUE), by = .(gene, event_type)],
                       gene ~ event_type, value.var = "present", fill = FALSE)
    # Optional: hide hybrid redundancy if both AFE/HFE or ALE/HLE exist for same gene
    for (pair in list(c("AFE","HFE"), c("ALE","HLE"))) {
      if (all(pair %in% names(gene_wide))) {
        gene_wide[[pair[2]]] <- ifelse(gene_wide[[pair[1]]] & gene_wide[[pair[2]]], FALSE, gene_wide[[pair[2]]])
      }
    }
    if (have_upset && ncol(gene_wide) > 2) {
      # gene_upset <- ComplexUpset::upset(
      #   gene_wide,
      #   intersect = setdiff(names(gene_wide), "gene"),
      #   name = "Gene intersections",
      #   set_sizes = ComplexUpset::upset_set_size(),
      #   base_annotations = list(
      #     'Counts' = ComplexUpset::intersection_size(text = list(size = 3))
      #   ),
      #   width_ratio = 0.15
      # )
      gene_upset <- ComplexUpset::upset(
        gene_wide,
        intersect = setdiff(names(gene_wide), "gene"),
        name = "Gene intersections",
        set_sizes = ComplexUpset::upset_set_size(),
        base_annotations = list(
          'Counts' = ComplexUpset::intersection_size(text = list(size = 3))
        ),
        width_ratio = 0.15,
        themes = list(
          # valid ggplot theme for elements upset will override
          `intersections_matrix` = theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
          )
        )
      )
    } else {
      gene_upset <- ggplot() + theme_void() + ggtitle("Gene UpSet (ComplexUpset not installed / not enough sets)")
    }
  } else {
    gene_upset <- ggplot() + theme_void() + ggtitle("Gene UpSet (no gene column found)")
  }

  domain_list_cols <- c("either_domains_list","inc_only_domains_list","exc_only_domains_list")
  dom_col_present <- domain_list_cols[domain_list_cols %in% names(DT)]

  if (length(dom_col_present)) {
    # Build a "domain" vector per row (prefers 'either' if present)
    if ("either_domains_list" %in% dom_col_present) {
      Dtmp <- DT[, .(event_type,
                     domain = unlist(either_domains_list, use.names = FALSE)), by = .I]
    } else {
      # fallback: union inc-only and exc-only
      Dtmp <- DT[, {
        incv <- if ("inc_only_domains_list" %in% dom_col_present) unlist(inc_only_domains_list) else character(0)
        exec <- if ("exc_only_domains_list" %in% dom_col_present) unlist(exc_only_domains_list) else character(0)
        list(domain = unique(c(incv, exec)))
      }, by = .(event_type, .I)]
    }
    Dtmp <- Dtmp[!is.na(domain) & nzchar(domain)]
    if (nrow(Dtmp)) {
      Dwide <- dcast(unique(Dtmp[, .(domain, event_type, present = TRUE)]),
                     domain ~ event_type, value.var = "present", fill = FALSE)
      if (have_upset && ncol(Dwide) > 2) {
        domain_upset <- ComplexUpset::upset(
          Dwide,
          intersect = setdiff(names(Dwide), "gene"),
          name = "Domain intersections",
          set_sizes = ComplexUpset::upset_set_size(),
          base_annotations = list(
            'Counts' = ComplexUpset::intersection_size(text = list(size = 3))
          ),
          width_ratio = 0.15,
          themes = list(
            # valid ggplot theme for elements upset will override
            `intersections_matrix` = theme(
              axis.title.x = element_blank(),
              axis.title.y = element_blank()
            )
          )
        )

        # domain_upset <- ComplexUpset::upset(
        #   Dwide,
        #   intersect = setdiff(names(Dwide), "domain"),
        #   name = "Domain intersections",
        #   set_sizes = ComplexUpset::upset_set_size(),
        #   base_annotations = list(
        #     'Counts' = ComplexUpset::intersection_size(text = list(size = 3))
        #   ),
        #   width_ratio = 0.15
        # )
      } else {
        domain_upset <- ggplot() + theme_void() + ggtitle("Domain UpSet (ComplexUpset not installed / not enough sets)")
      }
    } else {
      domain_upset <- ggplot() + theme_void() + ggtitle("Domain UpSet (no domains present)")
    }
  } else {
    domain_upset <- ggplot() + theme_void() + ggtitle("Domain UpSet (domain list columns not found)")
  }
  if (!is.na(gene_col)) {
    # Binary matrix: genes x event types
    M <- dcast(unique(DT[!is.na(get(gene_col)) & nzchar(get(gene_col)),
                         .(gene = as.character(get(gene_col)), event_type, present = TRUE)]),
               gene ~ event_type, value.var = "present", fill = FALSE)
    evt <- setdiff(names(M), "gene")
    if (length(evt) >= 2) {
      B <- as.matrix(M[, ..evt]) * 1L

      # pairwise Jaccard between event types (columns)
      jacc <- matrix(NA_real_, nrow = length(evt), ncol = length(evt),
                     dimnames = list(evt, evt))
      for (i in seq_along(evt)) for (j in seq_along(evt)) {
        a <- B[, i] == 1; b <- B[, j] == 1
        u <- sum(a | b)
        jacc[i, j] <- if (u == 0) NA_real_ else sum(a & b)/u
      }

      # --- robust melt (no Var1/Var2/Freq assumptions) ---
      JT <- as.data.table(jacc)
      JT[, et1 := rownames(jacc)]
      JT <- melt(JT, id.vars = "et1", variable.name = "et2", value.name = "jaccard")

      coord_heatmap <- ggplot(JT, aes(et1, et2, fill = jaccard)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "white", high = "deeppink4", na.value = "grey95",
                            name = "Jaccard") +
        labs(
          x = "",  # not element_blank()
          y = "",
          title = "Event-type coordination across genes"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank()
        )

    } else {
      coord_heatmap <- ggplot() + theme_void() + ggtitle("Coordination heatmap (need >1 event types)")
    }
  } else {
    coord_heatmap <- ggplot() + theme_void() + ggtitle("Coordination heatmap (no gene column found)")
  }

  ru <- relative_use_pre_post(pre_filter_hits, hits)

  top_row    <- wrap_plots(list(p1, p2), ncol = 2)
  top2_row <- wrap_plots(list(p3, ppi1, ppi2), ncol = 3)
  middle_row <- wrap_plots(list(gene_upset, ru$plot), ncol = 2, widths = c(2.5, .5))
  bottom_row <- wrap_plots(list(domain_upset, coord_heatmap), ncol = 2, widths = c(2.5, .5))
  combined   <- (top_row / top2_row / middle_row / bottom_row) + plot_layout(heights = c(1, 1, 1, 1))

  return(list(
    summaries = list(
      by_type         = by_type[],
      class_counts    = class_counts[],
      score_summary   = score_summary[],
      domain_prevalence = dom_prev[],
      relative_use = ru$summary[]
    ),
    plot = combined
  ))
}

#' @title Relative pre- vs post-filter event retention
#'
#' @description
#' Computes how many events of each type remain after filtering and
#' visualizes post/pre ratios.
#' Used within [integrated_event_summary()] to display the fraction
#' of retained events after alignment and domain processing.
#'
#' @param pre_filter_hits `data.frame` or `data.table` of events prior to filtering.
#' @param hits `data.frame` or `data.table` of events after filtering or processing.
#'
#' @return A named list with:
#' \describe{
#'   \item{`summary`}{Data.table with `n_pre`, `n_post`, and `relative_use` by event type.}
#'   \item{`plot`}{Bar plot showing relative retention per event type.}
#' }
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_y_continuous labs theme_minimal
#'   theme element_text expansion
#' @importFrom scales percent_format
#' @keywords internal
#' @noRd
relative_use_pre_post <- function(
    pre_filter_hits,
    hits
){
  PRE  <- as.data.table(pre_filter_hits)
  POST <- as.data.table(hits)

  # --- count unique events pre ---
  pre_counts <- unique(PRE[, .(event_type = event_type,
                               event_id   = event_id)])[
                                 , .(n_pre = .N), by = event_type]

  # --- count unique events post ---
  post_counts <- unique(POST[, .(event_type = event_type_inc,
                                 event_id   = event_id)])[
                                   , .(n_post = .N), by = event_type]

  # --- merge & compute relative use ---
  RU <- merge(pre_counts, post_counts, by = "event_type", all.x = TRUE)
  RU[is.na(n_post), n_post := 0L]
  RU[, relative_use := ifelse(n_pre > 0, n_post / n_pre, NA_real_)]

  # order by relative use (desc) for a nice plot
  RU <- RU[order(-relative_use)]
  RU[, event_type := factor(event_type, levels = event_type)]

  # --- plot ---
  p <- ggplot(RU, aes(event_type, relative_use)) +
    geom_col(fill = 'deeppink4') +
    geom_text(aes(label = paste0(n_post, "/", n_pre)),
              vjust = -0.35, size = 3.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "Relative use (post / pre)",
         title = "Relative Use by Event Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(summary = RU[], plot = p)
}

#' Gene Set Enrichment for Splicing-Linked Gene Lists
#'
#' Perform over-representation analysis for a foreground gene set, optionally
#' against a background universe, using GO, MSigDB, and Reactome categories.
#'
#' This is a convenience wrapper around `clusterProfiler`, `msigdbr`, and
#' optionally `ReactomePA`, producing a combined enrichment table and a quick
#' visualization of top significant terms.
#'
#' @param foreground Character vector of gene IDs (symbols or Ensembl IDs).
#' @param background Optional character vector of background genes
#'   (universe). If `NULL`, all genes present in annotation collections are used.
#' @param species Species for enrichment catalog (`"human"` or `"mouse"`).
#' @param gene_id_type Type of input gene IDs (`"symbol"` or `"ensembl"`).
#' @param sources Character vector selecting enrichment sources.
#'   Example options include:
#'   \itemize{
#'     \item `"GO:BP"`, `"GO:MF"`, `"GO:CC"`
#'     \item `"MSigDB:H"`, `"MSigDB:C2:CP:REACTOME"`
#'   }
#' @param min_size Minimum term size (default `10`).
#' @param max_size Maximum term size (default `2000`).
#' @param p_adjust_cutoff FDR cutoff for reporting significant terms.
#' @param simplify_go Whether to apply GO term redundancy reduction.
#' @param top_n_plot Number of terms to visualize in the quick plot.
#' @param plot_type `"dot"` (default) or `"bar"`.
#'
#' @return A list with:
#' \describe{
#'   \item{results_per_source}{List of enrichment result tables per source}
#'   \item{results_combined}{Combined enrichment table}
#'   \item{results_signif}{Filtered table by FDR cutoff}
#'   \item{plot}{A `ggplot2` object visualizing top terms}
#' }
#'
#' @details
#' Input genes are internally mapped to Entrez IDs. Enrichment tests are performed using:
#'
#' * `clusterProfiler::enrichGO`
#' * `clusterProfiler::enricher` (MSigDB)
#' * `ReactomePA::enrichPathway` (optional)
#'
#' @note
#' Requires these packages installed:
#' `clusterProfiler`, `msigdbr`, `data.table`, `AnnotationDbi`, `ggplot2`,
#' and `org.Hs.eg.db` or `org.Mm.eg.db`. For Reactome analysis you must also
#' install `ReactomePA`.
#'
#' @import data.table
#' @importFrom ggplot2 ggplot aes geom_point geom_col coord_flip scale_color_gradient
#'   labs theme_minimal theme element_text theme_void
#' @importFrom data.table data.table rbindlist as.data.table
#' @importFrom stats reorder
#' @examples
#' \donttest{
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
#'
#' annotation_df <- get_annotation(load = "test")
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' matched <- get_matched_events_chunked(res, annotation_df$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annotation_df$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#' enrichment <- get_enrichment(res$gene_id, bg$gene_id, species = 'human', 'ensembl', 'MSigDB:H')
#' }
#' @export
get_enrichment <- function(
    foreground,
    background        = NULL,           # universe; if NULL, use all genes seen in term2gene
    species           = c("human","mouse"),
    gene_id_type      = c("symbol","ensembl"),
    sources           = c("GO:BP","GO:MF","GO:CC","MSigDB:H","MSigDB:C2:CP:REACTOME"),
    min_size          = 10,
    max_size          = 2000,
    p_adjust_cutoff   = 0.05,
    simplify_go       = TRUE,           # remove GO redundancy via clusterProfiler::simplify
    top_n_plot        = 20,             # how many terms to show in the quick plot
    plot_type         = c("dot","bar")  # "dot" or "bar"
){
  suppressPackageStartupMessages({
    req <- c("clusterProfiler","msigdbr","AnnotationDbi","data.table","ggplot2")
    missing <- req[!vapply(req, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
    if (length(missing)) stop("Please install: ", paste(missing, collapse=", "))
    if (match.arg(species) == "human" && !requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("Please install org.Hs.eg.db")
    if (match.arg(species) == "mouse" && !requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop("Please install org.Mm.eg.db")
  })

  species   <- match.arg(species)
  gene_id_type <- match.arg(gene_id_type)
  plot_type <- match.arg(plot_type)

  ## ---- 1) ID mapping helper (to ENTREZ) ----
  map_to_entrez <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[nzchar(ids)]
    if (gene_id_type == "ensembl") ids <- sub("\\.\\d+$","", ids)
    org_pkg <- if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
    keytype <- switch(gene_id_type,
                      symbol  = "SYMBOL",
                      ensembl = "ENSEMBL",
                      entrez  = "ENTREZID")

    Org <- get(org_pkg)
    tab <- AnnotationDbi::select(Org, keys = ids, columns = "ENTREZID",
                                 keytype = keytype)
    unique(na.omit(tab$ENTREZID))
  }

  fg_entrez <- map_to_entrez(foreground)
  bg_entrez <- if (is.null(background)) NULL else map_to_entrez(background)

  ## ---- 2) Build MSigDB term2gene for selected categories ----
  make_msigdb_t2g <- function(msig_key) {
    # msig_key examples: "H", "C2:CP", "C2:CP:REACTOME", "C5:BP", etc.
    sp_name <- if (species == "human") "Homo sapiens" else "Mus musculus"
    mdbr <- msigdbr::msigdbr(species = sp_name)
    parts <- strsplit(msig_key, ":", fixed = TRUE)[[1]]
    if (length(parts) == 1) {
      md <- mdbr[mdbr$gs_cat == parts[1], ]
    } else if (length(parts) == 2) {
      md <- mdbr[mdbr$gs_cat == parts[1] & mdbr$gs_subcat == parts[2], ]
    } else {
      md <- mdbr[mdbr$gs_cat == parts[1] &
                   mdbr$gs_subcat == paste(parts[2], parts[3], sep = ":"),
      ]
    }
    if (!nrow(md)) return(NULL)
    # map SYMBOL to ENTREZ using msigdbr's ENTREZID if present; fall back to mapping
    if (!is.null(md$entrez_gene)) {
      t2g <- data.table(term = md$gs_name, gene = as.character(md$entrez_gene))
    } else if (!is.null(md$entrez_gene_id)) {
      t2g <- data.table(term = md$gs_name, gene = as.character(md$entrez_gene_id))
    } else {
      # convert from SYMBOL
      sym2ent <- map_to_entrez(md$gene_symbol)
      t2g <- data.table(term = md$gs_name, gene = sym2ent)
    }
    t2g[!is.na(gene)]
  }

  ## ---- 3) Run enrichments per source ----
  res_list <- list()

  # GO (BP/MF/CC)
  if (any(grepl("^GO:", sources))) {

    ontos <- sub("^GO:","", sources[grepl("^GO:", sources)], fixed = FALSE)
    ontos <- intersect(ontos, c("BP","MF","CC"))
    Org <- if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
    OrgPkg <- get(Org)
    for (ont in ontos) {
      ego <- clusterProfiler::enrichGO(
        gene          = fg_entrez,
        universe      = bg_entrez,
        OrgDb         = OrgPkg,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1,
        minGSSize     = min_size,
        maxGSSize     = max_size,
        readable      = TRUE
      )
      if (simplify_go && length(ego) > 0) {
        ego <- tryCatch(clusterProfiler::simplify(ego), error = function(e) ego)
      }
      if (length(ego) > 0 && nrow(ego@result)) {
        out <- as.data.table(ego@result)
        out[, source := paste0("GO:", ont)]
        res_list[[paste0("GO_", ont)]] <- out
      }
    }
  }

  # MSigDB categories
  msig_keys <- sub("^MSigDB:","", sources[grepl("^MSigDB:", sources)])
  if (length(msig_keys)) {
    for (k in msig_keys) {
      t2g <- make_msigdb_t2g(k)
      if (!is.null(t2g) && nrow(t2g)) {
        er <- clusterProfiler::enricher(
          gene          = fg_entrez,
          universe      = bg_entrez,
          TERM2GENE     = t2g,
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,
          minGSSize     = min_size,
          maxGSSize     = max_size
        )
        if (length(er) > 0 && nrow(er@result)) {
          out <- as.data.table(er@result)
          out[, source := paste0("MSigDB:", k)]
          res_list[[paste0("MSigDB_", gsub("[:]", "_", k))]] <- out
        }
      }
    }
  }

  # Merge results and filter by FDR
  combined <- if (length(res_list)) rbindlist(res_list, use.names = TRUE, fill = TRUE) else data.table()
  if (nrow(combined)) {
    # standardize columns we'll use in plotting
    if (!"p.adjust" %in% names(combined) && "qvalue" %in% names(combined)) {
      combined[, p.adjust := qvalue]
    }
    combined <- combined[!is.na(p.adjust)]
    combined[, GeneRatio_num := {
      # GeneRatio comes like "k/n"
      if ("GeneRatio" %in% names(combined)) {
        vapply(
          strsplit(GeneRatio, "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]),
          FUN.VALUE = numeric(1)
        )
      } else {
        Count / as.numeric(BgRatio)
      }
    }]
    combined_sig <- combined[p.adjust <= p_adjust_cutoff]
  } else {
    combined_sig <- combined
  }

  ## ---- 4) Quick plot (top_n_plot by FDR across all sources) ----
  if (nrow(combined_sig)) {
    top <- combined_sig[order(p.adjust, -Count)][seq_len(min(top_n_plot, .N))]
    top[, pretty := paste0(source, " | ", Description)]
    if (plot_type == "dot") {
      plt <- ggplot(top, aes(x = GeneRatio_num, y = reorder(pretty, GeneRatio_num))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "cadetblue4", high = "deeppink4", name = expression(-log[10]~FDR)) +
        labs(x = "Gene ratio", y = NULL, title = "Gene set enrichment (top by FDR)") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 9))
    } else {
      plt <- ggplot(top, aes(x = reorder(pretty, -log10(p.adjust)), y = -log10(p.adjust), fill = source)) +
        geom_col() +
        coord_flip() +
        labs(x = NULL, y = expression(-log[10]~FDR), title = "Gene set enrichment (top by FDR)") +
        theme_minimal()
    }
  } else {
    plt <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No enriched terms at selected cutoff")
  }

  list(
    results_per_source = res_list,   # each enrichResult@result as data.table
    results_combined   = combined,   # all results (unfiltered)
    results_signif     = combined_sig, # FDR-filtered
    plot               = plt
  )
}

#' Extract Domain-Altering Genes for Enrichment
#'
#' Return genes whose inclusion/exclusion isoforms show unique protein domain gain or loss.
#'
#' @param hits Data frame from `get_domains()`.
#'
#' @return Character vector of gene IDs.
#' @export
#'
#' @examples
#' \donttest{
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
#' enrichment <- get_enrichment(get_domain_gene_for_enrichment(hits_domain), bg$gene_id, species = 'human', 'ensembl', 'MSigDB:H')
#' }
get_domain_gene_for_enrichment <- function(hits) {
  enriched_genes <- hits$gene_id_inc[hits$diff_n > 0]
  return(enriched_genes)
}


#' Extract Protein-Interaction-Affected Genes for Enrichment
#'
#' Return genes whose splicing affects known protein-protein interactions.
#'
#' @param hits Data frame with PPI annotation columns.
#'
#' @return Character vector of gene IDs.
#' @export
#'
#' @examples
#' \donttest{
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
#' ppidm <- get_ppidm(test=TRUE)
#'
#' ppi <- get_isoform_interactions(protein_feature_total[ensembl_transcript_id %in%
#'                                                         annotation_df$annotations[, unique(transcript_id)]],
#'                                 ppidm, save = FALSE, load_dir = '/projectnb2/evolution/zwakefield/proteinImpacts/', init = TRUE)
#'
#'
#' hits_final <- get_ppi_switches(hits_domain, ppi)
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#' enrichment <- get_enrichment(get_ppi_gene_enrichment(hits_final), bg$gene_id, species = 'human', 'ensembl', 'MSigDB:H')
#' }
get_ppi_gene_enrichment <- function(hits) {
  ppi_genes <- hits$gene_id_inc[hits$n_ppi > 0]
  return(ppi_genes)
}



#' Extract Differential-Inclusion Genes for Enrichment
#'
#' Select genes whose splicing passes significance and effect-size thresholds.
#'
#' @param hits Data frame output from differential inclusion testing.
#' @param padj_threshold FDR cutoff (default `0.05`).
#' @param delta_psi_threshold Absolute deltaPSI threshold (default `0.1`).
#'
#' @return Character vector of gene IDs.
#' @export
#'
#' @examples
#' \donttest{
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
#'
#' annotation_df <- get_annotation(load = "test")
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' matched <- get_matched_events_chunked(res, annotation_df$annotations, chunk_size = 2000)
#' x_seq <- attach_sequences(matched, annotation_df$sequences)
#' pairs <- get_pairs(x_seq, source="multi")
#' seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
#' hits_domain <- get_domains(seq_compare, exon_features)
#'
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annotation_df$annotations,
#'                      protein_features = protein_feature_total)
#' enrichment <- get_enrichment(get_di_gene_enrichment(res, .05, .1), bg$gene_id, species = 'human', 'ensembl', 'MSigDB:H')
#' }
get_di_gene_enrichment <- function(hits, padj_threshold, delta_psi_threshold) {
  di_genes <- hits$gene_id[hits$padj < padj_threshold & abs(hits$delta_psi) > delta_psi_threshold]
}



















