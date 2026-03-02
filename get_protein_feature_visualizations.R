#' Parse a genomic span string ("start-end")
#'
#' Internal helper to parse span strings of the form `"start-end"` into a
#' length-2 numeric vector `(min, max)`. Returns `NULL` for missing, malformed,
#' or non-finite values.
#'
#' @param s Character scalar span string like `"115046409-115057425"`.
#' @param numeric Logical; if `TRUE`, parse as `numeric` (via `as.numeric()`).
#'   If `FALSE` (default), parse as integer (via `as.integer()`).
#'
#' @return `NULL` if `s` is invalid; otherwise a length-2 vector `c(xmin, xmax)`.
#'
#' @keywords internal
#' @noRd
.parse_span <- function(s, numeric = FALSE) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NULL)
  xy <- strsplit(s, "-", fixed = TRUE)[[1]]
  if (length(xy) != 2) return(NULL)
  if (numeric) {
    x1 <- suppressWarnings(as.numeric(xy[1]))
    x2 <- suppressWarnings(as.numeric(xy[2]))
  } else {
    x1 <- suppressWarnings(as.integer(xy[1]))
    x2 <- suppressWarnings(as.integer(xy[2]))
  }
  if (!is.finite(x1) || !is.finite(x2)) return(NULL)
  c(min(x1, x2), max(x1, x2))
}

#' Parse protein feature genomic coordinates encoded in the `name` field
#'
#' Internal helper to parse protein feature coordinates from strings formatted as
#' `"<label>;chr:start-end"` in the `name` column (e.g., `"PF00069;chr7:123-456"`).
#' Filters to the requested chromosome and adds `xmin`/`xmax` genomic interval
#' columns along with parsed helper fields.
#'
#' @param pf A data.frame/data.table of protein features with at least a `name`
#'   column. Typically also contains fields such as `feature_id`, `database`,
#'   and `ensembl_transcript_id`.
#' @param chr Character scalar chromosome name to keep (e.g., `"chr7"`).
#'
#' @return A `data.table` with parsed columns (`tmp`, `pf_chr`, `pf_start`,
#'   `pf_end`) and genomic interval columns `xmin`/`xmax`. Returns an empty
#'   `data.table` (0 rows) if no matches remain after filtering.
#'
#' @keywords internal
#' @noRd
.parse_pf_coords <- function(pf, chr) {
  pf <- data.table::copy(data.table::as.data.table(pf))
  if (!nrow(pf)) return(pf[0])
  
  pf[, c("tmp", "coord") := data.table::tstrsplit(name, ";", fixed = TRUE)]
  pf[, c("pf_chr", "range") := data.table::tstrsplit(coord, ":", fixed = TRUE)]
  pf[, c("pf_start", "pf_end") := data.table::tstrsplit(range, "-", fixed = TRUE)]
  pf[, `:=`(pf_start = as.integer(pf_start),
            pf_end   = as.integer(pf_end))]
  
  pf <- pf[pf_chr == chr]
  if (!nrow(pf)) return(pf[0])
  
  pf[, `:=`(
    xmin = pmin(pf_start, pf_end),
    xmax = pmax(pf_start, pf_end)
  )]
  
  pf
}

#' Build a genomic-to-compact exon map for a transcript
#'
#' Internal helper to construct an exon concatenation map where exons are placed
#' end-to-end in transcript 5'->3' order, producing "compact" coordinates that
#' remove introns. The mapping differs by strand: negative-strand exons are
#' ordered by decreasing genomic coordinate.
#'
#' @param exons_genomic A data.frame/data.table with exon genomic intervals. Must
#'   contain `xmin` and `xmax` (or fields coercible to them).
#' @param strand Character scalar `"+"` or `"-"`.
#'
#' @return A list with:
#' \describe{
#'   \item{map}{`data.table` with columns `gxmin`, `gxmax`, `cxmin`, `cxmax`.}
#'   \item{total_len}{Integer total exonic length (sum of exon lengths).}
#' }
#'
#' @keywords internal
#' @noRd
.build_exon_map <- function(exons_genomic, strand) {
  ex <- data.table::copy(data.table::as.data.table(exons_genomic))
  ex[, `:=`(
    gxmin = as.integer(pmin(xmin, xmax)),
    gxmax = as.integer(pmax(xmin, xmax))
  )]
  
  if (strand == "-") {
    ex <- ex[order(-gxmax, -gxmin)]
  } else {
    ex <- ex[order(gxmin, gxmax)]
  }
  
  ex[, exon_len := gxmax - gxmin + 1L]
  ex[, cxmin := 1L + data.table::shift(cumsum(exon_len), fill = 0L)]
  ex[, cxmax := cxmin + exon_len - 1L]
  
  list(
    map = ex[, .(gxmin, gxmax, cxmin, cxmax)],
    total_len = ex[, sum(exon_len, na.rm = TRUE)]
  )
}

#' Project genomic segments onto compact (introns-removed) coordinates
#'
#' Internal helper to convert genomic intervals into compact exon-concatenated
#' coordinates by intersecting segments with exon intervals and mapping overlap
#' positions into the transcript coordinate system.
#'
#' The function uses `data.table::foverlaps()` to intersect intervals and then
#' maps each overlap into `cxmin`/`cxmax`. For negative-strand transcripts, the
#' projection reverses within each exon to preserve 5'->3' directionality in
#' compact space.
#'
#' @param segments_dt A data.frame/data.table containing genomic intervals with
#'   columns `xmin` and `xmax`. Additional columns are preserved.
#' @param exon_map A data.frame/data.table produced by `.build_exon_map()` (i.e.,
#'   columns `gxmin`, `gxmax`, `cxmin`, `cxmax`).
#' @param strand Character scalar `"+"` or `"-"`.
#'
#' @return A `data.table` of overlapped segments with compact coordinate columns
#'   `cxmin` and `cxmax`. Returns `NULL` if there are no overlaps.
#'
#' @keywords internal
#' @noRd
.project_to_compact <- function(segments_dt, exon_map, strand) {
  if (is.null(segments_dt) || !nrow(segments_dt)) return(NULL)
  
  seg <- data.table::copy(data.table::as.data.table(segments_dt))
  exm <- data.table::copy(data.table::as.data.table(exon_map))
  
  # avoid collisions
  dup <- intersect(names(seg), c("cxmin", "cxmax"))
  if (length(dup)) seg[, (dup) := NULL]
  
  seg[, `:=`(
    sxmin = as.integer(pmin(xmin, xmax)),
    sxmax = as.integer(pmax(xmin, xmax))
  )]
  
  exm <- exm[, .(gxmin, gxmax, ex_cxmin = cxmin, ex_cxmax = cxmax)]
  data.table::setkey(exm, gxmin, gxmax)
  data.table::setkey(seg, sxmin, sxmax)
  
  ov <- data.table::foverlaps(
    seg, exm,
    by.x = c("sxmin", "sxmax"),
    by.y = c("gxmin", "gxmax"),
    type = "any",
    nomatch = 0L
  )
  if (!nrow(ov)) return(NULL)
  
  ov[, `:=`(
    gx1 = pmax(sxmin, gxmin),
    gx2 = pmin(sxmax, gxmax)
  )]
  
  if (strand == "-") {
    ov[, `:=`(
      c1 = ex_cxmin + (gxmax - gx1),
      c2 = ex_cxmin + (gxmax - gx2)
    )]
  } else {
    ov[, `:=`(
      c1 = ex_cxmin + (gx1 - gxmin),
      c2 = ex_cxmin + (gx2 - gxmin)
    )]
  }
  
  ov[, `:=`(
    cxmin = pmin(c1, c2),
    cxmax = pmax(c1, c2)
  )]
  
  ov[, c("sxmin","sxmax","gx1","gx2","c1","c2","ex_cxmin","ex_cxmax") := NULL]
  
  if (anyDuplicated(names(ov))) {
    dupn <- unique(names(ov)[duplicated(names(ov))])
    stop("Internal error: duplicate names after projection: ", paste(dupn, collapse = ", "))
  }
  
  ov
}

# ==========================
# Unified transcript prep
# ==========================

#' Prepare one transcript's plotting tracks (exons, introns, domains)
#'
#' Internal helper that builds the data structures required for plotting a single
#' transcript with exon structure and protein feature tracks. Supports two
#' coordinate modes:
#'
#' - `mode = "genomic"`: plot in genomic coordinates (`xmin`/`xmax`), including introns.
#' - `mode = "compact"`: plot in compact exonic coordinates (`cxmin`/`cxmax`), removing introns.
#'
#' Protein features are filtered to a transcript and optionally filtered by
#' `feature_db`, then clipped to exon intervals (genomic mode) or projected into
#' compact space (compact mode). Domains are stacked on separate rows to avoid
#' overlap; optionally combine repeated instances into shared tracks.
#'
#' @param transcript Character scalar Ensembl transcript ID.
#' @param gtf_df A data.frame/data.table of GTF-like annotations containing
#'   `transcript_id`, `type` (with `"exon"`), `exon_number`, `chr`, `strand`,
#'   `start`, `end`.
#' @param protein_features A data.frame/data.table of protein features with at
#'   least columns `ensembl_transcript_id`, `name`, `feature_id`, `database`.
#' @param feature_db Optional character vector of databases to retain (filters
#'   `protein_features$database`).
#' @param y_offset Numeric y-location of the transcript backbone.
#' @param mode Coordinate mode, `"genomic"` or `"compact"`.
#' @param domain_base_gap Numeric vertical gap between backbone and first domain row.
#' @param domain_track_step Numeric vertical spacing between domain rows.
#' @param domain_label_dy Numeric vertical offset for domain labels relative to
#'   domain rectangles.
#' @param combine_like_domains Logical; if `TRUE`, domains with identical labels
#'   share a track; otherwise tracks are per-instance.
#'
#' @return A list with components:
#' \describe{
#'   \item{meta}{List with `transcript`, `chr`, `strand`, `mode`.}
#'   \item{exon_map}{`data.table` exon map for compact mode; `NULL` for genomic.}
#'   \item{exons}{`data.table` exon rectangles (genomic or compact columns).}
#'   \item{introns}{`data.table` intron segments (genomic mode only) or `NULL`.}
#'   \item{domains}{`data.table` of domain rectangles/labels (or `NULL`).}
#'   \item{instance_spans}{`data.table` span per domain instance for connectors (or `NULL`).}
#'   \item{label}{`data.table` for transcript label placement.}
#' }
#'
#' @keywords internal
#' @noRd
.prepare_one_transcript <- function(transcript,
                                    gtf_df,
                                    protein_features,
                                    feature_db = NULL,
                                    y_offset = 0,
                                    mode = c("genomic", "compact"),
                                    domain_base_gap = 0.85,
                                    domain_track_step = 0.6,
                                    domain_label_dy = 0.20,
                                    combine_like_domains = FALSE) {
  
  mode <- match.arg(mode)
  
  gtf <- data.table::as.data.table(gtf_df)[
    !is.na(transcript_id) & transcript_id == transcript & type == "exon"
  ]
  if (!nrow(gtf)) stop(paste("Transcript not found:", transcript))
  
  gtf <- gtf[order(as.integer(exon_number))]
  chr <- as.character(gtf$chr[1])
  strand <- as.character(gtf$strand[1])
  
  exons_gen <- gtf[, .(
    xmin = as.integer(pmin(start, end)),
    xmax = as.integer(pmax(start, end))
  )]
  
  introns <- NULL
  exon_map <- NULL
  
  if (mode == "genomic") {
    exons <- exons_gen[, .(transcript = transcript, xmin, xmax, y = y_offset)]
    
    n_ex <- nrow(exons)
    if (n_ex > 1L) {
      introns <- data.table::data.table(
        transcript = transcript,
        x    = exons$xmax[seq_len(n_ex - 1L)],
        xend = exons$xmin[2:n_ex],
        y    = y_offset
      )
    }
  } else {
    exmap_obj <- .build_exon_map(exons_gen, strand)
    exon_map <- exmap_obj$map
    
    exons <- exon_map[, .(
      transcript = transcript,
      gxmin = gxmin, gxmax = gxmax,
      cxmin = cxmin, cxmax = cxmax,
      y = y_offset
    )]
  }
  
  # --- protein features ---
  pf <- data.table::as.data.table(protein_features)[ensembl_transcript_id == transcript]
  if (!is.null(feature_db)) pf <- pf[database %chin% feature_db]
  
  domains <- NULL
  instance_spans <- NULL
  
  if (nrow(pf)) {
    pf[, pf_instance_id := .I]
    pf <- .parse_pf_coords(pf, chr)
    
    if (nrow(pf)) {
      pf2 <- pf[, .(
        xmin = as.integer(xmin),
        xmax = as.integer(xmax),
        domain_id = as.character(feature_id),
        domain_name = as.character(tmp),
        database = as.character(database),
        pf_instance_id = as.integer(pf_instance_id)
      )]
      
      if (mode == "genomic") {
        # clip to exon bounds in genomic space using foverlaps
        ex_g <- exons[, .(xmin, xmax)]
        data.table::setkey(ex_g, xmin, xmax)
        data.table::setkey(pf2, xmin, xmax)
        
        clipped <- data.table::foverlaps(
          pf2, ex_g,
          by.x = c("xmin", "xmax"),
          by.y = c("xmin", "xmax"),
          type = "any",
          nomatch = 0L
        )[, .(
          xmin = pmax(i.xmin, xmin),
          xmax = pmin(i.xmax, xmax),
          domain_id, domain_name, database, pf_instance_id
        )]
        
        if (nrow(clipped)) domains <- clipped
      } else {
        # project to compact coords (exonic)
        clipped <- .project_to_compact(pf2, exon_map, strand)
        if (!is.null(clipped) && nrow(clipped)) domains <- clipped
      }
      
      if (!is.null(domains) && nrow(domains)) {
        domains[, label := paste0(domain_id, " [", database, "] (", domain_name, ")")]
        
        if (isTRUE(combine_like_domains)) {
          domains[, track := as.integer(factor(label, levels = unique(label)))]
        } else {
          domains[, track := as.integer(factor(pf_instance_id, levels = unique(pf_instance_id)))]
        }
        
        domains[, transcript := transcript]
        domains[, y := y_offset - domain_base_gap - (track - 1) * domain_track_step]
        domains[, label_y := y + (0.1 + domain_label_dy)]
        
        # instance spans (for connectors)
        if (mode == "genomic") {
          instance_spans <- domains[, .(
            xmin = min(xmin, na.rm = TRUE),
            xmax = max(xmax, na.rm = TRUE),
            y    = y[1]
          ), by = .(pf_instance_id)]
        } else {
          instance_spans <- domains[, .(
            cxmin = min(cxmin, na.rm = TRUE),
            cxmax = max(cxmax, na.rm = TRUE),
            y     = y[1]
          ), by = .(pf_instance_id)]
        }
      }
    }
  }
  
  list(
    meta = list(transcript = transcript, chr = chr, strand = strand, mode = mode),
    exon_map = exon_map,        # only for compact
    exons = exons,
    introns = introns,          # only for genomic
    domains = domains,
    instance_spans = instance_spans,
    label = data.table::data.table(transcript = transcript, chr = chr, strand = strand, y = y_offset + 0.55)
  )
}


#' Construct highlight rectangles for inclusion/exclusion spans
#'
#' Internal helper to convert event span annotations (e.g., `inc_inc`, `exc_inc`,
#' `inc_exc`, `exc_exc`) into per-transcript highlight rectangles. In genomic
#' mode, spans are used directly. In compact mode, spans are projected to compact
#' coordinates using each transcript's exon map and strand.
#'
#' The returned object includes `ymin`/`ymax` derived from transcript-specific
#' vertical bounds (`yr`) so each transcript's highlight band is independent.
#'
#' @param hits_row A single-row data.frame/data.table describing an event with at
#'   least `transcript_id_inc`, `transcript_id_exc`, and up to four span columns
#'   (`inc_inc`, `exc_inc`, `inc_exc`, `exc_exc`) containing `"start-end"` strings.
#' @param yr A data.table giving per-transcript vertical bounds with columns
#'   `transcript`, `y_top`, `y_bottom`.
#' @param tx_info Named list keyed by transcript ID, each element containing
#'   `exon_map` (from `.build_exon_map()`) and `strand`. Used only in compact mode.
#' @param mode Coordinate mode, `"genomic"` or `"compact"`.
#' @param col_blue Fill color for inclusion-associated spans (default `"blue"`).
#' @param col_red Fill color for exclusion-associated spans (default `"red"`).
#' @param alpha Alpha transparency for highlight rectangles.
#'
#' @return A `data.table` of highlight rectangles with columns:
#'   `transcript`, `xmin/xmax` (genomic) or `cxmin/cxmax` (compact),
#'   `fill`, `ymin`, `ymax`, `alpha`.
#'   Returns `NULL` if no valid spans are present.
#'
#' @keywords internal
#' @noRd
.make_highlights <- function(hits_row, yr, tx_info,
                             mode = c("genomic", "compact"),
                             col_blue = "blue",
                             col_red  = "red",
                             alpha = 0.18) {
  
  mode <- match.arg(mode)
  r <- as.list(hits_row)
  
  spans <- data.table::rbindlist(list(
    {
      tx <- r$transcript_id_inc
      a <- .parse_span(r$inc_inc)
      b <- .parse_span(r$exc_inc)
      data.table::rbindlist(list(
        if (!is.null(a)) data.table::data.table(transcript=tx, xmin=a[1], xmax=a[2], fill=col_blue),
        if (!is.null(b)) data.table::data.table(transcript=tx, xmin=b[1], xmax=b[2], fill=col_red)
      ), fill=TRUE)
    },
    {
      tx <- r$transcript_id_exc
      a <- .parse_span(r$inc_exc)
      b <- .parse_span(r$exc_exc)
      data.table::rbindlist(list(
        if (!is.null(a)) data.table::data.table(transcript=tx, xmin=a[1], xmax=a[2], fill=col_blue),
        if (!is.null(b)) data.table::data.table(transcript=tx, xmin=b[1], xmax=b[2], fill=col_red)
      ), fill=TRUE)
    }
  ), fill=TRUE)
  
  if (!nrow(spans)) return(NULL)
  
  if (mode == "compact") {
    # project each transcript's spans to compact coords using its exon_map + strand
    out <- data.table::rbindlist(lapply(split(spans, spans$transcript), function(d) {
      tx <- unique(d$transcript)
      info <- tx_info[[tx]]
      if (is.null(info)) return(NULL)
      
      proj <- .project_to_compact(d[, .(xmin, xmax, transcript, fill)], info$exon_map, info$strand)
      if (is.null(proj) || !nrow(proj)) return(NULL)
      
      # keep fill per-row (not just first)
      proj[, fill := d$fill[match(paste0(proj$xmin, "-", proj$xmax), paste0(d$xmin, "-", d$xmax))]]
      proj[, transcript := tx]
      proj
    }), fill = TRUE)
    
    if (is.null(out) || !nrow(out)) return(NULL)
    
    out <- merge(out, yr[, .(transcript, y_top, y_bottom)], by = "transcript", all.x = TRUE)
    out[, `:=`(ymin = y_bottom, ymax = y_top, alpha = alpha)]
    out[is.finite(ymin) & is.finite(ymax)]
  } else {
    # genomic highlights just use xmin/xmax directly
    out <- merge(spans, yr[, .(transcript, y_top, y_bottom)], by = "transcript", all.x = TRUE)
    out[, `:=`(ymin = y_bottom, ymax = y_top, alpha = alpha)]
    out[is.finite(ymin) & is.finite(ymax)]
  }
}

#' Plot two transcripts with exon structure and protein feature tracks
#'
#' Internal workhorse to visualize two transcripts, their exon structures, and
#' protein feature segments (e.g., Pfam/ELM/SEG/InterPro), optionally highlighting
#' event spans (e.g., inclusion/exclusion regions).
#'
#' Two coordinate systems are supported:
#' \describe{
#'   \item{`mode = "genomic"`}{Genomic x-axis (`xmin`/`xmax`) with introns drawn as
#'   connecting segments. If both transcripts are negative strand, the x-axis is
#'   reversed to read 5'->3' left-to-right.}
#'   \item{`mode = "compact"`}{Compact exon-concatenated x-axis (`cxmin`/`cxmax`)
#'   with introns removed.}
#' }
#'
#' Domain tracks are stacked vertically, and repeated domain instances can be
#' combined onto shared tracks via `combine_domains`.
#'
#' @param transcripts Character vector of length 2 with transcript IDs.
#' @param gtf_df A data.frame/data.table of GTF-like exon annotations.
#' @param protein_features A data.frame/data.table of protein features.
#' @param feature_db Optional character vector of feature databases to include.
#' @param wrap_width Integer width for wrapping long domain labels.
#' @param highlight_hits Optional table of event rows used for highlighting.
#' @param highlight_event_id Optional event ID to select from `highlight_hits`.
#' @param highlight_alpha Alpha for highlight bands.
#' @param highlight_box Logical; if `TRUE`, draw dashed vertical bounds around
#'   highlighted spans.
#' @param highlight_box_pad_frac Fraction of plotted x-range used to pad the
#'   highlight bounding box.
#' @param highlight_box_lwd Line width for highlight bounding lines.
#' @param combine_domains Logical; if `TRUE`, combine identical domain labels
#'   onto shared tracks.
#' @param mode Coordinate mode: `"genomic"` or `"compact"`.
#' @param domain_base_gap Vertical gap between transcript backbone and domain tracks.
#' @param domain_track_step Vertical spacing between domain rows.
#' @param domain_label_dy Vertical offset for domain labels.
#'
#' @return A `ggplot` object.
#'
#' @seealso \code{\link{plot_two_transcripts_with_domains_unified}} for the
#'   exported wrapper.
#'
#' @keywords internal
#' @importFrom ggplot2 ggplot geom_rect geom_segment geom_text geom_vline coord_cartesian
#' @importFrom ggplot2 scale_alpha_identity scale_x_reverse scale_x_continuous labs theme_minimal theme
#' @importFrom ggplot2 element_blank margin
#' @importFrom data.table as.data.table copy data.table rbindlist setkey foverlaps tstrsplit shift
#' @noRd
plot_two_transcripts_with_features <- function(transcripts,
                                               gtf_df,
                                               protein_features,
                                               feature_db = NULL,
                                               wrap_width = 55,
                                               highlight_hits = NULL,
                                               highlight_event_id = NULL,
                                               highlight_alpha = 0.30,
                                               highlight_box = TRUE,
                                               highlight_box_pad_frac = 0.006,
                                               highlight_box_lwd = 0.45,
                                               combine_domains = FALSE,
                                               mode = c("genomic", "compact"),
                                               domain_base_gap = 0.85,
                                               domain_track_step = 0.6,
                                               domain_label_dy = 0.20) {
  
  mode <- match.arg(mode)
  stopifnot(length(transcripts) == 2)
  
  # prep transcript 1
  y1 <- 0
  t1 <- .prepare_one_transcript(
    transcripts[1], gtf_df, protein_features, feature_db,
    y_offset = y1,
    mode = mode,
    domain_base_gap = domain_base_gap,
    domain_track_step = domain_track_step,
    domain_label_dy = domain_label_dy,
    combine_like_domains = combine_domains
  )
  
  # dynamic y2
  t1_bottom <- if (!is.null(t1$domains) && nrow(t1$domains)) min(t1$domains$y - 0.45, na.rm = TRUE) else y1 - 0.6
  y2 <- t1_bottom - 1.2
  
  t2 <- .prepare_one_transcript(
    transcripts[2], gtf_df, protein_features, feature_db,
    y_offset = y2,
    mode = mode,
    domain_base_gap = domain_base_gap,
    domain_track_step = domain_track_step,
    domain_label_dy = domain_label_dy,
    combine_like_domains = combine_domains
  )
  
  no_domains <- (is.null(t1$domains) || nrow(t1$domains) == 0L) &&
    (is.null(t2$domains) || nrow(t2$domains) == 0L)
  
  reverse_x <- FALSE
  if (mode == "genomic") {
    strands <- unique(c(t1$meta$strand, t2$meta$strand))
    reverse_x <- (length(strands) == 1L && strands[1] == "-")
  }
  
  exons   <- data.table::rbindlist(list(t1$exons, t2$exons), fill = TRUE)
  introns <- data.table::rbindlist(list(t1$introns, t2$introns), fill = TRUE)
  domains <- data.table::rbindlist(list(t1$domains, t2$domains), fill = TRUE)
  labels  <- data.table::rbindlist(list(t1$label, t2$label), fill = TRUE)
  instance_spans <- data.table::rbindlist(list(t1$instance_spans, t2$instance_spans), fill = TRUE)
  
  if (no_domains) {
    domains <- data.table::as.data.table(t1$domains)[0]  # preserves schema if t1$domains exists
    if (is.null(t1$domains) || !nrow(t1$domains)) {
      # fallback schema to avoid missing columns later
      domains <- data.table::data.table()
    }
    instance_spans <- data.table::data.table()
  } else {
    domains <- data.table::rbindlist(list(t1$domains, t2$domains), fill = TRUE)
    instance_spans <- data.table::rbindlist(list(t1$instance_spans, t2$instance_spans), fill = TRUE)
  }
  
  # yr per transcript (top/bottom) – same logic for both; depends on mode columns
  if (mode == "genomic") {
    yr <- data.table::rbindlist(list(
      exons[, .(transcript, y_top = y + 0.85, y_bottom = y - 0.35)],
      if (nrow(domains)) domains[, .(transcript, y_bottom = min(y - 0.35, na.rm = TRUE)), by = transcript] else NULL
    ), fill = TRUE)
  } else {
    yr <- data.table::rbindlist(list(
      exons[, .(transcript, y_top = y + 0.85, y_bottom = y - 0.35)],
      if (nrow(domains)) domains[, .(transcript, y_bottom = min(y - 0.35, na.rm = TRUE)), by = transcript] else NULL
    ), fill = TRUE)
  }
  yr <- yr[, .(y_top = max(y_top, na.rm = TRUE), y_bottom = min(y_bottom, na.rm = TRUE)), by = transcript]
  
  # build tx_info for compact projection
  tx_info <- setNames(
    list(
      list(exon_map = t1$exon_map, strand = t1$meta$strand),
      list(exon_map = t2$exon_map, strand = t2$meta$strand)
    ),
    transcripts
  )
  
  hl <- NULL
  event_type_label <- NULL
  if (!is.null(highlight_hits) && !is.null(highlight_event_id)) {
    row <- data.table::as.data.table(highlight_hits)[event_id == highlight_event_id][1]
    if (nrow(row)) {
      hl <- .make_highlights(row, yr, tx_info, mode = mode, alpha = highlight_alpha)
      event_type_label <- row$event_type_exc[1]
    }
  }
  
  # x column names based on mode
  x1 <- if (mode == "genomic") "xmin" else "cxmin"
  x2 <- if (mode == "genomic") "xmax" else "cxmax"
  
  # highlight bounding lines
  vlines <- NULL
  if (!is.null(hl) && nrow(hl) && isTRUE(highlight_box)) {
    xr_all <- range(exons[, c(get(x1), get(x2))], na.rm = TRUE)
    pad <- diff(xr_all) * highlight_box_pad_frac
    left  <- min(hl[[x1]], na.rm = TRUE) - pad
    right <- max(hl[[x2]], na.rm = TRUE) + pad
    vlines <- data.table::data.table(x = c(left, right))
  }
  
  # xlims and label position
  xr <- range(exons[, c(get(x1), get(x2))], na.rm = TRUE)
  dx <- diff(xr)
  xlim <- c(xr[1] - 0.03 * dx, xr[2] + 0.02 * dx)
  if (!reverse_x) {
    x_label <- xlim[1] + 0.005 * diff(xlim)  # near left edge (min x)
  } else {
    x_label <- xlim[2] - 0.005 * diff(xlim)  # near left edge after reverse (max x)
  }
  
  labels[, text := paste0(transcript, "  (", chr, ", strand ", strand, ")")]
  
  # wrap domain labels
  if (!no_domains) {
    domains[, label_wrapped := vapply(
      label,
      function(s) paste(strwrap(s, width = wrap_width), collapse = "\n"),
      character(1)
    )]
  } else {
    domains <- data.table::data.table()
    domains[, label_wrapped := character()]
  }
  
  # y limits
  y_exon_min <- min(exons$y - 0.25, na.rm = TRUE)
  y_exon_max <- max(exons$y + 0.80, na.rm = TRUE)
  y_dom_min <- if (!no_domains) min(domains$y - 0.35, na.rm = TRUE) else Inf
  y_dom_max <- if (!no_domains) max(domains$label_y + 0.35, na.rm = TRUE) else -Inf
  ylim <- c(min(y_exon_min, y_dom_min), max(y_exon_max, y_dom_max))
  
  event_label_x <- mean(xlim)
  event_label_y <- ylim[2] - 0.10
  
  # connectors (only really meaningful in genomic mode; but you can keep for compact too)
  show_instance_connectors <- TRUE
  instance_connector_dy <- 0.0
  instance_connector_lwd <- 0.9
  instance_connector_alpha <- 0.55
  
  p <- ggplot2::ggplot()
  
  # highlights
  if (!is.null(hl) && nrow(hl)) {
    p <- p + ggplot2::geom_rect(
      data = hl,
      ggplot2::aes(
        xmin = .data[[x1]], xmax = .data[[x2]],
        ymin = ymin, ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = hl$fill,
      alpha = hl$alpha,
      color = NA
    )
  }
  
  # dashed bounds
  if (!is.null(vlines) && nrow(vlines)) {
    p <- p + ggplot2::geom_vline(
      data = vlines,
      ggplot2::aes(xintercept = x),
      linetype = "dashed",
      linewidth = highlight_box_lwd,
      color = "grey35"
    )
  }
  
  # event label
  if (!is.null(event_type_label)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        x = event_label_x, y = event_label_y,
        label = paste0("Event type: ", event_type_label)
      ),
      hjust = 0, vjust = 1, fontface = "bold", size = 3.6
    )
  }
  
  # transcript labels
  p <- p + ggplot2::geom_text(
    data = labels,
    ggplot2::aes(x = x_label, y = y, label = text),
    hjust = 0, size = 3.6, fontface = "bold"
  )

  # transcript structure
  if (mode == "genomic") {
    p <- p +
      ggplot2::geom_segment(
        data = introns,
        ggplot2::aes(x = x, xend = xend, y = y, yend = y),
        linewidth = 0.45,
        color = "grey55"
      ) +
      ggplot2::geom_rect(
        data = exons,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - 0.20, ymax = y + 0.20),
        fill = "black",
        linewidth = 0.2
      )
  } else {
    # compact exon rectangles (no introns)
    # Optional: exon alternating alpha (like your compact version)
    exons[, exon_idx := seq_len(.N), by = transcript]
    exons[, exon_alpha := ifelse(exon_idx %% 2L == 0L, 1.0, 0.65)]
    exons[, `:=`(cxmin_eps = cxmin - 0.5, cxmax_eps = cxmax + 0.5)]
  p <- p +
    ggplot2::geom_rect(
      data = exons,
      ggplot2::aes(xmin = cxmin_eps, xmax = cxmax_eps, ymin = y - 0.20, ymax = y + 0.20, alpha = exon_alpha),
      fill = "black",
      linewidth = 0
    ) +
    ggplot2::scale_alpha_identity()
  }

  # instance connectors
  if (!no_domains && isTRUE(show_instance_connectors) && nrow(instance_spans)) {
    if (mode == "genomic") {
      p <- p +
        ggplot2::geom_segment(
          data = instance_spans,
          ggplot2::aes(x = xmin, xend = xmax, y = y + instance_connector_dy, yend = y + instance_connector_dy),
          inherit.aes = FALSE,
          linewidth = instance_connector_lwd,
          alpha = instance_connector_alpha
        )
    } else {
      p <- p +
        ggplot2::geom_segment(
          data = instance_spans,
          ggplot2::aes(x = cxmin, xend = cxmax, y = y + instance_connector_dy, yend = y + instance_connector_dy),
          inherit.aes = FALSE,
          linewidth = instance_connector_lwd,
          alpha = instance_connector_alpha
        )
    }
  }
  
  # protein features
  
  if (!no_domains) {
    if (mode == "compact") {
      domains[, `:=`(cxmin_eps = cxmin - 0.5, cxmax_eps = cxmax + 0.5)]
    }
    dom_xmin <- if (mode == "compact") "cxmin_eps" else x1
    dom_xmax <- if (mode == "compact") "cxmax_eps" else x2

    p <- p +
      ggplot2::geom_rect(
        data = domains,
        ggplot2::aes(
          xmin = .data[[dom_xmin]], xmax = .data[[dom_xmax]],
          ymin = y - 0.18, ymax = y + 0.18,
          fill = label
        ),
        alpha = 0.92,
        linewidth = 0
      ) 
  
    dom_label_dt <- if (nrow(domains)) {
      if (!reverse_x) {
        domains[, .SD[which.min(get(x1))], by = .(label, label_wrapped, y, label_y)]
      } else {
        domains[, .SD[which.max(get(x2))], by = .(label, label_wrapped, y, label_y)]
      }
    } else {
      domains
    }
    
    p <- p + ggplot2::geom_text(
      data = dom_label_dt,
      ggplot2::aes(
        x = if (reverse_x) .data[[x2]] else .data[[x1]],
        y = label_y,
        label = label_wrapped
      ),
      hjust = if (reverse_x) 1 else 0,
      size = 2.6,
      lineheight = 0.95
    )
  }

  p <- p + ggplot2::coord_cartesian(ylim = ylim, expand = FALSE, clip = "on")
  
  if (reverse_x) {
    p <- p + ggplot2::scale_x_reverse(limits = xlim)
  } else {
    p <- p + ggplot2::scale_x_continuous(limits = xlim)
  }
  p <- p + ggplot2::labs(
    x = if (mode == "genomic") "Genomic coordinate" else "Exonic coordinate (introns removed)",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(10, 12, 10, 12)
  )

p
}


#' Plot two transcripts with exon structure and protein feature tracks (unified view)
#'
#' Visualize **two Ensembl transcripts** side-by-side with their exon structures and
#' optional **protein feature/domain tracks** (e.g., InterPro, Pfam, ELM, SEG, SignalP),
#' using a single entrypoint. The plot can be rendered in either:
#'
#' - **Transcript view** (`view = "transcript"`): genomic x-axis with introns drawn and,
#'   when both transcripts are negative-strand, a strand-aware x-axis reversal so the
#'   display reads left-to-right in **5'→3'** direction.
#' - **Protein view** (`view = "protein"`): compact, intron-free x-axis where exons are
#'   concatenated end-to-end (exonic coordinates), making it easier to compare protein
#'   feature locations across isoforms without large genomic intron gaps.
#'
#' Optionally, event-specific genomic spans (e.g., inclusion/exclusion regions) can be
#' overlaid as translucent highlight bands per transcript.
#'
#' @details
#' **What gets drawn**
#' \itemize{
#'   \item Exons as black rectangles (alternating alpha in protein/compact view).
#'   \item Introns as grey connector segments (transcript/genomic view only).
#'   \item Protein features as stacked rectangles under each transcript (if present).
#'   \item Domain labels anchored to the left edge (or right when the x-axis is reversed).
#'   \item Optional highlighted spans with dashed bounding lines.
#' }
#'
#' **Protein/domain tracks**
#' Protein features are filtered to the two transcripts and optionally filtered by
#' \code{feature_db}. Features are clipped to exons in transcript/genomic view, and
#' projected into compact/exonic coordinates in protein view.
#'
#' \itemize{
#'   \item If \code{combine_domains = TRUE}, identical domain labels share a common
#'   vertical track to reduce redundancy.
#'   \item If \code{combine_domains = FALSE}, each feature instance is assigned its own
#'   track (potentially more vertical space, but preserves instance-level separation).
#' }
#'
#' **Event highlighting**
#' If \code{highlight_hits} and \code{highlight_event_id} are provided, event spans are
#' parsed from \code{inc_inc}, \code{exc_inc}, \code{inc_exc}, \code{exc_exc} columns.
#' In transcript/genomic view spans are used directly; in protein/compact view spans are
#' projected into compact coordinates per transcript using exon maps.
#'
#' @param ... Additional arguments forwarded to the internal workhorse
#'   \code{plot_two_transcripts_with_features()}. Common arguments include:
#'   \describe{
#'     \item{transcripts}{Character vector of length 2 of Ensembl transcript IDs.}
#'     \item{gtf_df}{GTF-like exon annotation table containing at least
#'       \code{transcript_id}, \code{type=="exon"}, \code{exon_number}, \code{chr},
#'       \code{strand}, \code{start}, \code{end}.}
#'     \item{protein_features}{Protein feature table with at least
#'       \code{ensembl_transcript_id}, \code{name}, \code{feature_id}, \code{database}.}
#'     \item{feature_db}{Optional character vector of databases to retain (e.g.,
#'       \code{c("interpro","pfam","elm","seg","signalp")}).}
#'     \item{wrap_width}{Integer; width for wrapping long domain labels.}
#'     \item{highlight_hits}{Optional data.frame/data.table of event rows used for highlighting.}
#'     \item{highlight_event_id}{Optional event ID to select from \code{highlight_hits}.}
#'     \item{highlight_alpha}{Alpha transparency for highlight bands.}
#'     \item{highlight_box}{Logical; draw dashed vertical bounds around highlighted spans.}
#'     \item{highlight_box_pad_frac}{Fraction of total x-range used to pad highlight bounds.}
#'     \item{highlight_box_lwd}{Line width for highlight bounding lines.}
#'     \item{combine_domains}{Logical; combine identical domain labels onto shared tracks.}
#'     \item{domain_base_gap}{Vertical gap between transcript backbone and first domain track.}
#'     \item{domain_track_step}{Vertical spacing between stacked domain tracks.}
#'     \item{domain_label_dy}{Vertical offset used to place domain labels.}
#'   }
#' @param view Character scalar selecting the visualization coordinate system:
#'   \code{"transcript"} (genomic/intron-aware) or \code{"protein"} (compact/exonic).
#'
#' @return A \code{ggplot} object, or \code{NULL} if no drawable content is available
#'   (e.g., when both transcripts have zero features and you have configured internal
#'   logic to early-return on missing domains).
#'
#' @section Expected input formats:
#' \subsection{protein_features \code{name} column}{
#' The internal parser expects feature genomic coordinates encoded in the \code{name}
#' field formatted as \code{"<label>;chr:start-end"} (e.g. \code{"PF00069;chr7:123-456"}).
#' }
#' \subsection{highlight_hits span columns}{
#' Span columns are expected as strings \code{"start-end"} using genomic coordinates.
#' Missing spans should be \code{NA_character_}.
#' }
#'
#' @section Highlight input (\code{custom_hits_domain}) example:
#' The \code{highlight_hits} object is expected to be a table (data.frame/data.table)
#' with at least one row per event. Use \code{highlight_event_id} to select the row to plot.
#' Required columns are \code{event_id}, \code{event_type_exc}, \code{transcript_id_inc},
#' \code{transcript_id_exc}. Span columns may include \code{inc_inc}, \code{exc_inc},
#' \code{inc_exc}, \code{exc_exc} and should be genomic coordinate strings of the form
#' \code{"start-end"} (or \code{NA_character_} when absent).
#'
#' \preformatted{
#' custom_hits_domain <- data.table::data.table(
#'   event_id = event:n,
#'   event_type_exc = event,
#'   transcript_id_inc = transcript_id,
#'   transcript_id_exc = transcript_id,
#'   inc_inc = inc_inc,
#'   inc_exc = inc_exc,
#'   exc_inc = exc_inc,
#'   exc_exc = exc_exc
#' )
#' }
#' @examples
#' \dontrun{
#' # Example highlight row (skipped exon / SE)
#' custom_hits_domain <- data.table::data.table(
#'   event_id = "SE:1",
#'   event_type_exc = "SE",
#'   transcript_id_inc = "ENST00000350763",
#'   transcript_id_exc = "ENST00000542877",
#'   inc_inc = "115046409-115057425",
#'   inc_exc = "115046409-115057425",
#'   exc_inc = "115048259-115048532",
#'   exc_exc = NA_character_
#' )
#'
#' # Transcript (genomic) view: introns included, strand-aware axis
#' p_tx <- plot_two_transcripts_with_domains_unified(
#'   transcripts = c("ENST00000350763", "ENST00000542877"),
#'   gtf_df = annotation_df$annotations,
#'   protein_features = protein_feature_total,
#'   feature_db = c("interpro", "pfam"),
#'   highlight_hits = custom_hits_domain,
#'   highlight_event_id = "SE:1",
#'   combine_domains = FALSE,
#'   view = "transcript"
#' )
#'
#' # Protein (compact) view: introns removed
#' p_prot <- plot_two_transcripts_with_domains_unified(
#'   transcripts = c("ENST00000350763", "ENST00000542877"),
#'   gtf_df = annotation_df$annotations,
#'   protein_features = protein_feature_total,
#'   feature_db = c("interpro", "pfam", "elm", "seg"),
#'   highlight_hits = custom_hits_domain,
#'   highlight_event_id = "SE:1",
#'   combine_domains = TRUE,
#'   view = "protein"
#' )
#' }
#'
#' @seealso \code{\link[ggplot2]{ggplot}} for rendering and theming.
#'
#' @importFrom ggplot2 ggplot geom_rect geom_segment geom_text geom_vline
#' @importFrom ggplot2 coord_cartesian scale_x_reverse scale_x_continuous
#' @importFrom ggplot2 scale_alpha_identity labs theme_minimal theme element_blank margin
#' @importFrom data.table as.data.table copy rbindlist data.table tstrsplit
#' @importFrom data.table setkey foverlaps shift
#'
#' @export
plot_two_transcripts_with_domains_unified <- function(...,
                                                      view = c("transcript", "protein")) {
  view <- match.arg(view)
  mode <- if (view == "transcript") "genomic" else "compact"
  plot_two_transcripts_with_features(..., mode = mode)
}

# 
p1 <- plot_two_transcripts_with_domains_unified(
  transcripts = c("ENST00000337907", "ENST00000400908"),
  gtf_df = annotation_df$annotations,
  protein_features = protein_feature_total,
  highlight_hits = custom_hits_domain,
  feature_db = c("interpro"),
  highlight_event_id = "SE:1",
  combine_domains = TRUE,
  view = "transcript"
)
