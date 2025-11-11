#' Quasi-binomial GLM with Cook's Distance Filtering
#'
#' Fits a quasi-binomial GLM to estimate deltaPSI between two conditions.
#' Outliers are removed based on Cook's distance before refitting.
#'
#' @param d A `data.frame` or `data.table` with columns `psi_adj`, `psi_raw`,
#'   `condition`, and `total`.
#' @param cooks_cutoff Numeric. Cook's distance cutoff, passed to
#'   [cutoff_num()].
#'
#' @return A `data.table` with columns `p.value`, `cooks_max`, `n`, `n_used`,
#'   `mean_psi_ctrl`, `mean_psi_case`, and `delta_psi`.
#'
#' @importFrom data.table as.data.table uniqueN setDTthreads getDTthreads
#' @importFrom stats glm anova cooks.distance
#' @keywords internal
.site_glm <- function(d, cooks_cutoff) {
  # avoid thread contention inside workers
  old <- data.table::getDTthreads(); data.table::setDTthreads(1L)
  on.exit(data.table::setDTthreads(old), add = TRUE)

  d <- as.data.table(d)
  d <- d[is.finite(total) & total > 0]
  if (nrow(d) < 2L) {
    return(data.table(p.value=NA_real_, cooks_max=NA_real_, n=as.integer(nrow(d)), n_used=NA_integer_,
                      mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_))
  }
  d[, condition := droplevels(condition)]
  if (data.table::uniqueN(d$condition) < 2L) {
    return(data.table(p.value=NA_real_, cooks_max=NA_real_, n=as.integer(nrow(d)), n_used=NA_integer_,
                      mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_))
  }

  full0 <- try(suppressWarnings(glm(psi_adj ~ condition,
                                    family = quasibinomial(),
                                    weights = total, data = d)), silent = TRUE)
  if (inherits(full0, "try-error")) {
    return(data.table(p.value=NA_real_, cooks_max=NA_real_, n=as.integer(nrow(d)), n_used=NA_integer_,
                      mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_))
  }

  thr <- if (nrow(d) > 10L) cutoff_num(nrow(d), cooks_cutoff) else Inf
  cd  <- try(suppressWarnings(cooks.distance(full0)), silent = TRUE)
  cd  <- if (inherits(cd, "try-error")) rep(NA_real_, nrow(d)) else as.numeric(cd)
  keep <- ifelse(is.na(cd), TRUE, cd <= thr)

  fit_once <- function(dd) {
    full <- try(suppressWarnings(glm(psi_adj ~ condition,
                                     family = quasibinomial(),
                                     weights = total, data = dd)), silent = TRUE)
    red  <- try(suppressWarnings(glm(psi_adj ~ 1,
                                     family = quasibinomial(),
                                     weights = total, data = dd)), silent = TRUE)
    p    <- tryCatch({
      aa <- anova(red, full, test = "F")
      if (nrow(aa) >= 2L) as.numeric(aa$`Pr(>F)`[2]) else NA_real_
    }, error = function(e) NA_real_)

    means <- dd[, .(psi_mean = mean(psi_raw, na.rm = TRUE)), by = condition]
    levs  <- levels(dd$condition)
    ctrlL <- levs[1]
    caseL <- levs[min(2L, length(levs))]
    ctrl_mean <- means[condition == ctrlL, psi_mean][1]
    case_mean <- means[condition == caseL, psi_mean][1]

    data.table(
      p.value       = p,
      cooks_max     = suppressWarnings(max(cd, na.rm = TRUE)),
      n             = as.integer(nrow(d)),
      n_used        = as.integer(nrow(dd)),
      mean_psi_ctrl = ifelse(is.finite(ctrl_mean), ctrl_mean, NA_real_),
      mean_psi_case = ifelse(is.finite(case_mean), case_mean, NA_real_),
      delta_psi     = ifelse(all(is.finite(c(ctrl_mean, case_mean))),
                             case_mean - ctrl_mean, NA_real_)
    )
  }

  if (any(!keep) && sum(keep) >= 2L && data.table::uniqueN(d$condition[keep]) >= 2L) {
    fit_once(d[keep])
  } else {
    fit_once(d)
  }
}


#' Parallel Quasi-binomial GLM Fitting Across Sites
#'
#' Fits quasi-binomial GLMs per splicing site (or event) in parallel,
#' splitting data into chunks to limit overhead. Each site-level fit
#' is performed by [`.site_glm()`].
#'
#' @param x A `data.frame` or `data.table` with columns:
#'   `site_id`, `condition`, `psi_adj`, `psi_raw`, and `total`.
#' @param cores Integer. Number of parallel workers (default 8).
#' @param chunk_size Integer. Approximate number of sites per chunk
#'   to send to each worker (default 2000).
#' @param progress Logical. Show progress bar (default `TRUE`).
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @param cooks_cutoff Numeric or character. Cook's distance cutoff,
#'   passed to [cutoff_num()].
#'
#' @return A `data.table` with one row per site and columns:
#'   `site_id`, `p.value`, `cooks_max`, `n`, `n_used`,
#'   `mean_psi_ctrl`, `mean_psi_case`, and `delta_psi`.
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom data.table as.data.table rbindlist setcolorder
#' @keywords internal
fit_sites_parallel <- function(x,
                               cores       = 8L,
                               chunk_size  = 2000L,
                               progress    = TRUE,
                               verbose     = TRUE,
                               cooks_cutoff) {

  if (verbose) cat("[PROCESSING] Fitting quasi-binomial GLMs per site...\n")

  x <- as.data.table(x)
  need <- c("site_id","condition","psi_adj","psi_raw","total")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  # Normalize lightweight columns once
  x[, condition := factor(condition)]
  x <- x[, ..need]

  # Build row-index list per site (very light to serialize)
  idx_list <- split(seq_len(nrow(x)), x$site_id)
  site_keys <- names(idx_list)

  # Chunk the site list to reduce task overhead
  n_sites  <- length(idx_list)
  n_chunks <- max(1L, ceiling(n_sites / chunk_size))
  chunk_bounds <- split(seq_len(n_sites), ceiling(seq_len(n_sites)/chunk_size))
  chunks <- lapply(chunk_bounds, function(ix) idx_list[ix])
  chunk_names <- lapply(chunk_bounds, function(ix) site_keys[ix])

  # Worker for a chunk (= many sites)
  chunk_worker <- function(ch, ch_names, cooks_cutoff) {
    # Evaluate each site in the chunk and stitch
    out <- vector("list", length(ch))
    for (j in seq_along(ch)) {
      res <- .site_glm(x[ch[[j]]], cooks_cutoff)
      # add site_id here to avoid another pass
      res[, site_id := ch_names[j]]
      out[[j]] <- res
    }
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }

  param <- BiocParallel::MulticoreParam(workers = cores, progressbar = progress)


  # Parallel over chunks with a progress bar
  res_chunks <- BiocParallel::bplapply(seq_along(chunks), function(k) {
    chunk_worker(chunks[[k]], chunk_names[[k]], cooks_cutoff)
  }, BPPARAM = param)

  RES <- data.table::rbindlist(res_chunks, use.names = TRUE, fill = TRUE)
  data.table::setcolorder(RES, c("site_id", setdiff(names(RES), "site_id")))
  if (verbose) cat("[DONE] Fitted", n_sites, "sites in", length(chunks), "chunks.\n")
  return(RES[])
}

#' Compute Cook's Distance Threshold
#'
#' Converts a Cook's distance cutoff specification into a numeric threshold.
#' Used internally by [`.site_glm()`] and [fit_sites_parallel()].
#'
#' @param n_rows Integer. Number of observations in the model.
#' @param cooks_cutoff Character or numeric. One of:
#'   \describe{
#'     \item{"Inf"}{No filtering.}
#'     \item{"4/n"}{Use the rule-of-thumb 4 / n threshold.}
#'     \item{numeric}{Explicit numeric value.}
#'   }
#'
#' @return Numeric scalar cutoff value.
#' @keywords internal
cutoff_num <- function(n_rows, cooks_cutoff) {
  if (identical(cooks_cutoff, "Inf")) return(Inf)
  if (identical(cooks_cutoff, "none")) return(Inf)
  if (identical(cooks_cutoff, "4/n")) return(4 / n_rows)
  as.numeric(cooks_cutoff)
}



#' Differential Inclusion Analysis from Hit Index Tables
#'
#' Performs per-site differential inclusion testing from a hit-index or
#' junction-form table. Each site is modeled with a quasi-binomial GLM
#' (`psi_adj ~ condition`) to estimate deltaPSI and significance, optionally
#' using parallel processing.
#'
#' @param DT A `data.frame` or `data.table` containing at least the columns
#'   `event_type`, `gene_id`, `chr`, `inc`, `exclusion_reads`,
#'   `inclusion_reads`, `condition`, and `sample`.
#' @param min_total_reads Integer. Minimum total reads per site/sample
#'   required for inclusion (default `10`).
#' @param minimum_proportion_containing_event Numeric in `[0,1]`. Minimum
#'   fraction of samples per condition that must contain the event
#'   (default `0.5`).
#' @param cooks_cutoff Character or numeric. Cook's distance cutoff:
#'   `"4/n"`, `"Inf"`, `"none"`, or a numeric value.
#' @param adjust_method Character. Multiple-testing correction method passed
#'   to [stats::p.adjust()] (default `"fdr"`).
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @param parallel_glm Logical. Use parallel fitting via
#'   [fit_sites_parallel()] (default `TRUE`).
#' @param chunk_size_glm Integer. Number of sites per parallel chunk
#'   (default `1000`).
#' @param cores_glm Integer. Number of parallel workers (default `1`).
#'
#' @return A `data.table` with one row per site containing:
#'   \itemize{
#'     \item Metadata columns (`site_id`, `event_type`, `event_id`, `gene_id`, ...)
#'     \item Sample counts (`n_samples`, `n_control`, `n_case`)
#'     \item Mean PSI per group (`mean_psi_ctrl`, `mean_psi_case`)
#'     \item deltaPSI (`delta_psi`)
#'     \item Raw and adjusted p-values (`p.value`, `padj`)
#'     \item Maximum Cook's distance (`cooks_max`)
#'   }
#'
#' @details
#' The function filters out low-coverage and low-presence events,
#' fills AFE/ALE sites with zero counts where necessary,
#' and then applies site-level GLMs. Parallelization uses
#' [BiocParallel] back-ends for reproducibility across platforms.
#'
#' @seealso [fit_sites_parallel()], [cutoff_num()], [stats::glm()],
#'   [BiocParallel::bplapply()]
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
#' head(res)
#'
#' @importFrom data.table as.data.table rbindlist setcolorder uniqueN %chin%
#' @importFrom stats glm anova cooks.distance p.adjust quasibinomial
#' @importFrom BiocParallel MulticoreParam bplapply
#' @export
get_differential_inclusion <- function(
    DT,
    min_total_reads = 10L,
    minimum_proportion_containing_event = 0.5,
    cooks_cutoff    = "Inf",
    adjust_method   = "fdr",
    verbose         = TRUE,
    parallel_glm = TRUE,
    chunk_size_glm = 1000,
    cores_glm = 1
){
  x <- data.table::as.data.table(DT)
  # --- guards ---
  need <- c("event_type","gene_id","chr","inc","inclusion_reads",
            "exclusion_reads","condition","sample")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))

  # types / baseline
  x[, `:=`(
    condition       = as.character(condition),
    inclusion_reads = as.numeric(inclusion_reads),
    exclusion_reads = as.numeric(exclusion_reads)
  )]
  lev  <- unique(x$condition)
  base <- if ("control" %chin% lev) "control" else sort(lev)[1]
  x[, condition := factor(condition, levels = c(base, setdiff(sort(lev), base)))]

  # site key (use inc only, like your earlier version)
  x[, site_id := paste(event_type, gene_id, chr, inc, exc, form, sep="|")]
  # genes = unique(gene_id)
  # events = unique(site_id)
  # event instances = dim(x)[1]
  if (verbose) print(paste0("[INFO] Input contains ", length(unique(x$gene_id)), " genes, ",
                            round(length(unique(x$event_id))), " events, and ",
                            round(dim(x)[1]/2), " event instances"))

  # total / psi

  x[, total := inclusion_reads + exclusion_reads]
  if (verbose) print(paste0("[PROCESSING/INFO] Filtering out low-coverage, removed ", sum(x$total < min_total_reads, na.rm = TRUE), " event instances, ",
                            sum(x$total >= min_total_reads, na.rm = TRUE), " event instances remain"))
  x <- x[total >= min_total_reads]

  if (!nrow(x)) return(x[0])

  x[, psi_raw := psi]
  x[, psi_adj := psi_raw]

  # ---------- (AFE/ALE) complete with zeros ----------
  if (any(x$event_type %chin% c("AFE","ALE"))) {
    more_types <- FALSE
    if (sum(x$event_type %in% c("AFE", "ALE")) != length(x$event_type)) {
      more_types <- TRUE
      non_type <- x[!(event_type %in% c("AFE", "ALE"))]
      x <- x[event_type %in% c("AFE", "ALE")]
    }
    preAFE <- dim(x[event_type == "AFE"])[1]
    preALE <- dim(x[event_type == "ALE"])[1]
    # minimal set of site columns
    site_cols <- c("gene_id","inc","exc","chr","strand","event_type", "event_id", "form")
    x <- do.call(rbind, lapply(c("AFE", "ALE"), function(terminal_type) {
      sites <- unique(x[event_type %chin% terminal_type, ..site_cols])
      smp   <- unique(x[event_type %chin% terminal_type, .(sample, condition, event_type)])
      # cross join sites x samples (within event_type)
      base_grid <- sites[smp, on="event_type", allow.cartesian=TRUE]
      keep_cols <- c("sample","condition", site_cols,
                     "inclusion_reads","exclusion_reads","total","psi_raw","psi_adj", "source_file")
      y <- base_grid[x[event_type == terminal_type, ..keep_cols], on=c("sample","condition", site_cols)]
      y <- x[event_type == terminal_type, ..keep_cols][
        base_grid,
        on = c("sample","condition", site_cols)
      ]
      y[, total_update := if (any(!is.na(total))) max(total, na.rm = TRUE) else 0,
        by = .(gene_id, sample, condition)]
      y[is.na(total), `:=` (total = total_update,
                            exclusion_reads = total_update
      )]
      y[is.na(psi_adj), `:=` (psi_adj = 0,
                              psi_raw = 0,
                              inclusion_reads = 0)]
      return(y)
    }))
    x[, site_id := paste(event_type, gene_id, chr, inc, exc, sep="|")]
    postAFE <- dim(x[event_type == "AFE"])[1]
    postALE <- dim(x[event_type == "ALE"])[1]
    if (verbose) print(paste0("[PROCESSING/INFO] Completing AFE/ALE with zeros per sample, filled in ", postAFE-preAFE, " AFE instances and ",
                              postALE-preALE, " ALE instances across samples, to a total of ", postAFE, " AFE instances and ", postALE, " ALE instances"))


    # ---------- drop genes all-zero within (sample,event_type) ----------
    all0 <- x[, .(all_zero = all((psi_adj %in% 0) | is.na(psi_adj))),
              by = .(gene_id, sample, event_type)]

    x <- x[all0[all_zero == FALSE], on=.(gene_id, sample, event_type)]
    if (verbose) print(paste0("[PROCESSING/INFO] Filtering genes that have no nonzero values in a given sample, removed ",
                              length(unique(all0[all_zero == TRUE, gene_id])), " genes from specific samples, ",
                              length(unique(x$gene_id)), " genes remain in total"))
    x <- x[, `:=` (
      psi = psi_adj,
      all_zero = NULL,
      total_update = NULL
    )]

    if (more_types == TRUE) {
      x <- rbind(x, non_type)
    }
  }


  # ---------- minimum presence per condition ----------
  present_dt <- x[total >= min_total_reads,
                  .(present = uniqueN(sample)),
                  by = .(condition, event_id)]
  ncond_dt <- x[, .(n_cond = uniqueN(sample)), by = condition]
  bad_pairs <- present_dt[ncond_dt, on="condition"][
    (present / n_cond) < minimum_proportion_containing_event,
    unique(event_id)
  ]
  x <- x[!(event_id %chin% bad_pairs)]
  if (verbose) print(paste0("[PROCESSING/INFO] Filtering events by presence per condition, removed ", length(unique(bad_pairs)), " events, ",
                            length(unique(x$event_id)), " events remain"))
  if (!nrow(x)) return(x[0])

  # ---------- drop non-changing events ----------
  site_multi <- x[, .(multi = uniqueN(psi_adj, na.rm = TRUE) > 1L), by = site_id]
  map_site_pair <- unique(x[, .(site_id, event_id, form)])
  pair_rule <- site_multi[map_site_pair, on = "site_id"][
    , .(keep_pair = any(multi)), by = event_id]
  bad_pairs <- pair_rule[keep_pair == FALSE, event_id]
  x <- x[!(event_id %chin% bad_pairs)]
  print(paste0("[PROCESSING/INFO] Removed ", length(unique(bad_pairs)), " nonchanging events, ",
               length(unique(x$event_id)), " events remain"))



  # ---------- keep with at least 2 rows and both phenotypes ----------
  two_levels <- x[, uniqueN(condition) > 1, by = .(event_id)]
  x <- x[two_levels[V1 == TRUE], on = .(event_id)]
  if (verbose) print(paste0("[PROCESSING/INFO] Filtering events not present in both conditions, removed ", length(unique(two_levels[V1 == FALSE, event_id])),
                            " events, ", length(unique(x$event_id)), " events remain"))

  if (parallel_glm == TRUE) {
    RES <- fit_sites_parallel(x, cores = cores_glm, chunk_size = chunk_size_glm, progress = verbose, verbose = verbose, cooks_cutoff = cooks_cutoff)
  } else {
    if (verbose) cat("[PROCESSING] Fitting quasi-binomial GLMs per site...\n")
    RES <- x[, {
      d <- .SD
      d <- d[is.finite(total) & total > 0]
      if (nrow(d) < 2) {
        list(p.value=NA_real_, cooks_max=NA_real_, n=.N, n_used=NA_integer_,
             mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_)
      } else {
        d[, condition := droplevels(condition)]
        if (data.table::uniqueN(d$condition) < 2) {
          list(p.value=NA_real_, cooks_max=NA_real_, n=.N, n_used=NA_integer_,
               mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_)
        } else {
          full0 <- try(suppressWarnings(glm(psi_adj ~ condition,
                                            family=quasibinomial(),
                                            weights=total, data=d)), silent=TRUE)
          if (inherits(full0,"try-error")) {
            list(p.value=NA_real_, cooks_max=NA_real_, n=.N, n_used=NA_integer_,
                 mean_psi_ctrl=NA_real_, mean_psi_case=NA_real_, delta_psi=NA_real_)
          } else {
            # Cook's only if >10 rows
            thr <- if (nrow(d) > 10) cutoff_num(nrow(d), cooks_cutoff) else Inf
            cd  <- try(suppressWarnings(cooks.distance(full0)), silent=TRUE)
            cd  <- if (inherits(cd,"try-error")) rep(NA_real_, nrow(d)) else as.numeric(cd)
            keep <- ifelse(is.na(cd), TRUE, cd <= thr)

            fit_once <- function(dd) {
              full <- try(suppressWarnings(glm(psi_adj ~ condition,
                                               family=quasibinomial(),
                                               weights=total, data=dd)), silent=TRUE)
              red  <- try(suppressWarnings(glm(psi_adj ~ 1,
                                               family=quasibinomial(),
                                               weights=total, data=dd)), silent=TRUE)
              p    <- tryCatch({
                aa <- anova(red, full, test="F")
                if (nrow(aa) >= 2) as.numeric(aa$`Pr(>F)`[2]) else NA_real_
              }, error=function(e) NA_real_)

              means <- dd[, .(psi_mean = mean(psi_raw, na.rm=TRUE)), by=condition]
              levs  <- levels(dd$condition)
              ctrlL <- levs[1]
              caseL <- levs[min(2L, length(levs))]
              ctrl_mean <- means[condition == ctrlL, psi_mean][1]
              case_mean <- means[condition == caseL, psi_mean][1]

              list(
                p.value        = p,
                cooks_max      = suppressWarnings(max(cd, na.rm=TRUE)),
                n              = nrow(d),
                n_used         = nrow(dd),
                mean_psi_ctrl  = ifelse(is.finite(ctrl_mean), ctrl_mean, NA_real_),
                mean_psi_case  = ifelse(is.finite(case_mean), case_mean, NA_real_),
                delta_psi      = ifelse(all(is.finite(c(ctrl_mean,case_mean))),
                                        case_mean - ctrl_mean, NA_real_)
              )
            }

            if (any(!keep) && sum(keep) >= 2 && data.table::uniqueN(d$condition[keep]) >= 2) {
              fit_once(d[keep])
            } else {
              fit_once(d)
            }
          }
        }
      }
    }, by = site_id, .SDcols = c("condition","psi_adj","psi_raw","total")]
  }

  RES[, padj := p.adjust(p.value, method = adjust_method)]

  meta_cols <- intersect(c("site_id","event_type","gene_id","chr","strand","inc","exc", "event_id", "form"), names(x))
  META <- unique(x[, ..meta_cols])
  out  <- META[RES, on="site_id"]

  CNT <- x[, .(n_samples = .N,
               n_control = sum(condition == levels(x$condition)[1]),
               n_case    = sum(condition == levels(x$condition)[2])),
           by = site_id]
  out <- CNT[out, on="site_id"]

  data.table::setcolorder(out, c("site_id","event_type","event_id", "gene_id","chr","strand","inc","exc",
                                 "n_samples","n_control","n_case",
                                 "mean_psi_ctrl","mean_psi_case","delta_psi",
                                 "p.value","padj","cooks_max", "form"))
  if (verbose) cat("[INFO] Done.\n")
  out[]
}

#' Filter event pairs by significance and deltaPSI thresholds
#'
#' Keeps all rows belonging to events where at least one isoform or site
#' passes adjusted p-value and deltaPSI significance criteria.
#'
#' @param DT A `data.frame` or `data.table` containing at least `event_id`,
#'   `padj`, and `delta_psi` columns.
#' @param padj_thr Numeric. Adjusted p-value threshold (default `0.05`).
#' @param dpsi_thr Numeric. Absolute deltaPSI threshold (default `0.1`).
#'
#' @return A `data.table` containing all rows from event pairs in which
#'   at least one row meets the significance criteria.
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
#' sig_di <- keep_sig_pairs(res)
#' @importFrom data.table as.data.table %chin%
#' @export
keep_sig_pairs <- function(DT, padj_thr = 0.05, dpsi_thr = 0.1) {
  x <- data.table::as.data.table(DT)

  # which rows pass?
  x[, pass := is.finite(padj) & is.finite(delta_psi) &
      (padj <= padj_thr) & (abs(delta_psi) >= dpsi_thr)]

  # keep all rows from pairs where any row passed
  keep_pairs <- x[, any(pass), by = event_id][V1 == TRUE, event_id]
  out <- x[event_id %chin% keep_pairs][, pass := NULL][]
  out[!is.na(delta_psi) & !is.na(p.value)]
}

#' Volcano plot for differential inclusion results
#'
#' Generates a volcano plot of deltaPSI vs. -log10(FDR) highlighting significant events.
#'
#' @param di A `data.frame` or `data.table` containing at least
#'   `delta_psi` and `padj` columns (optionally `event_type`).
#' @param padj_thr Numeric. Adjusted p-value threshold (default `0.05`).
#' @param dpsi_thr Numeric. Absolute deltaPSI threshold (default `0.1`).
#'
#' @return A `ggplot2` object showing differential inclusion significance.
#'
#' @details
#' Significant sites are colored in `deeppink4`; nonsignificant sites are
#' shown in light grey. Dashed and dotted lines indicate deltaPSI and FDR thresholds.
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
#' plot_di_volcano_dt(res)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_vline
#' @importFrom ggplot2 geom_hline labs theme_classic
#' @importFrom data.table as.data.table copy
#' @export
plot_di_volcano_dt <- function(di,
                               padj_thr = 0.05,
                               dpsi_thr = 0.10) {
  DT <- data.table::as.data.table(copy(di))

  # significance flag (data.table style)
  DT[, signif := is.finite(padj) & padj <= padj_thr &
       is.finite(delta_psi) & abs(delta_psi) >= dpsi_thr]

  # counts
  n_sig_total <- DT[ , sum(signif, na.rm = TRUE)]
  by_type <- if ("event_type" %in% names(DT)) {
    DT[ , .(n_sig = sum(signif, na.rm = TRUE)), by = event_type][order(-n_sig)]
  } else {
    data.table::data.table()
  }

  # plot (consistent palette: deeppink4 for sig, grey80 otherwise)
  p <- ggplot2::ggplot(DT, ggplot2::aes(x = delta_psi, y = -log10(padj))) +
    ggplot2::geom_point(ggplot2::aes(color = signif), alpha = 0.65, size = 1.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "deeppink4", `FALSE` = "grey80"),
                                guide = "none") +
    ggplot2::geom_vline(xintercept = c(-dpsi_thr, dpsi_thr), linetype = "dotted") +
    ggplot2::geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    ggplot2::labs(
      subtitle = paste0("Significant events (|DELTA PSI| > ",
                        dpsi_thr, ", FDR < ", padj_thr,
                        "): ", n_sig_total),
      x = expression(Delta~PSI),
      y = expression(-log[10]~FDR)
    ) +
    theme_classic()

  p
}


