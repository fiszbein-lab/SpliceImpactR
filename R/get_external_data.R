#' Standardize a differential inclusion (DI) result table
#'
#' Converts an arbitrary differential inclusion result table into the
#' standardized column format expected by SpliceImpactR functions.
#'
#' @param df A `data.frame` or `data.table` containing differential inclusion results.
#' @param colmap A named list mapping required fields in `df` to standard
#'   names (`gene_id`, `chr`, `strand`, `inc`, `exc`, `delta_psi`, `pvalue`,
#'   and optionally `event_type`).
#' @param default_event_type Character. Default `event_type` to assign if none
#'   provided (default `"SITE"`).
#' @param adjust_method Character. Multiple-testing correction method passed
#'   to [stats::p.adjust()] (default `"fdr"`).
#' @param add_chr_prefix Logical. Add `"chr"` prefix if absent (default `FALSE`).
#'
#' @return A standardized `data.table` with columns:
#'   `site_id`, `event_type`, `gene_id`, `chr`, `strand`, `inc`, `exc`,
#'   `delta_psi`, `p.value`, `padj`, and `form`.
#'
#' @details
#' This function provides a uniform interface for importing external DI results
#' (e.g. from rMATS, MAJIQ, or SUPPA2) so they can be compared or plotted
#' alongside SpliceImpactR outputs.
#'
#' @examples
#' df <- data.frame(
#'   gene_id = "ENSG00000280071",
#'   chr = "7",
#'   strand = "+",
#'   inc = "1940088-1940549",
#'   exc = "",
#'   delta_psi = 0.25,
#'   p.value = 0.01
#' )
#' di_std <- import_di_table(df)
#' head(di_std)

#'
#' @importFrom data.table as.data.table setcolorder
#' @importFrom stats p.adjust
#' @export
import_di_table <- function(df,
  colmap = list(
    gene_id   = "gene_id",
    chr       = "chr",
    strand    = "strand",
    inc       = "inc",
    exc       = "exc",
    delta_psi = "delta_psi",
    pvalue    = "p.value",   # name of the p-value column in input
    event_type= NULL         # if present in df, set this to its name
  ),
  default_event_type = "SITE",
  adjust_method = "fdr",
  add_chr_prefix = FALSE) {
    DT <- as.data.table(df)

    # resolve columns (allow user renames)
    need <- c("gene_id","chr","strand","inc","exc","delta_psi","pvalue")
    miss_map <- setdiff(need, names(colmap))
    if (length(miss_map)) stop("colmap missing mappings for: ", paste(miss_map, collapse=", "))

    # Check all required mapped columns present in df
    mapped_names <- vapply(need, function(nm) colmap[[nm]], character(1), USE.NAMES = FALSE)
    miss <- setdiff(mapped_names, names(DT))
    if (length(miss)) stop("Input is missing required columns: ", paste(miss, collapse=", "))

    # Build a small view with standard names
    out <- DT[, .(
      gene_id   = as.character(get(colmap$gene_id)),
      chr       = as.character(get(colmap$chr)),
      strand    = as.character(get(colmap$strand)),
      inc       = as.character(get(colmap$inc)),
      exc       = as.character(get(colmap$exc)),
      delta_psi = suppressWarnings(as.numeric(get(colmap$delta_psi))),
      p.value   = suppressWarnings(as.numeric(get(colmap$pvalue)))
    )]

    # Cleanups
    out[is.na(inc), inc := ""]
    out[is.na(exc), exc := ""]
    out[!strand %chin% c("+","-"), strand := "+"]

    if (isTRUE(add_chr_prefix)) {
      out[!grepl("^chr", chr, ignore.case = FALSE), chr := paste0("chr", chr)]
    }

    # event_type: from input if available, else default
    if (!is.null(colmap$event_type) && colmap$event_type %in% names(DT)) {
      out[, event_type := as.character(DT[[colmap$event_type]])]
    } else {
      out[, event_type := default_event_type]
    }

    # site_id = event_type|gene_id|chr|inc|exc  (strand can be appended if you prefer)
    out[, site_id := paste(event_type, gene_id, chr, inc, exc, sep="|")]

    # padj (multiple-testing correction on available p-values)
    out[, padj := p.adjust(p.value, method = adjust_method)]

    # form: this table is already per-site summary, so mark as SITE
    out[, form := "SITE"]

    # (optional) light, useful counts if present in input
    maybe_cols <- intersect(c("n","n_used","n_samples","n_control","n_case",
                              "mean_psi_ctrl","mean_psi_case"), names(DT))
    if (length(maybe_cols)) {
      out <- cbind(out, DT[, ..maybe_cols, with=FALSE])
    }

    # Tidy column order
    setcolorder(out, c("site_id","event_type","gene_id","chr","strand","inc","exc",
                       "delta_psi","p.value","padj","form",
                       setdiff(names(out),
                               c("site_id","event_type","gene_id","chr","strand","inc","exc",
                                 "delta_psi","p.value","padj","form"))))
    out[]
}



