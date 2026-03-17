#' @title Load bundled example inputs for documentation
#'
#' @description
#' Returns commonly used example objects (sample manifest, test annotations,
#' optional feature mappings, optional PPI table) so man-page examples can stay
#' focused on the documented function rather than setup boilerplate.
#'
#' @param include Character vector selecting which objects to return.
#'   Supported values are `"sample_frame"`, `"annotation_df"`,
#'   `"protein_feature_total"`, `"exon_features"`, `"ppi"`, and `"all"`.
#'   If `"all"` is present, it expands to
#'   `c("sample_frame", "annotation_df", "protein_feature_total", "exon_features")`.
#' @param biomaRt_databases Character vector passed to [get_protein_features()]
#'   when protein features are requested. Default is `"interpro"`.
#' @param test Logical passed to [get_protein_features()] (default `TRUE`).
#'
#' @return
#' Named list containing the requested objects.
#'
#' @examples
#' ex <- load_example_data(c("sample_frame", "annotation_df"))
#' sample_frame <- ex$sample_frame
#' annotation_df <- ex$annotation_df
#'
#' ex2 <- load_example_data(c("annotation_df", "exon_features"))
#' exon_features <- ex2$exon_features
#' print(exon_features)
#' @export
load_example_data <- function(
    include = c("sample_frame", "annotation_df"),
    biomaRt_databases = c("interpro"),
    test = TRUE
) {
  include <- unique(as.character(include))
  valid <- c(
    "sample_frame",
    "annotation_df",
    "protein_feature_total",
    "exon_features",
    "ppi",
    "all"
  )
  bad <- setdiff(include, valid)
  if (length(bad)) {
    stop(
      "load_example_data: unknown `include` value(s): ",
      paste(bad, collapse = ", ")
    )
  }

  if ("all" %in% include) {
    include <- union(
      include,
      c("sample_frame", "annotation_df", "protein_feature_total", "exon_features")
    )
    include <- setdiff(include, "all")
  }

  # ensure dependency closure
  if ("exon_features" %in% include) {
    include <- union(include, c("annotation_df", "protein_feature_total"))
  }
  if ("protein_feature_total" %in% include) {
    include <- union(include, "annotation_df")
  }

  out <- list()

  if ("sample_frame" %in% include) {
    ext <- system.file("extdata", package = "SpliceImpactR")
    out$sample_frame <- data.frame(
      path = c(
        file.path(ext, "rawData/control_S5/"),
        file.path(ext, "rawData/control_S6/"),
        file.path(ext, "rawData/control_S7/"),
        file.path(ext, "rawData/control_S8/"),
        file.path(ext, "rawData/case_S1/"),
        file.path(ext, "rawData/case_S2/"),
        file.path(ext, "rawData/case_S3/"),
        file.path(ext, "rawData/case_S4/")
      ),
      sample_name = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
      condition = c(
        "control", "control", "control", "control",
        "case", "case", "case", "case"
      ),
      stringsAsFactors = FALSE
    )
  }

  if ("annotation_df" %in% include) {
    out$annotation_df <- get_annotation(load = "test")
  }

  if ("protein_feature_total" %in% include) {
    pf <- get_protein_features(
      biomaRt_databases = biomaRt_databases,
      gtf_df = out$annotation_df$annotations,
      test = test
    )
    out$protein_feature_total <- get_comprehensive_annotations(list(pf))
  }

  if ("exon_features" %in% include) {
    out$exon_features <- get_exon_features(
      out$annotation_df$annotations,
      out$protein_feature_total
    )
  }

  if ("ppi" %in% include) {
    out$ppi <- get_ppi_interactions()
  }

  out
}
