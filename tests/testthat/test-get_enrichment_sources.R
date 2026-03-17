test_that("get_enrichment warns and returns empty for invalid source labels", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  warns <- character()
  out <- withCallingHandlers(
    get_enrichment(
      foreground = c("TP53", "EGFR", "MYC"),
      background = c("TP53", "EGFR", "MYC", "PTEN", "STAT3"),
      species = "human",
      gene_id_type = "symbol",
      sources = "MSyhigDB:H",
      min_size = 2
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("unknown source label\\(s\\) dropped", warns)))
  expect_true(any(grepl("no valid enrichment sources remain", warns)))

  expect_type(out, "list")
  expect_true(all(
    c("results_per_source", "results_combined", "results_signif", "plot") %in%
      names(out)
  ))
  expect_s3_class(out$results_combined, "data.table")
  expect_equal(nrow(out$results_combined), 0L)
  expect_equal(length(out$results_per_source), 0L)
})
