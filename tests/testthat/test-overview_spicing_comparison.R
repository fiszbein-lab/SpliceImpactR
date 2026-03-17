test_that("overview_spicing_comparison works on extdata HIT events", {

  skip_on_cran()

  sf <- data.frame(
    path        = c(
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/"),
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S1/")
    ),
    sample_name = c("S5", "S1"),
    condition   = c("control", "case"),
    stringsAsFactors = FALSE
  )

  # Get HIT events from our bundled rmats test files
  events <- get_rmats_hit(
    sf,
    event_types = c("AFE", "ALE"),
    keep_annotated_first_last = TRUE
  )


  # Run overview plot
  ov <- overview_spicing_comparison(
    events,
    sf,
    depth_norm = "exon_files",
    event_type = "AFE"
  )

  # Expectations
  expect_s3_class(ov, "patchwork")
  expect_false(is.null(ov))
})
