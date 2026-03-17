test_that("plot_di_volcano_dt generates a volcano plot", {
  sf <- data.frame(path = c(file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S6/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S7/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S8/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S1/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S2/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S3/"),
                            file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S4/")),
                   sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
                   condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
                   stringsAsFactors = FALSE)

  events <- get_rmats_hit(
    sf,
    event_types = c("AFE","ALE"),
    keep_annotated_first_last = TRUE
  )

  res <- get_differential_inclusion(
    events,
    min_total_reads = 10,
    verbose = FALSE,
    parallel_glm = FALSE
  )

  p <- plot_di_volcano_dt(res)

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})
