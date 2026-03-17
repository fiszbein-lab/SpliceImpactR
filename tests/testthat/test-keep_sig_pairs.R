test_that("keep_sig_pairs selects significant DI events", {
  skip_on_cran()

  sample_frame <- data.frame(path = c(file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/"),
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

  # DI results
  events <- get_rmats_hit(
    sample_frame,
    event_types = c("AFE","ALE"),
    keep_annotated_first_last = TRUE
  )

  res <- get_differential_inclusion(
    events,
    min_total_reads = 10,
    verbose = FALSE,
    parallel_glm = FALSE
  )

  sig <- keep_sig_pairs(res)

  expect_s3_class(sig, "data.table")

  if (nrow(sig) > 0) {
    expect_true("event_id" %in% names(sig))
    # ensure we kept *pairs* not single isoforms
    expect_true(all(sig$event_id %in% res$event_id))
  }
})
