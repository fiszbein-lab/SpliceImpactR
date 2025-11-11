test_that("get_differential_inclusion runs on extdata HIT events", {

  sf <- data.frame(path = c(check_extdata_dir('rawData/control_S5/'),
                                      check_extdata_dir('rawData/control_S6/'),
                                      check_extdata_dir('rawData/control_S7/'),
                                      check_extdata_dir('rawData/control_S8/'),
                                      check_extdata_dir('rawData/case_S1/'),
                                      check_extdata_dir('rawData/case_S2/'),
                                      check_extdata_dir('rawData/case_S3/'),
                                      check_extdata_dir('rawData/case_S4/')),
                             sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
                             condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
                             stringsAsFactors = FALSE)

  # Load test events
  events <- get_rmats_hit(
    sf,
    event_types = c("AFE", "ALE"),
    keep_annotated_first_last = TRUE
  )

  # DI test (parallel off for reproducible tiny test speed)
  res <- get_differential_inclusion(
    events,
    min_total_reads = 10,
    verbose = FALSE,
    parallel_glm = FALSE
  )

  expect_s3_class(res, "data.table")
  expect_gt(nrow(res), 0)

  # Key columns present
  req <- c(
    "site_id","event_type","event_id","gene_id","chr","strand","inc","exc",
    "n_samples","n_control","n_case",
    "mean_psi_ctrl","mean_psi_case","delta_psi",
    "p.value","padj","cooks_max","form"
  )
  expect_true(all(req %in% names(res)))

  # psi / p-values valid-ish
  expect_true(all(res$n_samples >= 2))
  expect_true(any(!is.na(res$padj)))
})
