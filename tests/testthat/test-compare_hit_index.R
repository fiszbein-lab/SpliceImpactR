test_that("compare_hit_index works on test extdata", {
  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/case_S1/')
    ),
    sample_name = c("S5","S1"),
    condition   = c("control","case"),
    stringsAsFactors = FALSE
  )

  cmp <- compare_hit_index(sf, condition_map = c(control = "control", test = "case"))

  expect_type(cmp, "list")
  expect_true("results" %in% names(cmp))
  expect_true("plot"    %in% names(cmp))

  res <- cmp$results
  expect_s3_class(res, "data.table")
  expect_true(nrow(res) > 0)
  expect_true(all(c("event_key","control","test","delta_HIT","fdr") %in% colnames(res)))
})
