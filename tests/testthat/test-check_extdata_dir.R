test_that("check_extdata_dir locates test data and errors on missing paths", {
  # This should exist in inst/extdata/rawData
  p <- check_extdata_dir("rawData/control_S5/")
  expect_true(dir.exists(p))

  # Bad path should error
  expect_error(check_extdata_dir("rawData/does_not_exist/"))
})
