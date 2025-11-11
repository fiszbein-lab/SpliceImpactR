test_that("get_comprehensive_annotations merges multiple feature sources", {
  ann <- get_annotation(load = "test")

  f1 <- get_protein_features("interpro", ann$annotations, test = TRUE)
  f2 <- get_protein_features("signalp", ann$annotations, test = TRUE)

  merged <- get_comprehensive_annotations(list(f1, f2))

  expect_s3_class(merged, "data.table")
  expect_gt(nrow(merged), 0)

  # Should contain rows from both feature sets
  expect_true(any(merged$database == "interpro"))
  expect_true(any(merged$database == "signalp"))
})
