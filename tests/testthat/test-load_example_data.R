test_that("load_example_data returns requested objects with dependencies", {
  ex_sf <- load_example_data("sample_frame")
  expect_true(is.data.frame(ex_sf$sample_frame))
  expect_true(all(c("path", "sample_name", "condition") %in% names(ex_sf$sample_frame)))
  expect_equal(nrow(ex_sf$sample_frame), 8L)

  ex_ann <- load_example_data("annotation_df")
  expect_true(is.list(ex_ann$annotation_df))
  expect_true(all(c("annotations", "sequences", "hybrids") %in% names(ex_ann$annotation_df)))

  ex_feat <- load_example_data(c("annotation_df", "exon_features"))
  expect_true("protein_feature_total" %in% names(ex_feat))
  expect_true("exon_features" %in% names(ex_feat))
  expect_gt(nrow(ex_feat$exon_features), 0L)
})

test_that("load_example_data rejects unknown include values", {
  expect_error(
    load_example_data("not_a_thing"),
    "unknown `include` value"
  )
})
