test_that("get_protein_features(test=TRUE) returns expected structure", {
  ann <- get_annotation(load = "test")
  feat <- get_protein_features(
    biomaRt_databases = c("interpro", "signalp"),
    gtf_df = ann$annotations,
    test = TRUE
  )

  expect_s3_class(feat, "data.table")
  expect_gt(nrow(feat), 0)

  req_cols <- c(
    "ensembl_transcript_id","ensembl_peptide_id",
    "database","feature_id","name","alt_name",
    "start","stop"
  )

  expect_true(all(req_cols %in% colnames(feat)))

  # Coordinates should be integers
  expect_true(is.integer(feat$start))
  expect_true(is.integer(feat$stop))

  # At least InterPro / signalp present
  expect_true(any(feat$database == "interpro"))
  expect_true(any(feat$database == "signalp"))
})
