test_that("get_background works end-to-end on extdata", {

  # --- Load annotation and protein features ---
  ann <- get_annotation(load = "test")

  interpro_features <- get_protein_features(
    c("interpro"), ann$annotations, timeout = 600, test = TRUE
  )
  protein_feats <- get_comprehensive_annotations(list(interpro_features))

  # ---- Sample frame (small for test speed) ----
  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/control_S6/'),
      check_extdata_dir('rawData/case_S1/'),
      check_extdata_dir('rawData/case_S2/')
    ),
    sample_name = c("S5", "S6", "S1", "S2"),
    condition   = c("control","control","case","case"),
    stringsAsFactors = FALSE
  )

  # ---- HIT-index background ----
  bg_hit <- get_background(
    source            = "hit_index",
    input             = sf,
    annotations       = ann$annotations,
    protein_features  = protein_feats,
    keep_annotated_first_last = TRUE,
    minOverlap        = 0.8
  )

  expect_s3_class(bg_hit, "data.table")
  expect_true(nrow(bg_hit) > 0)
  expect_true(all(c("gene_id","transcript_id","i.transcript_id") %in% colnames(bg_hit)))
  expect_true(all(c("domains_1","domains_2") %in% colnames(bg_hit)))
  expect_true(all(c("prot_aa_1","prot_aa_2") %in% colnames(bg_hit)))

  # ---- Annotated-only mode ----
  bg_ann <- get_background(
    source           = "annotated",
    input            = NULL,
    annotations      = ann$annotations,
    protein_features = protein_feats
  )
  expect_s3_class(bg_ann, "data.table")
  expect_true(nrow(bg_ann) > 0)

  # ---- User-given transcript test ----
  some_tx <- unique(ann$annotations$transcript_id)[1:20]

  bg_user <- get_background(
    source            = "user-given",
    input             = some_tx,
    annotations       = ann$annotations,
    protein_features  = protein_feats
  )

  expect_s3_class(bg_user, "data.table")
  expect_true(nrow(bg_user) >= 1)
  expect_true(all(c("transcript_id","i.transcript_id") %in% names(bg_user)))
})
