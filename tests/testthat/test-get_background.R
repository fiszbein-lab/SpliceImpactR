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
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/"),
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S6/"),
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S1/"),
      file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S2/")
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

test_that("get_domain_background keeps one-sided domain pairs", {
  background_pairs <- data.table::data.table(
    gene_id = "ENSG_TEST",
    transcript_id_1 = "TX1",
    transcript_id_2 = "TX2"
  )

  protein_features <- data.table::data.table(
    ensembl_transcript_id = "TX1",
    ensembl_peptide_id = "ENSP_TEST",
    database = "interpro",
    clean_name = "TestDomain",
    name = "TestDomain;chr1:100-200"
  )

  bg_domains <- SpliceImpactR:::get_domain_background(
    background = background_pairs,
    protein_features = protein_features,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(nrow(bg_domains) == 1L)
  expect_true(length(bg_domains$total_sd_domains[[1]]) == 1L)
  expect_true(length(bg_domains$domains_1[[1]]) == 1L)
  expect_true(length(bg_domains$domains_2[[1]]) == 0L)
})
