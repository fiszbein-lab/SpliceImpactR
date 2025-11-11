test_that("get_annotation(test) returns expected structure and fields", {

  ann <- get_annotation(load = "test")

  # Top-level structure
  expect_type(ann, "list")
  expect_named(ann, c("annotations", "sequences", "hybrids"), ignore.order = TRUE)

  # ---- annotations ----
  gtf <- ann$annotations
  expect_s3_class(gtf, "data.table")
  expect_true(nrow(gtf) > 0)

  # Required columns
  expect_true(all(c("gene_id", "transcript_id", "type", "chr", "start", "end") %in% colnames(gtf)))

  # exon fields added by pipeline
  expect_true(all(c("cds_has", "absolute_exon_position", "absolute_exon_class", "feature_length") %in% colnames(gtf)))

  # Frames columns may be partially NA, just check they exist
  expect_true(all(c("start_frame", "stop_frame") %in% colnames(gtf)))

  # Make sure at least one exon exists
  expect_gt(sum(gtf$type == "exon"), 0)

  # ---- sequences ----
  seqs <- ann$sequences
  expect_s3_class(seqs, "data.table")
  expect_true(nrow(seqs) > 0)
  expect_true(all(c("gene_id", "transcript_id", "protein_id") %in% colnames(seqs)))

  # Check sequences have string fields
  expect_true("transcript_seq" %in% colnames(seqs))
  expect_true("protein_seq" %in% colnames(seqs))
  expect_true(any(nchar(seqs$transcript_seq) > 0, na.rm = TRUE))
  expect_true(any(nchar(seqs$protein_seq) > 0, na.rm = TRUE))

  # ---- hybrids ----
  hyb <- ann$hybrids
  expect_type(hyb, "list")
  expect_named(hyb, c("first_hybrids", "last_hybrids"), ignore.order = TRUE)

  expect_s3_class(hyb$first_hybrids, "data.table")
  expect_s3_class(hyb$last_hybrids,  "data.table")

  # Hybrids tables may be empty for tiny test sets, so just check columns exist
  expect_true(all(c("gene_id", "transcript_id_terminal", "transcript_id_internal") %in% colnames(hyb$first_hybrids)))
  expect_true(all(c("gene_id", "transcript_id_terminal", "transcript_id_internal") %in% colnames(hyb$last_hybrids)))
})
