test_that("get_exon_features maps features to exons correctly in test mode", {
  ann <- get_annotation(load = "test")

  feat <- get_protein_features(
    biomaRt_databases = c("interpro","signalp"),
    gtf_df = ann$annotations,
    test = TRUE
  )
  merged <- get_comprehensive_annotations(list(feat))

  ex_feats <- get_exon_features(ann$annotations, merged)

  # Should return a data.table (may be empty for tiny test data)
  expect_s3_class(ex_feats, "data.table")

  if (nrow(ex_feats) > 0) {
    req_cols <- c(
      "gene_id","ensembl_transcript_id","ensembl_peptide_id",
      "exon_id","exon_number","strand",
      "database","feature_id","name","alt_name",
      "prot_start","prot_stop","exon_aa_start","exon_aa_end",
      "overlap_aa_start","overlap_aa_end","overlap_aa_len"
    )
    expect_true(all(req_cols %in% colnames(ex_feats)))

    expect_true(all(ex_feats$overlap_aa_len >= 1))
  }
})
