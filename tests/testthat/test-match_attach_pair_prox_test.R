test_that("terminal event matching → sequences → pairing → proximal shift runs on extdata", {

  ann <- get_annotation(load = "test")

  interpro_features <- get_protein_features(
    c("interpro"), ann$annotations, timeout = 600, test = TRUE
  )
  signalp_features <- get_protein_features(
    c("signalp"), ann$annotations, timeout = 600, test = TRUE
  )
  prot_feats <- get_comprehensive_annotations(list(signalp_features, interpro_features))

  exon_feats <- get_exon_features(ann$annotations, prot_feats)

  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/control_S6/'),
      check_extdata_dir('rawData/control_S7/'),
      check_extdata_dir('rawData/control_S8/'),
      check_extdata_dir('rawData/case_S1/'),
      check_extdata_dir('rawData/case_S2/'),
      check_extdata_dir('rawData/case_S3/'),
      check_extdata_dir('rawData/case_S4/')
    ),
    sample_name  = c("S5","S6","S7","S8","S1","S2","S3","S4"),
    condition    = c("control","control","control","control","case","case","case","case"),
    stringsAsFactors = FALSE
  )

  data <- get_rmats_hit(sf, event_types=c("ALE","AFE","MXE","SE","A3SS","A5SS","RI"))
  di   <- get_differential_inclusion(data, min_total_reads = 10, verbose = FALSE)
  di_sig <- keep_sig_pairs(di)

  # match to annotation
  matched <- get_matched_events_chunked(di_sig, ann$annotations, chunk_size = 2000)
  expect_s3_class(matched, "data.table")
  expect_true("transcript_id" %in% colnames(matched))

  # attach sequences
  with_seq <- attach_sequences(matched, ann$sequences)
  expect_s3_class(with_seq, "data.table")
  expect_true("transcript_seq" %in% colnames(with_seq))
  expect_true(nrow(with_seq) == nrow(matched))

  # pair events (HIT-mode)
  pairs <- get_pairs(with_seq, source = "multi")
  expect_s3_class(pairs, "data.table")

  if (nrow(pairs) > 0) {
    expect_true("event_id" %in% names(pairs))
  }

  # proximal shift calc (if user function exists)
  if ("get_proximal_shift_from_hits" %in% ls("package:SpliceImpactR")) {
    prox <- get_proximal_shift_from_hits(pairs)
    expect_s3_class(prox$data, "data.table")
    expect_s3_class(prox$plot, "ggplot")

  }
})
