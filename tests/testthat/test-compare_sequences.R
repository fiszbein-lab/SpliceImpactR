test_that("sequence comparison and length/identity visualizations run on extdata", {

  # --- Load annotation & protein features ---
  ann <- get_annotation(load = "test")

  interpro_features <- get_protein_features(
    c("interpro"), ann$annotations, timeout = 600, test = TRUE
  )
  signalp_features <- get_protein_features(
    c("signalp"), ann$annotations, timeout = 600, test = TRUE
  )

  prot_feats <- get_comprehensive_annotations(
    list(signalp_features, interpro_features)
  )
  exon_feats <- get_exon_features(ann$annotations, prot_feats)

  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/control_S6/'),
      check_extdata_dir('rawData/case_S1/'),
      check_extdata_dir('rawData/case_S2/')
    ),
    sample_name  = c("S5","S6","S1","S2"),
    condition    = c("control","control","case","case"),
    stringsAsFactors = FALSE
  )

  # --- HIT → DI → matched → sequences → pairs ---
  data <- get_rmats_hit(sf, event_types = c("ALE","AFE","MXE","SE","A3SS","A5SS","RI"))
  di   <- get_differential_inclusion(data, min_total_reads = 10, verbose = FALSE)
  di_sig <- keep_sig_pairs(di)

  matched   <- get_matched_events_chunked(di_sig, ann$annotations, chunk_size = 2000)
  with_seq  <- attach_sequences(matched, ann$sequences)
  pairs     <- get_pairs(with_seq, source = "multi")

  expect_s3_class(pairs, "data.table")

  # --- Sequence + Frame Comparison ---
  seq_compare <- compare_sequence_frame(pairs, ann$annotations)

  expect_s3_class(seq_compare, "data.table")
  expect_true(all(c("frame_call","summary_classification","prot_len_inc","prot_len_exc") %in% colnames(seq_compare)))
  expect_true(nrow(seq_compare) == nrow(pairs))

  # --- Alignment summary plot ---
  alignment_summary <- plot_alignment_summary(seq_compare)

  expect_s3_class(alignment_summary, "ggplot")

  # --- Length comparison plot ---
  length_output <- plot_length_comparison(seq_compare)

  # It’s a patchwork object, but it inherits ggplot parts
  expect_true(inherits(length_output, "patchwork") || inherits(length_output, "ggplot"))

})
