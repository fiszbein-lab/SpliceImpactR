test_that("integrated splicing → match → seq → compare → domains → enrichment runs on extdata", {

  # --- Load annotation + protein features ----
  annotation_df <- get_annotation(load = "test")

  interpro_features <- get_protein_features(
    c("interpro"), annotation_df$annotations,
    timeout = 600, test = TRUE
  )
  signalp_features <- get_protein_features(
    c("signalp"), annotation_df$annotations,
    timeout = 600, test = TRUE
  )

  protein_feature_total <- get_comprehensive_annotations(
    list(signalp_features, interpro_features)
  )

  exon_features <- get_exon_features(
    annotation_df$annotations, protein_feature_total
  )

  # --- Sample frame from extdata ----
  sample_frame <- data.frame(
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
    condition    = c("control","control","control","control",
                     "case","case","case","case"),
    stringsAsFactors = FALSE
  )

  # --- HIT index → differential splicing → filtering ----
  data <- get_rmats_hit(
    sample_frame,
    event_types = c("ALE","AFE","MXE","SE","A3SS","A5SS","RI")
  )

  res <- get_differential_inclusion(data, min_total_reads = 10)
  res_di <- keep_sig_pairs(res)

  # --- Map to annotation + sequences ----
  matched <- get_matched_events_chunked(
    res_di, annotation_df$annotations, chunk_size = 2000
  )
  expect_s3_class(matched, "data.table")

  hits_sequences <- attach_sequences(matched, annotation_df$sequences)
  expect_s3_class(hits_sequences, "data.table")
  expect_true("transcript_seq" %in% colnames(hits_sequences))

  # --- Isoform pairing + frame comparison ----
  pairs <- get_pairs(hits_sequences, source="multi")
  expect_s3_class(pairs, "data.table")

  seq_compare <- compare_sequence_frame(pairs, annotation_df$annotations)
  expect_s3_class(seq_compare, "data.table")
  expect_true(all(c("frame_call","summary_classification",
                    "prot_len_inc","prot_len_exc") %in% colnames(seq_compare)))

  alignment_summary <- plot_alignment_summary(seq_compare)
  expect_s3_class(alignment_summary, "ggplot")

  # --- Background (annotated) ----
  bg <- get_background(
    source = "annotated",
    annotations = annotation_df$annotations,
    protein_features = protein_feature_total
  )
  expect_s3_class(bg, "data.table")
  expect_true(nrow(bg) > 0)

  # --- Domain annotation ----
  hits_domain <- get_domains(seq_compare, exon_features)
  expect_s3_class(hits_domain, "data.table")
  expect_true(all(c("inc_only_n","exc_only_n","diff_n") %in% colnames(hits_domain)))

  # --- Hypergeometric enrichment ----
  enriched_domains <- enrich_domains_hypergeo(
    hits_domain, bg, db_filter = "interpro"
  )
  expect_s3_class(enriched_domains, "data.table")

  # --- Domain bar plot ----
  domain_plot <- plot_enriched_domains_counts(enriched_domains, top_n = 20)
  expect_s3_class(domain_plot, "ggplot")

  ppidm <- get_ppidm(test=TRUE)
  expect_s3_class(ppidm, "data.table")
  expect_true("I1" %in% names(ppidm))
  expect_true("I2" %in% names(ppidm))

  restrict_protein_features <- protein_feature_total[ensembl_transcript_id %in%
                                                       c(hits_domain$transcript_id_exc, hits_domain$transcript_id_exc)]
  ppi <- get_isoform_interactions(restrict_protein_features, ppidm, save = FALSE, load_dir = NULL, init = TRUE)
  expect_s3_class(ppi, "data.table")
  expect_true("IA" %in% names(ppi))
  expect_true("IB" %in% names(ppi))

  hits_final <- get_ppi_switches(hits_domain, ppi)
  expect_s3_class(hits_final, "data.table")
  expect_true("n_inc_ppi" %in% names(hits_final))
  expect_true("n_ppi" %in% names(hits_final))
  ppi_plot <- plot_ppi_summary(hits_final)
  expect_s3_class(ppi_plot, "ggplot")

})
