nmd_reference <- function() {
  # The selected control (the A5SS short form) is annotated as NMD.
  x <- rank_test_reference()
  x$annotations[transcript_id == "E_EXC", transcript_type := "nonsense_mediated_decay"]
  x
}

test_that("both matchers carry the biotype and the summary marks NMD pairs", {
  x <- nmd_reference()
  matched <- get_matched_events_chunked(x$events, x$annotations)
  expect_identical(matched[form == "EXC", transcript_type], "nonsense_mediated_decay")
  expect_identical(matched[form == "INC", transcript_type], "protein_coding")
  legacy <- get_pairs(attach_sequences(matched, x$sequences), source = "multi")
  orf <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)$pairs
  plain <- rank_test_reference()
  plain_frame <- compare_sequence_frame(
    get_ranked_pairs(plain$events, plain$annotations, plain$sequences, verbose = FALSE)$pairs,
    plain$annotations)
  for (pairs in list(legacy, orf)) {
    expect_identical(pairs$transcript_type_case, "protein_coding")
    expect_identical(pairs$transcript_type_control, "nonsense_mediated_decay")
    out <- compare_sequence_frame(pairs, x$annotations)
    expect_identical(out$summary_classification, "NMD")
    # NMD changes the summary only; the frame comparison is reported as before.
    expect_identical(out$frame_call, plain_frame$frame_call)
  }
  expect_false(identical(plain_frame$summary_classification, "NMD"))
})

test_that("the summary reads biotypes from the annotation when pairs lack them", {
  x <- nmd_reference()
  pairs <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)$pairs
  bare <- data.table::copy(pairs)[, c("transcript_type_case", "transcript_type_control") := NULL]
  out <- compare_sequence_frame(bare, x$annotations)
  expect_identical(out$transcript_type_control, "nonsense_mediated_decay")
  expect_identical(out$summary_classification, "NMD")
})

test_that("unmatched forms have no biotype", {
  x <- nmd_reference()
  x$events[form == "INC", inc := "5001-5100"]
  matched <- get_matched_events_chunked(x$events, x$annotations)
  expect_true(is.na(matched[form == "INC", transcript_type]))
  expect_identical(matched[form == "EXC", transcript_type], "nonsense_mediated_decay")
})

test_that("the integrated summary keeps NMD as its own class", {
  x <- nmd_reference()
  pf <- data.table::data.table(ensembl_transcript_id = x$sequences$transcript_id,
    ensembl_peptide_id = x$sequences$protein_id, database = "pfam", clean_name = "Test",
    name = "Test;chr1:128-160", feature_id = "PF_TEST", start = 10L, stop = 20L,
    chr = "chr1", strand = "+", alt_name = "Test")
  ppi <- data.table::data.table(geneA = character(), geneB = character(),
    DDI = logical(), DMI = logical(), ddi_for_A = character(), ddi_for_B = character(),
    dmi_for_A = character(), dmi_for_B = character())
  pairs <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)$pairs
  hits <- get_domains(compare_sequence_frame(pairs, x$annotations), get_exon_features(x$annotations, pf))
  summary <- integrated_event_summary(get_ppi_switches(hits, ppi, pf), x$events)
  counts <- summary$summaries$class_counts
  expect_true("NMD" %in% levels(counts$summary_classification))
  expect_identical(as.character(counts$summary_classification), "NMD")
})
