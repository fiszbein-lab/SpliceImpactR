empty_results_inputs <- function() {
  x <- rank_test_reference()
  pf <- data.table::data.table(ensembl_transcript_id = x$sequences$transcript_id,
    ensembl_peptide_id = x$sequences$protein_id, database = "pfam", clean_name = "Test",
    name = "Test;chr1:128-160", feature_id = "PF_TEST", start = 10L, stop = 20L,
    chr = "chr1", strand = "+", alt_name = "Test")
  ppi <- data.table::data.table(geneA = character(), geneB = character(),
    DDI = logical(), DMI = logical(), ddi_for_A = character(), ddi_for_B = character(),
    dmi_for_A = character(), dmi_for_B = character())
  list(x = x, pf = pf, ef = get_exon_features(x$annotations, pf), ppi = ppi,
       nonsig = data.table::copy(x$events)[, padj := 1])
}

empty_results_run <- function(inp, res, matching) {
  utils::capture.output(out <- get_splicing_impact(res = res,
    annotation_df = inp$x[c("annotations", "sequences")], protein_feature_total = inp$pf,
    exon_features = inp$ef, ppi = inp$ppi, matching = matching, verbose = FALSE,
    debug_steps = TRUE))
  out
}

column_schema <- function(d) vapply(d, function(v) class(v)[1], character(1))

renders <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(p)
  TRUE
}

test_that("an empty run returns the same columns and types as a run with pairs", {
  inp <- empty_results_inputs()
  for (matching in c("legacy", "orf")) {
    full <- empty_results_run(inp, inp$x$events, matching)
    empty <- empty_results_run(inp, inp$nonsig, matching)
    for (part in c("pairs", "seq_compare", "hits_domain", "hits_final")) {
      expect_gt(nrow(full[[part]]), 0L)
      expect_equal(nrow(empty[[part]]), 0L)
      expect_identical(column_schema(empty[[part]]), column_schema(full[[part]]),
                       info = paste(matching, part))
    }
  }
})

test_that("zero-row analysis steps need no reference resources", {
  inp <- empty_results_inputs()
  full <- empty_results_run(inp, inp$x$events, "legacy")
  seq0 <- compare_sequence_frame(full$pairs[0], NULL)
  dom0 <- get_domains(seq0, NULL)
  fin0 <- get_ppi_switches(dom0, NULL, NULL)
  expect_identical(column_schema(seq0), column_schema(full$seq_compare))
  expect_identical(column_schema(dom0), column_schema(full$hits_domain))
  expect_identical(column_schema(fin0), column_schema(full$hits_final))
})

test_that("summary helpers return empty results instead of failing", {
  inp <- empty_results_inputs()
  empty <- empty_results_run(inp, inp$nonsig, "legacy")

  expect_true(renders(plot_alignment_summary(empty$seq_compare)))
  expect_true(renders(plot_alignment_summary(empty$seq_compare, mode = "transcript")))

  summary <- integrated_event_summary(empty$hits_final, inp$nonsig)
  expect_true(renders(summary$plot))
  for (nm in c("by_type", "class_counts", "score_summary", "domain_prevalence")) {
    expect_equal(nrow(summary$summaries[[nm]]), 0L)
  }
  # The input events are still counted, with none retained.
  expect_identical(as.character(summary$summaries$relative_use$event_type), "A5SS")
  expect_identical(summary$summaries$relative_use$n_post, 0L)

  prox <- get_proximal_shift_from_hits(empty$pairs)
  expect_equal(nrow(prox$data), 0L)
  expect_identical(column_schema(prox$data), c(event_id = "character",
    event_type = "character", strand = "character", pos = "character",
    delta_psi_case = "numeric", neg = "character", delta_psi_control = "numeric",
    I = "integer", V1 = "character"))
  expect_true(renders(prox$plot))
})

test_that("summary helpers handle pairs that lack scores, coding pairs or terminal events", {
  inp <- empty_results_inputs()
  full <- empty_results_run(inp, inp$x$events, "legacy")
  # No pair has a protein alignment score.
  unscored <- data.table::copy(full$seq_compare)[, prot_pid := NA_real_]
  expect_true(renders(plot_alignment_summary(unscored)))
  # No protein-coding pair, so the PPI panels have nothing to show.
  noncoding <- data.table::copy(full$hits_final)[, pc_class := "onePC"]
  # ggplot warns that one pair cannot make a violin; that is not under test.
  expect_true(suppressWarnings(renders(integrated_event_summary(noncoding, inp$x$events)$plot)))
  # The only event is A5SS, so there is no AFE/ALE pair to classify.
  expect_equal(nrow(get_proximal_shift_from_hits(full$pairs)$data), 0L)
})
