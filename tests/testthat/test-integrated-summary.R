test_that("every integrated summary panel uses the same text size", {
  x <- rank_test_reference()
  pf <- data.table::data.table(ensembl_transcript_id = x$sequences$transcript_id,
    ensembl_peptide_id = x$sequences$protein_id, database = "pfam", clean_name = "Test",
    name = "Test;chr1:128-160", feature_id = "PF_TEST", start = 10L, stop = 20L,
    chr = "chr1", strand = "+", alt_name = "Test")
  ppi <- data.table::data.table(geneA = character(), geneB = character(),
    DDI = logical(), DMI = logical(), ddi_for_A = character(), ddi_for_B = character(),
    dmi_for_A = character(), dmi_for_B = character())
  pairs <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)$pairs
  hits <- get_ppi_switches(get_domains(compare_sequence_frame(pairs, x$annotations),
                                       get_exon_features(x$annotations, pf)), ppi, pf)
  # One event type, so the coordination heatmap is a placeholder panel.
  plot <- integrated_event_summary(hits, x$events)$plot
  panels <- function(p) {
    out <- list(p)
    if (inherits(p, "patchwork")) for (q in p$patches$plots) out <- c(out, panels(q))
    out
  }
  sizes <- vapply(panels(plot), function(p) {
    size <- p$theme$text$size
    if (is.null(size)) NA_real_ else as.numeric(size)
  }, numeric(1))
  expect_true(all(sizes == 13))
})
