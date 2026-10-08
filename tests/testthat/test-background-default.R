background_inputs <- function() {
  x <- rank_test_reference()
  tx <- x$sequences$transcript_id
  # A different domain on each transcript, so every within-gene pair differs.
  pf <- data.table::data.table(ensembl_transcript_id = tx, ensembl_peptide_id = paste0("P_", tx),
    database = "pfam", clean_name = paste0("D", seq_along(tx)),
    name = paste0("D", seq_along(tx), ";chr1:128-160"), feature_id = paste0("PF", seq_along(tx)),
    start = 10L, stop = 20L, chr = "chr1", strand = "+", alt_name = "D")
  list(annotations = x$annotations, protein_features = pf)
}

test_that("get_background defaults to the annotated background", {
  x <- background_inputs()
  bp <- BiocParallel::SerialParam()
  default <- get_background(annotations = x$annotations, protein_features = x$protein_features,
                            BPPARAM = bp)
  annotated <- get_background(source = "annotated", annotations = x$annotations,
                              protein_features = x$protein_features, BPPARAM = bp)
  expect_identical(default, annotated)
  expect_identical(nrow(default), 3L)
})

test_that("an input without a source is flagged instead of silently changing meaning", {
  x <- background_inputs()
  expect_warning(get_background(input = data.frame(path = "unused"), annotations = x$annotations,
                                protein_features = x$protein_features,
                                BPPARAM = BiocParallel::SerialParam()),
                 "`input` is ignored")
  expect_error(get_background(source = "hit_index", annotations = x$annotations,
                              protein_features = x$protein_features),
               "`input` is required")
})

test_that("a background missing foreground pairs names them and the annotated source", {
  fg <- data.table::data.table(event_id = c("E1", "E2"), transcript_id_case = c("T1", "T3"),
                               transcript_id_control = c("T2", "T4"),
                               either_domains_list = list("pfam;D1", "pfam;D2"))
  bg <- data.table::data.table(transcript_id_1 = "T1", transcript_id_2 = "T2",
                               total_sd_domains = list("pfam;D1"))
  expect_error(enrich_domains_hypergeo(fg, bg, min_fg_count = 1),
               "1 of 2 foreground transcript pairs are absent.*T3\\|T4.*source = \"annotated\"")
})
