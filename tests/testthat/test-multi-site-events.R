multi_site_events <- function() {
  # One AFE with a rising first exon and three falling ones; only the rising
  # site passes the thresholds itself.
  rank_events("AFE:1", "AFE", "SITE", c("101-200", "301-400", "501-600", "701-800"), "",
              c(0.6, -0.3, -0.2, -0.1))[, padj := c(0.001, 0.2, 0.3, 0.4)][]
}

test_that("keep_sig_pairs keeps a passing event's other sites and labels them", {
  quiet <- rank_events("SE:2", "SE", c("INC", "EXC"), c("101-200;301-400;501-600", "101-200;501-600"),
                       c("", "301-400"), c(0.05, -0.05))
  out <- keep_sig_pairs(data.table::rbindlist(list(multi_site_events(), quiet)))
  expect_identical(unique(out$event_id), "AFE:1")
  expect_identical(out[order(-delta_psi), site_significant], c(TRUE, FALSE, FALSE, FALSE))
})

test_that("multi-site events give one labelled comparison per site pair", {
  ev <- keep_sig_pairs(multi_site_events())
  ev[, `:=`(transcript_id = paste0("T", seq_len(.N)), exons = "", protein_id = NA_character_,
            transcript_seq = NA_character_, protein_seq = NA_character_)]
  pairs <- get_pairs(ev, source = "multi")
  expect_identical(pairs$n_event_comparisons, rep(3L, 3L))
  expect_identical(pairs$site_significant_case, rep(TRUE, 3L))
  expect_identical(pairs$site_significant_control, rep(FALSE, 3L))
  # A site without a transcript gives no row but is still one of the comparisons.
  ev[delta_psi == -0.1, transcript_id := NA_character_]
  pairs <- get_pairs(ev, source = "multi")
  expect_identical(pairs$n_event_comparisons, rep(3L, 2L))
  # A two-form event defines one comparison.
  two <- data.table::copy(ev[1:2])[, delta_psi := c(0.3, -0.3)]
  expect_identical(get_pairs(two, source = "multi")$n_event_comparisons, 1L)
})

test_that("both matchers report every comparison of a multi-site event", {
  structures <- list(T1 = matrix(c(101, 200, 901, 1000), byrow = TRUE, ncol = 2L),
                     T2 = matrix(c(301, 400, 901, 1000), byrow = TRUE, ncol = 2L),
                     T3 = matrix(c(501, 600, 901, 1000), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures)
  # The fourth site (701-800) has no transcript, so its comparison is unresolved.
  ev <- keep_sig_pairs(multi_site_events())
  out <- get_ranked_pairs(ev, x$annotations, x$sequences, verbose = FALSE)
  expect_equal(nrow(out$unmatched), 1L)
  expect_identical(out$pairs$n_event_comparisons, rep(3L, 2L))
  expect_identical(out$pairs$site_significant_case, rep(TRUE, 2L))
  expect_identical(out$pairs$site_significant_control, rep(FALSE, 2L))
  ann <- data.table::copy(x$annotations)[type == "exon", exon_number := seq_len(.N), by = transcript_id]
  legacy <- get_pairs(attach_sequences(get_matched_events_chunked(ev, ann), x$sequences), source = "multi")
  expect_identical(legacy$n_event_comparisons, rep(3L, 2L))
  expect_identical(legacy$site_significant_case, rep(TRUE, 2L))
  expect_identical(legacy$site_significant_control, rep(FALSE, 2L))
})
