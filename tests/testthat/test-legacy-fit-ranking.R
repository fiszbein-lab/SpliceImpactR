legacy_reference <- function(types) {
  structures <- list(
    # Starts at site A; site B is its second (internal) exon.
    T_pc = matrix(c(101, 200, 501, 600, 901, 1000), byrow = TRUE, ncol = 2L),
    # Starts at site B.
    T_alt = matrix(c(501, 600, 901, 1000), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures, cds = c(151, 950))
  x$annotations[type == "exon", exon_number := seq_len(.N), by = transcript_id]
  for (id in names(types)) x$annotations[transcript_id == id, transcript_type := types[[id]]]
  x
}

legacy_pairs <- function(x, events) {
  get_pairs(attach_sequences(get_matched_events_chunked(events, x$annotations), x$sequences), source = "multi")
}

test_that("the legacy matcher ranks structural fit before protein-coding status", {
  # Two alternative first exons: A rises, B falls.
  events <- rank_events("AFE:1", "AFE", "SITE", c("101-200", "501-600"), "", c(0.5, -0.5))
  x <- legacy_reference(c(T_pc = "protein_coding", T_alt = "retained_intron"))
  pairs <- legacy_pairs(x, events)
  # T_pc only contains site B as an internal exon; the transcript that starts
  # there represents site B, so the pair is not a transcript against itself.
  expect_identical(c(pairs$transcript_id_case, pairs$transcript_id_control), c("T_pc", "T_alt"))
})

test_that("protein-coding status still breaks ties between equal fits", {
  events <- rank_events("AFE:1", "AFE", "SITE", c("101-200", "501-600"), "", c(0.5, -0.5))
  x <- legacy_reference(c(T_pc = "protein_coding", T_alt = "retained_intron"))
  # A protein-coding copy of T_alt fits site B equally well and wins the tie.
  twin <- data.table::copy(x$annotations[transcript_id == "T_alt"])
  twin[, `:=`(transcript_id = "T_twin", exon_id = sub("T_alt", "T_twin", exon_id),
              transcript_type = "protein_coding")]
  x$annotations <- data.table::rbindlist(list(x$annotations, twin), fill = TRUE)
  x$sequences <- data.table::rbindlist(list(x$sequences,
    data.table::copy(x$sequences[transcript_id == "T_alt"])[, transcript_id := "T_twin"]))
  pairs <- legacy_pairs(x, events)
  expect_identical(pairs$transcript_id_control, "T_twin")
})
