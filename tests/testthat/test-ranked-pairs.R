test_that("joint ORF ranking chooses a closer coding background on both strands", {
  for (strand in c("+", "-")) {
    x <- rank_test_reference(strand)
    before <- lapply(x, data.table::copy)
    out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
    expect_identical(out$pairs$transcript_id_case, "Z_INC_close")
    expect_identical(out$pairs$transcript_id_control, "E_EXC")
    expect_equal(out$pairs$context_similarity, 1)
    expect_identical(out$pairs$similarity_source, "cds_sequences")
    expect_equal(out$pairs$n_candidate_pairs, 2L)
    expect_true(out$pairs$context_score_margin > 0)
    expect_true(all(out$rankings$context_similarity <= 1))
    expect_equal(out$events$n_candidates, c(2L, 1L))
    expect_equal(out$pairs$delta_psi_case, 0.3)
    expect_identical(x, before)
    shuffled <- get_ranked_pairs(x$events[2:1], x$annotations[nrow(x$annotations):1],
                                 x$sequences[3:1], verbose = FALSE)
    expect_equal(shuffled$pairs, out$pairs)
    old <- get_pairs(attach_sequences(get_matched_events_chunked(x$events, x$annotations), x$sequences), source = "multi")
    expect_identical(old$transcript_id_case, "A_INC_distant")
    expect_false(identical(old$transcript_id_case, out$pairs$transcript_id_case))
  }
})

test_that("context similarity compares exonic positions in the coding window", {
  pos <- SpliceImpactR:::.si_rank_positions
  ir <- function(start, end) IRanges::IRanges(start, end)
  a <- ir(c(101, 301), c(200, 400))
  b <- ir(c(101, 301), c(200, 350))
  expect_equal(pos(a, b, ir(101, 400))$similarity, 150 / 200)
  expect_equal(pos(a, b, ir(101, 400))$coverage_control, 1)
  expect_equal(pos(a, b, ir(101, 400), mask = ir(351, 400))$similarity, 1)
  expect_equal(pos(a, b, ir(101, 200))$similarity, 1)
  expect_true(is.na(pos(a, b, ir(1, 50))$similarity))
})

test_that("a frameshifting skip is not outranked by a frame-restoring isoform", {
  # Skipping the 100-nt exon shifts the frame, so the plain skip's annotated
  # CDS stops early. The rescue restores the frame with a 2-nt shorter acceptor
  # in a later exon; its CDS stays long but its structure differs.
  events <- rank_events("SE:1", "SE", c("INC", "EXC"),
                        c("101-200;301-400;501-700", "101-200;501-700"), c("", "301-400"), c(0.3, -0.3))
  structures <- list(
    INC = matrix(c(101, 200, 301, 400, 501, 700, 801, 1000), byrow = TRUE, ncol = 2L),
    P_skip = matrix(c(101, 200, 501, 700, 801, 1000), byrow = TRUE, ncol = 2L),
    A_rescue = matrix(c(101, 200, 501, 700, 803, 1000), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures,
    cds = list(INC = c(151, 950), P_skip = c(151, 560), A_rescue = c(151, 950)))
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_control, "P_skip")
  expect_equal(out$pairs$context_similarity, 1)
  expect_identical(out$pairs$context_window, "pair_cds_span")
  rescue <- out$rankings[transcript_id_control == "A_rescue"]
  expect_equal(rescue$context_similarity, 398 / 400)
})

test_that("an NMD-annotated frameshifting skip still counts as coding", {
  # Ensembl annotates the plain skip as nonsense-mediated decay because its
  # premature stop lies upstream of the last junction. It still has a CDS, so
  # its biotype must not hand the comparison to the frame-restoring isoform.
  events <- rank_events("SE:1", "SE", c("INC", "EXC"),
                        c("101-200;301-400;501-700", "101-200;501-700"), c("", "301-400"), c(0.3, -0.3))
  structures <- list(
    INC = matrix(c(101, 200, 301, 400, 501, 700, 801, 1000), byrow = TRUE, ncol = 2L),
    P_skip = matrix(c(101, 200, 501, 700, 801, 1000), byrow = TRUE, ncol = 2L),
    A_rescue = matrix(c(101, 200, 501, 700, 803, 1000), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures,
    cds = list(INC = c(151, 950), P_skip = c(151, 560), A_rescue = c(151, 950)))
  x$annotations[transcript_id == "P_skip", transcript_type := "nonsense_mediated_decay"]
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_control, "P_skip")
  expect_identical(out$pairs$coding_count, 2L)
  expect_identical(out$settings$protocol_version, 5L)
  # The biotype stays visible as a label.
  expect_identical(out$pairs$transcript_type_control, "nonsense_mediated_decay")
  expect_identical(out$pairs$transcript_type_case, "protein_coding")
  expect_identical(out$candidates[transcript_id == "P_skip", unique(transcript_type)], "nonsense_mediated_decay")

  # A protein_coding biotype without an annotated CDS is noncoding.
  x$annotations[transcript_id == "P_skip", `:=`(transcript_type = "protein_coding",
                                                cds_gen_start = NA_real_, cds_gen_stop = NA_real_)]
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_false(out$candidates[transcript_id == "P_skip", unique(coding)])
  expect_identical(out$pairs$transcript_id_control, "A_rescue")
})

test_that("eligibility enforces exon topology before ranking", {
  check <- function(ex, inc, exc = "", type = "SE", form = "INC", strand = "+") {
    dt <- data.table::data.table(start = ex[, 1L], end = ex[, 2L], exon_id = paste0("X", seq_len(nrow(ex))))
    SpliceImpactR:::.si_rank_compatible(dt, SpliceImpactR:::.si_rank_intervals(inc),
      SpliceImpactR:::.si_rank_intervals(exc), type, form, strand)
  }
  spliced <- matrix(c(101, 199, 301, 399, 501, 599), byrow = TRUE, ncol = 2L)
  fused <- matrix(c(101, 399, 501, 599), byrow = TRUE, ncol = 2L)
  skip <- spliced[c(1L, 3L), ]
  expect_true(check(spliced, "101-199;301-399;501-599")$eligible)
  expect_false(check(fused, "101-199;301-399;501-599")$eligible)
  expect_true(check(skip, "101-199;501-599", "301-399", form = "EXC")$eligible)
  expect_false(check(spliced, "101-199;501-599", "301-399", form = "EXC")$eligible)
  expect_true(check(spliced, "101-199;301-399;501-599", "701-799", type = "MXE")$eligible)
  expect_true(check(fused, "101-199;200-300;301-399", type = "RI")$eligible)
  expect_false(check(spliced, "101-199;200-300;301-399", type = "RI")$eligible)
  expect_true(check(spliced[1:2, ], "101-199;301-399", "200-300", type = "RI", form = "EXC")$eligible)
  expect_true(check(skip, "101-199;501-599", "301-399", form = "SITE")$eligible)
  expect_true(check(fused, "101-199;200-300;301-399", type = "RI", form = "SITE")$eligible)
  expect_true(check(spliced[1:2, ], "101-199;301-399", "200-300", type = "RI", form = "SITE")$eligible)
  expect_true(check(spliced, "101-199", type = "AFE")$eligible)
  expect_true(check(spliced, "501-599", type = "AFE", strand = "-")$eligible)
  expect_true(check(spliced, "101-199", type = "HLE", strand = "-")$eligible)
  expect_false(check(spliced, "301-399", type = "AFE")$eligible)
  expect_false(check(spliced, "101-199", type = "A3SS")$eligible)
  expect_true(check(spliced, "301-399", type = "A3SS")$eligible)
  expect_false(check(spliced, "301-400", type = "A5SS")$eligible)
  expect_true(check(spliced, "501-599", type = "A5SS", strand = "-")$eligible)
  expect_false(check(spliced, "101-199", type = "A5SS", strand = "-")$eligible)
})

test_that("relaxed matching keeps event boundaries exact but frees outer ends", {
  check <- function(ex, inc, exc = "", type = "SE", form = "INC", strand = "+") {
    dt <- data.table::data.table(start = ex[, 1L], end = ex[, 2L], exon_id = paste0("X", seq_len(nrow(ex))))
    run <- function(relaxed) SpliceImpactR:::.si_rank_compatible(dt, SpliceImpactR:::.si_rank_intervals(inc),
      SpliceImpactR:::.si_rank_intervals(exc), type, form, strand, relaxed = relaxed)$eligible
    c(exact = run(FALSE), relaxed = run(TRUE))
  }
  m <- function(...) matrix(c(...), byrow = TRUE, ncol = 2L)
  freed <- c(exact = FALSE, relaxed = TRUE)
  rejected <- c(exact = FALSE, relaxed = FALSE)
  se <- "101-199;301-399;501-599"
  expect_identical(check(m(151, 199, 301, 399, 501, 599), se), freed)              # alternative TSS
  expect_identical(check(m(101, 199, 301, 399, 501, 560), se), freed)              # alternative PAS
  expect_identical(check(m(101, 140, 161, 199, 301, 399, 501, 599), se), freed)    # flank split by an intron
  expect_identical(check(m(101, 195, 301, 399, 501, 599), se), rejected)           # donor moved
  expect_identical(check(m(101, 199, 301, 395, 501, 599), se), rejected)           # skipped exon changed
  expect_identical(check(m(101, 199, 501, 560), "101-199;501-599", "301-399", form = "EXC"), freed)
  expect_identical(check(m(151, 399), "101-199;200-300;301-399", type = "RI"), freed)
  expect_identical(check(m(200, 399), "101-199;200-300;301-399", type = "RI"), rejected)
  expect_identical(check(m(101, 199, 331, 450, 601, 700), "301-450", type = "A5SS"), freed)
  expect_identical(check(m(101, 199, 331, 445, 601, 700), "301-450", type = "A5SS"), rejected)
  expect_identical(check(m(101, 199, 301, 420, 601, 700), "301-450", type = "A3SS"), freed)
  expect_identical(check(m(101, 199, 301, 420, 601, 700), "301-450", type = "A5SS", strand = "-"), freed)
  expect_identical(check(m(121, 300, 401, 500), "101-300", type = "AFE"), freed)
  expect_identical(check(m(121, 300), "101-300", type = "AFE"), rejected)          # no donor
  expect_identical(check(m(11, 50, 121, 300, 401, 500), "101-300", type = "AFE"), rejected)
  expect_identical(check(m(101, 200, 401, 520), "401-500", type = "AFE", strand = "-"), freed)
  expect_identical(check(m(121, 300, 401, 500), "101-300", type = "ALE", strand = "-"), freed)
  expect_identical(check(m(101, 200, 405, 520), "401-500", type = "ALE"), rejected)
})

test_that("relaxed outer boundaries resolve events and lose ties to exact ones", {
  inc <- "101-199;301-399;501-599"
  events <- rank_events("SE:1", "SE", c("INC", "EXC"), c(inc, "101-199;501-599"),
                        c("", "301-399"), c(0.3, -0.3))
  # A_short differs from E_exact only by a transcription start in the 5' UTR.
  structures <- list(
    I_exact = matrix(c(101, 199, 301, 399, 501, 599), byrow = TRUE, ncol = 2L),
    E_exact = matrix(c(101, 199, 501, 599), byrow = TRUE, ncol = 2L),
    A_short = matrix(c(121, 199, 501, 599), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures, cds = c(151, 550))
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_control, "E_exact")
  expect_identical(out$pairs$structural_match_control, "exact")
  expect_equal(out$rankings[transcript_id_control == "A_short", n_inexact], 1L)
  expect_true(all(out$rankings$context_similarity == 1))
  expect_identical(out$candidates[transcript_id == "A_short" & form == "EXC", structural_match], "relaxed")

  # Without the exact transcript, the relaxed one represents the form.
  x <- rank_structure_reference(structures[c("I_exact", "A_short")], cds = c(151, 550))
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_control, "A_short")
  expect_identical(out$pairs$structural_match_control, "relaxed")
})

test_that("terminal sites sharing a splice site stay distinguishable", {
  events <- rank_events("AFE:1", "AFE", "SITE", c("101-300", "201-300"), "", c(0.3, -0.3))
  structures <- list(T_A = matrix(c(101, 300, 401, 500), byrow = TRUE, ncol = 2L),
                     T_B = matrix(c(201, 300, 401, 500), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures)
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(c(out$pairs$transcript_id_case, out$pairs$transcript_id_control), c("T_A", "T_B"))
  expect_identical(c(out$pairs$structural_match_case, out$pairs$structural_match_control), c("exact", "exact"))

  # Only a relaxed transcript start remains for the second site.
  structures$T_B <- matrix(c(221, 300, 401, 500), byrow = TRUE, ncol = 2L)
  x <- rank_structure_reference(structures)
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_control, "T_B")
  expect_identical(out$pairs$structural_match_control, "relaxed")
})

test_that("relaxed candidates must contain their form's alternative region", {
  events <- rank_events("A3SS:1", "A3SS", c("INC", "EXC"), c("301-450", "351-450"),
                        c("", "301-350"), c(0.3, -0.3))
  # X keeps the long acceptor on a separate short exon, so 321-350 is intronic.
  structures <- list(L = matrix(c(101, 200, 301, 450, 601, 700), byrow = TRUE, ncol = 2L),
                     S = matrix(c(101, 200, 351, 450, 601, 700), byrow = TRUE, ncol = 2L),
                     X = matrix(c(101, 200, 301, 320, 351, 450, 601, 700), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures)
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(c(out$pairs$transcript_id_case, out$pairs$transcript_id_control), c("L", "S"))
  expect_identical(out$pair_rejections[transcript_id_case == "X", reason], "event_region_not_covered")
  x <- rank_structure_reference(structures[c("S", "X")])
  out <- get_ranked_pairs(events, x$annotations, x$sequences, fallback = FALSE, verbose = FALSE)
  expect_equal(nrow(out$pairs), 0L)
  expect_identical(out$unmatched$reason, "event_region_not_covered")
  expect_identical(out$unmatched$fallback, "not_run")
  # The fallback keeps X as a labelled approximation: it covers 20 of the 50
  # alternative bases, which S lacks.
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(c(out$pairs$transcript_id_case, out$pairs$transcript_id_control), c("X", "S"))
  expect_identical(out$pairs$matching_tier, "approximate")
  expect_identical(out$pairs$structural_match_case, "approximate")
  expect_identical(out$pairs$structural_match_control, "exact")
  expect_equal(out$pairs$event_agreement, 0.4)
})

test_that("the approximate fallback uses legacy gates and needs an event difference", {
  approximate <- SpliceImpactR:::.si_rank_approximate
  ir <- function(start, end) IRanges::IRanges(start, end)
  exons <- ir(c(101, 301), c(200, 400))
  none <- IRanges::IRanges()
  expect_true(approximate(exons, ir(196, 295), none))                # 5 of 100 bases overlap
  expect_false(approximate(exons, ir(197, 296), none))               # 4 of 100
  expect_false(approximate(exons, ir(101, 200), ir(351, 450)))       # forbidden half covered
  expect_true(approximate(exons, ir(101, 200), ir(397, 496)))        # forbidden 4% covered

  events <- rank_events("SE:1", "SE", c("INC", "EXC"),
                        c("101-200;301-400;501-600", "101-200;501-600"), c("", "301-400"), c(0.3, -0.3))
  # No transcript has the skipped exon's exact boundaries.
  structures <- list(I = matrix(c(101, 200, 311, 400, 501, 600), byrow = TRUE, ncol = 2L),
                     E = matrix(c(101, 200, 501, 600), byrow = TRUE, ncol = 2L))
  x <- rank_structure_reference(structures)
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(c(out$pairs$transcript_id_case, out$pairs$transcript_id_control), c("I", "E"))
  expect_identical(out$pairs$matching_tier, "approximate")
  expect_equal(out$pairs$event_agreement, 0.9)
  # A pair without the event difference is never selected.
  x <- rank_structure_reference(structures["E"])
  out <- get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  expect_equal(nrow(out$pairs), 0L)
  expect_identical(out$unmatched$fallback, "missing_form_candidate")
})

test_that("coding and TSL priorities are explicit and missing sequence falls back", {
  x <- rank_test_reference()
  x$annotations[transcript_id == "Z_INC_close", transcript_support_level := "2"]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_case, "Z_INC_close")
  x$annotations[transcript_id == "Z_INC_close", transcript_support_level := "4"]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_case, "A_INC_distant")
  # Without an annotated CDS a transcript is noncoding, whatever its biotype.
  x$annotations[transcript_id == "A_INC_distant", c("cds_gen_start", "cds_gen_stop") := NA_integer_]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_case, "Z_INC_close")
  expect_false(out$candidates[transcript_id == "A_INC_distant", unique(coding)])
  x$sequences[, transcript_seq := NA_character_]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$similarity_source, "cds_annotation")
  expect_equal(out$pairs$context_similarity, 1)
  x$sequences[, transcript_seq := "ACGTA"]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_true(all(out$candidates[coding == TRUE, sequence_status] == "sequence_length_mismatch"))
  expect_identical(unique(out$candidates[coding == FALSE, sequence_status]), "no_cds_annotation")
})

test_that("a missing TSL column is treated as unknown support", {
  x <- rank_test_reference()
  with_tsl <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  annotations <- data.table::copy(x$annotations)[, transcript_support_level := NULL]
  out <- get_ranked_pairs(x$events, annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_case, with_tsl$pairs$transcript_id_case)
  expect_identical(out$pairs$transcript_id_control, with_tsl$pairs$transcript_id_control)
  expect_identical(with_tsl$pairs$tsl_tier, 0L)
  expect_identical(out$pairs$tsl_tier, 2L)
  expect_false("transcript_support_level" %in% names(annotations))
})

test_that("ranked matching does not require a transcript biotype column", {
  x <- rank_test_reference()
  with_type <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  annotations <- data.table::copy(x$annotations)[, transcript_type := NULL]
  out <- get_ranked_pairs(x$events, annotations, x$sequences, verbose = FALSE)
  labels <- c("transcript_type_case", "transcript_type_control")
  expect_equal(out$pairs[, !..labels], with_type$pairs[, !..labels])
  expect_equal(out$rankings, with_type$rankings)
  expect_identical(unlist(with_type$pairs[, ..labels], use.names = FALSE), c("protein_coding", "protein_coding"))
  expect_identical(unlist(out$pairs[, ..labels], use.names = FALSE), c(NA_character_, NA_character_))
})

test_that("candidate limits, duplicates and unmatched forms are never hidden", {
  x <- rank_test_reference()
  expect_message(get_ranked_pairs(x$events, x$annotations, x$sequences, max_candidates = 1L,
                                  verbose = FALSE), "max_candidates = 1")
  sq <- data.table::rbindlist(list(x$sequences, x$sequences[1L][, protein_seq := "OTHER"]))
  expect_error(get_ranked_pairs(x$events, x$annotations, sq), "Conflicting sequence")
  expect_error(get_ranked_pairs(data.table::rbindlist(list(x$events, x$events[1L])), x$annotations, x$sequences), "one differential row")
  x$events[form == "INC", inc := "101-500"]
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, fallback = FALSE, verbose = FALSE)
  expect_equal(nrow(out$events), 2L)
  expect_equal(nrow(out$pairs), 0L)
  expect_true(all(c("transcript_id_case", "matching_tier", "transcript_type_case") %in% names(out$pairs)))
  expect_identical(out$unmatched$reason, "missing_form_candidate")
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$pairs$matching_tier, "approximate")
  expect_identical(out$unmatched$reason, character())
  out <- get_ranked_pairs(x$events[1L], x$annotations, x$sequences, verbose = FALSE)
  expect_identical(out$events$matching_status, "no_opposite_direction_form")
})

test_that("ranked results retain the downstream table and S4 contracts", {
  x <- rank_test_reference()
  out <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  compare <- compare_sequence_frame(out$pairs, x$annotations)
  expect_equal(nrow(compare), 1L)
  expect_true("frame_call" %in% names(compare))
  obj <- as_splice_impact_result(res = x$events, res_di = x$events, matched = out$matched,
                                 hits_final = out$pairs)
  expect_true(methods::validObject(obj))
  expect_identical(as_dt_from_s4(obj, "paired_hits")$matching_pair_id, out$pairs$matching_pair_id)
  expect_equal(get_ranked_pairs(obj, x$annotations, x$sequences, verbose = FALSE)$pairs, out$pairs)
  paired <- data.table::copy(x$events)[, delta_psi := NULL]
  expect_equal(get_ranked_pairs(paired, x$annotations, x$sequences, source = "paired", verbose = FALSE)$pairs$transcript_id_case,
               "Z_INC_close")
})

test_that("ties, indistinguishable forms and multiple opposing forms are explicit", {
  x <- rank_test_reference()
  duplicate <- data.table::copy(x$annotations[transcript_id == "Z_INC_close"])
  duplicate[, `:=`(transcript_id = "Zzz_tied", exon_id = sub("Z_INC_close", "Zzz_tied", exon_id))]
  sq <- data.table::copy(x$sequences[transcript_id == "Z_INC_close"])[, transcript_id := "Zzz_tied"]
  ann <- data.table::rbindlist(list(x$annotations, duplicate), fill = TRUE)
  sequences <- data.table::rbindlist(list(x$sequences, sq))
  out <- get_ranked_pairs(x$events, ann, sequences, verbose = FALSE)
  expect_identical(out$pairs$transcript_id_case, "Z_INC_close")
  expect_equal(out$pairs$n_tied_best, 2L)
  expect_equal(out$pairs$context_score_margin, 0)
  ambiguous <- data.table::copy(x$events)[, `:=`(event_type = "GENERIC", inc = "101-199", exc = "")]
  out <- get_ranked_pairs(ambiguous, ann, sequences, verbose = FALSE)
  expect_equal(nrow(out$pairs), 0L)
  expect_true("both_forms_compatible" %in% out$pair_rejections$reason)
  expect_identical(out$unmatched$reason, "no_form_distinguishing_pair")

  extra_ann <- data.table::copy(x$annotations[transcript_id == "E_EXC"])
  extra_ann[, `:=`(transcript_id = "E_EXC2", exon_id = sub("E_EXC", "E_EXC2", exon_id))]
  extra_ann[type == "exon" & exon_number == 1L, `:=`(end = 196L, cds_gen_stop = 196L)]
  extra_seq <- data.table::copy(x$sequences[transcript_id == "E_EXC"])
  extra_seq[, `:=`(transcript_id = "E_EXC2", transcript_seq = substr(transcript_seq, 1L, 396L))]
  extra_event <- data.table::copy(x$events[form == "EXC"])[, `:=`(inc = "101-196", exc = "197-400")]
  out <- get_ranked_pairs(data.table::rbindlist(list(x$events, extra_event)),
    data.table::rbindlist(list(x$annotations, extra_ann)),
    data.table::rbindlist(list(x$sequences, extra_seq)), verbose = FALSE)
  expect_equal(nrow(out$pairs), 2L)
  expect_equal(nrow(out$matched), 4L)
  expect_false(anyDuplicated(out$pairs$matching_pair_id) > 0L)
})

test_that("the opt-in wrapper returns diagnostics and complete consequences", {
  x <- rank_test_reference()
  pf <- data.table::data.table(ensembl_transcript_id = x$sequences$transcript_id,
    ensembl_peptide_id = x$sequences$protein_id, database = "pfam", clean_name = "Test",
    name = "Test;chr1:128-160", feature_id = "PF_TEST", start = 10L, stop = 20L,
    chr = "chr1", strand = "+", alt_name = "Test")
  ef <- get_exon_features(x$annotations, pf)
  ppi <- data.table::data.table(geneA = character(), geneB = character(),
    DDI = logical(), DMI = logical(), ddi_for_A = character(), ddi_for_B = character(),
    dmi_for_A = character(), dmi_for_B = character())
  for (cls in c("data.table", "S4")) {
    out <- get_splicing_impact(res = x$events,
      annotation_df = x[c("annotations", "sequences")], protein_feature_total = pf,
      exon_features = ef, ppi = ppi, matching = "orf", return_class = cls, verbose = FALSE)
    if (cls == "S4") {
      expect_true(methods::validObject(out))
      expect_identical(out@metadata$matching, "orf")
      expect_equal(nrow(out@metadata$matching_diagnostics$rankings), 2L)
    } else {
      expect_identical(out$hits_final$transcript_id_case, "Z_INC_close")
      expect_identical(out$hits_final$transcript_type_case, "protein_coding")
      expect_equal(nrow(out$matching$rankings), 2L)
    }
  }
  # The wrapper forwards the fallback choice.
  unmatched <- data.table::copy(x$events)[form == "INC", inc := "101-500"]
  fit <- function(...) get_splicing_impact(res = unmatched,
    annotation_df = x[c("annotations", "sequences")], protein_feature_total = pf,
    exon_features = ef, ppi = ppi, matching = "orf", verbose = FALSE, ...)
  expect_equal(nrow(fit(matching_fallback = FALSE)$hits_final), 0L)
  expect_identical(fit()$hits_final$matching_tier, "approximate")
  expect_identical(formals(get_splicing_impact)$matching_fallback, formals(get_ranked_pairs)$fallback)

  # Re-running on an S4 result replaces earlier provenance instead of appending.
  run <- function(matching, data = NULL, res = NULL) get_splicing_impact(
    data = data, res = res, annotation_df = x[c("annotations", "sequences")],
    protein_feature_total = pf, exon_features = ef, ppi = ppi, matching = matching,
    return_class = "S4", verbose = FALSE)
  orf <- run("orf", data = run("legacy", res = x$events))
  expect_identical(anyDuplicated(names(orf@metadata)), 0L)
  expect_identical(orf@metadata$matching, "orf")
  expect_false(is.null(orf@metadata$matching_diagnostics))
  legacy <- run("legacy", data = orf)
  expect_identical(anyDuplicated(names(legacy@metadata)), 0L)
  expect_identical(legacy@metadata$matching, "legacy")
  expect_null(legacy@metadata$matching_diagnostics)
})

test_that("S4 input uses significance-filtered events when present", {
  x <- rank_test_reference()
  nonsig <- data.table::copy(x$events)[, `:=`(event_id = "E_nonsig", padj = 0.9, p.value = 0.5)]
  res <- data.table::rbindlist(list(x$events, nonsig))
  filtered <- get_ranked_pairs(keep_sig_pairs(as_splice_impact_result(res = res)),
                               x$annotations, x$sequences, verbose = FALSE)
  expect_identical(unique(filtered$pairs$event_id), "E")
  # Without res_di, every differential event is ranked.
  all_events <- get_ranked_pairs(as_splice_impact_result(res = res),
                                 x$annotations, x$sequences, verbose = FALSE)
  expect_setequal(unique(all_events$pairs$event_id), c("E", "E_nonsig"))
})

