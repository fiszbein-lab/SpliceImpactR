partner_row <- function(strand = "+", flank = c(400L, 450L), id = 1L) data.table::data.table(
  sample = "S1", condition = "case", source_file = "synthetic", ID = id, GeneID = "GA",
  chr = "chr1", strand = strand, event_type = "A5SS",
  longExonStart_0base = 200L, longExonEnd = 300L,
  shortES = if (strand == "+") 200L else 240L, shortEE = if (strand == "+") 260L else 300L,
  flankingES = flank[1L], flankingEE = flank[2L],
  IJC_SAMPLE_1 = 20L, SJC_SAMPLE_1 = 10L, IncLevel1 = 2/3, IncLevelDifference = 0.3,
  PValue = 0.001, FDR = 0.01
)

test_that("A3SS/A5SS forms keep their partner exon in genomic order", {
  # On the minus strand the partner lies at lower coordinates.
  raw <- get_rmats(partner_row("-", c(50L, 100L)))
  expect_identical(raw[form == "INC", inc], "51-100;201-300")
  expect_identical(raw[form == "EXC", inc], "51-100;241-300")
  expect_identical(raw[form == "EXC", exc], "201-240")
  post <- suppressMessages(get_rmats_post_di(partner_row("-", c(50L, 100L)), event_type = "A5SS"))
  expect_identical(post[form == "INC", inc], "51-100;201-300")
  # Without partner coordinates the variable exon stands alone.
  bare <- get_rmats(partner_row(flank = c(NA_integer_, NA_integer_)))
  expect_identical(bare[form == "INC", inc], "201-300")
})

test_that("events that differ only in their partner exon stay separate", {
  two <- data.table::rbindlist(list(partner_row(), partner_row(flank = c(600L, 650L), id = 2L)))
  raw <- get_rmats(data.table::copy(two))
  expect_equal(data.table::uniqueN(raw$event_id), 2L)
  expect_setequal(raw[form == "INC", inc], c("201-300;401-450", "201-300;601-650"))
  post <- suppressMessages(get_rmats_post_di(data.table::copy(two), event_type = "A5SS"))
  expect_equal(data.table::uniqueN(post$event_id), 2L)
})

test_that("a partner exon lets ranked matching check the A5SS junction", {
  check <- function(ex, inc) {
    dt <- data.table::data.table(start = ex[, 1L], end = ex[, 2L], exon_id = paste0("X", seq_len(nrow(ex))))
    SpliceImpactR:::.si_rank_compatible(dt, SpliceImpactR:::.si_rank_intervals(inc),
      SpliceImpactR:::.si_rank_intervals(""), "A5SS", "INC", "+")
  }
  m <- function(...) matrix(c(...), byrow = TRUE, ncol = 2L)
  right <- check(m(101, 199, 201, 300, 401, 450), "201-300;401-450")
  expect_true(right$eligible)
  expect_identical(right$context, "junction_chain")
  # The same donor spliced to another exon is not this event.
  expect_false(check(m(201, 300, 601, 650), "201-300;401-450")$eligible)
  between <- check(m(201, 300, 351, 380, 401, 450), "201-300;401-450")
  expect_identical(between$reason, "splice_junction_chain_mismatch")
  # A lone variable exon still verifies the splice site only.
  expect_identical(check(m(201, 300, 601, 650), "201-300")$context, "splice_site_only")
})

test_that("A5SS frame checks stay at the variable exon when the partner is listed", {
  for (strand in c("+", "-")) {
    x <- rank_test_reference(strand)
    partner <- if (strand == "+") "1001-1300" else "701-1000"
    with_partner <- data.table::copy(x$events)
    joined <- if (strand == "+") paste(with_partner$inc, partner, sep = ";") else
      paste(partner, with_partner$inc, sep = ";")
    with_partner[, inc := joined]
    frames <- lapply(list(x$events, with_partner), function(ev) {
      pairs <- get_ranked_pairs(ev, x$annotations, x$sequences, verbose = FALSE)$pairs
      compare_sequence_frame(pairs, x$annotations)
    })
    expect_identical(frames[[2L]]$structural_context_case, "junction_chain")
    expect_identical(frames[[2L]]$frame_check_exon1, frames[[1L]]$frame_check_exon1)
    expect_identical(frames[[2L]]$frame_check_exon2, frames[[1L]]$frame_check_exon2)
    expect_identical(frames[[2L]]$frame_call, frames[[1L]]$frame_call)
  }
})
