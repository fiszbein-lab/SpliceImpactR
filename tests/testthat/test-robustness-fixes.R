test_that("a form over max_candidates keeps its best-supported candidates, labelled", {
  x <- rank_test_reference()
  full <- get_ranked_pairs(x$events, x$annotations, x$sequences, verbose = FALSE)
  expect_identical(c(full$pairs$candidates_dropped_case, full$pairs$candidates_dropped_control), c(0L, 0L))
  expect_message(
    out <- get_ranked_pairs(x$events, x$annotations, x$sequences, max_candidates = 1L, verbose = FALSE),
    "1 comparison\\(s\\) had a form with more than max_candidates = 1")
  expect_equal(nrow(out$pairs), 1L)
  expect_identical(c(out$pairs$candidates_dropped_case, out$pairs$candidates_dropped_control), c(1L, 0L))
  # Candidate counts in `events` are taken before the cut.
  expect_equal(out$events$n_candidates, c(2L, 1L))
  expect_identical(out$settings$max_candidates, 1L)
})

test_that("the candidate cut keeps annotated CDS before transcript support", {
  # Skipped exon: the two inclusion transcripts share a structure; the coding
  # one has weaker support (TSL 4, a lower tier) than the non-coding one (TSL 1).
  # A support-first cut would keep the non-coding transcript.
  x <- rank_structure_reference(list(
    T_nc = matrix(c(101, 200, 301, 400, 501, 600), byrow = TRUE, ncol = 2L),
    T_cd = matrix(c(101, 200, 301, 400, 501, 600), byrow = TRUE, ncol = 2L),
    T_skip = matrix(c(101, 200, 501, 600), byrow = TRUE, ncol = 2L)),
    cds = list(T_nc = NULL, T_cd = c(151, 550), T_skip = c(151, 550)))
  x$annotations[transcript_id == "T_cd", transcript_support_level := "4"]
  events <- rank_events("SE:1", "SE", c("INC", "EXC"), c("101-200;301-400;501-600", "101-200;501-600"),
                        c("", "301-400"), c(0.5, -0.5))
  out <- suppressMessages(get_ranked_pairs(events, x$annotations, x$sequences, max_candidates = 1L,
                                           verbose = FALSE))
  expect_identical(out$pairs$transcript_id_case, "T_cd")
  expect_identical(out$pairs$candidates_dropped_case, 1L)
})

test_that("ORF input errors stop before matching and name the events", {
  x <- rank_test_reference()
  run <- function(events) get_ranked_pairs(events, x$annotations, x$sequences, verbose = FALSE)
  strand <- data.table::copy(x$events)[1L, strand := "*"]
  expect_error(run(strand), "strands \\(found \\*\\)\\. Affected: E\\.")
  gene <- data.table::copy(x$events)[, gene_id := NA_character_]
  expect_error(run(gene), "cannot be missing: gene_id\\. Affected: E\\.")
  interval <- data.table::copy(x$events)[1L, inc := "400-101"]
  expect_error(run(interval), "start-end spans.*Affected: E\\.")
  form <- data.table::copy(x$events)[1L, form := "OTHER"]
  expect_error(run(form), "found OTHER.*Affected: E\\.")
})

test_that("import_di_table rejects unknown strands instead of guessing +", {
  df <- data.frame(gene_id = "G1", chr = "chr1", strand = c("+", "*"), inc = c("1-10", "21-30"),
                   exc = "", delta_psi = c(0.2, -0.2), p.value = 0.01, event_id = c("A", "B"))
  expect_error(import_di_table(df), "strand must be \\+ or - \\(found \\*\\)\\. Affected: B\\.")
})

test_that("get_rmats_hit rejects unsupported event types", {
  sf <- data.frame(path = tempdir(), sample_name = "S1", condition = "case")
  expect_error(get_rmats_hit(sf, event_types = c("SE", "alt3")), "unsupported event_types: alt3")
})

test_that("old exact-name cache records are ignored and nothing is written outside the cache", {
  root <- tempfile("bfc-old-")
  work <- tempfile("cwd-")
  dir.create(work)
  on.exit(unlink(c(root, work), recursive = TRUE), add = TRUE)
  original_dir <- setwd(work)
  on.exit(setwd(original_dir), add = TRUE, after = FALSE)
  bfc <- SpliceImpactR:::.si_bfc(root)
  # Earlier versions created records whose bare file name resolves against
  # the working directory.
  old <- BiocFileCache::bfcnew(bfc, rname = "key", rtype = "local", ext = ".rds", fname = "exact")
  saveRDS("stale", unname(old))
  expect_null(SpliceImpactR:::.si_bfc_get_rds(bfc, "key"))
  path <- SpliceImpactR:::.si_bfc_put_rds(bfc, "key", "fresh")
  expect_identical(normalizePath(dirname(path)), normalizePath(BiocFileCache::bfccache(bfc)))
  expect_identical(SpliceImpactR:::.si_bfc_get_rds(bfc, "key"), "fresh")
  # The old file is left as it was and nothing new appears in the working folder.
  expect_identical(list.files(work), basename(unname(old)))
  expect_identical(readRDS(file.path(work, basename(unname(old)))), "stale")
})

test_that("cache fingerprints ignore the R version and the table class", {
  x <- list(a = 1:3, b = c("x", NA), d = data.frame(u = c(1.5, 2), v = c("p", "q")))
  # Fixed content gives a fixed hash: the R version recorded in the
  # serialization header is not hashed.
  expect_identical(SpliceImpactR:::.si_md5_object(x), "26be85d9f069573dec0d31e511f10c99")
  y <- x
  y$d <- data.table::as.data.table(y$d)
  expect_identical(SpliceImpactR:::.si_md5_object(y), SpliceImpactR:::.si_md5_object(x))
})

test_that("get_pairs names conflicting event IDs in input order", {
  dt <- data.table::data.table(
    event_id = rep(c("E1", "E3", "E2"), each = 2L), gene_id = c("G", "G", "G", "H", "G", "G"),
    chr = "chr1", strand = c("+", "+", "+", "+", "+", "-"), event_type = "SE",
    form = c("INC", "EXC"), transcript_id = c("T1", "T2"), exons = "X", protein_id = "P",
    inc = c("101-200", "301-400"), exc = "", delta_psi = c(0.3, -0.3), p.value = 0.001,
    padj = 0.01, transcript_seq = "ATG", protein_seq = "M")
  expect_error(get_pairs(dt, source = "multi"), "Conflicting IDs: E3, E2$")
  expect_equal(nrow(get_pairs(dt[event_id == "E1"], source = "multi")), 1L)
})

test_that("legacy matched output keeps event columns before the matched transcript", {
  x <- rank_test_reference()
  x$annotations[type == "exon", exon_number := seq_len(.N), by = transcript_id]
  m <- get_matched_events_chunked(x$events, x$annotations, verbose = FALSE)
  expect_identical(names(m), c(names(x$events), "transcript_id", "exons", "transcript_type"))
})

test_that("legacy matcher reports chunks as messages unless verbose = FALSE", {
  x <- rank_test_reference()
  x$annotations[type == "exon", exon_number := seq_len(.N), by = transcript_id]
  expect_message(get_matched_events_chunked(x$events, x$annotations), "^Chunk 1/1: rows 1\\.\\.")
  expect_silent(get_matched_events_chunked(x$events, x$annotations, verbose = FALSE))
})
