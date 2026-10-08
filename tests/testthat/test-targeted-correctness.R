test_that("overlapping genes retain their own transcripts across chunk sizes", {
  ann <- data.table::data.table(
    type = rep(c("transcript", "exon"), 2), row_uid = 1:4,
    chr = "chr1", start = 101L, end = 200L, strand = "+",
    gene_id = rep(c("GA", "GB"), each = 2),
    gene_name = rep(c("GA", "GB"), each = 2),
    transcript_id = rep(c("TA", "TB"), each = 2),
    transcript_name = rep(c("TA", "TB"), each = 2),
    transcript_type = "protein_coding",
    protein_id = rep(c("PA", "PB"), each = 2),
    exon_id = c(NA_character_, "EA", NA_character_, "EB"),
    exon_number = rep(c(NA_integer_, 1L), 2),
    transcript_support_level = c("5", "5", "1", "1")
  )
  events <- data.table::data.table(
    event_id = c("a", "b"), event_type = "AFE", form = "SITE",
    gene_id = c("GA", "GB"), chr = "chr1", strand = "+",
    inc = "101-200", exc = "", delta_psi = c(0.3, -0.3),
    p.value = 0.001, padj = 0.01,
    n_samples = 4L, n_control = 2L, n_case = 2L
  )

  for (chunk_size in c(1L, 2L, 10L)) {
    for (rows in list(1:2, 2:1)) {
      matched <- get_matched_events_chunked(
        data.table::copy(events[rows]), ann, chunk_size = chunk_size
      )
      matched <- matched[order(event_id)]
      expect_identical(matched$gene_id, c("GA", "GB"))
      expect_identical(matched$transcript_id, c("TA", "TB"))
    }
  }
})

test_that("legacy matching uses TSL only to break structural ties", {
  transcript <- function(id, tsl, exons) {
    data.table::rbindlist(list(
      data.table::data.table(type = "transcript", transcript_id = id,
        start = min(exons[, 1L]), end = max(exons[, 2L]),
        exon_id = NA_character_, exon_number = NA_integer_),
      data.table::data.table(type = "exon", transcript_id = id,
        start = exons[, 1L], end = exons[, 2L],
        exon_id = paste0(id, "_", seq_len(nrow(exons))),
        exon_number = seq_len(nrow(exons)))
    ))[, transcript_support_level := tsl][]
  }
  reference <- function(...) {
    data.table::rbindlist(list(...))[, `:=`(
      row_uid = .I, chr = "chr1", strand = "+", gene_id = "G", gene_name = "G",
      transcript_name = transcript_id, transcript_type = "protein_coding",
      protein_id = paste0("P_", transcript_id)
    )][]
  }
  # A5SS: INC is the long exon; EXC excludes the alternative 5' extension.
  events <- data.table::data.table(
    event_id = "A5SS:1", event_type = "A5SS", form = c("INC", "EXC"),
    gene_id = "G", chr = "chr1", strand = "+",
    inc = c("1000-1200", "1000-1150"), exc = c("", "1151-1200"),
    delta_psi = c(0.4, -0.4), p.value = 0.001, padj = 0.01
  )
  pick <- function(ann) {
    matched <- get_matched_events_chunked(data.table::copy(events), ann)
    stats::setNames(matched$transcript_id, matched$form)
  }
  short <- rbind(c(100L, 200L), c(1000L, 1150L), c(2000L, 2100L))
  long <- rbind(c(100L, 200L), c(1000L, 1200L), c(2000L, 2100L))

  # A better-supported transcript lacking the extension must not take the INC form.
  fit <- pick(reference(transcript("T_long", "2", long),
                        transcript("T_short", "1", short)))
  expect_identical(fit[["INC"]], "T_long")
  expect_identical(fit[["EXC"]], "T_short")

  # Among structurally identical candidates, the better TSL wins regardless of ID order.
  tie <- pick(reference(transcript("T_a", "2", long), transcript("T_b", "1", long),
                        transcript("T_short", "1", short)))
  expect_identical(tie[["INC"]], "T_b")
  expect_identical(tie[["EXC"]], "T_short")

  # Ensembl GTF values carry a suffix; the leading level is used.
  ensembl <- pick(reference(
    transcript("T_a", "2 (assigned to previous version 3)", long),
    transcript("T_b", "1 (assigned to previous version 5)", long),
    transcript("T_short", "1", short)
  ))
  expect_identical(ensembl[["INC"]], "T_b")

  # Without a TSL column (e.g. RefSeq or custom GTFs), matching still runs.
  untiered <- reference(transcript("T_a", "2", long), transcript("T_b", "1", long),
                        transcript("T_short", "1", short))
  untiered[, transcript_support_level := NULL]
  none <- pick(untiered)
  expect_identical(none[["INC"]], "T_a")
  expect_identical(none[["EXC"]], "T_short")
  expect_false("transcript_support_level" %in% names(untiered))
})

test_that("TSL values are ranked 1-5 with unknown support last", {
  rank <- SpliceImpactR:::.si_tsl_rank
  expect_identical(
    rank(c("1", "1 (assigned to previous version 5)", " 5", "NA", NA, "6", "tsl2", "")),
    c(1L, 1L, 5L, 6L, 6L, 6L, 6L, 6L)
  )
  expect_identical(rank(c(4L, NA)), c(4L, 6L))
  expect_identical(rank(factor(c("3", "2"))), c(3L, 2L))
  expect_identical(rank(character()), integer())
})

test_that("HIT background and overview accept normalized sample directories", {
  root <- tempfile("hit-paths-")
  sample_dir <- file.path(root, "sample_1")
  dir.create(sample_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  exon <- data.table::data.table(
    gene = "GA", exon = "chr1:101-200", ID = "internal",
    nFE = 0L, nLE = 0L, nUP = 20L, nDOWN = 20L, HITindex = 0.25
  )
  data.table::fwrite(exon, file.path(sample_dir, "sample_1.exon"), sep = "\t")

  for (path in c(normalizePath(sample_dir), paste0(sample_dir, "/"))) {
    samples <- data.frame(path = path, sample_name = "alias", condition = "case")
    background <- SpliceImpactR:::read_background(samples)
    expect_identical(background$gene, "GA")
    expect_identical(background$exon, "chr1:101-200")
    overview <- SpliceImpactR:::.getHITindex(samples)
    expect_identical(overview$HITindex, 0.25)
    expect_identical(overview$sample, "alias")
    expect_identical(overview$condition, "case")
  }
})

test_that("manual features at identical coordinates keep their identity", {
  ann <- data.table::data.table(
    type = c("transcript", "exon"), transcript_id = "T1", chr = "chr1", strand = "+",
    cds_has = c(NA, TRUE), cds_rel_start = c(NA_integer_, 1L),
    cds_rel_stop = c(NA_integer_, 300L), cds_gen_start = c(NA_integer_, 1001L),
    cds_gen_stop = c(NA_integer_, 1300L)
  )
  manual <- data.table::data.table(
    ensembl_transcript_id = "T1",
    name = c("Signal peptide", "Transmembrane helix", "Motif A", "Motif B", "Motif B"),
    database = c("signalp", "tmhmm", "user", "user", "user"),
    feature_id = NA_character_, start = 1L, stop = 20L
  )
  mapped <- get_manual_features(manual, ann)
  # Only the exact duplicate collapses; features from other sources or with
  # other labels at the same span are retained.
  expect_equal(nrow(mapped), 4L)
  expect_setequal(sub(";.*$", "", mapped$name),
                  c("Signal peptide", "Transmembrane helix", "Motif A", "Motif B"))
  expect_setequal(mapped$database, c("signalp", "tmhmm", "user"))
})

test_that("manual features given only a peptide ID are placed; unplaced rows are reported", {
  ann <- data.table::data.table(
    type = c("transcript", "exon", "transcript", "exon"),
    transcript_id = c("T1", "T1", "T2", "T2"), protein_id = c(NA, "P1", NA, NA),
    chr = "chr1", strand = "+", cds_has = c(NA, TRUE, NA, FALSE),
    cds_rel_start = c(NA, 1L, NA, NA), cds_rel_stop = c(NA, 300L, NA, NA),
    cds_gen_start = c(NA, 1001L, NA, NA), cds_gen_stop = c(NA, 1300L, NA, NA)
  )
  feature <- function(...) data.table::data.table(..., name = "Motif", start = 1L, stop = 20L)
  by_transcript <- suppressMessages(get_manual_features(feature(ensembl_transcript_id = "T1"), ann))
  expect_message(by_peptide <- get_manual_features(feature(ensembl_peptide_id = "P1"), ann),
                 "^Placed 1 of 1 manual feature row\\(s\\)\\.\n$")
  expect_identical(by_peptide$ensembl_peptide_id, "P1")
  # A peptide ID gives the same placement as its transcript ID.
  expect_identical(by_peptide[, !"ensembl_peptide_id"], by_transcript[, !"ensembl_peptide_id"])

  rows <- feature(ensembl_transcript_id = c("T1", NA, "T9", "T2"),
                  ensembl_peptide_id = c(NA, "P9", NA, NA))
  expect_message(mapped <- get_manual_features(rows, ann), paste0(
    "Placed 1 of 4 manual feature row\\(s\\)\\. Not in the annotation: 2 \\(P9, T9\\)\\. ",
    "No annotated CDS: 1 \\(T2\\)\\."))
  expect_identical(mapped$ensembl_transcript_id, "T1")
})

test_that("manual features follow CDS direction within and across coding exons", {
  for (strand in c("+", "-")) {
    ann <- data.table::data.table(
      type = c("transcript", "exon", "exon"), transcript_id = "T1",
      chr = "chr1", strand = strand, cds_has = c(NA, TRUE, TRUE),
      cds_rel_start = c(NA_integer_, 1L, 151L),
      cds_rel_stop = c(NA_integer_, 150L, 300L),
      cds_gen_start = if (strand == "+") c(NA_integer_, 1001L, 1301L) else c(NA_integer_, 1151L, 801L),
      cds_gen_stop = if (strand == "+") c(NA_integer_, 1150L, 1450L) else c(NA_integer_, 1300L, 950L)
    )
    manual <- data.table::data.table(
      ensembl_transcript_id = "T1", name = c("first", "spanning", "last"),
      feature_id = c("F1", "F2", "F3"),
      start = c(1L, 46L, 91L), stop = c(10L, 55L, 100L)
    )
    mapped <- get_manual_features(manual, ann)
    mapped <- mapped[match(c("F1", "F2", "F3"), feature_id)]
    expect_equal(nrow(mapped), 3L)
    # start/stop remain amino-acid coordinates; name stores genomic endpoints
    # in transcript direction, so the negative-strand endpoints descend.
    expect_equal(mapped$start, c(1L, 46L, 91L))
    expect_equal(mapped$stop, c(10L, 55L, 100L))
    expected_spans <- if (strand == "+") {
      c("1001-1030", "1136-1315", "1421-1450")
    } else {
      c("1300-1271", "1165-936", "830-801")
    }
    expect_identical(mapped$name, paste0(c("first", "spanning", "last"), ";chr1:", expected_spans))
  }
})
