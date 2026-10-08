test_that("feature caches respond to coordinate and sequence changes", {
  gtf <- data.table::data.table(transcript_id = "T", start = 10L, end = 20L)
  changed <- data.table::copy(gtf)[, start := 11L]
  fingerprint <- SpliceImpactR:::.si_pf_fingerprint_gtf
  expect_identical(fingerprint(gtf), fingerprint(data.table::copy(gtf)))
  expect_false(identical(fingerprint(gtf), fingerprint(changed)))
  seq <- data.table::data.table(transcript_id = "T", protein_seq = "MAA")
  changed <- data.table::copy(seq)[, protein_seq := "MAV"]
  fingerprint <- SpliceImpactR:::.si_pf_fingerprint_sequences
  expect_false(identical(fingerprint(seq), fingerprint(changed)))
  expect_false(identical(fingerprint(list(sequences = seq)),
                         fingerprint(list(sequences = changed))))
})

test_that("BioMart batches bound transcript IDs rather than whole chromosomes", {
  ann <- data.table::data.table(type = "transcript", transcript_type = "protein_coding",
                               transcript_id = paste0("T", 1:7, ".1"), chr = "chr1")
  groups <- SpliceImpactR:::split_into_bits(ann, 3L)
  expect_equal(lengths(groups), c(`1` = 3L, `2` = 3L, `3` = 1L))
  expect_setequal(unlist(groups, use.names = FALSE), paste0("T", 1:7))
  expect_error(SpliceImpactR:::split_into_bits(ann, 0L), "positive integer")
})

test_that("compressed HIT companion files and requested hybrid types work", {
  root <- tempfile("compressed-hit-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ann <- data.table::data.table(gene = "G", exon = "chr1:101-200", ID = "first",
                               nUP = 10, nDOWN = 30, nFE = 1, nLE = 0)
  data.table::fwrite(ann, file.path(root, "sample.exon.gz"), sep = "\t")
  psi <- data.table::copy(ann)[, `:=`(HFEPSI = 0.5, `sumR-L` = 40)]
  data.table::fwrite(psi, file.path(root, "sample.HFEPSI.gz"), sep = "\t")
  expect_identical(SpliceImpactR:::.read_exon_files(file.path(root, "sample."))$gene, "G")
  out <- get_rmats_hit(data.frame(path = root, sample_name = "s", condition = "case"),
                       event_types = "HFE", keep_annotated_first_last = FALSE)
  expect_identical(out$event_type, "HFE")
  expect_identical(out$inc, "101-200")
  expect_error(get_rmats_hit(data.frame(), event_types = "unknown"), "supported")
})

test_that("plain and gzipped copies of one input are read once", {
  root <- tempfile("duplicate-inputs-")
  plain <- file.path(root, "plain")
  both <- file.path(root, "both")
  dir.create(plain, recursive = TRUE)
  dir.create(both)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ann <- data.table::data.table(gene = "G", exon = c("chr1:101-200", "chr1:301-400"),
                               ID = "first", nUP = 10, nDOWN = 30, nFE = 1, nLE = 0)
  psi <- data.table::copy(ann)[, `:=`(AFEPSI = c(0.25, 0.75), `sumR-L` = 40)]
  rmats <- system.file("extdata", "rawData", "case_S1", "SE.MATS.JCEC.txt",
                       package = "SpliceImpactR")
  for (dir in c(plain, both)) {
    data.table::fwrite(ann, file.path(dir, "sample.exon"), sep = "\t")
    data.table::fwrite(psi, file.path(dir, "sample.AFEPSI"), sep = "\t")
    file.copy(rmats, file.path(dir, "SE.MATS.JCEC.txt"))
  }
  data.table::fwrite(psi, file.path(both, "sample.AFEPSI.gz"), sep = "\t")
  R.utils::gzip(rmats, destname = file.path(both, "SE.MATS.JCEC.txt.gz"), remove = FALSE)

  expect_identical(basename(SpliceImpactR:::.find_hitindex_files(both)), "sample.AFEPSI")
  samples <- function(dir) data.frame(path = dir, sample_name = "s", condition = "case")
  hit_plain <- get_rmats_hit(samples(plain), event_types = "AFE")
  hit_both <- get_rmats_hit(samples(both), event_types = "AFE")
  expect_equal(hit_both$psi, c(0.25, 0.75))
  columns <- setdiff(names(hit_plain), "source_file")
  expect_equal(hit_both[, columns, with = FALSE], hit_plain[, columns, with = FALSE])
  expect_equal(nrow(load_rmats(samples(both), use = "JCEC", event_types = "SE")),
               nrow(load_rmats(samples(plain), use = "JCEC", event_types = "SE")))
})

test_that("size factors use explicit metadata and handle normalized paths", {
  sf <- data.frame(sample_name = c("a", "b"), sizeFactor = c(0.5, 2))
  expect_equal(SpliceImpactR:::.get_size_factors_from_exons(sf, "user-given"), sf)
  expect_equal(SpliceImpactR:::.get_size_factors_from_exons(sf, "user_given"), sf)
  sf$sizeFactor[1L] <- 0
  expect_error(SpliceImpactR:::.get_size_factors_from_exons(sf, "user-given"), "positive")
})

test_that("domain enrichment counts unique pairs and respects supplied columns", {
  fg <- data.table::data.table(event_id = c("E1", "E2"), transcript_id_case = "a",
                              transcript_id_control = "b", custom = list("db;D", "db;D"))
  bg <- data.table::data.table(transcript_id_1 = c("a", "c", "e", "g"),
                              transcript_id_2 = c("b", "d", "f", "h"),
                              custom_bg = list("db;D", "db;D", "db;X", "db;X"))
  original <- data.table::copy(fg)
  out <- enrich_domains_hypergeo(fg, bg, domain_col_fg = "custom",
                                 domain_col_bg = "custom_bg", min_fg_count = 1)
  expect_identical(out$domain_id, "D")
  expect_equal(out[, .(k, K, M, B)], data.table::data.table(k = 1, K = 1, M = 2, B = 4))
  expect_equal(out$pval, 0.5)
  expect_equal(out$OR, 5)
  expect_identical(fg, original)
  empty <- enrich_domains_hypergeo(fg[0], bg, "custom", "custom_bg")
  expect_equal(nrow(empty), 0L)
  expect_identical(names(empty), names(out))
  fg[, transcript_id_case := "absent"]
  expect_error(enrich_domains_hypergeo(fg, bg, "custom", "custom_bg"), "absent")
})

test_that("domain enrichment counts only pairs with a tested domain change", {
  fg <- data.table::data.table(event_id = "E1", transcript_id_case = "a",
                              transcript_id_control = "b",
                              custom = list(c("db;D", "other;Y")))
  bg <- data.table::data.table(
    transcript_id_1 = c("a", "c", "e", "g", "i"),
    transcript_id_2 = c("b", "d", "f", "h", "j"),
    # i|j has no remaining change, e.g. domains differing only by overlap.
    custom_bg = list(c("db;D", "other;Y"), "db;D", "db;X", "other;Y", character())
  )
  all_db <- enrich_domains_hypergeo(fg, bg, "custom", "custom_bg", min_fg_count = 1)
  expect_setequal(all_db$domain_id, c("D", "Y"))
  expect_true(all(all_db$K == 1L & all_db$B == 4L & all_db$M == 2L))
  expect_equal(all_db$pval, rep(0.5, 2))

  # With db_filter, pairs without a change in that database leave both populations.
  one_db <- enrich_domains_hypergeo(fg, bg, "custom", "custom_bg", min_fg_count = 1,
                                    db_filter = "db")
  expect_equal(one_db[, .(domain_id, k, K, M, B)],
               data.table::data.table(domain_id = "D", k = 1L, K = 1L, M = 2L, B = 3L))
  expect_equal(one_db$pval, stats::phyper(0, 2, 1, 1, lower.tail = FALSE))
})

test_that("identical sequence shortcuts return the actual matrix score", {
  matrix <- pwalign::nucleotideSubstitutionMatrix(match = 3, mismatch = -1, baseOnly = TRUE)
  expect_equal(SpliceImpactR:::.align_dna("ACGT", "ACGT", matrix)$score, 12)
  env <- new.env()
  utils::data("BLOSUM62", package = "Biostrings", envir = env)
  expect_equal(SpliceImpactR:::.align_aa("MW", "MW", env$BLOSUM62)$score,
               as.numeric(pwalign::pairwiseAlignment(Biostrings::AAString("MW"),
                 Biostrings::AAString("MW"), substitutionMatrix = env$BLOSUM62, scoreOnly = TRUE)))
})

test_that("legacy matching retains unmatched forms without creating absent isoforms", {
  x <- rank_test_reference()
  x$events[form == "INC", inc := "5001-6000"]
  for (chunk in c(1L, 20L)) {
    matched <- get_matched_events_chunked(x$events, x$annotations, chunk_size = chunk)
    expect_equal(nrow(matched), 2L)
    expect_true(is.na(matched[form == "INC", transcript_id]))
    pairs <- get_pairs(attach_sequences(matched, x$sequences), source = "multi")
    expect_equal(nrow(pairs), 0L)
  }
})

test_that("the feature gene universe holds genes with a feature-changing pair", {
  bg <- data.table::data.table(
    gene_id = c("unchanged", "unannotated", "changed", "overlap_only"),
    transcript_id_1 = c("A", "C", "E", "G"), transcript_id_2 = c("B", "D", "F", "H"))
  pf <- data.table::data.table(ensembl_transcript_id = c("A", "B", "E", "G", "H"),
    ensembl_peptide_id = paste0("P", c("A", "B", "E", "G", "H")), database = "db",
    clean_name = "D",
    name = c(rep("D;chr1:101-200", 4L), "D;chr1:101-190"))
  result <- SpliceImpactR:::get_domain_background(bg, pf, BPPARAM = BiocParallel::SerialParam())
  # Identical or overlap-equivalent features can never yield a domain foreground.
  expect_identical(attr(result, "feature_gene_universe"), "changed")
  expect_identical(attr(result, "feature_gene_universe"),
                   sort(unique(result[lengths(total_sd_domains) > 0L, gene_id])))
})

test_that("alignment and length plots export the plot they return", {
  hits <- data.table::data.table(prot_pid = c(50, 70, 90), event_type = "SE",
    summary_classification = "protein_coding", prot_len_case = c(110, 160, 220),
    prot_len_control = c(100, 130, 170), prot_len_diff = c(10, 30, 50))
  root <- tempfile("plot-export-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  alignment <- file.path(root, "alignment.pdf")
  lengths <- file.path(root, "lengths.pdf")
  expect_s3_class(plot_alignment_summary(hits, output_file = alignment), "ggplot")
  expect_s3_class(plot_length_comparison(hits, output_file = lengths), "ggplot")
  expect_true(all(file.info(c(alignment, lengths))$size > 1000L))
})
