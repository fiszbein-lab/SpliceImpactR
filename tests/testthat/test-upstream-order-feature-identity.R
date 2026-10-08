test_that("upstream shifts follow transcript order on both strands", {
  for (strand in c("+", "-")) {
    starts <- if (strand == "+") c(100L, 300L, 500L) else c(500L, 300L, 100L)
    ann <- data.table::data.table(
      type = "exon", transcript_id = rep(c("T1", "T2"), each = 3),
      exon_id = rep(c("Z_UPSTREAM", "A_EVENT", "M_DOWNSTREAM"), 2),
      strand = strand, cds_gen_start = rep(starts, 2),
      cds_gen_stop = rep(starts + 98L, 2),
      start_frame = c(0L, 0L, 0L, 1L, 0L, 0L),
      stop_frame = c(2L, 2L, 2L, 0L, 2L, 2L),
      absolute_exon_position = rep(1:3, 2),
      coding_exon_position = rep(1:3, 2),
      coding_exon_class = rep(c("first", "internal", "last"), 2)
    )
    index <- SpliceImpactR:::build_coding_index(ann)
    before <- data.table::copy(index)
    expect_true(SpliceImpactR:::.incoming_upstream_shift(
      index, "T1", "T2", "A_EVENT", "A_EVENT"
    ))
    expect_false(SpliceImpactR:::.incoming_upstream_shift(
      index, "T1", "T2", "Z_UPSTREAM", "Z_UPSTREAM"
    ))
    expect_identical(index, before)

    # A downstream frame difference must not be reported as an upstream shift.
    ann[transcript_id == "T2" & exon_id == "Z_UPSTREAM",
        `:=`(start_frame = 0L, stop_frame = 2L)]
    ann[transcript_id == "T2" & exon_id == "M_DOWNSTREAM",
        `:=`(start_frame = 1L, stop_frame = 0L)]
    index <- SpliceImpactR:::build_coding_index(ann)
    expect_false(SpliceImpactR:::.incoming_upstream_shift(
      index, "T1", "T2", "A_EVENT", "A_EVENT"
    ))
  }
})

test_that("co-located distinct features survive while exact duplicates collapse", {
  ann <- data.table::data.table(
    type = c("transcript", "exon"), transcript_id = "T1",
    chr = "chr1", strand = "+", cds_has = c(NA, TRUE),
    cds_rel_start = c(NA_integer_, 1L), cds_rel_stop = c(NA_integer_, 300L),
    cds_gen_start = c(NA_integer_, 1001L), cds_gen_stop = c(NA_integer_, 1300L)
  )
  biomart <- data.table::data.table(
    ensembl_transcript_id = "T1", ensembl_peptide_id = "P1",
    interpro = "IPR_TEST", interpro_description = "Test domain",
    interpro_short_description = "Test", interpro_start = 1L, interpro_end = 10L,
    pfam = c("PF_A", "PF_B", "PF_A"), pfam_start = 1L, pfam_end = 10L
  )

  # Substitute only the external service boundary; run real conversion/mapping.
  get_features <- SpliceImpactR::get_protein_features
  fixture_env <- new.env(parent = environment(get_features))
  fixture_env$get_biomart_protein_features <- function(...) data.table::copy(biomart)
  environment(get_features) <- fixture_env
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (combine in c(FALSE, TRUE)) {
    result <- get_features(c("interpro", "pfam"), ann,
                           use_cache = FALSE, combine_overlaps = combine)
    expect_equal(nrow(result), 3L)
    expect_setequal(paste(result$database, result$feature_id),
                    c("interpro IPR_TEST", "pfam PF_A", "pfam PF_B"))
    expect_equal(result$start, rep(1L, 3))
    expect_equal(result$stop, rep(10L, 3))
    expect_true(all(grepl(";chr1:1001-1030$", result$name)))
    expect_identical(names(result), c(
      "ensembl_transcript_id", "start", "stop", "chr", "strand",
      "feature_id", "clean_name", "alt_name", "database",
      "ensembl_peptide_id", "method", "name"
    ))
  }

  # Database identity must also be preserved when feature labels coincide.
  biomart[, `:=`(interpro = "shared_label", pfam = "shared_label")]
  result <- get_features(c("interpro", "pfam"), ann, use_cache = FALSE)
  expect_equal(nrow(result), 2L)
  expect_setequal(result$database, c("interpro", "pfam"))

  # A pre-fix cache entry must not hide the new results. Emulate cache storage
  # in memory so this regression does not need a network or a user cache.
  legacy_key <- paste0(
    "protein_features/v", utils::packageVersion("SpliceImpactR"),
    "/species-hsapiens_gene_ensembl/release-109/db-interpro,pfam",
    "/combine-FALSE/gtf-", SpliceImpactR:::.si_pf_fingerprint_gtf(ann),
    "/seq-no_elm.rds"
  )
  cache <- new.env(parent = emptyenv())
  assign(legacy_key, data.table::copy(result[1]), envir = cache)
  fixture_env$.si_bfc <- function(...) cache
  fixture_env$.si_bfc_get_rds <- function(bfc, key) get0(key, envir = bfc, inherits = FALSE)
  fixture_env$.si_bfc_put_rds <- function(bfc, key, value) assign(key, value, envir = bfc)
  requests <- 0L
  fixture_env$get_biomart_protein_features <- function(...) {
    requests <<- requests + 1L
    data.table::copy(biomart)
  }
  refreshed <- get_features(c("interpro", "pfam"), ann)
  cached <- get_features(c("interpro", "pfam"), ann)
  expect_equal(nrow(refreshed), 2L)
  expect_identical(cached, refreshed)
  expect_identical(requests, 1L)
  expect_equal(nrow(get(legacy_key, envir = cache)), 1L)
})
