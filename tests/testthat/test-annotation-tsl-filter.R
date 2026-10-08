.annotation_tsl_reference <- function(root, include_tsl = TRUE) {
  dir.create(root, recursive = TRUE)
  levels <- c("1", "2", "3", "4", "5", NA_character_, "NA", "unknown")
  ids <- paste0("ENST_TSL_", seq_along(levels))
  ids <- c(ids, "ENST_TSL_INCOMPLETE")
  levels <- c(levels, "1")
  lines <- paste(
    "chr1", "test", "gene", "1", "2000", ".", "+", ".",
    'gene_id "ENSG_TSL"; gene_type "protein_coding"; gene_name "TSL"; tag "basic";',
    sep = "\t"
  )
  for (i in seq_along(ids)) {
    attrs <- paste0(
      'gene_id "ENSG_TSL"; gene_type "protein_coding"; gene_name "TSL"; ',
      'transcript_id "', ids[[i]], '"; transcript_type "protein_coding"; ',
      'transcript_name "', ids[[i]], '"; protein_id "ENSP_TSL_', i, '"; ',
      'exon_id "ENSE_TSL_', i, '"; exon_number "1"; ',
      'tag "', if (i == length(ids)) "cds_start_NF" else "basic", '";'
    )
    if (isTRUE(include_tsl) && !is.na(levels[[i]])) {
      attrs <- paste0(attrs, ' transcript_support_level "', levels[[i]], '";')
    }
    start <- 100L * i
    for (type in c("transcript", "exon", "CDS")) {
      lines <- c(lines, paste(
        "chr1", "test", type, start, start + 8L, ".", "+",
        if (type == "CDS") "0" else ".", attrs, sep = "\t"
      ))
    }
  }
  gtf <- file.path(root, "annotation.gtf")
  transcripts <- file.path(root, "transcripts.fa")
  proteins <- file.path(root, "proteins.fa")
  writeLines(lines, gtf)
  writeLines(as.vector(rbind(paste0(">", ids, "|CDS:1-9"), "ATGGCTGCT")), transcripts)
  writeLines(as.vector(rbind(paste0(">ENSP_TSL_", seq_along(ids), "|", ids), "MAA")), proteins)
  list(
    base_dir = file.path(root, "cache"), gtf_path = gtf,
    transcript_path = transcripts, translation_path = proteins
  )
}

test_that("NULL TSL filtering retains all support values and isolates processed caches", {
  root <- tempfile("annotation-tsl-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- .annotation_tsl_reference(root)
  tx_ids <- function(x) sort(x$annotations[type == "transcript", transcript_id])

  filtered <- do.call(get_annotation, c(list(load = "path"), files))
  expect_identical(tx_ids(filtered), paste0("ENST_TSL_", 1:3))
  expect_equal(nrow(filtered$annotations[type == "gene"]), 1L)
  expect_error(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = NULL), "tsl-none")

  all_tsl <- do.call(get_annotation, c(list(load = "path", filter_tsl = NULL), files))
  expect_identical(tx_ids(all_tsl), paste0("ENST_TSL_", 1:8))
  expect_equal(nrow(all_tsl$annotations[type == "exon"]), 8L)
  expect_equal(nrow(all_tsl$annotations[type == "gene"]), 1L)
  expect_true(is.na(all_tsl$annotations[
    type == "transcript" & transcript_id == "ENST_TSL_6", transcript_support_level
  ]))
  expect_identical(all_tsl$annotations[
    type == "transcript" & transcript_id == "ENST_TSL_8", transcript_support_level
  ], "unknown")
  expect_true(all(!is.na(all_tsl$sequences[
    transcript_id %in% paste0("ENST_TSL_", 1:8), transcript_seq
  ])))

  cached_all <- get_annotation(load = "cached", base_dir = files$base_dir, filter_tsl = NULL)
  cached_default <- get_annotation(load = "cached", base_dir = files$base_dir)
  expect_equal(cached_all, all_tsl)
  expect_equal(cached_default, filtered)
  expect_equal(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = character()), filtered)

  numeric_tsl <- do.call(get_annotation, c(list(load = "path", filter_tsl = 1:5), files))
  expect_identical(tx_ids(numeric_tsl), paste0("ENST_TSL_", 1:5))
  expect_equal(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = 1:5), numeric_tsl)
  expect_equal(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = NULL), all_tsl)

  cache_paths <- BiocFileCache::bfcinfo(SpliceImpactR:::.si_bfc(files$base_dir))$rpath
  expect_true(all(startsWith(cache_paths, paste0(normalizePath(files$base_dir), "/"))))
  expect_true(all(file.exists(cache_paths)))
  load_elsewhere <- function() {
    original_dir <- getwd()
    on.exit(setwd(original_dir), add = TRUE)
    setwd(root)
    expect_equal(get_annotation(load = "cached", base_dir = files$base_dir,
                                filter_tsl = NULL), all_tsl)
    expect_equal(get_annotation(load = "cached", base_dir = files$base_dir), filtered)
  }
  load_elsewhere()

  # Substitute only asset acquisition; import, preparation and caching are real.
  load_link <- get_annotation
  fixture_env <- new.env(parent = environment(load_link))
  fixture_env$.si_prepare_assets <- function(base_dir, species, release, mode) {
    list(bfc = SpliceImpactR:::.si_bfc(base_dir), paths = list(
      gtf_gz = files$gtf_path, txfa_gz = files$transcript_path,
      aafa_gz = files$translation_path
    ))
  }
  environment(load_link) <- fixture_env
  expect_equal(load_link(load = "link", base_dir = files$base_dir,
                         filter_tsl = NULL), all_tsl)
})

test_that("TSL filtering can be disabled for GTFs without a support-level column", {
  root <- tempfile("annotation-no-tsl-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- .annotation_tsl_reference(root, include_tsl = FALSE)
  out <- do.call(get_annotation, c(list(load = "path", filter_tsl = NULL), files))
  expect_false("transcript_support_level" %in% names(out$annotations))
  expect_setequal(out$annotations[type == "transcript", transcript_id],
                  paste0("ENST_TSL_", 1:8))
  # The unfiltered reference remains usable for transcript matching.
  matched <- get_matched_events_chunked(data.table::data.table(
    event_id = "AFE:1", event_type = "AFE", form = "SITE", gene_id = "ENSG_TSL",
    chr = "chr1", strand = "+", inc = "100-108", exc = ""
  ), out$annotations)
  expect_identical(matched$transcript_id, "ENST_TSL_1")
  expect_equal(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = NULL), out)
  expect_error(get_annotation(load = "cached", base_dir = files$base_dir,
                              filter_tsl = "unknown"), "arg")
})
