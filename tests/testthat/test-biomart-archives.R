test_that("BioMart connections use the release's e<release>.ensembl.org host", {
  calls <- list()
  local_mocked_bindings(useEnsembl = function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    "mart"
  }, .package = "biomaRt")
  connect <- SpliceImpactR:::.si_use_ensembl_mart
  expect_identical(connect("hsapiens_gene_ensembl", version = 111), "mart")
  expect_identical(calls[[1]]$host, "https://e111.ensembl.org")
  expect_null(calls[[1]]$version)
  expect_null(calls[[1]]$mirror)
  # An explicit host is used as given.
  connect("hsapiens_gene_ensembl", version = 111, host = "https://may2025.archive.ensembl.org")
  expect_identical(calls[[2]]$host, "https://may2025.archive.ensembl.org")
  # Releases without BioMart stop before connecting.
  expect_error(connect("hsapiens_gene_ensembl", version = 104), "105 to 116")
  expect_error(connect("hsapiens_gene_ensembl", version = 117), "105 to 116")
  expect_length(calls, 2L)
  # Connection failures name the host.
  local_mocked_bindings(useEnsembl = function(...) stop("Timeout was reached"), .package = "biomaRt")
  expect_error(connect("hsapiens_gene_ensembl", version = 111),
               "e111\\.ensembl\\.org.*Timeout was reached.*temporarily unavailable")
})

test_that("biomaRt before 2.70 reaches every release, including 116", {
  hosts <- character()
  local_mocked_bindings(useEnsembl = function(biomart, dataset, host, ...) {
    hosts <<- c(hosts, host)
    "mart"
  }, .package = "biomaRt")
  connect <- SpliceImpactR:::.si_use_ensembl_mart
  local_mocked_bindings(.si_biomart_handles_archive_hosts = function() FALSE)
  # e116 avoids old biomaRt sending the newest archive's queries to www.
  expect_identical(connect("hsapiens_gene_ensembl", version = 116), "mart")
  expect_identical(hosts, "https://e116.ensembl.org")
  # The 116 archive host itself cannot work there; the error names e116.
  expect_error(connect("hsapiens_gene_ensembl", version = 111,
                       host = "https://JUN2026.archive.ensembl.org/"), "e116\\.ensembl\\.org")
  expect_length(hosts, 1L)
  # Other archive hosts fail once old biomaRt cannot get Ensembl's archive list.
  local_mocked_bindings(useEnsembl = function(...) stop("Unable to contact any Ensembl mirror"),
                        .package = "biomaRt")
  archive <- "https://jan2024.archive.ensembl.org"
  expect_error(connect("hsapiens_gene_ensembl", version = 111, host = archive),
               "Unable to contact any Ensembl mirror.*biomaRt before 2.70.*e<release>")
  local_mocked_bindings(.si_biomart_handles_archive_hosts = function() TRUE)
  expect_error(connect("hsapiens_gene_ensembl", version = 111, host = archive),
               "temporarily unavailable")
  local_mocked_bindings(useEnsembl = function(...) "mart", .package = "biomaRt")
  expect_identical(connect("hsapiens_gene_ensembl", version = 111,
                           host = "https://jun2026.archive.ensembl.org"), "mart")
})

test_that("get_protein_features queries Ensembl 111 by default", {
  ann <- get_annotation(load = "test")$annotations
  tx <- ann[type == "exon" & cds_has == TRUE, transcript_id][1]
  hosts <- character()
  local_mocked_bindings(
    useEnsembl = function(biomart, dataset, host, ...) {
      hosts <<- c(hosts, host)
      "mart"
    },
    getBM = function(attributes, ...) {
      data.frame(ensembl_transcript_id = tx, ensembl_peptide_id = "P1",
                 pfam = "PF_TEST", pfam_start = 1L, pfam_end = 10L)[, attributes]
    },
    .package = "biomaRt"
  )
  expect_identical(formals(get_protein_features)$release, 111)
  out <- get_protein_features("pfam", ann, use_cache = FALSE)
  expect_identical(hosts, "https://e111.ensembl.org")
  expect_identical(unique(out$feature_id), "PF_TEST")
  get_protein_features("pfam", ann, use_cache = FALSE,
                       ensembl_host = "https://sep2025.archive.ensembl.org")
  expect_identical(hosts[2], "https://sep2025.archive.ensembl.org")
  # Ensembl retired its mirrors; the argument is ignored with a warning.
  expect_warning(get_protein_features("pfam", ann, use_cache = FALSE, ensembl_mirror = "useast"),
                 "retired")
  expect_identical(hosts[3], "https://e111.ensembl.org")
  expect_error(get_protein_features("pfam", ann, use_cache = FALSE, ensembl_host = NA_character_),
               "ensembl_host")
  # The ELM route maps UniProt IDs through the same archive.
  local_mocked_bindings(useEnsembl = function(...) stop("offline"), .package = "biomaRt")
  expect_error(SpliceImpactR:::get_linear_motifs(ann, NULL, release = 111,
                                                 ensembl_host = "https://may2024.archive.ensembl.org"),
               "may2024\\.archive\\.ensembl\\.org.*offline")
})

test_that("the protein-feature cache key records an explicit host", {
  ann <- get_annotation(load = "test")$annotations
  key <- function(host = NULL) SpliceImpactR:::.si_pf_cache_key(
    "pfam", ann, NULL, "hsapiens_gene_ensembl", 111L, FALSE, ensembl_host = host)
  expect_false(grepl("host-", key()))
  expect_false(identical(key(), key("https://sep2025.archive.ensembl.org")))
  expect_identical(key("https://sep2025.archive.ensembl.org"), key("https://sep2025.archive.ensembl.org/"))
})
