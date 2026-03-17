test_that("gene enrichment selectors work for data.table and S4 inputs", {
  di <- data.table::data.table(
    event_id = c("e1", "e2", "e3"),
    form = c("INC", "INC", "INC"),
    gene_id = c("g1", "g2", ""),
    chr = c("chr1", "chr1", "chr1"),
    strand = c("+", "+", "+"),
    inc = c("10-20", "20-30", "30-40"),
    exc = c("", "", ""),
    padj = c(0.01, 0.2, 0.01),
    delta_psi = c(0.2, 0.3, 0.05)
  )

  hits_final_like <- data.table::data.table(
    event_id = c("e1", "e2", "e3", "e4"),
    gene_id = c("g1", "g2", NA_character_, ""),
    chr = c("chr1", "chr1", "chr1", "chr1"),
    strand = c("+", "+", "+", "+"),
    inc_case = c("10-20", "20-30", "30-40", "40-50"),
    inc_control = c("11-21", "21-31", "31-41", "41-51"),
    exc_case = c("", "", "", ""),
    exc_control = c("", "", "", ""),
    diff_n = c(1, 0, 2, 1),
    n_ppi = c(0, 3, 1, 1)
  )

  # data.table inputs
  expect_setequal(get_di_gene_enrichment(di, 0.05, 0.1), "g1")
  expect_setequal(get_domain_gene_for_enrichment(hits_final_like), "g1")
  expect_setequal(get_ppi_gene_enrichment(hits_final_like), "g2")

  # S4 input
  obj <- as_splice_impact_result(
    res = di,
    hits_final = hits_final_like
  )

  expect_setequal(get_di_gene_enrichment(obj, 0.05, 0.1), "g1")
  expect_setequal(get_domain_gene_for_enrichment(obj), "g1")
  expect_setequal(get_ppi_gene_enrichment(obj), "g2")

  # unified wrapper: x input
  expect_setequal(get_gene_enrichment("di", x = di, padj_threshold = 0.05, delta_psi_threshold = 0.1), "g1")
  expect_setequal(get_gene_enrichment("domain", x = hits_final_like), "g1")
  expect_setequal(get_gene_enrichment("ppi", x = hits_final_like), "g2")

  # unified wrapper: explicit res/hits
  expect_setequal(get_gene_enrichment("di", res = di, padj_threshold = 0.05, delta_psi_threshold = 0.1), "g1")
  expect_setequal(get_gene_enrichment("domain", hits = hits_final_like), "g1")
  expect_setequal(get_gene_enrichment("ppi", hits = hits_final_like), "g2")

  # unified wrapper: S4
  expect_setequal(get_gene_enrichment("di", x = obj, padj_threshold = 0.05, delta_psi_threshold = 0.1), "g1")
  expect_setequal(get_gene_enrichment("domain", x = obj), "g1")
  expect_setequal(get_gene_enrichment("ppi", x = obj), "g2")
})
