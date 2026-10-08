.di_pairing_fixture <- function() {
  data.table::data.table(
    event_id = rep(c("up", "down"), each = 2), gene_id = "GA",
    chr = "chr1", strand = "+", event_type = "SE",
    form = rep(c("INC", "EXC"), 2),
    transcript_id = rep(c("TI", "TE"), 2), exons = rep(c("EI", "EE"), 2),
    protein_id = rep(c("PI", "PE"), 2),
    inc = rep(c("101-200", "301-400"), 2), exc = "",
    delta_psi = c(0.3, -0.3, -0.4, 0.4), p.value = 0.001, padj = 0.01,
    n_samples = 4L, n_control = 2L, n_case = 2L,
    transcript_seq = "ATG", protein_seq = "M"
  )
}

test_that("paired mode shares multi schema and follows differential direction", {
  input <- .di_pairing_fixture()
  before <- data.table::copy(input)
  paired <- get_pairs(input, source = "paired")
  multi <- get_pairs(input, source = "multi")
  expect_identical(paired, multi)
  expect_identical(paired[event_id == "up", transcript_id_case], "TI")
  expect_identical(paired[event_id == "down", transcript_id_case], "TE")
  expect_identical(paired[event_id == "down", form_case], "EXC")
  expect_identical(input, before)

  minimal <- data.table::data.table(event_id = "E1", form = c("INC", "EXC"))
  expect_equal(nrow(get_pairs(minimal, "paired")), 1L)
  positive <- data.table::copy(input)
  positive[, delta_psi := abs(delta_psi)]
  for (mode in c("paired", "multi")) {
    empty <- get_pairs(positive, mode)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), names(paired))
    expect_equal(nrow(get_pairs(input[0], mode)), 0L)
  }
  invalid <- data.table::copy(input)
  invalid$gene_id[1] <- "GB"
  expect_error(get_pairs(invalid, "paired"), "one gene")
  expect_error(get_pairs(invalid, "multi"), "one gene")
  expect_equal(nrow(get_pairs(input[form == "INC"], "paired")), 0L)
})

test_that("generic DI preserves mapped grouping, form and adjusted p-values", {
  input <- .di_pairing_fixture()
  custom <- data.table::copy(input)
  data.table::setnames(custom, c("event_id", "form", "padj"), c("group", "direction", "FDR"))
  custom[, c("n_samples", "n_control", "n_case") := NULL]
  mapping <- list(gene_id = "gene_id", chr = "chr", strand = "strand", inc = "inc", exc = "exc",
                  delta_psi = "delta_psi", pvalue = "p.value",
                  event_id = "group", form = "direction", padj = "FDR")
  parsed <- import_di_table(custom, colmap = mapping)
  expect_identical(parsed$event_id, input$event_id)
  expect_identical(parsed$form, input$form)
  expect_identical(parsed$padj, input$padj)
  expect_identical(parsed$event_type, input$event_type)
  expect_true(all(is.na(parsed$n_samples)))
  expect_true(all(is.na(parsed$n_control)))
  expect_true(all(is.na(parsed$n_case)))
  for (column in c("transcript_id", "exons", "protein_id", "transcript_seq", "protein_seq")) {
    parsed[, (column) := input[[column]]]
  }
  expect_equal(nrow(get_pairs(parsed, "paired")), 2L)
  expect_equal(nrow(get_pairs(parsed, "multi")), 2L)
  native_names <- import_di_table(input)
  expect_identical(native_names$event_id, input$event_id)
  expect_identical(native_names$form, input$form)
  expect_identical(native_names$padj, input$padj)

  mapping$event_id <- "missing"
  expect_error(import_di_table(custom, colmap = mapping), "missing mapped column")
  without_group <- data.table::copy(input[1:2])
  without_group[, c("event_id", "form") := NULL]
  ungrouped <- import_di_table(without_group)
  expect_equal(data.table::uniqueN(ungrouped$event_id), 2L)
  expect_true(all(ungrouped$form == "SITE"))
})

test_that("pairing and zero-significance results support S4", {
  input <- .di_pairing_fixture()
  obj <- as_splice_impact_result(matched = input)
  result <- get_pairs(obj, "paired")
  expect_true(methods::validObject(result))
  expect_equal(nrow(as_dt_from_s4(result, "paired_hits")), 2L)
  positive <- data.table::copy(input)
  positive[, delta_psi := abs(delta_psi)]
  empty <- get_pairs(as_splice_impact_result(matched = positive), "paired")
  expect_true(methods::validObject(empty))
  expect_equal(length(empty@paired_hits), 0L)

  nonsignificant <- data.table::copy(input)
  nonsignificant[, padj := 1]
  for (di in list(nonsignificant, nonsignificant[0])) {
    before <- data.table::copy(di)
    out <- get_splicing_impact(res = di, verbose = FALSE, debug_steps = TRUE)
    expect_identical(di, before)
    expect_null(out$data)
    for (part in c("res", "matched", "hits_sequences", "pairs", "seq_compare", "hits_domain", "hits_final")) {
      expect_equal(nrow(out[[part]]), 0L)
    }
    expect_true(all(c("event_id", "chr", "inc_case", "inc_control") %in% names(out$hits_final)))
    out_s4 <- get_splicing_impact(res = di, return_class = "S4", verbose = FALSE)
    expect_true(methods::validObject(out_s4))
    expect_equal(length(out_s4@paired_hits), 0L)
    expect_false(out_s4@metadata$has_hits)
    again <- get_pairs(out_s4, "paired")
    expect_true(methods::validObject(again))
    expect_equal(length(again@paired_hits), 0L)
    supplied_s4 <- get_splicing_impact(
      data = as_splice_impact_result(res = di), return_class = "S4", verbose = FALSE
    )
    expect_true(methods::validObject(supplied_s4))
  }
})
