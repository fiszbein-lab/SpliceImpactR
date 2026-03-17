test_that("compare_transcript_pairs matches `matched` schema and cleans invalid IDs", {
  ann <- get_annotation(load = "test")
  ann_dt <- data.table::as.data.table(ann$annotations)
  tx <- unique(
    ann_dt[type == "exon" & !is.na(transcript_id) & nzchar(transcript_id),
           transcript_id]
  )
  expect_true(length(tx) >= 2)

  expected_cols <- c(
    "event_id", "event_type", "form", "gene_id", "chr", "strand",
    "inc", "exc", "delta_psi", "p.value", "padj",
    "n_samples", "n_control", "n_case", "transcript_id", "exons"
  )

  pairs <- data.frame(
    transcript1 = c(tx[1], tx[2], "", "NOT_A_TX"),
    transcript2 = c(tx[2], tx[1], tx[1], tx[2]),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(compare_transcript_pairs(pairs, ann$annotations))

  expect_s3_class(out, "data.table")
  expect_identical(names(out), expected_cols)
  expect_setequal(unique(out$form), c("INC", "EXC"))
  expect_true(all(out[, .N, by = event_id]$N == 2L))
  expect_false(any(is.na(out$gene_id) | out$gene_id == ""))
  expect_false(any(is.na(out$transcript_id) | out$transcript_id == ""))

  out_empty <- suppressMessages(compare_transcript_pairs(
    data.frame(
      transcript1 = c("", "NOPE"),
      transcript2 = c("", "MISSING"),
      stringsAsFactors = FALSE
    ),
    ann$annotations
  ))
  expect_s3_class(out_empty, "data.table")
  expect_identical(names(out_empty), expected_cols)
  expect_equal(nrow(out_empty), 0L)
})
