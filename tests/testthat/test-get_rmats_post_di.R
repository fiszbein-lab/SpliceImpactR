test_that("get_rmats_post_di works on a single rMATS row", {
  df <- data.table(
    ID = 1L,
    GeneID = "ENSG00000182871",
    geneSymbol = "COL18A1",
    chr = "chr21",
    strand = "+",
    longExonStart_0base = 45505834L,
    longExonEnd = 45505966L,
    shortES = 45505837L,
    shortEE = 45505966L,
    flankingES = 45505357L,
    flankingEE = 45505431L,
    ID.2 = 2L,
    IJC_SAMPLE_1 = "4,1,0",
    SJC_SAMPLE_1 = "9,12,3",
    IJC_SAMPLE_2 = "0,4,5",
    SJC_SAMPLE_2 = "11,15,15",
    IncFormLen = 52L,
    SkipFormLen = 49L,
    PValue = 0.6967562,
    FDR = 1,
    IncLevel1 = "0.295,0.073,0.0",
    IncLevel2 = "0.0,0.201,0.239",
    IncLevelDifference = -0.024,
    stringsAsFactors = FALSE
  )

  # run
  res <- get_rmats_post_di(df, event_type = "A3SS")

  # structure checks
  expect_s3_class(res, "data.table")
  expect_true(all(c("event_id","event_type","form","gene_id","chr","strand",
                    "inc","exc","delta_psi","p.value","padj") %in% names(res)))

  # should return INC + EXC = 2 rows
  expect_equal(nrow(res), 2)

  # event_id stability
  expect_equal(length(unique(res$event_id)), 1)

  # forms
  expect_setequal(unique(res$form), c("INC","EXC"))

  # delta psi direction (A3SS logic)
  expect_true(res$form[res$delta_psi > 0] == "EXC")
  expect_true(res$form[res$delta_psi < 0] == "INC")

  # numeric checks
  expect_equal(res$p.value[1], 0.6967562, tolerance = 1e-6)
  expect_equal(res$padj[1], 1, tolerance = 1e-6)

  # coordinate sanity (not NA and string format)
  expect_false(any(is.na(res$inc)))
  expect_true(is.character(res$inc))
})
