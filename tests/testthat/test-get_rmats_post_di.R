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

test_that("post-DI files of one event type receive distinct event IDs", {
  root <- tempfile("post-di-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  event <- function(gene, start) data.table::data.table(
    ID = 1L, GeneID = paste0(gene, ".1"), geneSymbol = gene, chr = "chr1", strand = "+",
    longExonStart_0base = start, longExonEnd = start + 100L, shortES = start + 10L,
    shortEE = start + 100L, flankingES = start - 500L, flankingEE = start - 400L,
    IJC_SAMPLE_1 = "1", SJC_SAMPLE_1 = "1", IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "1",
    IncFormLen = 1L, SkipFormLen = 1L, PValue = 0.001, FDR = 0.01,
    IncLevel1 = "0.9", IncLevel2 = "0.1", IncLevelDifference = 0.8
  )
  files <- file.path(root, c("part1.A3SS.MATS.JC.txt", "part2.A3SS.MATS.JC.txt"))
  data.table::fwrite(event("GA", 1000L), files[1], sep = "\t")
  data.table::fwrite(event("GB", 5000L), files[2], sep = "\t")
  input <- function(path, grp2 = "KO") {
    data.frame(path = path, grp1 = "WT", grp2 = grp2, event_type = "A3SS")
  }

  di <- get_rmats_post_di(input(files))
  expect_setequal(unique(di$event_id), c("A3SS:1", "A3SS:2"))
  expect_true(all(di[, data.table::uniqueN(gene_id), by = event_id]$V1 == 1L))
  pairs <- get_pairs(di, source = "paired")
  expect_setequal(pairs$gene_id, c("GA", "GB"))
  # A file listed twice contributes its events once.
  expect_equal(nrow(get_rmats_post_di(input(files[c(1L, 1L)]))), 2L)
  # One DI table represents one comparison.
  expect_error(get_rmats_post_di(input(files, grp2 = c("KO1", "KO2"))), "one comparison")
})

rmats_case_row <- function() data.table::data.table(
  ID = 1L, GeneID = "GA.1", geneSymbol = "GA", chr = "chr1", strand = "+",
  longExonStart_0base = 1000L, longExonEnd = 1100L, shortES = 1010L, shortEE = 1100L,
  flankingES = 500L, flankingEE = 600L, IJC_SAMPLE_1 = "9", SJC_SAMPLE_1 = "1",
  IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "9", IncFormLen = 1L, SkipFormLen = 1L,
  PValue = 0.001, FDR = 0.01, IncLevel1 = "0.9", IncLevel2 = "0.1", IncLevelDifference = 0.8
)

test_that("case_group orients delta PSI as case minus control", {
  row <- rmats_case_row()
  expect_message(by_b1 <- get_rmats_post_di(row, event_type = "A3SS"),
                 "group 1 \\(--b1\\) is the case")
  expect_equal(by_b1[form == "INC", delta_psi], 0.8)
  expect_message(by_b2 <- get_rmats_post_di(row, event_type = "A3SS", case_group = 2),
                 "group 2 \\(--b2\\) is the case")
  expect_equal(by_b2[form == "INC", delta_psi], -0.8)
  expect_equal(by_b2[form == "EXC", delta_psi], 0.8)
  # Only the sign changes, so the case role moves to the other form.
  expect_equal(by_b2[, !"delta_psi"], by_b1[, !"delta_psi"])
  expect_error(get_rmats_post_di(row, event_type = "A3SS", case_group = 3), "case_group")
  expect_error(get_rmats_post_di(row, event_type = "A3SS", case_group = "KO"), "grp1 and grp2")
})

test_that("a case_group label is matched in each file's grp1/grp2", {
  root <- tempfile("post-di-case-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- file.path(root, c("A3SS.MATS.JC.txt", "A5SS.MATS.JC.txt"))
  for (f in files) data.table::fwrite(rmats_case_row(), f, sep = "\t")
  # These two runs listed the groups in opposite orders.
  input <- data.frame(path = files, grp1 = c("WT", "KO"), grp2 = c("KO", "WT"),
                      event_type = c("A3SS", "A5SS"))
  expect_message(di <- get_rmats_post_di(input, case_group = "KO"), "by file")
  expect_equal(di[event_type == "A3SS" & form == "INC", delta_psi], -0.8)
  expect_equal(di[event_type == "A5SS" & form == "INC", delta_psi], 0.8)
  expect_message(get_rmats_post_di(input[1L, ], case_group = "KO"), 'group 2 \\(--b2, "KO"\\)')
  expect_error(get_rmats_post_di(input, case_group = "HET"), "exactly one")
})

test_that("standardized DI tables keep their orientation", {
  di <- data.table::data.table(site_id = "s1", event_type = "SE", event_id = "SE:1",
    gene_id = "G", chr = "chr1", strand = "+", inc = "1-10", exc = "", n_samples = 4L,
    n_control = 2L, n_case = 2L, mean_psi_ctrl = 0.2, mean_psi_case = 0.6,
    delta_psi = 0.4, p.value = 0.01, padj = 0.05, cooks_max = 1, form = "INC",
    n = 4L, n_used = 4L)
  expect_equal(get_rmats_post_di(di), di)
  expect_error(get_rmats_post_di(di, case_group = 2), "already in SpliceImpactR format")
})
