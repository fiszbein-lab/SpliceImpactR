test_that("rMATS raw and post-DI converters use closed one-based intervals", {
  for (strand in c("+", "-")) {
    for (event_type in c("SE", "MXE", "A3SS", "A5SS", "RI")) {
      rm <- data.table::data.table(
        sample = "S1", condition = "case", source_file = "synthetic",
        ID = 1L, GeneID = "GA", chr = "chr1", strand = strand,
        event_type = event_type,
        upstreamES = 100L, upstreamEE = 150L,
        downstreamES = 400L, downstreamEE = 450L,
        exonStart_0base = 200L, exonEnd = 201L,
        longExonStart_0base = 200L, longExonEnd = 300L,
        shortES = if (strand == "+") 200L else 201L,
        shortEE = if (strand == "+") 299L else 300L,
        flankingES = 400L, flankingEE = 450L,
        IJC_SAMPLE_1 = 20L, SJC_SAMPLE_1 = 10L, IncLevel1 = 2/3,
        IncLevelDifference = 0.3, PValue = 0.001, FDR = 0.01
      )
      rm[, `:=`(`1stExonStart_0base` = 200L, `1stExonEnd` = 250L,
                `2ndExonStart_0base` = 300L, `2ndExonEnd` = 350L)]
      raw <- get_rmats(data.table::copy(rm))
      post <- get_rmats_post_di(data.table::copy(rm), event_type = event_type)
      for (result in list(raw, post)) {
        inc <- result[form == "INC"]
        exc <- result[form == "EXC"]
        expect_equal(nrow(result), 2L)
        if (event_type == "SE") {
          expect_identical(inc$inc, "101-150;201-201;401-450")
          expect_identical(exc$exc, "201-201")
        } else if (event_type == "MXE") {
          exon1 <- if (strand == "+") "201-250" else "301-350"
          exon2 <- if (strand == "+") "301-350" else "201-250"
          expect_identical(inc$inc, paste0("101-150;", exon1, ";401-450"))
          expect_identical(inc$exc, exon2)
          expect_identical(exc$inc, paste0("101-150;", exon2, ";401-450"))
          expect_identical(exc$exc, exon1)
        } else if (event_type %in% c("A3SS", "A5SS")) {
          # Each form keeps its partner (flanking) exon, in genomic order.
          expect_identical(inc$inc, "201-300;401-450")
          expect_identical(exc$inc, paste0(if (strand == "+") "201-299" else "202-300", ";401-450"))
          expect_identical(exc$exc, if (strand == "+") "300-300" else "201-201")
        } else {
          expect_identical(inc$inc, "101-150;151-400;401-450")
          expect_identical(exc$inc, "101-150;401-450")
          expect_identical(exc$exc, "151-400")
        }
      }
    }
  }
})

test_that("event lengths use the selected transcript and retain unmatched rows", {
  ann <- data.table::data.table(
    type = c("exon", "exon", "exon", "CDS"),
    transcript_id = c("T1", "T1", "T2", "T1"), exon_id = "E1",
    cds_len = c(90L, 90L, 60L, 90L), feature_length = c(100L, 100L, 100L, 90L)
  )
  hits <- data.table::data.table(
    event_type = "SE", transcript_id_case = c("T1", "T2", "missing", "T1"),
    transcript_id_control = c("T2", "T1", "T1", "T1"),
    exons_case = c("E1;E1", "E1", "E1", ""), exons_control = "E1"
  )
  case <- SpliceImpactR:::.sum_exon_lengths(hits, ann, "exons_case")
  control <- SpliceImpactR:::.sum_exon_lengths(hits, ann, "exons_control")
  expect_equal(case$row_id, 1:4)
  expect_equal(case$cds_len, c(90L, 60L, NA_integer_, NA_integer_))
  expect_equal(case$exon_len, c(100L, 100L, NA_integer_, NA_integer_))
  expect_equal(control$cds_len, c(60L, 90L, 90L, 90L))
  partial <- data.table::copy(hits[1])
  partial[, exons_case := "E1;missing"]
  expect_true(is.na(SpliceImpactR:::.sum_exon_lengths(partial, ann, "exons_case")$cds_len))
  expect_equal(nrow(SpliceImpactR:::.sum_exon_lengths(hits[0], ann, "exons_case")), 0L)
  conflicting <- data.table::copy(ann)
  conflicting$cds_len[2] <- 99L
  expect_error(SpliceImpactR:::.sum_exon_lengths(hits, conflicting, "exons_case"), "conflicting lengths")
})

test_that("PPI evidence is attributed to the corresponding gene endpoint", {
  ppi <- data.table::data.table(
    geneA = c("GA", "GB", "GA", "GA"), geneB = c("GB", "GA", "GA", "GB"),
    DDI = c(TRUE, TRUE, TRUE, NA), DMI = c(TRUE, TRUE, TRUE, NA),
    DDI_A = list("PF_A", "PF_B", "PF_SELF", "PF_A"),
    DDI_B = list("PF_B", "PF_A", "PF_A", "PF_B"),
    DMI_A = list("PF_A", "LIG_B", "PF_SELF", "PF_A"),
    DMI_B = list("LIG_B", "PF_A", "LIG_A", "LIG_B")
  )
  before <- data.table::copy(ppi)
  own <- SpliceImpactR:::mark_changing_partners_split(ppi, "GA", "PF_A", character())
  expect_equal(own$DDI_changed_case, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(own$DMI_changed_case, c(TRUE, TRUE, FALSE, FALSE))
  partner_only <- SpliceImpactR:::mark_changing_partners_split(
    ppi, "GA", "PF_B", character(), changed_motif_case = "LIG_B"
  )
  expect_false(any(partner_only$interaction_changed_case))
  motif <- SpliceImpactR:::mark_changing_partners_split(
    ppi, "GB", character(), character(), changed_motif_control = "LIG_B"
  )
  expect_equal(motif$DMI_changed_control, c(TRUE, TRUE, FALSE))
  expect_false(any(motif$DDI_changed_control))
  expect_identical(ppi, before)
  empty <- SpliceImpactR:::mark_changing_partners_split(ppi, "absent", "PF_A", character())
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("partner_gene", "interaction_changed_case") %in% names(empty)))
})
