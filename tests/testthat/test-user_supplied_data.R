test_that("get_user_data works on example", {
  df <- data.frame(
    event_id = rep("A3SS:1", 8),
    event_type = "A3SS",
    form = rep(c("INC","EXC"), each=4),
    gene_id="ENSG00000158286",
    chr="chrX",
    strand="-",
    inc=c(rep("149608626-149608834",4), rep("149608626-149608829",4)),
    exc=c(rep("",4), rep("149608830-149608834",4)),
    inclusion_reads=c(30,32,29,31, 2,3,4,3),
    exclusion_reads=c(1,1,2,1, 28,27,26,30),
    sample=c("S1","S2","S3","S4","S1","S2","S3","S4"),
    condition=rep(c("case","case","control","control"),2)
  )
  out <- get_user_data(df)
  expect_true("psi" %in% names(out))
  expect_true(all(out$event_type == "A3SS"))
})
test_that("get_user_data_post_di expands SITE events", {
  example_user_data <- data.frame(
    event_id = rep("ENSG00000158286:AFE", 8),
    event_type = "AFE",
    gene_id = "ENSG00000158286",
    chr = "chrX",
    strand = "-",
    form = rep(c("SITE"), each = 8),
    inc = c(
      rep("149608626-149608834", 4),
      rep("149608626-149608829", 4)
    ),
    exc = c(
      rep("", 4),
      rep("", 4)
    ),
    inclusion_reads = c(30, 28, 25, 32,  2, 3, 4, 3),
    exclusion_reads = c(1, 2, 1, 1, 28, 27, 26, 30),
    sample = c("S1","S2","S3","S4","S1","S2","S3","S4"),
    condition = rep(c("case","case","control","control"), 2),
    stringsAsFactors = FALSE
  )

  out <- get_user_data_post_di(example_user_data)
  expect_equal(nrow(out), 16) # duplicated +1 / -1
  expect_setequal(out$delta_psi, rep(c(1,-1), each = 8))
})

test_that("get_user_data_post_di handles INC/EXC", {
  example_user_data <- data.frame(
    event_id = rep("A3SS:1", 8),
    event_type = "A3SS",
    gene_id = "ENSG00000158286",
    chr = "chrX",
    strand = "-",
    form = rep(c("INC","EXC"), each = 4),
    inc = c(
      rep("149608626-149608834", 4),
      rep("149608626-149608829", 4)
    ),
    exc = c(
      rep("", 4),
      rep("149608830-149608834", 4)
    ),
    inclusion_reads = c(30, 28, 25, 32,  2, 3, 4, 3),
    exclusion_reads = c(1, 2, 1, 1, 28, 27, 26, 30),
    sample = c("S1","S2","S3","S4","S1","S2","S3","S4"),
    condition = rep(c("case","case","control","control"), 2),
    stringsAsFactors = FALSE
  )
  out <- get_user_data_post_di(example_user_data)
  expect_setequal(out$delta_psi, c(1,-1))
})

test_that("get_user_data_post_di requires event IDs and INC+EXC or SITE forms", {
  ev <- data.frame(event_id = "E1", event_type = "SE", gene_id = "G1", chr = "chr1",
    strand = "+", form = c("INC", "EXC"),
    inc = c("100-200;300-400;500-600", "100-200;500-600"), exc = c("", "300-400"),
    delta_psi = c(0.3, -0.3))
  expect_equal(nrow(get_user_data_post_di(ev)), 2L)
  # Generated IDs would give each form its own event, so none are generated.
  expect_error(get_user_data_post_di(ev[, names(ev) != "event_id"]),
               "Missing required columns: event_id")
  expect_error(get_user_data_post_di(transform(ev, event_id = c("E1", ""))),
               "missing or empty in 1 row\\(s\\): 2$")
  expect_error(get_user_data_post_di(transform(ev, event_id = c(NA, "E1"))),
               "missing or empty in 1 row\\(s\\): 1$")
  expect_error(get_user_data_post_di(rbind(ev, transform(ev[1, ], event_id = "E2"))),
               "do not have INC\\+EXC or SITE: E2$")
  expect_error(get_user_data_post_di(rbind(ev, transform(ev[1, ], form = "SITE"))),
               "both SITE and INC/EXC forms: E1$")
})

test_that("get_user_data names events with invalid forms", {
  df <- data.frame(event_id = c("E1", "E1", "E2", "E2"), form = c("INC", "EXC", "INC", "INC"),
    gene_id = "G1", chr = "chr1", strand = "+", inc = "100-200", exc = "",
    inclusion_reads = 10, exclusion_reads = 5, sample = c("S1", "S2"), condition = "case")
  expect_error(get_user_data(df), "do not have INC\\+EXC or SITE: E2$")
})
