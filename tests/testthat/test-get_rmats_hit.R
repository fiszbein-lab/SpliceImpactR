test_that("get_rmats_hit runs on test extdata and returns expected structure", {
  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/case_S1/')
    ),
    sample_name = c("S5", "S1"),
    condition = c("control", "case"),
    stringsAsFactors = FALSE
  )

  # restrict event types to minimize runtime for test mode
  data <- get_rmats_hit(
    sf,
    event_types = c("AFE", "ALE"), # just HIT-based events
    keep_annotated_first_last = TRUE
  )

  expect_s3_class(data, "data.table")
  expect_gt(nrow(data), 0)

  req <- c("event_id","event_type","form","gene_id","chr","strand",
           "inc","exc","inclusion_reads","exclusion_reads","psi",
           "sample","condition","source_file")
  expect_true(all(req %in% colnames(data)))

  # event types limited to AFE/ALE
  expect_true(all(data$event_type %in% c("AFE","ALE")))
})

test_that("get_rmats_hit includes rmats when requested", {
  sf <- data.frame(
    path = c(
      check_extdata_dir('rawData/control_S5/'),
      check_extdata_dir('rawData/case_S1/')
    ),
    sample_name = c("S5", "S1"),
    condition = c("control", "case"),
    stringsAsFactors = FALSE
  )

  data <- get_rmats_hit(
    sf,
    event_types = c("SE","AFE"),  # SE triggers rmats branch + AFE triggers HIT
    keep_annotated_first_last = FALSE
  )

  expect_true(any(data$event_type %in% "SE"))
  expect_true(any(data$event_type %in% "AFE"))
})

test_that("get_rmats_hit errors when sample_frame is missing required columns", {
  bad <- data.frame(path = ".", stringsAsFactors = FALSE)  # missing condition & sample_name
  expect_error(get_rmats_hit(bad))
})
