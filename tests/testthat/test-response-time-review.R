test_that("response_time_review summarizes valid timing rows", {
  toy <- load_mfrmr_data("example_core")
  toy$ResponseTime <- 8 + (seq_len(nrow(toy)) %% 9) +
    as.numeric(toy$Score)
  toy$ResponseTime[1] <- 1
  toy$ResponseTime[2] <- 40

  rt <- response_time_review(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    time = "ResponseTime",
    rapid_quantile = 0.10,
    slow_quantile = 0.90
  )

  expect_s3_class(rt, "mfrm_response_time_review")
  expect_true(all(c(
    "overview", "thresholds", "observations", "person_summary",
    "facet_summary", "score_summary", "flags", "notes"
  ) %in% names(rt)))
  expect_equal(rt$overview$ValidRows, nrow(toy))
  expect_true(all(c("RapidFlag", "SlowFlag", "LogTime")
                  %in% names(rt$observations)))
  expect_true(all(c("Facet", "Level", "MedianTime", "RapidRate", "SlowRate")
                  %in% names(rt$facet_summary)))
  expect_gt(nrow(rt$score_summary), 0L)

  sx <- summary(rt)
  expect_s3_class(sx, "summary.mfrm_response_time_review")
  expect_true("top_rapid_persons" %in% names(sx))
})

test_that("response_time_review rejects unusable timing inputs", {
  toy <- load_mfrmr_data("example_core")
  toy$ResponseTime <- NA_real_

  expect_error(
    response_time_review(toy, "Person", c("Rater"), "ResponseTime"),
    "No rows contain finite response times"
  )
  expect_error(
    response_time_review(toy, "Person", c("MissingFacet"), "ResponseTime"),
    "facets"
  )
  expect_error(
    response_time_review(toy, "Person", c("Rater"), "ResponseTime",
                         rapid_quantile = 0.9, slow_quantile = 0.1),
    "rapid_quantile"
  )
})

test_that("plot_response_time_review returns reusable plot data", {
  toy <- load_mfrmr_data("example_core")
  toy$ResponseTime <- 10 + (seq_len(nrow(toy)) %% 6) +
    as.numeric(toy$Score)
  rt <- response_time_review(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    time = "ResponseTime"
  )

  for (tp in c("distribution", "person", "facet", "score")) {
    p <- plot_response_time_review(rt, type = tp, draw = FALSE)
    expect_s3_class(p, "mfrm_plot_data")
    expect_identical(p$name, "response_time_review")
    expect_true(is.data.frame(p$data$table))
    expect_true(is.data.frame(p$data$thresholds))
    expect_identical(p$data$type, tp)
  }

  p_method <- plot(rt, type = "person", draw = FALSE)
  expect_s3_class(p_method, "mfrm_plot_data")
  components <- plot_data_components(p_method)
  expect_true(any(components$Component == "table" &
                    components$Role == "primary_data"))
  expect_true(any(components$Component == "thresholds" &
                    components$Role == "settings"))
  rt_table <- plot_data(p_method, component = "table")
  expect_s3_class(rt_table, "data.frame")
  expect_true(all(c("Person", "MedianTime", "RapidRate", "SlowRate")
                  %in% names(rt_table)))
  expect_error(
    plot_response_time_review(
      response_time_review(toy, "Person", time = "ResponseTime"),
      type = "score",
      draw = FALSE
    ),
    "Supply `score`"
  )
})
