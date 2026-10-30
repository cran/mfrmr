test_that("fit progress is display-only and restores the enclosing context", {
  events <- list()
  record <- function(type, value) events[[length(events) + 1L]] <<- list(type, value)
  local_mocked_bindings(
    cli_progress_bar = function(...) { record("start", list(...)); "test-fit" },
    cli_progress_update = function(...) record("update", list(...)),
    cli_progress_done = function(...) record("done", list(...)), .package = "mfrmr")
  withr::local_options(list(mfrmr._fit_progress = "outer-context"))
  d <- load_mfrmr_data("example_operational")
  args <- list(data = d, person = "Person", facets = c("Rater", "Criterion"),
    score = "Score", attach_diagnostics = FALSE)
  plain <- suppressWarnings(do.call(fit_mfrm, args))
  expect_length(events, 0L)
  shown <- suppressWarnings(do.call(fit_mfrm, c(args, list(progress = TRUE))))
  expect_identical(shown$opt$par, plain$opt$par)
  expect_identical(shown$opt$value, plain$opt$value)
  expect_identical(shown$config$replay_inputs, plain$config$replay_inputs)
  expect_false("progress" %in% names(shown$config$replay_inputs))
  expect_identical(getOption("mfrmr._fit_progress"), "outer-context")
  expect_identical(events[[1]][[1]], "start")
  expect_true(any(vapply(events, function(x) x[[1]] == "update", FALSE)))
  expect_identical(tail(events, 1)[[1]][[2]]$result, "clear")
  for (bad in list(NA, NULL, "yes", 1, c(TRUE, FALSE)))
    expect_error(do.call(fit_mfrm, c(args, list(progress = bad))), "progress")
  expect_error(fit_mfrm(d[FALSE, ], "Person", c("Rater", "Criterion"), "Score",
    progress = TRUE), "zero rows")
  expect_identical(getOption("mfrmr._fit_progress"), "outer-context")
  expect_identical(tail(events, 1)[[1]][[2]]$result, "clear")
})

test_that("two-family progress is cleared on interruption", {
  closed <- FALSE
  local_mocked_bindings(cli_progress_bar = function(...) "test-interrupt",
    cli_progress_done = function(...) closed <<- TRUE, .package = "mfrmr")
  withr::local_options(list(mfrmr._fit_progress = NULL))
  local_mocked_bindings(mfrm_fit_product_slopes = function(...) {
    stop(structure(list(message = "interrupted", call = NULL),
      class = c("interrupt", "condition")))
  }, .package = "mfrmr")
  d <- data.frame(Person = "p", Task = "t", Rater = "r", Score = 1)
  interrupted <- tryCatch(fit_mfrm(d, "Person", c("Task", "Rater"), "Score",
    model = "GPCM", slope_facet = c("Task", "Rater"), progress = TRUE),
    interrupt = identity)
  expect_s3_class(interrupted, "interrupt")
  expect_true(closed)
  expect_null(getOption("mfrmr._fit_progress"))
})
