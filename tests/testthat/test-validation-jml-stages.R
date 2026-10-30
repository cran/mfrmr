jml_script <- test_path("..", "..", "inst", "validation", "mfrm-wide-map-jml-stages-20261001.R")
skip_if_not(file.exists(jml_script), "Repository-only JML execution helper")
jml_env <- new.env(parent = globalenv())
sys.source(test_path("..", "..", "inst", "validation", "mfrm-wide-map-mml-stages-20261001.R"), jml_env)
sys.source(jml_script, jml_env)
jml_fixture <- function(order = NULL) {
  root <- test_path("..", "..", "validation-results", "mfrm-facet-structure-pilot-20261001")
  p <- file.path(root, "fits", paste0("20-base_two_facet-JML", order, ".rds"))
  skip_if_not(file.exists(p), "Repository-only retained JML pilot")
  fit <- readRDS(p)$fit
  d <- readRDS(file.path(root, "inputs", "20-base_two_facet.rds"))$data
  args <- list(data = d, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
    model = "GPCM", method = "JML", step_facet = "Criterion", slope_facet = "Criterion",
    category_policy = "preserve", rating_min = 0, rating_max = 2, maxit = 400L)
  if (is.null(order)) { args$optimizer <- "BFGS"; args$reltol <- 1e-9 } else {
    args$jml_correction_order <- as.integer(order); args$jml_correction_sampling <- "fixed_rosters"
  }
  list(fit = fit, spec = list(ConditionId = "retained-control-flow", Arm = paste0("JML", order),
    Replicate = 1L, InputId = "same-input", args = args, facet = "Rater", purpose = "workflow_witness_only"))
}
jml_mock <- function(...) {
  bindings <- list(...); old <- mget(names(bindings), envir = jml_env)
  list2env(bindings, jml_env); withr::defer(list2env(old, jml_env), envir = parent.frame())
}

test_that("JML admits explicit procedures without MML controls or implicit order selection", {
  for (order in list(NULL, 2L, 4L)) expect_true(jml_env$wide_jml_validate(jml_fixture(order)$spec))
  bad <- jml_fixture(2L)$spec; bad$args$optimizer <- "BFGS"
  expect_error(jml_env$wide_jml_validate(bad))
  bad <- jml_fixture()$spec; bad$args$quad_points <- 31L
  expect_error(jml_env$wide_jml_validate(bad))
  bad <- jml_fixture(2L)$spec; bad$args$jml_correction_order <- 3L
  expect_error(jml_env$wide_jml_validate(bad))
})

test_that("JML resumes source review without repeating its public fit", {
  x <- jml_fixture(); calls <- 0L; interrupted <- TRUE
  jml_mock(wide_jml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) { calls <<- calls + 1L; expect_identical(args, x$spec$args); x$fit },
    wide_jml_source = function(fit) {
      if (interrupted) signalCondition(structure(list(message = "interrupted"), class = c("interrupt", "condition")))
      list(ready = TRUE)
    })
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_identical(tryCatch(jml_env$wide_jml_run(x$spec, out), interrupt = function(e) "interrupted"), "interrupted")
  expect_identical(readRDS(file.path(out, "status.rds"))$status, "interrupted")
  interrupted <- FALSE; z <- jml_env$wide_jml_run(x$spec, out)
  expect_identical(calls, 1L); expect_null(z$root); expect_false(z$formal_structural_intervals)
  jml_env$wide_mml_fit <- jml_env$wide_jml_source <- function(...) stop("No recomputation")
  expect_identical(jml_env$wide_jml_run(x$spec, out), z)
  expect_true(all(jml_env$wide_jml_records(x$spec, z)$Available))
  bad <- x$spec; bad$InputId <- "different"
  expect_error(jml_env$wide_jml_records(bad, z), "exact job")
  expect_error(jml_env$wide_jml_run(bad, out), "data/settings/source changed")
})

test_that("corrected roots and covariance remain separate from source admission and truth intervals", {
  x <- jml_fixture(2L); x$fit$jml_adjustment$covariance$available <- FALSE
  x$fit$jml_adjustment$point$status <- "root_with_unresolved_starts"
  jml_mock(wide_jml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) x$fit, wide_jml_source = function(fit) list(ready = TRUE))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- jml_env$wide_jml_run(x$spec, out)
  expect_identical(z$root$status, "root_with_unresolved_starts")
  expect_false(z$root_covariance$available)
  expect_identical(z$root$theta, x$fit$jml_adjustment$point$theta)
  expect_true(any(is.infinite(z$root$theta)))
  records <- jml_env$wide_jml_records(x$spec, z)
  expect_true(all(records$Available)); expect_true(all(is.na(records$SE)))
  plan <- jml_env$wide_jml_output_plan(x$spec)
  expect_true(all(plan$Eligibility[plan$Output == "interval"] == "unsupported"))
  expect_identical(unique(records$Output), "point")
})

test_that("a JML source error is settled once and keeps its initial structural points", {
  x <- jml_fixture(4L); calls <- 0L
  jml_mock(wide_jml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) { calls <<- calls + 1L; x$fit },
    wide_jml_source = function(...) stop("unresolved source"))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- jml_env$wide_jml_run(x$spec, out)
  expect_identical(calls, 1L)
  records <- jml_env$wide_jml_records(x$spec, z)
  expect_true(all(records$Status == "error" & !records$Available))
  expect_true(all(is.finite(records$Estimate)))
  expect_identical(jml_env$wide_jml_run(x$spec, out), z)
  expect_identical(calls, 1L)
})
