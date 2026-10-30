multi_path <- test_path("..", "..", "inst", "validation", "mfrm-wide-map-multi-output-20261001.R")
skip_if_not(file.exists(multi_path), "Repository-only multi-output execution")
multi <- new.env(parent = globalenv())
for (name in c("mml-stages", "jml-stages", "output-summary", "multi-output"))
  sys.source(test_path("..", "..", "inst", "validation", paste0("mfrm-wide-map-", name, "-20261001.R")), multi)
multi_fixture <- function() {
  path <- test_path("..", "..", "validation-results", "native-location-branches-20261001", "fixed-gpcm-manifest.rds")
  skip_if_not(file.exists(path), "Retained witness arguments")
  args <- readRDS(path)$args; args$quad_points <- NULL
  list(fit = readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds")),
    spec = list(ConditionId = "control-flow", Arm = "MML", Replicate = 1L, InputId = "same",
      args = args, level = .95, purpose = "workflow_witness_only", outputs = list(
        raters = list(kind = "location", facet = "Rater"),
        criteria = list(kind = "location", facet = "Criterion"), slopes = list(kind = "slopes"))))
}
multi_mock <- function(...) {
  bindings <- list(...); old <- mget(names(bindings), envir = multi)
  list2env(bindings, multi); withr::defer(list2env(old, multi), envir = parent.frame())
}
multi_controls <- function(x, events, interrupt = function() FALSE, fail = function(q) FALSE) {
  list(wide_multi_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) {
      q <- args$quad_points; events$fits <- c(events$fits, q)
      if (fail(q)) stop("higher-order fit error")
      f <- x$fit; f$config$estimation_control$quad_points <- q; f
    }, wide_mml_score = function(fit) {
      q <- fit$config$estimation_control$quad_points; events$scores <- c(events$scores, q)
      list(ready = q >= 61L, order = q)
    }, wide_mml_scoring_retry = function(x) !x$ready,
    wide_multi_consume = function(fit, consumer, level) {
      key <- if (consumer$kind == "slopes") "slopes" else consumer$facet
      q <- fit$config$estimation_control$quad_points
      if (key == "Criterion" && q == 31L && interrupt())
        signalCondition(structure(list(message = "stop between consumers"), class = c("interrupt", "condition")))
      events$consumers <- c(events$consumers, paste(q, key))
      list(order = q, pass = q >= switch(key, Rater = 31L, Criterion = 61L, slopes = 121L))
    }, wide_multi_interval_table = function(value, spec, consumer, fit) {
      t <- multi$wide_multi_targets(spec, consumer)
      t$Estimate <- value$order; t$SE <- 1; t$Lower <- value$order - 1; t$Upper <- value$order + 1
      t$Available <- value$pass; t$Reason <- if (value$pass) "" else "integration review"
      t[!t$Available, c("SE", "Lower", "Upper")] <- NA_real_; t
    }, wide_multi_retry = function(value, consumer) !value$pass)
}

test_that("three consumers share fits and keep distinct first-passing stages", {
  x <- multi_fixture(); events <- new.env()
  do.call(multi_mock, multi_controls(x, events))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- multi$wide_multi_run(x$spec, out)
  expect_identical(events$fits, c(31L, 61L, 121L))
  expect_identical(events$scores, c(31L, 61L))
  expect_identical(events$consumers, c("31 Rater", "31 Criterion", "31 slopes", "61 Criterion", "61 slopes", "121 slopes"))
  expect_match(z$point$stage, "^q61:"); expect_match(z$intervals$raters$stage, "^q31:")
  expect_match(z$intervals$criteria$stage, "^q61:"); expect_match(z$intervals$slopes$stage, "^q121:")
  r <- multi$wide_multi_records(x$spec, z)
  expect_equal(nrow(r), 18L); expect_false(anyDuplicated(multi$wide_output_key(r)) > 0)
  expect_true(all(r$Available)); expect_equal(unique(r$Estimate[r$Output == "interval"]), c(31, 61, 121))
  multi$wide_mml_fit <- multi$wide_mml_score <- multi$wide_multi_consume <- function(...) stop("No repeated computation")
  expect_identical(multi$wide_multi_run(x$spec, out), z)
  bad <- x$spec; bad$outputs$slopes <- NULL
  expect_error(multi$wide_multi_run(bad, out), "output plan changed")
  expect_error(multi$wide_multi_records(bad, z), "exact shared job")
})

test_that("interruption between consumers preserves completed fits and intervals", {
  x <- multi_fixture(); events <- new.env(); interrupt <- TRUE
  do.call(multi_mock, multi_controls(x, events, interrupt = function() interrupt))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_identical(tryCatch(multi$wide_multi_run(x$spec, out), interrupt = function(e) "interrupted"), "interrupted")
  expect_identical(readRDS(file.path(out, "status.rds"))$status, "interrupted")
  path <- file.path(out, "q31-interval-raters.rds"); hash <- tools::md5sum(path)
  interrupt <- FALSE; z <- multi$wide_multi_run(x$spec, out)
  expect_identical(events$fits, c(31L, 61L, 121L))
  expect_identical(events$scores, c(31L, 61L)); expect_identical(tools::md5sum(path), hash)
  expect_equal(sum(events$consumers == "31 Rater"), 1L)
  expect_true(all(multi$wide_multi_records(x$spec, z)$Available))
})

test_that("a later fit error leaves an earlier qualified interval intact", {
  x <- multi_fixture(); events <- new.env()
  do.call(multi_mock, multi_controls(x, events, fail = function(q) q == 61L))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- multi$wide_multi_run(x$spec, out); r <- multi$wide_multi_records(x$spec, z)
  expect_identical(events$fits, c(31L, 61L))
  expect_match(z$intervals$raters$stage, "^q31:")
  expect_identical(z$intervals$criteria$error, "higher-order fit error")
  expect_equal(sum(r$Available), 3L); expect_equal(sum(r$Status == "error"), 15L)
})

test_that("owners, locations and log slopes cannot collide or be silently compared across scales", {
  x <- multi_fixture(); p <- multi$wide_multi_plan(x$spec)
  expect_equal(length(unique(p$Target)), 9L)
  expect_setequal(p$Coordinate, c("location", "log_slope"))
  x$spec$outputs$duplicate <- x$spec$outputs$raters
  expect_error(multi$wide_multi_validate(x$spec), "Duplicate target")
  # Equal numerical reference values are insufficient when the scale differs.
  p <- p[p$Output == "point", ]; p$Reference <- "generating_parameter"; p$ReferenceValue <- 0
  b <- p; b$Arm <- "B"; b$ScaleReference <- "different-scale"; plan <- rbind(p, b)
  r <- plan[c("ConditionId", "Arm", "Replicate", "Target", "Output")]
  r$Status <- "returned"; r$Available <- TRUE; r$Estimate <- 0; r$SE <- r$Lower <- r$Upper <- NA_real_
  r$SourceStage <- "fixture"; r$Reason <- ""
  expect_error(multi$wide_paired_point_summary(plan, r, "MML", "B"), "scale references")
  plan$ScaleReference <- rep(p$ScaleReference, 2)
  expect_true(all(multi$wide_paired_point_summary(plan, r, "MML", "B")$Complete))
})

test_that("bad retained-stage arguments cannot be imported", {
  x <- multi_fixture(); args <- x$spec$args; args$quad_points <- 61L
  bad <- list(`31` = list(args = args, fitting = list(value = x$fit, error = ""), provenance = "fixture"))
  expect_error(multi$wide_multi_run(x$spec, tempfile(), bad))
})

test_that("a consumer error does not suppress other outputs or cause its own refit", {
  x <- multi_fixture(); events <- new.env(); controls <- multi_controls(x, events)
  original <- controls$wide_multi_consume
  controls$wide_multi_consume <- function(fit, consumer, level) {
    if (identical(consumer$facet, "Criterion")) stop("consumer error")
    original(fit, consumer, level)
  }
  do.call(multi_mock, controls)
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- multi$wide_multi_run(x$spec, out); r <- multi$wide_multi_records(x$spec, z)
  expect_identical(events$fits, c(31L, 61L, 121L))
  expect_match(z$intervals$criteria$stage, "^q31:")
  expect_identical(z$intervals$criteria$error, "consumer error")
  expect_equal(sum(r$Available), 15L); expect_equal(sum(r$Status == "error"), 3L)
})

test_that("JML shares one fit and one source check across location and slope points", {
  root <- test_path("..", "..", "validation-results", "mfrm-wide-map-jml-stages-20261001", "JML2")
  skip_if_not(file.exists(file.path(root, "manifest.rds")))
  x <- readRDS(file.path(root, "manifest.rds"))$payload
  spec <- x$spec; spec$facet <- NULL; spec$level <- .95
  spec$outputs <- multi_fixture()$spec$outputs
  fit_calls <- score_calls <- 0L
  multi_mock(wide_multi_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) { fit_calls <<- fit_calls + 1L; expect_null(args$quad_points); x$initial$fit },
    wide_jml_source = function(fit) { score_calls <<- score_calls + 1L; list(ready = TRUE) },
    wide_multi_consume = function(...) stop("No JML structural CI procedure"))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- multi$wide_multi_run(spec, out); r <- multi$wide_multi_records(spec, z)
  expect_identical(fit_calls, 1L); expect_identical(score_calls, 1L)
  expect_equal(length(z$intervals), 0L)
  expect_true(all(z$plan$Eligibility[z$plan$Output == "interval"] == "unsupported"))
  expect_equal(nrow(r), 7L); expect_true(all(r$Available))
  expect_true(all(is.na(r$SE) & is.na(r$Lower)))
})
