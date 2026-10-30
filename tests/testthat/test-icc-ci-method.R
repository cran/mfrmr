# ICC intervals must use the jointly refitted ratio and retain failed attempts.
skip_if_not_installed("lme4")

.icc_data <- load_mfrmr_data("example_core")
.icc_call <- function(...) {
  compute_facet_icc(.icc_data, facets = c("Rater", "Criterion"),
                    score = "Score", person = "Person", ...)
}

test_that("point-only results remain usable and profile intervals are refused", {
  x <- .icc_call()
  expect_s3_class(x, "mfrm_facet_icc")
  expect_true(all(is.na(x$ICC_CI_Lower) & is.na(x$ICC_CI_Upper)))
  expect_true(all(x$ICC_CI_Status == "Not requested"))
  attr(x, "icc_ci") <- NULL
  expect_output(print(x), "Not requested")
  expect_error(.icc_call(ci_method = "profile"), "withdrawn.*ICC confidence")
  expect_error(.icc_call(ci_method = "nonsense"), "arg")
  expect_error(analyze_hierarchical_structure(
    .icc_data, c("Rater", "Criterion"), ci_method = "profile"
  ), "withdrawn")
  expect_error(suppressWarnings(analyze_hierarchical_structure(
    .icc_data, c("Rater", "Criterion"), icc_ci_method = "profile"
  )), "withdrawn")
})

test_that("bootstrap controls reject invalid values instead of truncating", {
  for (level in list(0, 1.2, NA_real_, c(.9, .95), .95 + 1i)) {
    expect_error(.icc_call(ci_level = level), "ci_level")
  }
  for (reps in list(0, 1, 2.5, NA_real_, c(5, 10), Inf, 2^32)) {
    expect_error(.icc_call(ci_method = "boot", ci_boot_reps = reps), "ci_boot_reps")
  }
  for (ncpu in list(0, 1.5, NA_real_)) {
    expect_error(.icc_call(ci_method = "boot", ci_boot_ncpus = ncpu), "ci_boot_ncpus")
  }
  for (seed in list(-1, 1.5, NA_real_, c(1, 2), 2^32)) {
    expect_error(.icc_call(ci_method = "boot", ci_boot_seed = seed), "ci_boot_seed")
  }
})

test_that("bootstrap matches independently calculated joint variance ratios", {
  # A small count keeps this regression test quick; it is not a reporting recommendation.
  n <- 24L
  data <- .icc_data
  data$Score <- as.numeric(data$Score)
  for (column in c("Person", "Rater", "Criterion")) data[[column]] <- factor(data[[column]])
  fit <- lme4::lmer(Score ~ 1 + (1 | Person) + (1 | Rater) + (1 | Criterion),
                    data = data, REML = TRUE)
  reference <- suppressWarnings(lme4::bootMer(
    fit, nsim = n, seed = 2026L, type = "parametric", use.u = FALSE,
    FUN = function(f) {
      vc <- as.data.frame(lme4::VarCorr(f))
      stats::setNames(vc$vcov / sum(vc$vcov), vc$grp)
    }
  ))
  x <- .icc_call(ci_method = "boot", ci_boot_reps = n, ci_boot_seed = 2026L)
  details <- attr(x, "icc_ci")$bootstrap
  expect_equal(details$draws, reference$t[, x$Facet, drop = FALSE])
  expect_equal(rowSums(details$draws), rep(1, n), tolerance = 1e-10)
  expect_true(all(x$ICC_CI_NRequested == n))
  expect_true(all(x$ICC_CI_NReps == sum(details$usable)))
  expect_true(all(x$ICC_CI_NReps + x$ICC_CI_NUnavailable == n))
  expect_length(details$converged, n)
  expect_length(details$singular, n)
  if (all(x$ICC_CI_Status == "Available")) {
    bounds <- apply(reference$t, 2, quantile, probs = c(.025, .975))
    expect_equal(x$ICC_CI_Lower, unname(round(bounds[1, x$Facet], 4)))
    expect_equal(x$ICC_CI_Upper, unname(round(bounds[2, x$Facet], 4)))
  } else {
    expect_true(all(is.na(x$ICC_CI_Lower) & is.na(x$ICC_CI_Upper)))
  }
  y <- .icc_call(ci_method = "boot", ci_boot_reps = n, ci_boot_seed = 2026L)
  expect_identical(attr(x, "icc_ci")$bootstrap, attr(y, "icc_ci")$bootstrap)
})

# Mock only the expensive lme4 bootstrap, preserving the public fit and result path.
# Indicators follow the function evaluated by bootMer: all component ICCs, then
# convergence and singularity. This exercises deliberate missing/error outcomes.
.icc_boot_result <- function(x, FUN, nsim, ...) {
  initial <- FUN(x)
  t <- matrix(rep(initial, each = nsim), nrow = nsim)
  list(t = t)
}

test_that("an unavailable draw withholds intervals and retains lme4 error messages", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(...) {
    b <- .icc_boot_result(...)
    b$t[2, ] <- NA_real_
    attr(b, "bootFail") <- 1L
    attr(b, "boot.fail.msgs") <- table("refit failed")
    attr(b, "boot.all.msgs") <- list(`factory-error` = table("refit failed"))
    warning("some bootstrap runs failed (1/5)")
    b
  }, .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 5)
  expect_true(all(is.na(x$ICC_CI_Lower) & is.na(x$ICC_CI_Upper)))
  expect_true(all(x$ICC_CI_Status == "Incomplete bootstrap"))
  expect_true(all(x$ICC_CI_NRequested == 5 & x$ICC_CI_NReps == 4 &
                    x$ICC_CI_NUnavailable == 1))
  d <- attr(x, "icc_ci")$bootstrap
  expect_equal(nrow(d$draws), 5)
  expect_true(all(is.na(d$draws[2, ])))
  expect_equal(d$n_errors, 1L)
  expect_equal(names(d$failure_messages), "refit failed")
  expect_match(d$warnings, "some bootstrap runs failed")
  expect_output(print(x), "Incomplete bootstrap")
  expect_output(summary(x), "Incomplete bootstrap")
})

test_that("nonconverged draws remain visible but withhold intervals", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(...) {
    b <- .icc_boot_result(...)
    b$t[2, ncol(b$t) - 1L] <- 0
    b
  }, .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 5)
  expect_true(all(x$ICC_CI_NReps == 4))
  expect_true(all(is.finite(attr(x, "icc_ci")$bootstrap$draws)))
  expect_true(all(is.na(x$ICC_CI_Lower)))
})

test_that("converged boundary draws enter percentile intervals", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(x, FUN, nsim, ...) {
    b <- .icc_boot_result(x, FUN, nsim)
    # Two full decompositions, with a zero component in the first draw.
    b$t[, 1:4] <- rbind(c(0, .2, .3, .5), c(.4, .1, .1, .4))
    b$t[, 5:6] <- cbind(c(1, 1), c(1, 0))
    attr(b, "boot.all.msgs") <- list(`factory-message` = table("boundary (singular) fit"))
    b
  }, .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 2, ci_level = .8)
  expect_true(all(x$ICC_CI_Status == "Available"))
  expect_equal(x$ICC_CI_Lower, c(.04, .11, .12, .41))
  expect_equal(x$ICC_CI_Upper, c(.36, .19, .28, .49))
  expect_equal(attr(x, "icc_ci")$bootstrap$singular, c(TRUE, FALSE))
})

test_that("warnings and aborted bootstraps remain visible without interval claims", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(...) {
    b <- .icc_boot_result(...)
    attr(b, "boot.all.msgs") <- list(`factory-warning` = table("Hessian warning"))
    b
  }, .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 5)
  expect_true(all(x$ICC_CI_Status == "Bootstrap warnings require review"))
  expect_true(all(x$ICC_CI_NReps == 5 & is.na(x$ICC_CI_Lower)))
  expect_equal(names(attr(x, "icc_ci")$bootstrap$messages$`factory-warning`),
               "Hessian warning")
})

test_that("an aborted bootstrap reports unknown replicate counts and its error", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(...) stop("simulation failed"),
                                .package = "mfrmr")
  expect_message(x <- .icc_call(ci_method = "boot", ci_boot_reps = 5), "simulation failed")
  expect_true(all(x$ICC_CI_Status == "Bootstrap failed"))
  expect_true(all(x$ICC_CI_NRequested == 5))
  expect_true(all(is.na(x$ICC_CI_NReps) & is.na(x$ICC_CI_NUnavailable)))
  expect_true(all(is.na(x$ICC_CI_Lower)))
  expect_identical(attr(x, "icc_ci")$error, "simulation failed")
})

test_that("original convergence diagnostics stop interval computation", {
  original_lmer <- lme4::lmer
  testthat::local_mocked_bindings(.mfrm_icc_lmer = function(...) {
    f <- original_lmer(...)
    f@optinfo$conv$lme4 <- list(code = -1L, messages = "gradient failure")
    f
  }, .mfrm_boot_mer = function(...) stop("must not run"), .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 5)
  expect_true(all(x$ICC_CI_Status == "Original fit requires review"))
  expect_true(all(x$ICC_CI_NReps == 0 & x$ICC_CI_NUnavailable == 5))
  expect_false(attr(x, "icc_ci")$fit$converged)
  expect_identical(attr(x, "icc_ci")$fit$convergence$lme4$messages, "gradient failure")
})

test_that("hierarchical wrapper retains CI results and saved intervals require rerunning", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = .icc_boot_result, .package = "mfrmr")
  h <- analyze_hierarchical_structure(.icc_data, c("Rater", "Criterion"),
                                      ci_method = "boot", ci_boot_reps = 5,
                                      igraph_layout = FALSE)
  expect_true(all(h$icc$ICC_CI_Status == "Available"))
  expect_equal(attr(h$icc, "icc_ci")$calculation_version, 2L)
  old <- h$icc
  attr(old, "icc_ci")$calculation_version <- 1L
  expect_error(print(old), "recomputed")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(h, type = "icc"))
  # Keep annotations out of the figure and available in the returned notes.
  title <- NULL
  testthat::local_mocked_bindings(barplot = function(height, main, ...) {
    title <<- main
    graphics::plot.new()
    seq_along(height)
  }, arrows = function(...) NULL, .package = "mfrmr")
  before <- h
  expect_no_error(p <- plot(h, type = "icc"))
  expect_identical(title, "")
  expect_identical(h, before)
  expect_identical(p$icc, h$icc)
  expect_match(plot_data(p, "notes")$Text, "95% parametric bootstrap CI")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(p, path)
  expect_identical(plot_data(readRDS(path), "notes"), plot_data(p, "notes"))
  if (requireNamespace("callr", quietly = TRUE)) {
    # A source-loaded test must not reopen its result with an older installed release.
    source_path <- getNamespaceInfo(asNamespace("mfrmr"), "path")
    if (!file.exists(file.path(source_path, "R", "api-hierarchical-audit.R"))) source_path <- NULL
    reopened <- callr::r(function(file, source_path) {
      if (!is.null(source_path)) pkgload::load_all(source_path, quiet = TRUE)
      mfrmr::plot_data(readRDS(file), "notes")
    }, args = list(path, source_path), libpath = .libPaths())
    expect_identical(reopened, plot_data(p, "notes"))
  }
  expect_no_error(crossed <- plot(p, type = "crosstab"))
  expect_null(attr(crossed, "mfrmr_notes"))
  partial <- h
  partial$icc$ICC_CI_Lower[1] <- partial$icc$ICC_CI_Upper[1] <- NA_real_
  partial$icc$ICC_CI_Status[1] <- "Incomplete bootstrap"
  expect_no_error(p <- plot(partial, type = "icc"))
  expect_match(plot_data(p, "notes")$Text, "95% parametric bootstrap CI")
  expect_match(plot_data(p, "notes")$Text, "Intervals unavailable: Incomplete bootstrap")
  for (method in c("profile", "boot")) {
    old <- h$icc
    old$ICC_CI_Method <- method
    attr(old, "icc_ci") <- NULL
    expect_error(print(old), "recomputed")
    expect_error(summary(old), "recomputed")
    h$icc <- old
    expect_error(summary(h), "recomputed")
    expect_error(plot(h, type = "icc"), "recomputed")
  }
  # A current failed result explicitly says why the plot has no intervals.
  h$icc <- .icc_call(ci_method = "boot", ci_boot_reps = 5)
  h$icc$ICC_CI_Lower <- h$icc$ICC_CI_Upper <- NA_real_
  h$icc$ICC_CI_Status <- "Incomplete bootstrap"
  expect_no_error(p <- plot(h, type = "icc"))
  expect_identical(title, "")
  expect_match(plot_data(p, "notes")$Text, "Intervals unavailable: Incomplete bootstrap")
})

test_that("constant-response bootstrap refits cannot contribute numerical variance residue", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(x, FUN, nsim, ...) {
    constant <- suppressMessages(suppressWarnings(lme4::refit(x, rep(5, stats::nobs(x)))))
    draw <- FUN(constant)
    list(t = matrix(rep(draw, each = nsim), nrow = nsim))
  }, .package = "mfrmr")
  x <- .icc_call(ci_method = "boot", ci_boot_reps = 2)
  expect_true(all(x$ICC_CI_Status == "Incomplete bootstrap"))
  expect_true(all(x$ICC_CI_NReps == 0 & x$ICC_CI_NUnavailable == 2))
  expect_true(all(is.na(attr(x, "icc_ci")$bootstrap$draws)))
})
