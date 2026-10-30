# Isolate refit orchestration using a retained calibration. Actual subgroup
# numerical results are replayed separately in the validation record.
dff_refit_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      fit <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))$RSM$fit
      cached <<- list(fit = fit, diag = diagnose_mfrm(fit, residual_pca = "none"),
        data = load_mfrmr_data("example_core"))
    }
    cached
  }
})

dff_refit_mock_success <- function(x, env = parent.frame()) {
  local_mocked_bindings(
    fit_mfrm = function(data, ...) {
      sub <- x$fit
      sub$prep$data <- data
      sub$.test_group <- as.character(data$Group[1])
      sub
    },
    diagnose_mfrm = function(fit, ...) {
      out <- x$diag
      take <- out$measures$Facet != "Rater" |
        out$measures$Level %in% fit$prep$data$Rater
      out$measures <- out$measures[take, ]
      i <- out$measures$Facet == "Rater"
      out$measures$Estimate[i & out$measures$Level == "R02"] <-
        if (fit$.test_group == sort(unique(x$data$Group))[1]) .5 else -.5
      out$measures$Estimate[i & out$measures$Level == "R04"] <- .1
      out$measures$SE[i] <- .2
      out
    },
    resolve_dff_subgroup_precision = function(...) list(method = "MML",
      converged = TRUE, precision_tier = "model_based", supports_formal = TRUE),
    summarize_dff_group_linkage = function(...) list(status = "linked",
      ets_eligible = TRUE, anchored_levels = 6L, detail = "Test linking support."),
    .package = "mfrmr", .env = env
  )
}

test_that("failed or skipped DFF refits never borrow baseline estimates", {
  x <- dff_refit_fixture()
  calls <- 0L
  local_mocked_bindings(fit_mfrm = function(...) {
    calls <<- calls + 1L
    stop("forced subgroup failure")
  }, .package = "mfrmr")
  groups <- unique(x$data$Group)
  data <- x$data[x$data$Group != groups[1] |
    seq_len(nrow(x$data)) %in% which(x$data$Group == groups[1])[1:3], ]
  out <- analyze_dff(x$fit, x$diag, "Rater", "Group", data, method = "refit")
  expect_identical(calls, 1L)
  expect_equal(nrow(out$dif_table), 4L)
  expect_true(all(is.na(out$dif_table$Estimate1)))
  expect_true(all(is.na(out$dif_table$Estimate2)))
  expect_true(all(is.na(out$dif_table$Contrast)))
  expect_true(all(out$dif_table$Classification == "Unavailable"))
  expect_true(all(out$dif_table$ComparisonStatus == "unavailable"))
  expect_true(all(grepl("forced subgroup failure", out$dif_table$ComparisonReason)))
  expect_equal(out$summary$Count[out$summary$Classification == "Unavailable"], 4)
  expect_equal(dif_report(out)$counts[["Unavailable"]], 4)
  expect_match(dif_report(out)$narrative, "unavailable")
  expect_error(plot(out, draw = FALSE), "No finite contrast")
})

test_that("DFF keeps missing target cells and counts each target's retained rows", {
  x <- dff_refit_fixture()
  dff_refit_mock_success(x)
  a <- sort(unique(x$data$Group))[1]
  b <- sort(unique(x$data$Group))[2]
  x$data <- x$data[!(x$data$Group == b & x$data$Rater == "R03"), ]
  thin <- which(x$data$Group == a & x$data$Rater == "R02")
  x$data <- x$data[-thin[-(1:3)], ]
  out <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data,
    method = "refit", focal = b, min_obs = 10)
  rows <- out$dif_table
  expect_equal(nrow(rows), 4L)
  i <- match("R02", rows$Level)
  expect_equal(rows$N_Group1[i], 3L)
  expect_gt(rows$N_Subgroup1[i], rows$N_Group1[i])
  expect_true(rows$sparse[i])
  expect_equal(rows$Contrast[i], 1)
  expect_false(rows$ConditionalRefitScreenEligible[i])
  expect_match(rows$ConditionalScreenReason[i], "sparse")
  j <- match("R03", rows$Level)
  expect_equal(rows$N_Group2[j], 0L)
  expect_true(is.na(rows$Estimate2[j]))
  expect_true(is.na(rows$Contrast[j]))
  expect_match(rows$ComparisonReason[j], "absent")
  expect_identical(rows$Classification[j], "Unavailable")
  expect_equal(sum(out$summary$Count), nrow(rows))
  expect_identical(rows$ContrastDirection[i], "larger_location_in_group1")
  expect_identical(rows$ContrastDirection[rows$Level == "R04"], "equal_locations")
  pairs <- plot_data(plot(out, draw = FALSE), "pairs")
  expect_identical(pairs$ComparisonReason, rows$ComparisonReason)
  expect_identical(pairs$PlotStatus[j], "Contrast unavailable")
})

test_that("a fixed value in one group does not hide failure of the other refit", {
  x <- dff_refit_fixture()
  dff_refit_mock_success(x)
  original <- get("fit_mfrm", asNamespace("mfrmr"))
  local_mocked_bindings(fit_mfrm = function(data, ...) {
    if (data$Group[1] == sort(unique(x$data$Group))[2]) stop("second group failed")
    original(data, ...)
  }, .package = "mfrmr")
  out <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data, method = "refit")
  expect_true(out$dif_table$Fixed1[out$dif_table$Level == "R01"])
  expect_true(all(out$dif_table$Classification == "Unavailable"))
  expect_true(all(out$dif_table$ComparisonStatus == "unavailable"))
  expect_true(all(grepl("second group failed", out$dif_table$ComparisonReason)))
  expect_equal(dif_report(out)$counts[["Unavailable"]], 4L)
})

test_that("DFF screening requires complete precision and disjoint Persons", {
  x <- dff_refit_fixture()
  dff_refit_mock_success(x)
  ordinary <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data,
    method = "refit", min_obs = 1)
  i <- match("R02", ordinary$dif_table$Level)
  row <- ordinary$dif_table[i, ]
  expect_true(row$ConditionalRefitScreenEligible)
  expect_equal(row$SE, sqrt(.08))
  expect_equal(row$df, 2 * (row$N_Group1 - 1))
  expect_equal(row$p_value, 2 * pt(abs(row$Contrast) / sqrt(.08),
    df = row$df, lower.tail = FALSE))
  # This is a location contrast, not a whole response-curve direction.
  positive <- x$fit
  positive$config$positive_facets <- "Rater"
  positive$config$facet_signs[["Rater"]] <- 1
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) NULL,
    .package = "mfrmr")
  same <- analyze_dff(positive, x$diag, "Rater", "Group", x$data,
    method = "refit", min_obs = 1)
  expect_identical(same$dif_table$ContrastDirection, ordinary$dif_table$ContrastDirection)

  groups <- sort(unique(x$data$Group))
  p <- unique(x$data$Person[x$data$Group == groups[1]])[1]
  q <- unique(x$data$Person[x$data$Group == groups[2]])[1]
  x$data$Person[x$data$Person == q] <- p
  shared <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data,
    method = "refit", min_obs = 1)
  expect_true(all(shared$dif_table$NPersons_Overlap == 1L))
  expect_false(any(shared$dif_table$ConditionalRefitScreenEligible))
  expect_true(all(is.na(shared$dif_table$SE)))
  expect_match(shared$dif_table$ConditionalScreenReason[i], "shared Persons")
  expect_equal(shared$dif_table$Contrast, ordinary$dif_table$Contrast)
})

test_that("missing DFF SEs do not retain conditional-screen eligibility", {
  x <- dff_refit_fixture()
  dff_refit_mock_success(x)
  original <- get("diagnose_mfrm", asNamespace("mfrmr"))
  local_mocked_bindings(diagnose_mfrm = function(...) {
    out <- original(...)
    out$measures$SE[out$measures$Facet == "Rater"] <- NA_real_
    out
  }, .package = "mfrmr")
  out <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data, method = "refit")
  expect_false(any(out$dif_table$ConditionalRefitScreenEligible))
  expect_true(all(is.na(out$dif_table$SE)))
  expect_true(all(is.finite(out$dif_table$Contrast[out$dif_table$Level != "R01"])))
  expect_match(out$dif_table$ConditionalScreenReason[out$dif_table$Level == "R02"], "SE")
  old <- out
  old$config$refit_output_version <- NULL
  expect_error(summary(old), "Rerun.*subgroup")
  expect_error(dif_report(old), "Rerun.*subgroup")
  expect_error(plot(old, draw = FALSE), "Rerun.*subgroup")
})

test_that("DFF screens preserve representable SEs across extreme scales", {
  x <- dff_refit_fixture()
  dff_refit_mock_success(x)
  original <- get("diagnose_mfrm", asNamespace("mfrmr"))
  se_scale <- 1
  local_mocked_bindings(diagnose_mfrm = function(fit, ...) {
    out <- original(fit, ...)
    out$measures$SE[out$measures$Facet == "Rater"] <-
      se_scale * if (fit$.test_group == sort(unique(x$data$Group))[1]) .2 else .3
    out
  }, .package = "mfrmr")
  for (s in c(1e-200, 1e200)) {
    se_scale <- s
    out <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data,
      method = "refit", min_obs = 1)
    row <- out$dif_table[out$dif_table$Level == "R02", ]
    expect_true(row$ConditionalRefitScreenEligible)
    expect_equal(row$SE / s, sqrt(.13), tolerance = 1e-14)
    expected_df <- .13^2 / (.04^2 / (row$N_Group1 - 1) +
      .09^2 / (row$N_Group2 - 1))
    expect_equal(row$df, expected_df, tolerance = 1e-14)
    expect_true(is.finite(row$t))
    expect_true(is.finite(row$p_value))
  }
  se_scale <- 1e-309
  out <- analyze_dff(x$fit, x$diag, "Rater", "Group", x$data,
    method = "refit", min_obs = 1)
  row <- out$dif_table[out$dif_table$Level == "R02", ]
  expect_equal(row$Contrast, 1)
  expect_false(row$ConditionalRefitScreenEligible)
  expect_true(is.na(row$SE))
  expect_match(row$ConditionalScreenReason, "representable")
})
