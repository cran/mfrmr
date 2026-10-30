test_that("G-study omission is explicit and D-study retains its source-row accounting", {
  skip_if_not_installed("lme4")
  fit <- make_toy_fit(maxit = 150)
  data <- as.data.frame(fit$prep$data)
  baseline <- mfrm_generalizability(fit, data = data)
  variance <- setNames(baseline$variance_components$Variance, baseline$variance_components$Source)
  expect_equal(baseline$coefficients$G, unname(variance["Person"] /
    (variance["Person"] + variance["Residual"])), tolerance = 1e-14)
  expect_equal(baseline$coefficients$Phi, unname(variance["Person"] / sum(variance)), tolerance = 1e-14)
  n <- nrow(data)
  extra <- data[1:3, ]
  extra$Score[c(1, 3)] <- NA_real_
  extra$Rater[c(2, 3)] <- NA
  input <- rbind(data, extra)
  input$Unselected <- NA_real_
  expect_error(mfrm_generalizability(fit, data = input), "3 row.*missing")
  out <- mfrm_generalizability(fit, data = input, missing = "omit")
  expect_identical(out$variance_components, baseline$variance_components)
  expect_identical(out$coefficients[c("G", "Phi")], baseline$coefficients[c("G", "Phi")])
  expect_identical(out$data_usage$source, "Supplied data")
  expect_identical(out$data_usage$missing, "omit")
  expect_equal(out$data_usage$counts, c(InputRows = n + 3L, UsedRows = n, ExcludedRows = 3L))
  expect_equal(out$data_usage$excluded_rows, n + 1:3)
  expect_equal(out$data_usage$missing_cells,
    data.frame(InputRow = n + c(2L, 3L, 1L, 3L), Column = c("Rater", "Rater", "Score", "Score")))
  expect_equal(unname(unlist(out$coefficients[c("InputRows", "UsedRows", "ExcludedRows")])),
               c(n + 3L, n, 3L))
  expect_output(print(out), "3 excluded")
  expect_output(print(out), "no missing values were imputed")

  ds <- mfrm_d_study(out, data.frame(Rater = c(2, 4), Criterion = 3))
  expect_equal(ds$InputRows, rep(n + 3L, 2))
  expect_equal(ds$UsedRows, rep(n, 2))
  expect_equal(ds$ExcludedRows, c(3L, 3L))
  expect_identical(attr(ds, "data_usage"), out$data_usage)
  expect_identical(attr(ds[1, c("G", "Phi"), drop = FALSE], "data_usage"), out$data_usage)
  expect_output(print(ds[1, c("G", "Phi"), drop = FALSE]), "3 excluded")
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv))
  utils::write.csv(ds, csv, row.names = FALSE)
  expect_equal(utils::read.csv(csv)$ExcludedRows, ds$ExcludedRows)
  expect_equal(utils::read.csv(csv)$GStudyDataSource, rep("Supplied data", 2))
  direct <- mfrm_d_study(fit, data = input, missing = "omit")
  expect_identical(attr(direct, "data_usage"), out$data_usage)
  expect_error(mfrm_d_study(fit, data = input), "missing.*omit")

  stored <- mfrm_generalizability(fit)
  expect_identical(stored$data_usage$source, "Stored fitted rows")
  expect_equal(stored$data_usage$counts, c(InputRows = n, UsedRows = n, ExcludedRows = 0L))
  expect_equal(nrow(stored$data_usage$missing_cells), 0L)
  expect_output(print(stored), "earlier MFRM filtering is not included")
  legacy <- baseline
  legacy$data_usage <- NULL
  legacy$coefficients[c("InputRows", "UsedRows", "ExcludedRows", "GStudyDataSource")] <- NULL
  legacy_ds <- mfrm_d_study(legacy)
  expect_true(all(is.na(legacy_ds[c("InputRows", "UsedRows", "ExcludedRows")])))
  expect_true(is.na(legacy_ds$GStudyDataSource))
  expect_output(print(legacy_ds), "accounting is unavailable")
  expect_output(print(legacy), "accounting is unavailable")
})

test_that("G-study refuses malformed scores instead of silently discarding them", {
  skip_if_not_installed("lme4")
  fit <- make_toy_fit(maxit = 150)
  data <- as.data.frame(fit$prep$data)
  invalid <- data
  invalid$Score <- as.character(invalid$Score)
  invalid$Score[1] <- "not recorded"
  expect_error(mfrm_generalizability(fit, invalid), "nonnumeric or infinite")
  expect_error(mfrm_generalizability(fit, invalid, missing = "omit"), "nonnumeric or infinite")
  invalid$Score[1] <- " "
  expect_error(mfrm_generalizability(fit, invalid, missing = "omit"), "use NA")
  invalid <- data
  invalid$Score[1] <- Inf
  expect_error(mfrm_generalizability(fit, invalid, missing = "omit"), "infinite")
  invalid$Score <- as.Date("2026-01-01") + seq_len(nrow(data))
  expect_error(mfrm_generalizability(fit, invalid), "must be numeric")
  invalid$Score <- rep(NA_real_, nrow(data))
  expect_error(mfrm_generalizability(fit, invalid, missing = "omit"), "No complete rows")
  invalid <- data
  invalid$Rater <- as.character(invalid$Rater)
  invalid$Rater[1] <- " "
  expect_error(mfrm_generalizability(fit, invalid, missing = "omit"), "blank labels")
  expect_error(mfrm_generalizability(fit, reml = "TRUE"), "TRUE or FALSE")
  expect_error(mfrm_generalizability(fit, data, object_facet = "Score"), "reserved")
  invalid <- data
  invalid$Rater <- as.numeric(invalid$Rater)
  invalid$Rater[1] <- NaN
  expect_error(mfrm_generalizability(fit, invalid), "1 row.*missing")

  # Numeric labels retain their values rather than becoming factor level codes.
  baseline <- mfrm_generalizability(fit, data)
  data$Score <- factor(data$Score, levels = rev(sort(unique(data$Score))))
  labelled <- mfrm_generalizability(fit, data)
  expect_identical(labelled$variance_components, baseline$variance_components)
  # Literal facet names must not be interpreted as formula expressions.
  names(data)[names(data) == "Rater"] <- "Rater group"
  named <- mfrm_generalizability(fit, data, random_facets = c("Rater group", "Criterion"))
  expect_equal(named$coefficients[c("G", "Phi")], baseline$coefficients[c("G", "Phi")])
  expect_true("Rater group" %in% named$variance_components$Source)
  named_d <- mfrm_d_study(named)
  expect_true("n_Rater group" %in% names(named_d))
  expect_equal(named_d$Phi, mfrm_d_study(baseline)$Phi)
  expect_equal(mfrm_d_study(named, list(`n_Rater group` = 2, Criterion = 4))$Phi,
    mfrm_d_study(baseline, list(Rater = 2, n_Criterion = 4))$Phi)
})

test_that("D-study count aliases are resolved jointly without reusing columns", {
  gt <- structure(list(
    variance_components = data.frame(Source = c("Person", "Rater", "n_Rater", "Residual"),
      Variance = c(1, .4, .2, .8)),
    design = list(object_facet = "Person", random_facets = c("Rater", "n_Rater"),
      observed_levels = c(Rater = 2L, n_Rater = 4L),
      calculation_version = 2L, identification_status = "identified")),
    class = c("mfrm_generalizability", "list"))
  raw <- mfrm_d_study(gt, data.frame(Rater = 2, n_Rater = 4))
  prefixed <- mfrm_d_study(gt, data.frame(n_Rater = 2, n_n_Rater = 4))
  mixed <- mfrm_d_study(gt, list(n_n_Rater = 4, Rater = 2))
  expect_identical(raw, prefixed)
  expect_identical(raw, mixed)
  expect_identical(raw, mfrm_d_study(gt))
  expect_equal(raw$RelativeErrorVariance, .8 / (2 * 4))
  expect_equal(raw$AbsoluteErrorVariance, .4 / 2 + .2 / 4 + .8 / (2 * 4))
  expect_equal(raw$G, 1 / (1 + .1), tolerance = 1e-14)
  expect_equal(raw$Phi, 1 / (1 + .35), tolerance = 1e-14)
  expect_error(mfrm_d_study(gt, data.frame(n_Rater = 2)), "exactly one")
  expect_error(mfrm_d_study(gt, data.frame(Rater = 2, n_Rater = 4, n_n_Rater = 8)), "exactly one")
  expect_error(mfrm_d_study(gt, data.frame(n_Rater = 2, Unused = 4)), "one-to-one")
  expect_error(mfrm_d_study(gt,
    data.frame(Rater = 2, Rater = 3, n_Rater = 4, check.names = FALSE)), "unique")
  expect_error(mfrm_d_study(gt, list(Rater = 2)), "missing count")
  # A longer alias chain and mixed naming must also preserve each count.
  gt$design$random_facets <- c("n_n_Rater", "Rater", "n_Rater")
  gt$variance_components <- rbind(gt$variance_components,
    data.frame(Source = "n_n_Rater", Variance = .3))
  chained <- mfrm_d_study(gt, list(n_Rater = 2, n_n_Rater = 4, n_n_n_Rater = 5))
  expect_equal(unname(unlist(chained[c("n_Rater", "n_n_Rater", "n_n_n_Rater")])), c(2, 4, 5))
  expect_equal(chained$RelativeErrorVariance, .8 / (2 * 4 * 5))
})

test_that("D-study coefficient values agree with their reference statuses before rounding", {
  target <- .79996
  gt <- structure(list(
    variance_components = data.frame(Source = c("Person", "Rater", "Residual"),
      Variance = c(1, .1, 1 / target - 1)),
    design = list(object_facet = "Person", random_facets = "Rater",
      calculation_version = 2L, identification_status = "identified")),
    class = c("mfrm_generalizability", "list"))
  d <- mfrm_d_study(gt, list(Rater = 1))
  expect_equal(d$G, target, tolerance = 1e-14)
  expect_lt(d$G, .8)
  expect_identical(d$GStatus, "at_or_above_0.70_reference")
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  write.csv(d, path, row.names = FALSE)
  expect_lt(read.csv(path)$G, .8)
  expect_equal(plot_data(plot(d, draw = FALSE))$table$G, target, tolerance = 1e-14)
})

test_that("main-effects boundary points and errors have independent availability", {
  g <- structure(list(variance_components = data.frame(
    Source = c("Person", "Rater", "Residual"), Variance = c(0, .4, .8)),
    design = list(object_facet = "Person", random_facets = "Rater",
      calculation_version = 2L, identification_status = "boundary_or_singular_fit")),
    class = c("mfrm_generalizability", "list"))
  d <- mfrm_d_study(g, list(Rater = c(1, 2)), residual_scaling = "sensitivity")
  expect_equal(d$G, rep(0, 6))
  expect_equal(d$Phi, rep(0, 6))
  expect_equal(d$RelativeErrorVariance, c(.8, .8, .8, .4, .4, .8))
  expect_equal(d$AbsoluteErrorVariance, c(1.2, 1.2, 1.2, .6, .6, 1))
  expect_true(all(d$GStatus == "identification_warning"))
  expect_true(all(d$GReason == "" & d$PhiReason == ""))
  g$variance_components$Variance[3] <- 0
  d <- mfrm_d_study(g, list(Rater = 2))
  expect_true(is.na(d$G))
  expect_identical(d$GReason, "zero_total_variance")
  expect_equal(d$Phi, 0)
  expect_equal(d$RelativeErrorVariance, 0)
  g$variance_components$Variance[1] <- NA_real_
  d <- mfrm_d_study(g, list(Rater = 2))
  expect_true(is.na(d$G) && is.na(d$Phi))
  expect_identical(d$PhiReason, "object_variance_missing_or_invalid")
  expect_equal(d$AbsoluteErrorVariance, .2)
  expect_equal(d$ResidualDivisor, 2)
})

test_that("changing the object facet retains persons among the default random conditions", {
  skip_if_not_installed("lme4")
  fit <- make_toy_fit(maxit = 150)
  g <- mfrm_generalizability(fit, object_facet = "Rater")
  explicit <- mfrm_generalizability(fit, object_facet = "Rater",
    random_facets = c("Person", "Criterion"))
  expect_identical(g, explicit)
  expect_identical(g$design$random_facets, c("Person", "Criterion"))
  expect_identical(g$design$object_facet, "Rater")
  d <- mfrm_d_study(g, list(Person = 10, Criterion = 3))
  expect_equal(d$ObjectVariance,
    g$variance_components$Variance[g$variance_components$Source == "Rater"])
  expect_true("n_Person" %in% names(d))
})

test_that("main-effects projections avoid overflow and disclose unrepresentable quantities", {
  g <- structure(list(variance_components = data.frame(
    Source = c("Person", "Rater", "Task", "Residual"), Variance = c(1e308, 1e308, 0, 1e308)),
    design = list(object_facet = "Person", random_facets = c("Rater", "Task"),
      calculation_version = 2L, identification_status = "identified")),
    class = c("mfrm_generalizability", "list"))
  d <- mfrm_d_study(g, list(Rater = 1, Task = 1))
  expect_equal(d$G, .5)
  expect_equal(d$Phi, 1/3)
  expect_true(is.na(d$AbsoluteErrorVariance))
  expect_identical(d$AbsoluteErrorReason, "absolute_error_exceeds_numerical_range")
  expect_identical(d$PhiReason, "")
  g$variance_components$Variance <- c(1, 0, 0, 1e308)
  huge <- mfrm_d_study(g, list(Rater = 1e200, Task = 1e200))
  expect_equal(huge$RelativeErrorVariance / 1e-92, 1, tolerance = 1e-14)
  expect_true(is.na(huge$ResidualDivisor))
  expect_identical(huge$ResidualDivisorReason, "divisor_exceeds_numerical_range")
  g$variance_components$Variance <- c(1, 0, 0, 1e-300)
  tiny <- mfrm_d_study(g, list(Rater = 1e200, Task = 1e200))
  expect_true(is.na(tiny$G))
  expect_identical(tiny$GReason, "relative_error_below_numerical_range")
  expect_true(is.na(tiny$RelativeErrorVariance))
  g$variance_components$Variance <- c(1e-300, 0, 0, 1e308)
  underflow <- mfrm_d_study(g, list(Rater = 1, Task = 1))
  expect_true(is.na(underflow$G) && is.na(underflow$Phi))
  expect_identical(underflow$GReason, "coefficient_below_numerical_range")
  expect_equal(underflow$RelativeErrorVariance, 1e308)
})

test_that("main-effects kernel information matches independent dense covariance derivatives", {
  d <- expand.grid(Person = 1:5, Rater = 1:3, Task = 1:2)
  d <- d[-c(2, 9, 21), ]
  terms <- names(d)
  n <- nrow(d)
  kernels <- c(lapply(d, function(x) outer(x, x, `==`)), list(Residual = diag(n)))
  for (reml in c(FALSE, TRUE)) {
    h <- if (reml) diag(n)-matrix(1/n, n, n) else diag(n)
    projected <- lapply(kernels, function(k) h %*% k %*% h)
    reference <- outer(seq_along(kernels), seq_along(kernels),
      Vectorize(function(i, j) sum(projected[[i]]*projected[[j]])))
    dimnames(reference) <- list(names(kernels), names(kernels))
    info <- mfrmr:::.mfrmr_gt_geometry(d, terms, reml)
    expect_equal(info$kernel_gram, reference, tolerance = 1e-12)
    expect_equal(info$rank, 4)
    expect_equal(ncol(info$null_space), 0L)
  }
})

test_that("confounded main effects retain identified sums but refuse arbitrary allocations", {
  skip_if_not_installed("lme4")
  set.seed(20261020)
  data <- expand.grid(Person = paste0("P", 1:16), Rater = paste0("R", 1:4), Repeat = 1:2)
  data$Task <- paste0("T", as.integer(factor(data$Rater)))
  data$Score <- rnorm(16)[match(data$Person, unique(data$Person))] +
    2*rnorm(4)[match(data$Rater, unique(data$Rater))] + rnorm(nrow(data))
  fit <- structure(list(config = list(facet_names = c("Rater", "Task")),
    prep = list(data = data)), class = "mfrm_fit")
  for (reml in c(FALSE, TRUE)) {
    g <- mfrm_generalizability(fit, reml = reml)
    expect_identical(g$design$identification_status, "components_not_identified")
    expect_identical(g$design$component_identification$rank, 3L)
    expect_identical(g$variance_components$Identified,
      !g$variance_components$Source %in% c("Rater", "Task"))
    grid <- data.frame(Rater = c(1, 2, 1, 2), Task = c(1, 1, 2, 2))
    d <- mfrm_d_study(g, grid)
    expect_true(all(is.finite(d$G)))
    expect_true(all(is.finite(d$Phi[c(1, 4)])))
    expect_true(all(is.na(d$Phi[c(2, 3)])))
    expect_equal(d$PhiReason, c("", "absolute_error_not_identified",
      "absolute_error_not_identified", ""))
    expect_equal(g$coefficients$Phi, d$Phi[1])
    # Move along the likelihood's null direction. Identified projections
    # must be invariant to this arbitrary split of the same total variance.
    other <- g
    at <- match(c("Rater", "Task"), other$variance_components$Source)
    total <- sum(other$variance_components$Variance[at])
    other$variance_components$Variance[at] <- c(.25, .75)*total
    expect_equal(mfrm_d_study(other, grid)$Phi, d$Phi)
    expect_equal(mfrm_d_study(other, grid)$G, d$G)
    before <- serialize(g, NULL)
    p <- plot(d, draw = FALSE)
    expect_identical(plot_data(p)$table$Phi, d$Phi)
    expect_identical(serialize(g, NULL), before)
    file <- tempfile(fileext = ".rds")
    saveRDS(g, file)
    expect_identical(mfrm_d_study(readRDS(file), grid), d)
    unlink(file)
  }
})
