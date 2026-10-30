response_testlet_fixture <- function(v = .64, ps = 1.2, missing = FALSE) {
  d <- expand.grid(Repeat = 1:2, Block = c("x", "y"), Person = c("A", "B", "C"))
  d$Score <- c(0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0)
  if (missing) { extra <- d[1, ]; extra$Person <- "missing"; extra$Score <- NA; d <- rbind(d, extra) }
  inp <- mfrm_testlet_data(d, "Person", "Score", "Block", character(), 0:1, "omit")
  structure(list(input = inp, parameters = c(.2, v, ps^2),
    calibration = list(steps = .2, variance = v, person_variance = ps^2),
    calibration_table = data.frame(Parameter = "Step", Facet = "Score", Level = "1", Estimate = .2),
    checks = data.frame(NumericalReady = TRUE, InformationPositive = TRUE),
    settings = list(quad_points = 61L, missing = "omit", method = "Fixed reference calibration")), class = "mfrm_testlet")
}

response_random_fixture <- function(sd = .7, ps = 1.2) {
  d <- expand.grid(Person = c("A", "B", "C"), Rater = c("R1", "R2"))
  d$Score <- c(0, 1, 2, 1, 0, 2)
  inp <- mfrm_random_rater_data(d, "Person", "Rater", character(), "Score", 0:2)
  structure(list(input = inp, calibration = list(beta = numeric(), steps = c(-.6, .6),
    rater_sd = sd, person_sd = ps, person_variance = ps^2),
    calibration_table = data.frame(Parameter = "Step", Facet = "Score", Level = c("1", "2"), Estimate = c(-.6, .6)),
    raters = data.frame(Rater = c("R1", "R2"), Estimate = 0),
    checks = data.frame(NumericalReady = TRUE, InformationPositive = TRUE),
    settings = list(quad_points = 61L, missing = "fail", method = "Fixed reference calibration")), class = "mfrm_random_rater")
}

test_that("nested testlet predictions agree with independent continuous integration", {
  f <- response_testlet_fixture()
  diag <- mfrm_response_diagnostics(f, rows = 1:4, group_by = "Person")
  # For A, x has 00 and y has 10; the replicate shares A's x effect.
  integral <- function(power, variance_part = FALSE) integrate(function(theta) vapply(theta, function(t) {
    inner <- function(g) plogis(t + .8 * g - .2)
    x <- integrate(function(g) {
      p <- inner(g)
      dnorm(g) * (1 - p)^2 * if (variance_part) p * (1 - p) else p^power
    }, -Inf, Inf, rel.tol = 1e-9)$value
    y <- integrate(function(g) { p <- inner(g); dnorm(g) * p * (1 - p) },
      -Inf, Inf, rel.tol = 1e-9)$value
    x * y * dnorm(t, sd = 1.2)
  }, numeric(1)), -Inf, Inf, rel.tol = 1e-8)$value
  mass <- integral(0); mean <- integral(1) / mass
  expect_equal(diag$probabilities[1, 2], mean, tolerance = 1e-7)
  within <- integral(0, TRUE) / mass; between <- integral(2) / mass - mean^2
  expect_gt(between, .01)
  expect_equal(diag$rows$PredictiveVariance[1], within + between, tolerance = 1e-7)
  expect_true(all(diag$rows$Status == "available_conditional"))
  expect_equal(diag$measures$Infit, sum(diag$rows$SquaredResidual) / sum(diag$rows$PredictiveVariance))
  expect_equal(diag$measures$Outfit, mean(diag$rows$SquaredResidual / diag$rows$PredictiveVariance))
  expect_false(any(c("Flag", "ZSTD", "p_value") %in% names(diag$measures)))
})

test_that("row selection retains evidence, category origins and missing or failed rows", {
  f <- response_testlet_fixture(missing = TRUE)
  all <- mfrm_response_diagnostics(f, group_by = "Person")
  selected <- mfrm_response_diagnostics(f, rows = c(3, 1, 13), group_by = "Person")
  expect_equal(selected$probabilities, all$probabilities[c(3, 1, 13), ])
  expect_equal(selected$rows$Status, c("available_conditional", "available_conditional", "missing_score"))
  expect_true(is.na(selected$measures$Infit[2]))
  selected_results <- mfrm_results(f, response_diagnostics = selected)
  overview <- selected_results$tables$response_overview
  expect_equal(unname(unlist(overview[1:7])), c(13, 3, 2, 0, 1, 0, 10))
  expect_identical(overview$Status, "review")
  expect_identical(selected_results$status$Status[selected_results$status$Section == "response_diagnostics"], "review")
  # Reoriginating the categories changes means, not residual variances.
  shifted <- f; shifted$input$score_levels <- 5:6
  shifted$input$data$Score <- shifted$input$data$Score + 5
  shifted$input$assigned_data$Score <- shifted$input$assigned_data$Score + 5
  moved <- mfrm_response_diagnostics(shifted)
  expect_equal(moved$rows$ExpectedScore, all$rows$ExpectedScore + 5)
  expect_equal(moved$rows$PredictiveVariance, all$rows$PredictiveVariance)
  expect_error(mfrm_response_diagnostics(f, rows = c(1, 1)), "distinct")
  expect_error(mfrm_response_diagnostics(f, group_by = "Score"), "identifier")
  expect_error(mfrm_response_diagnostics(f, quad_points = 242), "quad_points")
  broken <- f; broken$checks$NumericalReady <- FALSE
  expect_error(mfrm_response_diagnostics(broken), "checks")
  local_mocked_bindings(mfrm_testlet_response_probabilities = function(...) stop("unresolved reference"))
  bad <- mfrm_response_diagnostics(f, group_by = "Person")
  expect_true(all(bad$rows$Status[1:12] == "unavailable"))
  expect_match(bad$rows$Reason[1], "unresolved reference")
  expect_true(all(is.na(bad$measures$Infit)))
  result <- mfrm_results(f, response_diagnostics = bad)
  expect_identical(result$status$Status[result$status$Section == "response_diagnostics"], "not_available")
  expect_equal(result$tables$response_overview$Unresolved, 12L)
  expect_equal(result$tables$response_overview$Missing, 1L)
  expect_identical(mfrm_report(result)$tables$response_overview, result$tables$response_overview)
  expect_identical(mfrm_report(result)$first_screen$Status[4L], "unavailable")
})

test_that("zero variance, partial failures and unconverged quadrature cannot make valid summaries", {
  f <- response_testlet_fixture(missing = TRUE)
  local_mocked_bindings(mfrm_testlet_response_probabilities = function(input, par, rows, order) {
    p <- matrix(.5, length(rows), 2); p[1, ] <- c(1, 0)
    list(probabilities = p, normalization_error = rep(0, length(rows)))
  })
  d <- mfrm_response_diagnostics(f, rows = c(1, 2, 13), group_by = "Person")
  expect_identical(d$rows$Status, c("unavailable", "available_conditional", "missing_score"))
  expect_match(d$rows$Reason[1], "variance")
  expect_true(all(is.na(d$measures$Infit)))
  expect_identical(d$measures$Available, c(1L, 0L))
  local_mocked_bindings(mfrm_testlet_response_probabilities = function(input, par, rows, order) {
    p <- if (order == 61) c(.3, .7) else c(.4, .6)
    list(probabilities = matrix(rep(p, each = length(rows)), length(rows)), normalization_error = rep(0, length(rows)))
  })
  unresolved <- mfrm_response_diagnostics(f, rows = 1, group_by = "Person")
  expect_identical(unresolved$rows$Status, "unavailable")
  expect_match(unresolved$rows$Reason, "integration is unresolved")
})

test_that("zero effects reduce to independently integrated Rasch predictions", {
  f <- response_testlet_fixture(v = 0, ps = .7)
  diag <- mfrm_response_diagnostics(f, rows = 1)
  density <- function(t) { p <- plogis(t - .2); p * (1 - p)^3 * dnorm(t, sd = .7) }
  reference <- integrate(function(t) plogis(t - .2) * density(t), -Inf, Inf)$value /
    integrate(density, -Inf, Inf)$value
  expect_equal(diag$probabilities[1, 2], reference, tolerance = 1e-7)
  zero <- mfrm_response_diagnostics(response_testlet_fixture(v = 0, ps = 0))
  expect_equal(zero$rows$ExpectedScore, rep(plogis(-.2), 12))
  expect_true(all(zero$rows$Status == "available_conditional"))
  expect_identical(mfrm_response_diagnostic_overview(zero)$Status, "available")
  skip_if_not_installed("RTMB", "2.0")
  r <- response_random_fixture(sd = 0, ps = .7)
  calc <- mfrm_response_diagnostics(r, rows = 1)
  prob <- function(t) { w <- exp(cbind(0, t + .6, 2 * t)); w / rowSums(w) }
  den <- integrate(function(t) { p <- prob(t); p[, 1] * p[, 2] * dnorm(t, sd = .7) }, -15, 15)$value
  ref <- vapply(1:3, function(k) integrate(function(t) {
    p <- prob(t); p[, 1] * p[, 2] * p[, k] * dnorm(t, sd = .7)
  }, -15, 15)$value / den, numeric(1))
  expect_equal(as.numeric(calc$probabilities), ref, tolerance = 1e-7)
  expect_lt(abs(calc$rows$NormalizationError), 1e-8)
})

test_that("category-specific Laplace integrals include the hypothetical shared replicate", {
  skip_if_not_installed("RTMB", "2.0")
  f <- response_random_fixture(); input <- f$input
  # Independent likelihood for a fixed rater vector; no RTMB operations here.
  direct <- function(z, k) {
    ll <- vapply(1:3, function(p) {
      ids <- which(input$person == p)
      integrate(function(t) vapply(t, function(a) {
        probs <- lapply(ids, function(i) {
          eta <- a - .7 * z[input$rater[i]]
          w <- exp(c(0, eta + .6, 2 * eta)); w / sum(w)
        })
        prod(vapply(seq_along(ids), function(j) probs[[j]][input$y[ids[j]] + 1], numeric(1))) *
          if (p == 1) probs[[1]][k + 1] else 1
      }, numeric(1)) * dnorm(t, sd = 1.2), -12, 12, rel.tol = 1e-9)$value
    }, numeric(1))
    -sum(log(ll)) - sum(dnorm(z, log = TRUE))
  }
  obj <- mfrm_random_rater_objective(input, 61, fixed_sd = .7,
    fixed_person_sd = 1.2, replicate_row = 1, random = FALSE)
  for (k in 0:2) expect_equal(as.numeric(obj$fn(c(-.6, .6, .2, -.3, k, as.numeric(1:2 <= k)))),
    direct(c(.2, -.3), k), tolerance = 1e-8)
  laplace <- vapply(0:2, function(k) {
    opt <- optim(c(0, 0), direct, k = k, method = "BFGS", control = list(reltol = 1e-12))
    H <- optimHess(opt$par, function(z) direct(z, k))
    -opt$value - as.numeric(determinant(H, logarithm = TRUE)$modulus) / 2 + log(2 * pi)
  }, numeric(1))
  reference <- exp(laplace - max(laplace)); reference <- reference / sum(reference)
  diag <- mfrm_response_diagnostics(f, rows = 1, group_by = "Rater")
  expect_equal(as.numeric(diag$probabilities), reference, tolerance = 2e-6)
  changed <- f; changed$input$y[c(2, 3, 5, 6)] <- 2
  changed$input$data$Score <- changed$input$y; changed$input$assigned_data <- changed$input$data
  expect_gt(max(abs(mfrm_response_diagnostics(changed, rows = 1)$probabilities - diag$probabilities)), .001)
  expect_true(is.finite(diag$rows$NormalizationError))
  all <- mfrm_response_diagnostics(f, rows = c(1, 4), group_by = "Person")
  expect_equal(as.numeric(all$probabilities[1, ]), as.numeric(diag$probabilities))
  expect_equal(as.numeric(all$probabilities[2, ]),
    as.numeric(mfrm_response_diagnostics(f, rows = 4)$probabilities))
})

test_that("older saved response results restore consistent display metadata without recalculation", {
  fit <- response_testlet_fixture(missing = TRUE)
  d <- mfrm_response_diagnostics(fit, rows = c(1, 13), group_by = "Person")
  current <- mfrm_results(fit, response_diagnostics = d)
  old <- current
  old$tables$response_overview <- NULL
  old$table_index <- mfrm_results_table_index(old$tables)
  at <- old$status$Section == "response_diagnostics"
  old$status$Status[at] <- "available"
  old$status$Detail[at] <- "Saved diagnostics."
  old$tables$section_status <- old$status
  before <- serialize(old, NULL)
  local_mocked_bindings(mfrm_response_diagnostics = function(...) stop("no new diagnostics"),
    mfrm_testlet_response_probabilities = function(...) stop("no integration"),
    fit_mfrm_testlet = function(...) stop("no fitting"), .package = "mfrmr")
  expect_identical(summary(old)$status, current$status)
  report <- mfrm_report(old)
  expect_identical(report$tables$response_overview, current$tables$response_overview)
  expect_identical(report$tables$section_status, current$status)
  expect_identical(report$source$fit, old$fit)
  expect_identical(report$source$diagnostics, old$diagnostics)
  expect_identical(mfrm_results_restore_response_overview(report$source), report$source)
  html <- mfrm_results_html(old)
  expect_identical(html$summary$status, current$status)
  expect_match(html$html, "response_overview")
  unlink(html$path)
  folder <- withr::local_tempdir()
  exported <- export_mfrm_results(old, output_dir = folder,
    include = c("summary", "tables", "report", "replay"), acknowledge_sensitive = TRUE)
  for (key in c("summary_status", "table_section_status")) {
    path <- exported$written_files$Path[exported$written_files$Component == key]
    expect_length(path, 1L)
    expect_identical(read.csv(path)$Status, current$status$Status)
  }
  path <- exported$written_files$Path[exported$written_files$Component == "results_rds"]
  restored <- readRDS(path)
  expect_identical(restored$tables$response_overview, current$tables$response_overview)
  expect_identical(restored$diagnostics, d)
  replay <- exported$written_files$Path[exported$written_files$Component == "replay_code"]
  env <- new.env(parent = globalenv())
  withr::with_dir(folder, sys.source(replay, envir = env))
  expect_identical(env$res, restored)
  expect_identical(serialize(old, NULL), before)
  no_diagnostics <- mfrm_results(fit)
  expect_identical(mfrm_results_restore_response_overview(no_diagnostics), no_diagnostics)
})

test_that("stored diagnostic plots and reports preserve limits and exact source identity", {
  f <- response_testlet_fixture(missing = TRUE)
  d <- mfrm_response_diagnostics(f, group_by = "Person")
  res <- mfrm_results(f, diagnostics = d)
  expect_equal(res$tables$response_measures, d$measures)
  expect_equal(res$tables$response_residuals, d$rows)
  expect_equal(res$tables$response_probabilities$Probability, as.vector(d$probabilities))
  expect_match(mfrm_report(res)$markdown, "no calibrated fit cutoffs")
  expect_identical(mfrm_report(res)$first_screen$Status[4L], "review")
  different <- f; different$input$assigned_data$Score[1] <- 1
  expect_error(mfrm_results(different, diagnostics = d), "exact source roster")
  different <- f; different$input$data$Score[1] <- 1
  expect_error(mfrm_results(different, diagnostics = d), "exact source roster")
  expect_error(mfrm_results(f, diagnostics = list()), "mfrm_response_diagnostics")
  path <- tempfile(); saveRDS(res, path); expect_identical(readRDS(path), res)
  for (style in c("paired", "scatter")) {
    payload <- plot(res, type = "response_diagnostics", style = style, draw = FALSE,
      show_title = FALSE, show_notes = FALSE, palette = "mono")
    expect_equal(plot_data(payload)$table, d$measures)
    expect_match(plot_data(payload)$alt_text, "calibration fixed")
    expect_false(plot_data(payload)$display$show_title)
    before <- dev.cur(); plot(d, draw = FALSE); expect_identical(dev.cur(), before)
    grDevices::pdf(tempfile(fileext = ".pdf")); expect_silent(plot(d, style = style)); grDevices::dev.off()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- as_ggplot(payload); expect_s3_class(p, "ggplot")
      expect_silent(ggplot2::ggplotGrob(p))
      expect_true(is.null(p$labels$title)); expect_true(is.null(p$labels$caption))
      expect_match(attr(p, "mfrmr_alt_text"), "no expectation-one")
    }
  }
  expect_error(plot(d, cutoff = 1.5), "empty")
  local_mocked_bindings(mfrm_testlet_response_probabilities = function(...) stop("must not recompute"))
  expect_s3_class(mfrm_report(mfrm_results(f, diagnostics = d)), "mfrm_report")
})
