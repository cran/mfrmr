mvdc_fixture <- function() {
  fixture <- new.env()
  eval(parse(testthat::test_path("test-multivariate-gtheory.R"))[[1L]], fixture)
  fixture$mvgt_fixture()$data
}

test_that("paired intervals match independent quadratic forms and numerical gradients", {
  full <- mvdc_fixture()
  grid <- data.frame(Raters = c(2, 3, 4), Tasks = c(6, 4, 3))
  for (incomplete in c(FALSE, TRUE)) {
    data <- if (incomplete) full[-c(1, 8, 23), ] else full
    data[1:3] <- lapply(data[1:3], as.character)
    g <- mfrm_g_study(data, c("Content", "Organization"),
      method = if (incomplete) "minque0" else "anova")
    d <- mfrm_d_study(g, grid)
    result <- mfrm_multivariate_d_compare(d, score = "Content", assumption = "normal")
    metrics <- c("G", "Phi", "RelativeSEM", "AbsoluteSEM")
    N <- nrow(data); H <- diag(N) - 1/N
    subsets <- list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
    K <- lapply(subsets, function(s) Reduce(`*`, lapply(data[s], function(x) outer(x, x, `==`) * 1)))
    B <- lapply(K, function(k) H %*% k %*% H)
    inverse <- solve(crossprod(vapply(B, as.vector, numeric(N^2))))
    A <- lapply(1:7, function(i) Reduce(`+`, Map(`*`, B, inverse[i, ])))
    theta <- vapply(g$components, function(a) a[1, 1], numeric(1))
    V <- Reduce(`+`, Map(`*`, K, theta))
    AV <- lapply(A, function(a) a %*% V)
    C <- outer(1:7, 1:7, Vectorize(function(i, j) 2 * sum(AV[[i]] * t(AV[[j]]))))
    contrast <- function(g) {
      tab <- mfrm_d_study(g, grid)$coefficients
      values <- as.matrix(tab[tab$Score == "Content", metrics])
      as.vector(t(sweep(values[-1, ], 2, values[1, ], "-")))
    }
    J <- vapply(1:7, function(j) {
      plus <- minus <- g
      plus$components[[j]][1, 1] <- plus$components[[j]][1, 1] + 1e-6
      minus$components[[j]][1, 1] <- minus$components[[j]][1, 1] - 1e-6
      (contrast(plus) - contrast(minus)) / 2e-6
    }, numeric(8))
    expect_equal(result$comparisons$Difference, contrast(g), tolerance = 1e-12)
    expect_equal(unname(result$covariance), J %*% C %*% t(J), tolerance = 1e-8)
    expect_equal(result$comparisons$SE^2, unname(diag(result$covariance)), tolerance = 1e-12)
    expect_equal(result$comparisons$Upper - result$comparisons$Lower,
      2 * qnorm(.975) * result$comparisons$SE)
    if (!incomplete) {
      minque <- mfrm_g_study(data, c("Content", "Organization"), method = "minque0")
      alternate <- mfrm_multivariate_d_compare(mfrm_d_study(minque, grid),
        score = "Content", assumption = "normal")
      expect_equal(alternate$comparisons, result$comparisons, tolerance = 1e-10)
    }
    reverse <- mfrm_multivariate_d_compare(d, reference = 2, score = "Content", assumption = "normal")
    expect_equal(reverse$comparisons$Difference[1:4], -result$comparisons$Difference[1:4])
    expect_equal(reverse$comparisons$Lower[1:4], -result$comparisons$Upper[1:4])
    expect_equal(summary(result), result$comparisons)
    expect_output(print(result), "normal random effects required")
  }
})

test_that("composite scaling, facet names and replay preserve paired intervals", {
  data <- mvdc_fixture()
  g <- mfrm_g_study(data, c("Content", "Organization"),
    facets = c(Assessor = "Rater", Occasion = "Task"))
  grid <- data.frame(Assessor = c(2, 3), Occasion = c(6, 4))
  w <- c(Content = .5, Organization = .5)
  d <- mfrm_d_study(g, grid, cbind(Equal = w, Scaled = -3 * w))
  a <- mfrm_multivariate_d_compare(d, composite = "Equal", assumption = "normal")
  b <- mfrm_multivariate_d_compare(d, composite = "Scaled", assumption = "normal")
  for (column in c("Difference", "SE", "Lower", "Upper")) {
    expect_equal(b$comparisons[[column]], a$comparisons[[column]] * c(1, 1, 3, 3))
  }
  expect_equal(b$covariance, a$covariance * outer(c(1, 1, 3, 3), c(1, 1, 3, 3)))
  data$Total <- rowSums(data[c("Content", "Organization")]) / 2
  scalar <- mfrm_multivariate_d_compare(mfrm_d_study(
    mfrm_g_study(data, "Total", facets = c(Assessor = "Rater", Occasion = "Task")), grid),
    assumption = "normal")
  expect_equal(scalar$comparisons$SE, a$comparisons$SE)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(a, path)
  expect_identical(summary(readRDS(path)), summary(a))
  expect_identical(plot_data(plot(readRDS(path), draw = FALSE)), plot_data(plot(a, draw = FALSE)))
  expect_error(mfrm_multivariate_d_compare(d, assumption = "normal"), "Select one")
  expect_error(mfrm_multivariate_d_compare(d, composite = "Missing", assumption = "normal"), "present")
})

test_that("boundary and identical-plan outcomes remain explicit", {
  g <- mfrm_g_study(mvdc_fixture(), c("Content", "Organization"))
  g$components$Person[,] <- 0
  d <- mfrm_d_study(g, data.frame(Raters = c(2, 2, 3), Tasks = c(6, 6, 4)))
  result <- mfrm_multivariate_d_compare(d, score = "Content", assumption = "normal")
  expect_equal(result$comparisons$Difference[1:4], rep(0, 4))
  expect_equal(result$comparisons$SE[1:4], rep(0, 4))
  expect_true(all(is.na(result$comparisons$SE[5:6])))
  expect_true(all(grepl("boundary", result$comparisons$Status[5:6])))
  expect_true(all(is.finite(result$comparisons$SE[7:8])))
  plans <- plot_data(plot(result, draw = FALSE))
  expect_identical(plans$view, "plans")
  expect_equal(plans$series$Value, rep(0, 6))
  expect_equal(nrow(plans$unavailable), 0L)
  expect_equal(nrow(plans$interval_unavailable), 2L)
  payload <- plot_data(plot(result, view = "differences", draw = FALSE))
  expect_equal(nrow(payload$unavailable), 2L)
  expect_true(all(payload$unavailable$Difference == 0))
  expect_true(all(is.na(payload$unavailable$Lower)))
  g$components$Person[,] <- -1
  missing <- mfrm_multivariate_d_compare(mfrm_d_study(g, d$design_grid),
    score = "Content", assumption = "normal")
  expect_true(all(is.na(missing$comparisons$Difference[missing$comparisons$Metric %in% c("G", "Phi")])))
  expect_true(all(grepl("Point projection", missing$comparisons$Status[1:2])))
  unavailable <- plot_data(plot(missing, draw = FALSE))
  expect_true(all(is.na(unavailable$series$Value)))
  expect_equal(nrow(unavailable$unavailable), 6L)
})

test_that("saved comparisons retain source row accounting and unnormalized weights", {
  data <- mvdc_fixture()
  data$Content[1] <- NA_real_
  g <- mfrm_g_study(data, c("Content", "Organization"),
    facets = c(Assessor = "Rater", Occasion = "Task"), method = "minque0", missing = "omit")
  d <- mfrm_d_study(g, data.frame(Assessor = c(2, 3), Occasion = c(6, 4)),
    c(Content = 2, Organization = -1))
  before <- serialize(d, NULL)
  a <- mfrm_multivariate_d_compare(d, assumption = "normal")
  expect_identical(a$source_design$counts, c(Person = 4L, Assessor = 3L, Occasion = 4L))
  expect_identical(a$source_rows, c(InputRows = 48L, UsedRows = 47L, ExcludedRows = 1L))
  expect_identical(a$source_design$complete, FALSE)
  expect_null(a$source_design$levels)
  expect_null(a$data)
  output <- capture.output(print(a))
  expect_true(all(capture.output(print(c(Content = 2, Organization = -1))) %in% output))
  expect_match(paste(output, collapse = "\n"), "Source rows: 47 used of 48 ; 1 explicitly omitted.", fixed = TRUE)
  expect_match(paste(output, collapse = "\n"), "future counts do not reproduce its observed assignments", fixed = TRUE)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(a, path)
  expect_identical(capture.output(print(readRDS(path))), output)
  expect_identical(summary(readRDS(path)), a$comparisons)
  old <- a; old[c("source_design", "source_rows")] <- NULL
  expect_output(print(old), "Source design not recorded")
  expect_identical(summary(old), summary(a))
  expect_identical(plot_data(plot(old, draw = FALSE)), plot_data(plot(a, draw = FALSE)))
  expect_identical(serialize(d, NULL), before)
})

test_that("unsupported assumptions and malformed comparison requests fail clearly", {
  g <- mfrm_g_study(mvdc_fixture(), "Content")
  d <- mfrm_d_study(g, data.frame(Raters = c(2, 3), Tasks = c(6, 4)))
  expect_error(mfrm_multivariate_d_compare(d), "assumption")
  expect_error(mfrm_multivariate_d_compare(d, assumption = "robust"), "assumption")
  for (reference in list(0, 3, NA_real_, 1.1, 1+1i, matrix(1))) {
    expect_error(mfrm_multivariate_d_compare(d, reference = reference, assumption = "normal"), "reference")
  }
  for (level in list(0, 1, NA_real_, c(.8, .9), .9+1i, matrix(.9))) {
    expect_error(mfrm_multivariate_d_compare(d, level = level, assumption = "normal"), "level")
  }
  expect_error(mfrm_multivariate_d_compare(d, score = "Content", composite = "Equal",
    assumption = "normal"), "either")
  one <- mvdc_fixture(); one <- one[one$Rater == one$Rater[1], ]
  d1 <- mfrm_d_study(mfrm_g_study(one, "Content", rater = NULL),
    data.frame(Tasks = c(2, 4)))
  expect_error(mfrm_multivariate_d_compare(d1, assumption = "normal"), "two common")
  d$gstudy$data <- NULL
  expect_error(mfrm_multivariate_d_compare(d, assumption = "normal"), "retain")
})

test_that("comparison plots retain interval endpoints and reject generic conversion", {
  g <- mfrm_g_study(mvdc_fixture(), "Content")
  d <- mfrm_d_study(g, data.frame(Raters = c(2, 3, 4), Tasks = c(6, 4, 3)))
  result <- mfrm_multivariate_d_compare(d, reference = 2, assumption = "normal")
  before <- serialize(result, NULL)
  for (type in c("coefficients", "sem")) {
    payload <- plot_data(plot(result, type = type, draw = FALSE))
    wanted <- result$comparisons[result$comparisons$Metric %in% payload$metrics, ]
    expect_equal(payload$table[names(wanted)], wanted)
    expect_equal(payload$series$Value, unlist(d$coefficients[payload$metrics], use.names = FALSE))
    expect_equal(payload$series$Scenario[payload$series$IsReference], c(2L, 2L))
    expect_false(any(c("Lower", "Upper", "SE") %in% names(payload$series)))
    expect_match(payload$title, "Raters = 3, Tasks = 4", fixed = TRUE)
    difference <- plot_data(plot(result, type = type, view = "differences", draw = FALSE))
    expect_equal(difference$table[names(wanted)], wanted)
    expect_match(difference$subtitle, "Normal random effects")
  }
  expect_identical(serialize(result, NULL), before)
  expect_error(plot(result, view = "invalid", draw = FALSE), "arg")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_error(as_ggplot(result), "plot_data")
    expect_error(as_ggplot(plot(result, draw = FALSE)), "plot_data")
  }
})

test_that("comparison plots draw and replay on a standard landscape device", {
  g <- mfrm_g_study(mvdc_fixture(), c("Content", "Organization"))
  d <- mfrm_d_study(g,
    data.frame(Raters = c(2, 3, 4), Tasks = c(6, 4, 3)),
    weights = c(Content = .5, Organization = .5))
  result <- mfrm_multivariate_d_compare(d, assumption = "normal")
  on_device <- function(width, height, code) {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path, width = width, height = height)
    on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
    grDevices::dev.control(displaylist = "enable")
    force(code)
  }
  for (view in c("plans", "differences")) for (type in c("coefficients", "sem")) {
    recorded <- on_device(7, 7, {
      old <- graphics::par(c("mfrow", "mar", "oma", "cex", "mex"))
      plot(result, view = view, type = type)
      expect_identical(graphics::par(names(old)), old)
      grDevices::recordPlot()
    })
    on_device(700 / 96, 433 / 96, {
      expect_silent(grDevices::replayPlot(recorded))
      expect_silent(plot(result, view = view, type = type))
    })
  }
})
