mvgt_nested_fixture <- function() {
  n <- c(P = 8, T = 4, R = 3)
  data <- expand.grid(Person = paste0("P", 1:n[1]), Task = paste0("T", 1:n[2]),
    Rater = paste0("R", 1:n[3]), stringsAsFactors = FALSE)
  N <- nrow(data)
  projection <- function(formula) {
    X <- model.matrix(formula, data)
    qr <- qr(X); Q <- qr.Q(qr)[, seq_len(qr$rank), drop = FALSE]
    tcrossprod(Q)
  }
  # Independent fixed-effects projection spaces for the five orthogonal strata.
  P <- projection(~ Person); T <- projection(~ Task)
  RT <- projection(~ Task:Rater); PT <- projection(~ Person * Task)
  grand <- matrix(1/N, N, N)
  strata <- list(P - grand, T - grand, RT - T, PT - P - T + grand,
    diag(N) - PT - RT + T)
  df <- vapply(strata, function(a) round(sum(diag(a))), numeric(1))
  components <- lapply(list(c(.9, .1, .1, .7), c(.3, .05, .05, .2),
    c(.2, -.03, -.03, .25), c(.4, .1, .1, .3), c(.5, .15, .15, .6)),
    matrix, nrow = 2, dimnames = list(c("Content", "Organization"), c("Content", "Organization")))
  names(components) <- c("Person", "Task", "Rater(Task)", "Person:Task", "Residual")
  c <- components
  M <- list(n[3]*n[2]*c[[1]] + n[3]*c[[4]] + c[[5]],
    n[1]*n[3]*c[[2]] + n[1]*c[[3]] + n[3]*c[[4]] + c[[5]],
    n[1]*c[[3]] + c[[5]], n[3]*c[[4]] + c[[5]], c[[5]])
  y <- matrix(0, N, 2)
  for (i in 1:5) {
    basis <- eigen(strata[[i]], symmetric = TRUE)$vectors[, 1:2]
    y <- y + basis %*% chol(df[i] * M[[i]])
  }
  data[c("Content", "Organization")] <- sweep(y, 2, c(10, 20), "+")
  list(data = data, components = c, mean_products = setNames(M, names(c)), df = df)
}

test_that("balanced nested ANOVA and MINQUE recover independent QR components", {
  f <- mvgt_nested_fixture()
  for (method in c("anova", "minque0")) {
    g <- mfrm_g_study(f$data, c("Content", "Organization"),
      nesting = c(Rater = "Task"), method = method)
    expect_equal(g$components, f$components, tolerance = 1e-10)
    expect_identical(g$design$calculation_version, 3L)
    expect_identical(g$design$counts, c(Person = 8L, Rater = 12L, Task = 4L))
    expect_identical(g$design$child_counts, setNames(rep(3L, 4), paste0("T", 1:4)))
    expect_true(g$design$complete && g$design$balanced)
    expect_equal(g$design$observed_fraction, 1)
    if (method == "anova") {
      expect_equal(g$mean_products, f$mean_products, tolerance = 1e-10)
      expect_equal(g$degrees_of_freedom$DF, f$df)
    }
    expect_output(print(g), "total across parents")
    expect_equal(mfrm_d_study(g)$design_grid, data.frame(Raters = 3L, Tasks = 4L))
  }
})

test_that("large nested plan products remain finite in projections and intervals", {
  g <- mfrm_g_study(mvgt_nested_fixture()$data,
    c("Content", "Organization"), nesting = c(Rater = "Task"))
  grid <- data.frame(Raters = c(30000L, 50000L), Tasks = c(30000L, 50000L))
  d <- expect_warning(mfrm_d_study(g, grid), NA)
  reference <- mfrm_d_study(g, as.data.frame(lapply(grid, as.double)))
  expect_identical(d$design_grid, grid)
  expect_identical(d$covariances, reference$covariances)
  expect_equal(d$coefficients, reference$coefficients)
  expect_equal(d$coefficients$RelativeErrorVariance[3L],
    .4 / 50000 + .5 / 2500000000, tolerance = 1e-14)
  # This mapping is shared by delta and REML-profile coefficient calculations.
  expect_equal(mfrmr:::.mfrm_mvgt_profile_error(d, 2L, "Phi"),
    setNames(c(0, 1/50000, 1/2500000000, 1/50000, 1/2500000000), names(g$components)))
  for (method in c("delta", "approximate_f", "mls")) {
    metric <- if (method == "mls") "AbsoluteSEM" else "Phi"
    a <- expect_warning(confint(d, parm = metric, score = "Content",
      method = method, assumption = "normal"), NA)
    b <- confint(reference, parm = metric, score = "Content", method = method, assumption = "normal")
    expect_equal(a, b)
    expect_true(all(is.finite(a$Lower) & is.finite(a$Upper)))
  }
})

test_that("nested incomplete and unequal sources agree with independent kernels", {
  full <- mvgt_nested_fixture()$data
  unequal <- full[!(full$Task == "T1" & full$Rater == "R3"), ]
  for (data in list(full[-c(1, 10, 70), ], unequal)) {
    g <- mfrm_g_study(data, c("Content", "Organization"),
      nesting = c(Rater = "Task"), method = "minque0")
    groups <- list(data$Person, data$Task, interaction(data$Task, data$Rater),
      interaction(data$Person, data$Task), seq_len(nrow(data)))
    kernels <- lapply(groups, function(x) outer(x, x, `==`) * 1)
    H <- diag(nrow(data)) - 1/nrow(data)
    B <- lapply(kernels, function(k) H %*% k %*% H)
    S <- crossprod(vapply(B, as.vector, numeric(nrow(data)^2)))
    y <- as.matrix(data[c("Content", "Organization")])
    Q <- t(vapply(B, function(b) as.vector(crossprod(y, b %*% y)), numeric(4)))
    expect_equal(unname(t(vapply(g$components, as.vector, numeric(4)))), solve(S, Q), tolerance = 1e-10)
    # Expected quadratic products must recover known population components.
    theta <- c(.9, .3, .2, .4, .5)
    V <- Reduce(`+`, Map(`*`, kernels, theta))
    expected_Q <- vapply(B, function(b) sum(b * V), numeric(1))
    expect_equal(as.vector(solve(S, expected_Q)), theta, tolerance = 1e-10)
    expect_error(mfrm_g_study(data, "Content", nesting = c(Rater = "Task")), "balanced")
    expect_error(mfrm_d_study(g), "explicit.*design_grid")
    expect_s3_class(mfrm_d_study(g, data.frame(Raters = 2, Tasks = 6)), "mfrm_multivariate_d_study")
  }
  expect_true(g$design$complete)
  expect_false(g$design$balanced)
  expect_equal(g$design$potential_cells, nrow(unequal))
  expect_equal(unname(g$design$child_counts), c(2, 3, 3, 3))
  # Omission preserves the original input row and excludes the whole score vector.
  bad <- full; bad$Content[1] <- NA
  g <- mfrm_g_study(bad, c("Content", "Organization"), nesting = c(Rater = "Task"),
    method = "minque0", missing = "omit")
  expect_identical(g$data_usage$excluded_rows, 1L)
  expect_equal(g$design$potential_cells, nrow(full))
})

test_that("nested identities, facet order and weighted projections are invariant", {
  data <- mvgt_nested_fixture()$data
  g <- mfrm_g_study(data, c("Content", "Organization"), nesting = c(Rater = "Task"))
  unique_ids <- data; unique_ids$Rater <- paste(data$Task, data$Rater, sep = ".")
  unique_g <- mfrm_g_study(unique_ids[nrow(data):1, ], c("Content", "Organization"),
    nesting = c(Rater = "Task"))
  expect_equal(unique_g$components, g$components)
  expect_equal(unique_g$design$child_counts, g$design$child_counts)
  named <- mfrm_g_study(data, c("Content", "Organization"),
    facets = c(Station = "Task", Judge = "Rater"), nesting = c(Judge = "Station"))
  expect_equal(unname(named$components), unname(g$components))
  grid <- data.frame(Raters = c(2, 4), Tasks = c(3, 6))
  weights <- c(Content = .6, Organization = .4)
  d <- mfrm_d_study(g, grid, weights)
  output <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(output, "Person = 8, Rater = 12, Task = 4", fixed = TRUE)
  expect_match(output, "observed child count is the total across parents", fixed = TRUE)
  expect_match(output, "Raters is the count per Task", fixed = TRUE)
  named_d <- mfrm_d_study(named, data.frame(Judge = grid$Raters, Station = grid$Tasks), weights)
  expect_equal(named_d$coefficients[-c(2, 3)], d$coefficients[-c(2, 3)])
  expect_equal(named_d$covariances, d$covariances)
  data$Total <- .6 * data$Content + .4 * data$Organization
  scalar <- mfrm_d_study(mfrm_g_study(data, "Total",
    nesting = c(Rater = "Task")), grid)$coefficients
  measures <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance", "G", "Phi", "RelativeSEM", "AbsoluteSEM")
  expect_equal(d$coefficients[d$coefficients$Kind == "Composite", measures], scalar[measures], ignore_attr = TRUE)
  expected_relative <- g$components$`Person:Task` / 3 + g$components$Residual / 6
  expected_absolute <- expected_relative + g$components$Task / 3 + g$components$`Rater(Task)` / 6
  expect_equal(d$covariances[[1]]$RelativeError, expected_relative)
  expect_equal(d$covariances[[1]]$AbsoluteError, expected_absolute)
  expect_error(mfrm_multivariate_d_compare(d, assumption = "normal"), "nested-design intervals")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(g, path)
  expect_equal(mfrm_d_study(readRDS(path), grid, weights), d)
})

test_that("nested D-study agrees with published random-facet divisors", {
  g <- mfrm_g_study(mvgt_nested_fixture()$data, "Content", nesting = c(Rater = "Task"))
  # Brennan (1992), NCME module pp. 31-32, equations 14-16/Table 3.
  # This checks the published component-to-projection calculation, not raw-data recovery.
  for (j in 1:5) g$components[[j]][, ] <- c(.25, .06, .02, .15, .16)[j]
  d <- mfrm_d_study(g, data.frame(Raters = 2, Tasks = 3))$coefficients
  expect_equal(d$UniverseVariance, .25)
  expect_equal(d$RelativeErrorVariance, .15/3 + .16/6)
  expect_equal(d$AbsoluteErrorVariance, .1)
  expect_equal(d$G, .25/(.25 + .15/3 + .16/6))
  expect_equal(d$Phi, .25/.35)
})

test_that("nested design refuses unsupported or confounded specifications", {
  data <- mvgt_nested_fixture()$data
  for (nesting in list(c(Rater = "Person"), c(Rater = "Rater"), "Task", character(),
    c(Rater = "Task", Task = "Rater"), c(Unknown = "Task"))) {
    expect_error(mfrm_g_study(data, "Content", nesting = nesting), "nesting")
  }
  expect_error(mfrm_g_study(data, "Content", rater = NULL,
    nesting = c(Rater = "Task")), "nesting")
  one <- data[data$Rater == "R1", ]; one$Rater <- paste0("R_", one$Task)
  expect_error(mfrm_g_study(one, "Content", nesting = c(Rater = "Task")), "Cannot separate")
  expect_error(mfrm_g_study(one, "Content", nesting = c(Rater = "Task"),
    method = "minque0"), "Cannot separate")
  duplicate <- rbind(data, data[1, ])
  expect_error(mfrm_g_study(duplicate, "Content", nesting = c(Rater = "Task")), "Duplicate")
  g <- mfrm_g_study(data, "Content", nesting = c(Rater = "Task"))
  g$design$child_counts[1] <- 9L
  expect_error(mfrm_d_study(g), "child counts")
  g$design$nesting <- NULL
  expect_error(mfrm_d_study(g), "nesting specification")
})

test_that("nested plots and ggplot retain per-parent counts and point values", {
  g <- mfrm_g_study(mvgt_nested_fixture()$data, "Content", nesting = c(Rater = "Task"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 4), Tasks = c(3, 6)))
  payload <- plot_data(plot(d, x_var = "Raters", draw = FALSE))
  expect_match(payload$x_label, "per Task", fixed = TRUE)
  expect_identical(payload$nesting, c(Rater = "Task"))
  expect_match(payload$subtitle, "nested designs: Rater within Task", fixed = TRUE)
  other <- plot_data(plot(d, x_var = "Tasks", draw = FALSE))
  expect_true(all(grepl("Raters per Task", other$legend$label, fixed = TRUE)))
  expect_equal(payload$series$Value[payload$series$Metric == "G"], d$coefficients$G)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(d, x_var = "Raters")
    expect_match(gg$labels$x, "per Task", fixed = TRUE)
    expect_match(paste(plot_data(gg, "notes")$Text, collapse = "\n"), "nested designs", fixed = TRUE)
    expect_equal(gg$data$Value, payload$series$Value)
  }
})

# GENOVA 3.1 cc.manual, first D-study; urGENOVA 2.1 cc.manual, balanced
# p x (i:r) runs. Executed on Windows 2026-10-23. External r maps to Task,
# i:r to Rater(Task); SEVC uses df+2. Tolerances cover printed rounding.
test_that("executed GENOVA nested components and urGENOVA intervals agree", {
  x <- read.csv(test_path("fixtures", "genova-31-nested.csv"))
  g <- mfrm_g_study(x, "Score", nesting = c(Rater = "Task"))
  tab <- attr(vcov(g, assumption = "normal"), "summary")
  expect_equal(tab$Source, c("Person", "Task", "Rater(Task)", "Person:Task", "Residual"))
  expect_lt(max(abs(tab$Estimate - c(.47315, .32515, .64753, .55957, 2.38025))), 5.1e-6)
  expect_lt(max(abs(tab$SE - c(.38558, .43799, .37941, .37663, .36949))), 5.1e-6)
  ci <- confint(g, assumption = "normal", level = .8)
  expect_lt(max(abs(ci$Lower - c(.0606, 0, .2987, .1746, 1.9738))), 5.1e-5)
  expect_lt(max(abs(ci$Upper - c(1.4657, 5.4144, 1.6724, 1.3202, 2.9581))), 5.1e-5)
  d <- mfrm_d_study(g, data.frame(Raters = 1:4, Tasks = 3))$coefficients
  expect_lt(max(abs(d$G - c(.32562, .44790, .51199, .55144))), 5.1e-6)
  expect_lt(max(abs(d$Phi - c(.26622, .37177, .42839, .46370))), 5.1e-6)
  expect_lt(max(abs(d$RelativeSEM - c(.98992, .76370, .67156, .62038))), 5.1e-6)
  expect_lt(max(abs(d$AbsoluteSEM - c(1.14200, .89417, .79456, .73974))), 5.1e-6)
})
