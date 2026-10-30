dstudy_crossed_example <- function() {
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed,
    c("Content", "Organization"))
  mfrm_d_study(g, data.frame(Raters = c(2, 4), Tasks = c(6, 12)),
    weights = c(Content = .5, Organization = .5))
}

test_that("method selection preserves defaults and rejects mismatched targets", {
  data <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  d <- mfrm_d_study(mfrm_g_study(data, c("V", "W"), rater = NULL))
  expect_identical(confint(d, assumption = "normal"),
    confint(d, assumption = "normal", method = "exact_f"))
  expect_identical(confint(d, parm = "Phi", assumption = "normal", draws = 4000),
    confint(d, parm = "Phi", assumption = "normal", method = "gci", draws = 4000))
  expect_error(confint(d, assumption = "normal", method = "gci"), "Unsupported")
  expect_error(confint(d, parm = "Phi", assumption = "normal", method = "exact_f"), "Unsupported")
})

test_that("explicit two-facet methods preserve estimates and source context", {
  skip_if_not_installed("lpSolve")
  d <- dstudy_crossed_example(); original <- d
  ci <- confint(d, parm = "Phi", assumption = "normal",
    method = c("approximate_f", "mean_square_region"))
  expect_s3_class(ci, "mfrm_dstudy_intervals")
  expect_identical(ci$MethodCode, rep(c("approximate_f", "mean_square_region"), each = 2))
  expect_equal(ci$Estimate, rep(subset(d$coefficients, Kind == "Composite")$Phi, 2))
  expect_true(all(ci$Available)); expect_true(all(grepl("Experimental", ci$Validation)))
  expect_true(all(grepl("undercoverage documented",
    ci$Validation[ci$MethodCode == "approximate_f"], fixed = TRUE)))
  expect_true(all(grepl("conservative joint coverage under the normal model",
    ci$Validation[ci$MethodCode == "mean_square_region"], fixed = TRUE)))
  expect_identical(d, original)
  expect_equal(attr(ci, "weights"), c(Content = .5, Organization = .5))
  expect_length(attr(ci, "diagnostics"), 4)
  expect_identical(confint(d, assumption = "normal"),
    confint(d, assumption = "normal", method = "mls"))
  expect_error(confint(d, assumption = "normal", method = "exact_f"), "Unsupported")
  expect_error(confint(d, assumption = "normal", method = c("auto", "approximate_f")), "Unsupported")
  expect_error(confint(d, assumption = "normal", method = rep("approximate_f", 2)), "distinct")
  expect_error(confint(d, assumption = "normal", method = "approximate_f", reference = 1), "reference")
  expect_error(confint(d, method = "approximate_f"), "normal")
  for (bad in list(NA_character_, character(), 1, matrix("approximate_f")))
    expect_error(confint(d, assumption = "normal", method = bad), "method")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(ci, path); expect_identical(readRDS(path), ci)
  d$gstudy$data <- NULL
  expect_identical(confint(d, parm = "Phi", assumption = "normal",
    method = c("approximate_f", "mean_square_region")), ci)
  scaled <- mfrm_d_study(d$gstudy, d$design_grid, -10*d$weights)
  other <- confint(scaled, parm = "Phi", assumption = "normal",
    method = c("approximate_f", "mean_square_region"))
  expect_equal(other$Lower, ci$Lower, tolerance = 1e-8)
  expect_equal(other$Upper, ci$Upper, tolerance = 1e-8)
})

test_that("approximate F reproduces published G and independently derived Phi", {
  # Wong & McGraw Table 3. Independent expected-MS design, not package point code.
  n <- 40; r <- 8; t <- 3
  L <- rbind(c(r*t,0,0,t,r,0,1), c(0,n*t,0,t,0,n,1),
    c(0,0,n*r,0,r,n,1), c(0,0,0,t,0,0,1), c(0,0,0,0,r,0,1),
    c(0,0,0,0,0,n,1), c(0,0,0,0,0,0,1))
  df <- c(39,7,2,273,78,14,546); map <- list(A = solve(L), df = df)
  ms <- c(439.957,68.574,89.752,336.301,118.415,14.348,666.152)/df
  errors <- list(c(0,0,0,1/r,1/t,0,1/(r*t)), c(0,0,0,1,1,0,1),
    c(0,1/r,1/t,1/r,1/t,1/(r*t),1/(r*t)))
  expected <- list(c(.76376405,.92501745), c(.16331788,.36333347), c(.50391444,.89063001))
  for (i in seq_along(errors)) {
    z <- .mfrm_mvgt_approximate_f(ms, map, errors[[i]], .95)
    expect_equal(z$bounds, expected[[i]], tolerance = 1e-7)
    for (scale in c(1e-200, 1e200)) expect_equal(
      .mfrm_mvgt_approximate_f(ms*scale, map, errors[[i]], .95)$bounds, z$bounds)
    wide <- .mfrm_mvgt_approximate_f(ms, map, errors[[i]], .99)$bounds
    expect_true(wide[1] <= z$bounds[1] && wide[2] >= z$bounds[2])
  }
})

test_that("negative components and a missing optional solver stay explicit", {
  x <- mfrmr_example_gtheory_crossed
  small <- subset(x, Person %in% head(unique(Person), 10) &
    Rater %in% head(unique(Rater), 2) & Task %in% head(unique(Task), 4))
  g <- mfrm_g_study(small, c("Content", "Organization"))
  d <- mfrm_d_study(g, weights = c(Content = .5, Organization = .5))
  ci <- confint(d, parm = "Phi", assumption = "normal", method = "approximate_f")
  expect_true(all(ci$NegativeComponent)); expect_true(all(ci$Available))
  expect_true(any(grepl("experimental extension", attr(ci, "notes"))))
  # Isolate dependency lookup without altering the loaded package namespace.
  f <- .mfrm_mvgt_two_facet_intervals
  environment(f) <- list2env(list(requireNamespace = function(...) FALSE), parent = environment(f))
  expect_error(f(d, NULL, NULL, "Phi", .95, "mean_square_region"), "optional lpSolve")
})

test_that("F intervals invert inequalities when future-count divisors reverse", {
  d <- dstudy_crossed_example(); p <- .mfrm_mvgt_anova_projection(d$gstudy, weights = d$weights)
  # Signed weights and small sample df are stressed with a deterministic MS grid.
  map <- .mfrm_mvgt_anova_interval_design(list(counts = c(Person = 10, Rater = 2, Task = 4)))
  grid <- seq(0, 1, length.out = 1001)
  for (ms in list(c(1,1,1,.6,.7,.3,1), c(.1,1,1,.6,.7,.3,1), p$ms)) {
    for (future in list(c(1,1), c(2,8), c(20,40))) {
      e <- c(0,0,0,1/future[1],1/future[2],0,1/prod(future))
      z <- .mfrm_mvgt_approximate_f(ms, map, e, .95)
      if (!is.finite(z$denominator_df) || z$denominator_df <= 0 ||
          z$status == "Invalid F moments") next
      q <- qf(c(.025, .975), map$df[1], z$denominator_df)
      if (any(!is.finite(q) | q <= 0)) {
        expect_identical(z$status, "Unrepresentable F quantile")
        expect_true(all(is.na(z$bounds)))
        next
      }
      a <- map$A[1,]; delta <- drop(e %*% map$A)
      den <- vapply(grid, function(rho) sum((-a[-1]+rho*(a[-1]+delta[-1]))/a[1]*ms[-1]), numeric(1))
      F <- (1-grid)*ms[1]/den
      expected <- den > 0 & F >= qf(.025, map$df[1], z$denominator_df) &
        F <= qf(.975, map$df[1], z$denominator_df)
      returned <- if (all(is.finite(z$bounds))) grid >= z$bounds[1] & grid <= z$bounds[2] else rep(FALSE, length(grid))
      expect_identical(returned, expected)
    }
  }
})

test_that("mean-square region includes every known feasible parameter vector", {
  skip_if_not_installed("lpSolve")
  map <- .mfrm_mvgt_anova_interval_design(list(counts = c(Person = 10, Rater = 2, Task = 4)))
  theta <- c(1,.2,.3,.3,.4,.1,1); ms <- drop(solve(map$A) %*% theta)
  for (level in c(.8,.95,.99)) {
    setup <- .mfrm_mvgt_ms_region_setup(map, level)
    expect_equal((1-2*setup$tail)^7, level)
    for (e in list(c(0,0,0,.5,.25,0,.125), c(0,.5,.25,.5,.25,.125,.125))) {
      z <- .mfrm_mvgt_ms_region(ms, setup, e)
      truth <- theta[1]/(theta[1] + sum(e*theta))
      expect_true(z$bounds[1] <= truth && z$bounds[2] >= truth)
      expect_true(all(z$residual <= 1e-7))
      expect_equal(.mfrm_mvgt_ms_region(ms*1e150, setup, e)$bounds, z$bounds, tolerance = 1e-8)
      bad <- ms; bad[1] <- 1e-12
      failed <- .mfrm_mvgt_ms_region(bad, setup, e)
      expect_true(failed$status %in% c("Empty mean-square region", "Linear-program solver failure"))
      expect_true(all(is.na(failed$bounds)))
      expect_identical(.mfrm_mvgt_ms_region(rep(0,7),setup,e)$status, "Degenerate mean squares")
    }
  }
})

test_that("forest plots retain unavailable rows, saved data and annotations", {
  skip_if_not_installed("lpSolve")
  ci <- confint(dstudy_crossed_example(), assumption = "normal",
    method = c("approximate_f", "mean_square_region"))
  ci$Lower[2] <- ci$Upper[2] <- ci$Estimate[2] <- NA_real_
  p <- plot(ci, draw = FALSE, preset = "monochrome")
  expect_equal(nrow(plot_data(p, "table")), 4)
  expect_equal(nrow(plot_data(p, "unavailable")), 1)
  expect_identical(plot_data(p, "intervals"), ci)
  expect_true(all(plot_data(p, "table")$Colour == "black"))
  tab <- plot_data(p, "table")
  expect_equal(length(unique(tab$Y)), 4)
  expect_identical(plot_data(p, "layout"), "combined")
  expect_equal(tab$Lower, ci$Lower); expect_equal(tab$Upper, ci$Upper)
  expect_equal(tab$Estimate, ci$Estimate)
  expect_true(any(grepl("Experimental", plot_data(p, "notes")$Text)))
  expect_error(plot(ci, draw = NA), "draw")
  expect_error(plot(ci, ncol = 0), "ncol")
  expect_error(plot(ci, layout = "invalid"), "arg")
  path <- tempfile(fileext = ".pdf"); on.exit(unlink(path), add = TRUE)
  grDevices::pdf(path, width = 8, height = 4)
  tryCatch(expect_identical(plot(ci), plot(ci, draw = FALSE)), finally = grDevices::dev.off())
  skip_if_not_installed("ggplot2")
  gg <- as_ggplot(p); built <- ggplot2::ggplot_build(gg)
  expect_equal(nrow(built$data[[1]]), 4) # blank geometry preserves missing rows
  expect_equal(nrow(built$data[[3]]), 3)
  expect_equal(nrow(built$data[[4]]), 3)
  expect_equal(nrow(built$data[[5]]), 1)
  expect_equal(nrow(built$layout$layout), 1)
  expect_null(gg$labels$title); expect_null(gg$labels$caption)
  panels <- ggplot2::ggplot_build(as_ggplot(plot(ci, layout = "panels", draw = FALSE)))
  expect_equal(nrow(panels$layout$layout), 2)
  expect_equal(nrow(panels$data[[2]]), 3)
  ci$Lower[] <- ci$Upper[] <- ci$Estimate[] <- NA_real_
  empty <- ggplot2::ggplot_build(as_ggplot(plot(ci, draw = FALSE)))
  expect_equal(nrow(empty$data[[1]]), 4)
  expect_equal(nrow(empty$data[[3]]), 0)
  expect_equal(nrow(empty$data[[4]]), 0)
  expect_equal(nrow(empty$data[[5]]), 4)
})
