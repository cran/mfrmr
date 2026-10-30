delta_source <- function(type = "crossed", incomplete = FALSE, method = "anova") {
  x <- mfrmr_example_gtheory_crossed
  x <- x[x$Person %in% head(unique(x$Person), 8) &
    x$Rater %in% head(unique(x$Rater), 3) & x$Task %in% head(unique(x$Task), 4), ]
  if (type == "one") x <- aggregate(x[c("Content", "Organization")], x[c("Person", "Task")], mean)
  if (incomplete) x <- x[-c(1, 8, 23), ]
  mfrm_g_study(x, c("Content", "Organization"), method = method,
    rater = if (type == "one") NULL else "Rater",
    nesting = if (type == "nested") c(Rater = "Task") else NULL)
}

test_that("crossed coefficient SE agrees with independent MS covariance and ratio formula", {
  g <- delta_source()
  d <- mfrm_d_study(g, data.frame(Raters = c(2, 5), Tasks = c(3, 9)))
  n <- 8; r <- 3; t <- 4
  L <- rbind(c(r*t,0,0,t,r,0,1), c(0,n*t,0,t,0,n,1),
    c(0,0,n*r,0,r,n,1), c(0,0,0,t,0,0,1), c(0,0,0,0,r,0,1),
    c(0,0,0,0,0,n,1), c(0,0,0,0,0,0,1))
  A <- solve(L)
  ms <- vapply(g$mean_products, function(m) m[1,1], numeric(1))
  theta <- drop(A %*% ms)
  C <- A %*% diag(2*ms^2/(g$degrees_of_freedom$DF+2)) %*% t(A)
  for (metric in c("G", "Phi")) {
    ci <- confint(d, parm = metric, method = "delta", assumption = "normal")
    for (i in 1:2) {
      r <- d$design_grid$Raters[i]; t <- d$design_grid$Tasks[i]
      e <- c(0,0,0,1/r,1/t,0,1/(r*t))
      if (metric == "Phi") e <- e+c(0,1/r,1/t,0,0,1/(r*t),0)
      rho <- theta[1]/(theta[1]+sum(e*theta))
      # Almehrizi (2021), equations 24--25, in universe/error coordinates.
      variance <- ((1-rho)^2*C[1,1]+rho^2*drop(t(e)%*%C%*%e)-
        2*rho*(1-rho)*sum(C[1,]*e))/(theta[1]+sum(e*theta))^2
      expect_equal(ci$Estimate[i], unname(rho), tolerance = 1e-12)
      expect_equal(ci$SE[i]^2, unname(variance), tolerance = 1e-12)
      expect_equal(c(ci$Lower[i], ci$Upper[i]), unname(rho+c(-1,1)*qnorm(.975)*sqrt(variance)))
    }
  }
})

test_that("delta retains each design's covariance and signed composite gradient", {
  for (type in c("one", "crossed", "nested")) for (incomplete in c(FALSE, TRUE)) {
    g <- delta_source(type, incomplete, if (incomplete) "minque0" else "anova")
    grid <- if (type == "one") data.frame(Tasks = c(2, 7)) else
      data.frame(Raters = c(2, 5), Tasks = c(3, 7))
    w <- c(Content = .6, Organization = -.2)
    d <- mfrm_d_study(g, grid, w)
    before <- serialize(d, NULL)
    theta <- vapply(g$components, function(m) drop(t(w)%*%m%*%w), numeric(1))
    C <- vcov(g, weights = w, assumption = "normal")
    for (metric in c("G", "Phi")) {
      ci <- confint(d, parm = metric, method = "delta", assumption = "normal")
      for (i in 1:2) {
        m <- grid$Tasks[i]; r <- if (type != "one") grid$Raters[i]
        e <- switch(type, one = c(0,0,1/m), crossed = c(0,0,0,1/r,1/m,0,1/(r*m)),
          nested = c(0,0,0,1/m,1/(r*m)))
        if (metric == "Phi") e <- e + switch(type, one = c(0,1/m,0),
          crossed = c(0,1/r,1/m,0,0,1/(r*m),0), nested = c(0,1/m,1/(r*m),0,0))
        ratio <- function(x) x[1]/(x[1]+sum(e*x))
        h <- max(abs(theta))*1e-5
        gradient <- vapply(seq_along(theta), function(j) {
          step <- numeric(length(theta)); step[j] <- h
          (ratio(theta+step)-ratio(theta-step))/(2*h)
        }, numeric(1))
        variance <- drop(t(gradient)%*%C%*%gradient)
        expect_equal(ci$SamplingVariance[i], variance, tolerance = 1e-8)
      }
      scaled <- confint(mfrm_d_study(g, grid, -1e60*w),
        parm = metric, method = "delta", assumption = "normal")
      expect_equal(scaled$SE, ci$SE, tolerance = 1e-10)
      expect_equal(scaled$Lower, ci$Lower, tolerance = 1e-10)
      wide <- confint(d, parm = metric, method = "delta", level = .99, assumption = "normal")
      expect_true(all(wide$Lower <= ci$Lower & wide$Upper >= ci$Upper))
    }
    expect_identical(serialize(d, NULL), before)
    if (!incomplete) {
      b <- delta_source(type, FALSE, "minque0")
      ci <- confint(mfrm_d_study(b, grid, w), method = "delta", assumption = "normal")
      a <- confint(d, method = "delta", assumption = "normal")
      expect_equal(ci$SE, a$SE, tolerance = 1e-9)
    }
  }
})

test_that("positive delta SEs survive numerically unresolved confidence widths", {
  g <- delta_source("one")
  d <- mfrm_d_study(g, data.frame(Tasks = c(2L, 7L)))
  for (metric in c("G", "Phi")) {
    ci <- confint(d, parm = metric, method = "delta", assumption = "normal", level = 1e-18)
    expect_true(all(is.finite(ci$Estimate) & is.finite(ci$SE) & ci$SE > 0))
    expect_true(all(ci$SEStatus == "Available"))
    expect_false(any(ci$Available))
    expect_true(all(is.na(ci$Lower) & is.na(ci$Upper)))
    expect_true(all(ci$Status == "Delta interval below numeric resolution"))
    resolved <- confint(d, parm = metric, method = "delta", assumption = "normal", level = 1e-6)
    expect_identical(ci$Estimate, resolved$Estimate)
    expect_identical(ci$SE, resolved$SE)
    expect_true(all(resolved$Available & resolved$Upper > resolved$Lower))
  }
})

test_that("delta exposes boundary, undefined and negative-variance states", {
  x <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  g <- mfrm_g_study(x, "V", rater = NULL)
  n <- length(unique(x$Person)); m <- length(unique(x$Task))
  # Consistent saved ANOVA summaries give exact boundary cases.
  make <- function(ms) {
    a <- g; a$data <- NULL
    for (j in 1:3) a$mean_products[[j]][1,1] <- ms[j]
    values <- c((ms[1]-ms[3])/m, (ms[2]-ms[3])/n, ms[3])
    for (j in 1:3) a$components[[j]][1,1] <- values[j]
    mfrm_d_study(a, data.frame(Tasks = c(1, 6, 12)))
  }
  d <- make(c(1,2,1))
  zero_person <- confint(d, method = "delta", assumption = "normal")
  expect_true(all(zero_person$AtBoundary & zero_person$Available))
  expect_true(all(zero_person$Estimate == 0 & zero_person$SE > 0))
  expect_true(all(zero_person$Lower < 0 & zero_person$OutsideParameterSpace))
  for (ms in list(c(.5,2,1), c(0,0,0))) {
    ci <- confint(make(ms), method = "delta", assumption = "normal")
    expect_true(all(is.na(ci$Estimate) & is.na(ci$SE) & !ci$Available))
    expect_true(all(ci$Status == "Point projection unavailable"))
  }
  ci <- confint(make(c(2,1,0)), method = "delta", assumption = "normal")
  expect_true(all(ci$SE == 0 & !ci$Available & ci$AtBoundary))
  expect_true(all(ci$SEStatus == "Zero first-order variance"))
  # Test propagation independently of the moment-system algebra checked elsewhere.
  original <- vcov.mfrm_multivariate_gstudy
  local_mocked_bindings(vcov.mfrm_multivariate_gstudy = function(...) {
    result <- original(...); info <- attr(result, "diagnostics")
    info$normalized_covariance <- -diag(nrow(result)); info$positive_semidefinite <- FALSE
    attr(result, "diagnostics") <- info; result
  })
  ci <- confint(d, method = "delta", assumption = "normal")
  expect_true(all(ci$SamplingVariance < 0 & is.na(ci$SE) & !ci$Available))
  expect_true(all(ci$SEStatus == "Negative sampling-variance estimate"))
  expect_false(any(ci$SamplingCovariancePSD))
})

test_that("delta limits outside zero and one survive summaries, plots and saved reuse", {
  g <- delta_source("one")
  d <- mfrm_d_study(g, data.frame(Tasks = c(1, 4, 12)))
  ci <- confint(d, method = "delta", level = .9999, assumption = "normal")
  expect_true(any(ci$OutsideParameterSpace))
  s <- summary(d, intervals = ci)
  expect_identical(s$intervals[[1]], ci)
  expect_match(paste(capture.output(print(s)), collapse = "\n"), "SE")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(list(d = d, ci = ci), path); saved <- readRDS(path)
  local_mocked_bindings(confint.mfrm_multivariate_d_study = function(...) stop("Unexpected recalculation"))
  expect_identical(summary(saved$d, intervals = saved$ci), s)
  forest <- plot(saved$ci, draw = FALSE)
  expect_lte(forest$data$xrange[1], min(ci$Lower))
  expect_gte(forest$data$xrange[2], max(ci$Upper))
  curve <- plot(saved$d, metric = "G", intervals = saved$ci, draw = FALSE)
  expect_equal(curve$data$series$Lower, ci$Lower)
  limits <- .mfrm_mvds_y_limits(curve$data$series, TRUE, numeric())
  expect_lte(limits[1], min(ci$Lower)); expect_gte(limits[2], max(ci$Upper))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot_build(as_ggplot(curve))
    expect_equal(p$layout$panel_scales_y[[1]]$get_limits(), limits)
    expect_true(all(is.finite(p$data[[1]]$ymin)))
    expect_equal(range(p$data[[1]]$ymax), range(ci$Upper))
    expect_no_error(ggplot2::ggplot_build(as_ggplot(forest)))
  }
  changed <- ci; changed$Upper[1] <- changed$Lower[1]-1
  expect_error(summary(d, intervals = changed), "endpoints")
  expect_error(plot(d, metric = "G", intervals = changed, draw = FALSE), "endpoints")
  changed <- ci; changed$MethodCode[] <- "exact_f"
  expect_error(summary(d, intervals = changed), "endpoints")
})

test_that("delta is explicit, composes with other methods and checks source identity", {
  g <- delta_source("one")
  d <- mfrm_d_study(g, data.frame(Tasks = c(2, 8)))
  a <- confint(d, method = "delta", assumption = "normal")
  both <- confint(d, method = c("exact_f", "delta"), assumption = "normal")
  expect_identical(attr(both, "method_results")$delta, a)
  expect_no_error(summary(d, intervals = both))
  expect_identical(confint(d, assumption = "normal"), confint(d, method = "exact_f", assumption = "normal"))
  expect_no_error(plot(d, metric = "G", intervals = TRUE, method = "delta", assumption = "normal", draw = FALSE))
  expect_error(confint(d, method = "delta"), "normal")
  expect_error(confint(d, parm = "RelativeSEM", method = "delta", assumption = "normal"), "Unsupported")
  expect_error(confint(d, reference = 1, method = "delta", assumption = "normal"), "reference")
  g <- delta_source("crossed", TRUE, "minque0")
  d <- mfrm_d_study(g, data.frame(Raters = 3, Tasks = 4))
  ci <- confint(d, parm = "Phi", method = "delta", assumption = "normal")
  expect_no_error(summary(d, intervals = ci))
  expect_error(confint(d, assumption = "normal"), "explicitly select")
  bad <- ci; attr(bad, "source_data")$Content[1] <- -999
  expect_error(summary(d, intervals = bad), "source ratings")
  expect_error(plot(d, metric = "Phi", intervals = bad, draw = FALSE), "source ratings")
})
