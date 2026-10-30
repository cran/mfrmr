# Target definitions and selected-interval accounting; no fitting here.
contrast <- new.env(parent = globalenv())
for (name in c("mml-stages-20261001", "jml-stages-20261001", "output-summary-20261001",
    "multi-output-20261001", "recovery-20261002", "interval-recovery-20261002")) {
  path <- test_path("..", "..", "inst", "validation", paste0("mfrm-wide-map-", name, ".R"))
  skip_if_not(file.exists(path), "Repository-only interval recovery")
  sys.source(path, contrast)
}
contrast_fixture <- function(arm = "PCM-MML-free") {
  root <- test_path("..", "..", "validation-results", "mfrm-wide-map-recovery-20261002")
  skip_if_not(file.exists(file.path(root, "witness-manifest.rds")), "Retained PCM input/fit")
  m <- readRDS(file.path(root, "witness-manifest.rds"))$payload
  q <- if (m$specs[[arm]]$args$method == "MML") 31L else 0L
  list(input = m$input, spec = contrast$wide_interval_spec(m$input, m$specs[[arm]]),
    fit = readRDS(file.path(root, arm, paste0("q",q,"-fit.rds")))$payload$value)
}

test_that("declared contrasts map by level label and retain their coefficient identity", {
  x <- contrast_fixture(); c <- x$spec$outputs$rater
  p <- contrast$wide_interval_plan(x$input, x$spec)
  expect_equal(p$ReferenceValue[p$Output == "interval"], c(-.3,-.5))
  t <- contrast$wide_multi_points(x$fit, x$spec)
  f <- x$fit$facets$others
  expect_equal(t$Estimate[1], unname(f$Estimate[f$Facet == "Rater" & f$Level == "1"]-
    f$Estimate[f$Facet == "Rater" & f$Level == "2"]))
  reverse <- c; reverse$contrasts <- c$contrasts[, 3:1, drop = FALSE]
  expect_identical(contrast$wide_multi_targets(x$spec, c), contrast$wide_multi_targets(x$spec, reverse))
  changed <- c; changed$contrasts <- -c$contrasts
  expect_false(identical(contrast$wide_multi_targets(x$spec, c)$TargetDefinition,
    contrast$wide_multi_targets(x$spec, changed)$TargetDefinition))
  bad <- x$spec; bad$outputs$rater$contrasts[1, 1] <- 2
  expect_error(contrast$wide_multi_validate(bad))
  bad <- x$spec; colnames(bad$outputs$rater$contrasts)[1] <- "unknown"
  expect_error(contrast$wide_multi_validate(bad))
  a <- p[p$Output == "point", ]; b <- a; b$Arm <- "different"
  b$TargetDefinition <- "different coefficients"
  plan <- rbind(a,b); r <- plan[c("ConditionId","Arm","Replicate","Target","Output")]
  r$Status <- "returned"; r$Available <- TRUE; r$Estimate <- plan$ReferenceValue
  r$SE <- r$Lower <- r$Upper <- NA_real_; r$SourceStage <- "fixture"; r$Reason <- ""
  expect_error(contrast$wide_paired_point_summary(plan,r,a$Arm[1],"different"), "contrast definitions")
})

test_that("fixed normal PCM and JML retain their different interval procedures", {
  x <- contrast_fixture("PCM-MML-fixed")
  expect_true(contrast$wide_mml_validate(contrast$wide_multi_view(x$spec)))
  p <- contrast$wide_interval_plan(x$input, x$spec)
  expect_true(all(p$Eligibility == "eligible"))
  expect_true(all(p$Procedure[p$Output == "interval"] == "fixed_population_model_existing"))
  for (arm in c("PCM-JML","GPCM-JML","GPCM-JML2","GPCM-JML4")) {
    z <- contrast_fixture(arm); plan <- contrast$wide_interval_plan(z$input, z$spec)
    expect_true(all(plan$Eligibility[plan$Output == "interval"] == "unsupported"))
    expect_true(all(plan$Eligibility[plan$Output == "point"] == "eligible"))
  }
  bad <- x$spec; bad$args$population_formula <- ~group
  expect_error(contrast$wide_multi_validate(bad))
})

test_that("legacy Status is normalized without editing public interval objects", {
  x <- contrast_fixture("PCM-MML-fixed"); c <- x$spec$outputs$rater
  ci <- structure(list(fit = x$fit, contrasts = c$contrasts,
    settings = list(facet = "Rater", method = "model", level = .95),
    table = data.frame(Target = "1-2", Estimate = .2, SE = .3, Lower = -.4, Upper = .8,
      Status = "available")), class = "mfrm_facet_intervals")
  original <- ci
  z <- contrast$wide_multi_interval_table(ci, x$spec, c, x$fit)
  expect_true(z$Available); expect_equal(z$SE,.3)
  expect_identical(ci, original)
  expect_false(contrast$wide_multi_retry(ci,c))
  ci$table$Status <- "nonpositive_variance"
  z <- contrast$wide_multi_interval_table(ci, x$spec, c, x$fit)
  expect_false(z$Available); expect_true(all(is.na(z[c("SE","Lower","Upper")])))
  expect_match(z$Reason,"nonpositive_variance")
  ci$contrasts <- -ci$contrasts
  expect_error(contrast$wide_multi_interval_table(ci, x$spec, c, x$fit))
})

test_that("selected contrasts use the public full covariance and their own estimate", {
  x <- contrast_fixture(); C <- x$spec$outputs$rater$contrasts
  # A covariance fixture with a nonzero off-diagonal makes the independence
  # shortcut wrong. The public consumer must apply its constraint Jacobian.
  p <- length(x$fit$opt$par); V <- diag(p)
  sizes <- mfrmr:::build_param_sizes(x$fit$config)
  slice <- mfrmr:::build_param_slices(sizes)$Rater
  V[slice[1],slice[2]] <- V[slice[2],slice[1]] <- .4
  local_mocked_bindings(mfrm_native_location_inference = function(fit) list(
    covariance = V, information = list(), check = list(eligible = TRUE, review = "fixture"),
    checks = data.frame(Check = "fixture", Passed = TRUE)), .package = "mfrmr")
  ci <- suppressWarnings(contrast$wide_multi_consume(x$fit, x$spec$outputs$rater, .95))
  jac <- C %*% mfrmr:::constraint_jacobian(x$fit$config$facet_specs$Rater)
  expected <- drop(jac %*% V[slice,slice] %*% t(jac))
  expect_equal(ci$table$SE^2, expected)
  expect_equal(expected, 1.2)
  t <- contrast$wide_multi_interval_table(ci, x$spec, x$spec$outputs$rater, x$fit)
  expect_true(t$Available)
  expect_equal(t$Estimate, contrast$wide_multi_points(x$fit,x$spec)$Estimate[1])
  ci$settings$level <- .9
  expect_error(contrast$wide_multi_interval_table(ci,x$spec,x$spec$outputs$rater,x$fit))
})

test_that("interval truth joins preserve unavailability and unsupported denominators", {
  x <- contrast_fixture(); plan <- contrast$wide_multi_plan(x$spec)
  point_table <- contrast$wide_multi_points(x$fit,x$spec)
  selected <- list(spec_hash = contrast$wide_mml_hash(x$spec), plan = plan,
    point = list(error = "", stage = "q61:fixture", value = list(ready = TRUE), table = point_table),
    intervals = list())
  for (id in names(x$spec$outputs)) {
    t <- contrast$wide_multi_targets(x$spec,x$spec$outputs[[id]])
    t$Estimate <- if(id == "rater") -.3 else -.5
    t$SE <- .1; t$Lower <- t$Estimate-.2; t$Upper <- t$Estimate+.2
    t$Available <- id == "rater"; t$Reason <- if(id == "rater") "" else "numerical review"
    t[!t$Available,c("SE","Lower","Upper")] <- NA_real_
    selected$intervals[[id]] <- list(error = "", stage = "q31:fixture", table = t)
  }
  z <- contrast$wide_interval_recovery(x$input,x$spec,selected)
  interval <- z$rows[z$rows$Output == "interval", ]
  expect_identical(interval$SourceStage,rep("q31:fixture",2))
  expect_equal(interval$ReferenceValue,c(-.3,-.5))
  expect_identical(interval$Available,c(TRUE,FALSE))
  expect_true(all(z$summary$Complete))
  expect_equal(z$summary$ReturnedAndCovered[z$summary$Output == "interval"],c(0,1))
  changed <- x$spec; changed$outputs$rater$contrasts <- -changed$outputs$rater$contrasts
  expect_error(contrast$wide_interval_recovery(x$input,changed,selected),"exact shared job")
})
