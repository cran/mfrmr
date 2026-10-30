test_that("two-facet intervals retain points when their quantiles lose numerical separation", {
  skip_if_not_installed("lpSolve")
  for (nested in c(FALSE, TRUE)) {
    x <- if(nested)mfrmr_example_gtheory_nested else mfrmr_example_gtheory_crossed
    g <- mfrm_g_study(x,"Content",nesting=if(nested)c(Rater="Task"))
    d <- mfrm_d_study(g,data.frame(Raters=c(3,3,6),Tasks=c(6,6,12)))
    for (method in c("approximate_f","mean_square_region")) {
      levels <- if(method=="approximate_f")c(1e-16,1e-100) else c(1e-100,1e-300)
      for (level in levels) {
        ci <- confint(d,method=method,assumption="normal",level=level)
        expect_false(any(ci$Available))
        expect_true(all(is.na(ci$Lower)&is.na(ci$Upper)))
        expect_identical(ci$Estimate,d$coefficients$G)
        expect_true(all(grepl("below numeric resolution",ci$Status,fixed=TRUE)))
        for (diagnostic in attr(ci,"diagnostics")) {
          if (method=="approximate_f") {
            expect_length(diagnostic$f_quantiles,2)
            expect_lte(diagnostic$quantile_log_gap,diagnostic$quantile_resolution)
          } else {
            expect_true(all(diagnostic$relative_width<=diagnostic$region_resolution))
            expect_true(all(is.na(diagnostic$solver_status)))
          }
        }
      }
    }
    for (metric in c("G","Phi","RelativeSEM","AbsoluteSEM")) {
      ci <- confint(d,parm=metric,method="mean_square_region",assumption="normal",
        level=1e-100,reference=1)
      expect_false(any(ci$Available))
      expect_true(all(grepl("Reference interval unavailable",ci$Status,fixed=TRUE)))
      expect_equal(ci$Estimate,d$coefficients[[metric]][-1]-d$coefficients[[metric]][1])
    }
    # A tiny joint level need not collapse its much wider per-stratum regions.
    ci <- confint(d,method="mean_square_region",assumption="normal",level=1e-16)
    expect_true(all(ci$Available & ci$Lower<ci$Upper))
    expect_true(all(vapply(attr(ci,"diagnostics"),function(z)
      all(z$mean_square_bounds[,2]>z$mean_square_bounds[,1]),logical(1))))
    ci <- confint(d,method="approximate_f",assumption="normal",level=1e-6)
    expect_true(all(ci$Available & ci$Lower<ci$Upper))
  }
})

test_that("region checks distinguish degeneracy, numerical loss and true infeasibility", {
  skip_if_not_installed("lpSolve")
  L <- rbind(c(6,0,0,2,1),c(0,20,10,2,1),c(0,0,10,0,1),
    c(0,0,0,2,1),c(0,0,0,0,1))
  map <- list(A=solve(L),df=c(9,2,3,18,27))
  setup <- .mfrm_mvgt_ms_region_setup(map,.95)
  e <- c(0,0,0,1/3,1/6)
  ms <- drop(L%*%c(.6,.2,.3,.1,.8))
  expect_identical(.mfrm_mvgt_ms_region(rep(0,5),setup,e)$status,"Degenerate mean squares")
  tiny <- ms; tiny[1] <- 1e-300
  expect_identical(.mfrm_mvgt_ms_region(tiny,setup,e)$status,"Empty mean-square region")
  # Normalization must not turn positive mean squares into exact zero constraints.
  lost <- ms; lost[1] <- 1e-300; lost[2] <- 1e300
  ans <- .mfrm_mvgt_ms_region(lost,setup,e)
  expect_identical(ans$status,"Mean-square region below numeric resolution")
  expect_true(all(is.na(ans$solver_status)))
  expect_true(all(is.na(ans$bounds)))
})

test_that("region SEM transformations report positive widths lost by square roots", {
  skip_if_not_installed("lpSolve")
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed,"Content")
  d <- mfrm_d_study(g,data.frame(Raters=c(3,3),Tasks=c(6,6)))
  candidate <- c(1,1+.Machine$double.eps)
  expect_gt(diff(candidate),0)
  expect_equal(sqrt(candidate[1]),sqrt(candidate[2]),tolerance=0)
  local_mocked_bindings(.mfrm_mvgt_ms_region=function(...)
    list(bounds=candidate,status="Computed conservative interval"),.package="mfrmr")
  ci <- confint(d,parm="RelativeSEM",method="mean_square_region",assumption="normal")
  expect_false(any(ci$Available))
  expect_identical(ci$Status,rep("Unrepresentable SEM endpoints",2))
  expect_identical(ci$Estimate,d$coefficients$RelativeSEM)
  expect_identical(attr(ci,"diagnostics")[[1]]$bounds,candidate)
  change <- confint(d,parm="RelativeSEM",method="mean_square_region",
    assumption="normal",reference=1)
  expect_false(change$Available)
  expect_match(change$Status,"Reference interval unavailable")
})
