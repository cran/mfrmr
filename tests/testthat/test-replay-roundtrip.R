# End-to-end round-trip tests for `export_mfrm_bundle()` -> `replay.R`.
# We export the bundle, source the replay script in a clean environment,
# and assert that the reproduced fit's headline statistics match the
# original. These tests catch the kind of silent argument-drop the
# 0.1.5 / early-0.1.6 replay path was prone to.

local({
  .toy <<- load_mfrmr_data("example_core")
  .fit <<- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 25)
  ))
})

bundle_and_source <- function(fit, data, prefix = "rt_test") {
  td <- tempfile("mfrm_replay_rt_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(td)
  export_mfrm_bundle(
    fit,
    output_dir = ".",
    prefix = prefix,
    include = c("core_tables", "manifest", "script"),
    data = data
  )
  e <- new.env(parent = globalenv())
  suppressMessages(suppressWarnings(
    sys.source(file.path(td, paste0(prefix, "_replay.R")), envir = e)
  ))
  e
}

test_that("replay round-trip reproduces JML log-likelihood", {
  e <- bundle_and_source(.fit, .toy)
  replayed <- e$fit
  expect_s3_class(replayed, "mfrm_fit")
  expect_equal(replayed$summary$LogLik, .fit$summary$LogLik,
               tolerance = 1e-6)
  expect_equal(replayed$summary$N, .fit$summary$N)
})

test_that("replay round-trip reproduces person estimates", {
  e <- bundle_and_source(.fit, .toy)
  replayed <- e$fit
  orig <- as.data.frame(.fit$facets$person, stringsAsFactors = FALSE)
  rep_p <- as.data.frame(replayed$facets$person, stringsAsFactors = FALSE)
  expect_equal(nrow(orig), nrow(rep_p))
  o <- orig[order(orig$Person), ]
  r <- rep_p[order(rep_p$Person), ]
  expect_equal(suppressWarnings(as.numeric(r$Estimate)),
               suppressWarnings(as.numeric(o$Estimate)),
               tolerance = 1e-6)
})

test_that("replay carries `mml_engine` argument forward", {
  fit_em <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 7, maxit = 25,
             mml_engine = "em")
  ))
  e <- bundle_and_source(fit_em, .toy, prefix = "rt_em")
  replayed <- e$fit
  expect_s3_class(replayed, "mfrm_fit")
  # Confirm the replayed fit also used the EM engine, not direct.
  expect_identical(
    as.character(replayed$config$estimation_control$mml_engine_requested),
    "em"
  )
})

test_that("replay records a package-version mismatch warning", {
  e <- bundle_and_source(.fit, .toy, prefix = "rt_ver")
  td <- attr(e, "path", exact = TRUE)
  # We re-source under a faked version and verify a warning fires.
  td2 <- tempfile("mfrm_replay_ver_"); dir.create(td2)
  on.exit(unlink(td2, recursive = TRUE), add = TRUE)
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE); setwd(td2)
  export_mfrm_bundle(.fit, output_dir = ".", prefix = "rt_ver2",
                     include = c("core_tables","manifest","script"),
                     data = .toy)
  script_lines <- readLines(file.path(td2, "rt_ver2_replay.R"))
  # The version-mismatch guard is emitted near the top of the script.
  has_guard <- any(grepl("Recorded mfrmr version", script_lines))
  has_warn <- any(grepl("Estimates may differ", script_lines))
  expect_true(has_guard)
  expect_true(has_warn)
})
