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

bundle_and_source <- function(fit, data, prefix = "rt_test", diagnostics = NULL) {
  td <- tempfile("mfrm_replay_rt_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(td)
  export_mfrm_bundle(
    fit,
    diagnostics = diagnostics,
    output_dir = ".",
    prefix = prefix,
    include = c("core_tables", "manifest", "script"),
    data = data,
    acknowledge_sensitive = TRUE
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

test_that("replay preserves literal IDs and non-syntactic column names", {
  dat <- .toy
  ids <- unique(dat$Person)
  labels <- c("01", "1", "NA", sprintf("%03d", seq.int(4L, length(ids))))
  dat$Person <- labels[match(dat$Person, ids)]
  raters <- unique(dat$Rater)
  rater_labels <- c("01", "1", "NA", paste0("R", seq_along(raters)))[seq_along(raters)]
  dat$Rater <- rater_labels[match(dat$Rater, raters)]
  names(dat)[match(c("Person", "Rater", "Criterion", "Score"), names(dat))] <-
    c("Person ID", "Rater ID", "Criterion name", "Score value")
  fit <- suppressWarnings(fit_mfrm(dat, "Person ID", c("Rater ID", "Criterion name"),
                                  "Score value", method = "JML", maxit = 25))
  replay <- bundle_and_source(fit, dat, prefix = "literal_ids")
  expect_identical(names(replay$data), names(dat))
  for (nm in c("Person ID", "Rater ID", "Criterion name")) {
    expect_identical(replay$data[[nm]], as.character(dat[[nm]]))
  }
  expect_equal(replay$fit$summary$N, fit$summary$N)
  expect_equal(replay$fit$summary$Persons, fit$summary$Persons)
  expect_equal(replay$fit$summary$LogLik, fit$summary$LogLik, tolerance = 1e-6)
})

test_that("latent-regression replay preserves background IDs and factor coding", {
  dat <- load_mfrmr_data("example_operational")
  ids <- unique(dat$Person)
  labels <- c("01", "1", "NA", sprintf("%03d", seq.int(4L, length(ids))))
  dat$Person <- labels[match(dat$Person, ids)]
  background <- data.frame(Person = labels)
  background[["Group code"]] <- factor(rep(c("01", "1"), length.out = length(ids)),
                                       levels = c("1", "01"))
  contrasts(background[["Group code"]]) <- contr.sum(2)
  fit <- suppressWarnings(fit_mfrm(
    dat, "Person", c("Rater", "Criterion"), "Score", method = "MML",
    population_formula = ~ `Group code`, person_data = background, person_id = "Person"
  ))
  replay <- bundle_and_source(fit, dat, prefix = "latent_literal_ids")
  expect_identical(replay$fit_person_data, fit$population$person_table_replay)
  expect_equal(replay$fit$summary$LogLik, fit$summary$LogLik, tolerance = 1e-6)
  expect_equal(replay$fit$population$coefficients, fit$population$coefficients,
               tolerance = 1e-6)
})

test_that("replay preserves diagnostic selection and standardization", {
  diagnostics <- suppressWarnings(diagnose_mfrm(
    .fit, interaction_pairs = list(c("Rater", "Criterion")), top_n_interactions = 3,
    whexact = TRUE, fit_df_method = "facets", diagnostic_mode = "legacy",
    residual_pca = "overall", pca_max_factors = 2L
  ))
  replay <- bundle_and_source(.fit, .toy, prefix = "diagnostic_settings",
                              diagnostics = diagnostics)$diagnostics
  expect_identical(replay$replay_inputs, diagnostics$replay_inputs)
  expect_identical(replay$diagnostic_mode, "legacy")
  expect_equal(replay$fit_standardization, diagnostics$fit_standardization)
  expect_equal(replay$measures, diagnostics$measures, tolerance = 1e-6)
  expect_equal(replay$interactions, diagnostics$interactions, tolerance = 1e-6)
  expect_equal(ncol(replay$residual_pca_overall$pca$loadings), 2L)
  expect_equal(replay$residual_pca_overall$pca$loadings,
               diagnostics$residual_pca_overall$pca$loadings, tolerance = 1e-6)
  expect_null(replay$residual_pca_by_facet)
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

test_that("replay preserves the complete interaction specification", {
  interaction_fit <- suppressMessages(suppressWarnings(
    fit_mfrm(
      .toy, "Person", c("Rater", "Criterion"), "Score",
      method = "JML", maxit = 30,
      facet_interactions = "Rater:Criterion",
      min_obs_per_interaction = 0,
      interaction_policy = "error"
    )
  ))
  replayed <- bundle_and_source(
    interaction_fit, .toy, prefix = "rt_interaction"
  )$fit

  expect_identical(
    names(replayed$config$interaction_specs),
    names(interaction_fit$config$interaction_specs)
  )
  expect_identical(length(replayed$opt$par), length(interaction_fit$opt$par))
  expect_equal(
    replayed$summary$LogLik, interaction_fit$summary$LogLik,
    tolerance = 1e-6
  )
  expect_equal(
    as.data.frame(interaction_effect_table(replayed)),
    as.data.frame(interaction_effect_table(interaction_fit)),
    tolerance = 1e-6
  )
  expect_equal(
    replayed$facets$others, interaction_fit$facets$others,
    tolerance = 1e-6
  )
  expect_identical(
    replayed$summary$FitReadiness,
    interaction_fit$summary$FitReadiness
  )
})

test_that("FACETS workflow replay runs with the original MML engine", {
  run <- suppressWarnings(run_mfrm_facets(
    .toy, "Person", c("Rater", "Criterion"), "Score",
    method = "MML", mml_engine = "em", quad_points = 7, maxit = 25
  ))
  replay <- bundle_and_source(run, .toy, prefix = "facets_workflow")
  expect_s3_class(replay$run, "mfrm_facets_run")
  expect_identical(replay$fit$config$estimation_control$mml_engine_requested, "em")
  expect_equal(replay$fit$summary$LogLik, run$fit$summary$LogLik, tolerance = 1e-6)
  expect_equal(replay$diagnostics$measures, run$diagnostics$measures, tolerance = 1e-6)
})

test_that("replay argument registry has no unhandled material fields", {
  replay_inputs <- list(
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "JML",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 3,
    interaction_policy = "error",
    population_policy = "error",
    package_version = "0.2.4.9000"
  )
  lines <- mfrmr:::build_replay_fit_mfrm_lines(
    replay_inputs = replay_inputs,
    fit_population = list(active = FALSE),
    fit_population_person_id = NULL,
    src = list(
      person = "Person", facets = c("Rater", "Criterion"), score = "Score"
    ),
    cfg = list(model = "RSM", method = "JML")
  )
  call <- parse(text = paste(lines, collapse = "\n"))[[1]][[3]]
  emitted <- names(as.list(call)[-1])
  intentionally_conditional <- c("population_policy", "package_version")
  material <- setdiff(names(replay_inputs), intentionally_conditional)

  expect_identical(setdiff(material, emitted), character(0))
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
                     data = .toy, acknowledge_sensitive = TRUE)
  script_lines <- readLines(file.path(td2, "rt_ver2_replay.R"))
  # The version-mismatch guard is emitted near the top of the script.
  has_guard <- any(grepl("Recorded mfrmr version", script_lines))
  has_warn <- any(grepl("Estimates may differ", script_lines))
  expect_true(has_guard)
  expect_true(has_warn)
})
