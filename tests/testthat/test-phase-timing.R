# --------------------------------------------------------------------------
# test-phase-timing.R
# Internal validation timing must remain opt-in and decision-nonintervening.
# --------------------------------------------------------------------------

test_that("phase timing is opt-in diagnostic metadata", {
  old_opt <- options(mfrmr.phase_timing = NULL)
  on.exit(options(old_opt), add = TRUE)

  d <- mfrmr:::sample_mfrm_data(seed = 42050)

  options(mfrmr.phase_timing = FALSE)
  fit_regular <- suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "JML",
    maxit = 20
  ))
  expect_null(fit_regular$config$phase_timing)

  options(mfrmr.phase_timing = TRUE)
  fit_timed <- suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "JML",
    maxit = 20
  ))

  timing <- fit_timed$config$phase_timing
  expected_phases <- c(
    "prepare_data",
    "resolve_structure",
    "build_design_config",
    "category_support_audit",
    "estimability_audit",
    "data_review_and_guards",
    "optimization",
    "post_optimization_audits",
    "parameter_and_person_tables",
    "person_boundary_audit",
    "structural_recession_audit",
    "joint_recession_audit",
    "gpcm_slope_boundary_audit",
    "gpcm_joint_boundary_audit",
    "boundary_application",
    "readiness_assembly",
    "remaining_output_tables",
    "mfrm_estimate_total"
  )

  expect_s3_class(timing, "data.frame")
  expect_identical(names(timing), c(
    "Order", "Phase", "ElapsedSeconds", "Clock", "Scope", "DecisionUse"
  ))
  expect_identical(timing$Order, seq_along(expected_phases))
  expect_identical(timing$Phase, expected_phases)
  expect_true(all(is.finite(timing$ElapsedSeconds)))
  expect_true(all(timing$ElapsedSeconds >= 0))
  expect_true(all(timing$Clock == "proc.time.elapsed"))
  expect_true(all(timing$DecisionUse == "diagnostic_only"))
  expect_gte(
    timing$ElapsedSeconds[timing$Phase == "mfrm_estimate_total"],
    max(timing$ElapsedSeconds[timing$Phase != "mfrm_estimate_total"])
  )
})

test_that("phase timing does not change estimates, audits, or readiness", {
  old_opt <- options(mfrmr.phase_timing = NULL)
  on.exit(options(old_opt), add = TRUE)

  d <- mfrmr:::sample_mfrm_data(seed = 42051)
  fit_once <- function(timed) {
    options(mfrmr.phase_timing = timed)
    suppressWarnings(fit_mfrm(
      d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      model = "PCM",
      step_facet = "Criterion",
      method = "JML",
      maxit = 20
    ))
  }

  fit_regular <- fit_once(FALSE)
  fit_timed <- fit_once(TRUE)

  normalize_elapsed <- function(x) {
    x$config$phase_timing <- NULL
    if (is.data.frame(x$opt$optimizer_polish$Stages) &&
        "ElapsedSeconds" %in% names(x$opt$optimizer_polish$Stages)) {
      x$opt$optimizer_polish$Stages$ElapsedSeconds <- 0
    }
    if (is.data.frame(x$config$estimation_control$optimizer_polish$Stages) &&
        "ElapsedSeconds" %in%
          names(x$config$estimation_control$optimizer_polish$Stages)) {
      x$config$estimation_control$optimizer_polish$Stages$ElapsedSeconds <- 0
    }
    x
  }
  fit_regular_normalized <- normalize_elapsed(fit_regular)
  fit_timed_normalized <- normalize_elapsed(fit_timed)

  expect_identical(fit_regular$summary, fit_timed$summary)
  expect_identical(fit_regular$facets, fit_timed$facets)
  expect_identical(fit_regular$interactions, fit_timed$interactions)
  expect_identical(fit_regular$steps, fit_timed$steps)
  expect_identical(fit_regular$slopes, fit_timed$slopes)
  expect_identical(fit_regular$readiness, fit_timed$readiness)
  expect_identical(fit_regular$data_review, fit_timed$data_review)
  expect_identical(fit_regular_normalized$opt, fit_timed_normalized$opt)
  expect_identical(fit_regular_normalized$config, fit_timed_normalized$config)
})
