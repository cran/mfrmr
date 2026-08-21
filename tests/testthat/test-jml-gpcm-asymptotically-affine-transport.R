fit_jml_gpcm_affine_transport_fixture <- function(
    kind = c("positive", "negative")) {
  kind <- match.arg(kind)
  if (identical(kind, "positive")) {
    data <- expand.grid(
      Person = c("P1", "P2"),
      Criterion = c("C1", "C2"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data$Score <- c(0L, 1L, 1L, 0L)
    data <- data[rep(seq_len(nrow(data)), each = 3L), , drop = FALSE]
  } else {
    data <- expand.grid(
      Person = c("P1", "P2"),
      Criterion = c("C1", "C2", "C3"),
      Score = 0:1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data <- data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE]
  }
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    rating_min = 0,
    rating_max = 1,
    maxit = 300,
    reltol = 1e-12
  ))
}

gpcm_affine_transport_log_probability <- function(utility, score) {
  maximum <- max(utility)
  utility[score + 1L] -
    (maximum + log(sum(exp(utility - maximum))))
}

fit_jml_gpcm_affine_slope_transport_fixture <- function() {
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- c(0L, 1L, 1L, 0L)
  data <- data[rep(seq_len(nrow(data)), each = 3L), , drop = FALSE]
  anchors <- data.frame(
    Facet = "Person",
    Level = c("P1", "P2"),
    Anchor = c(-1, 1),
    stringsAsFactors = FALSE
  )
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    anchors = anchors,
    rating_min = 0,
    rating_max = 1,
    maxit = 300,
    reltol = 1e-12
  ))
}

test_that("positive production certificates transport to vanishing residuals", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_affine_transport_fixture("positive")
  audit <-
    fit$config$boundary_audit$gpcm_asymptotically_affine_transport

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-asymptotically-affine-transport-0.2.3-v1"
  )
  expect_identical(
    audit$state, "certified_vanishing_residual_curved_neighborhood"
  )
  expect_true(audit$evaluated)
  expect_gt(audit$source_path_count, 0L)
  expect_gt(audit$transported_path_count, 0L)
  expect_true(all(audit$path_registry$TransportCertified))
  expect_true(all(
    audit$path_registry$SourceFamily == "joint_ordered_pair"
  ))
  expect_true(audit$same_boundary_limit_certified)
  expect_true(audit$positive_transport_theorem_applied)
  expect_true(audit$vanishing_residual_curved_neighborhood_classified)
  expect_false(audit$bounded_nonvanishing_residual_classified)
  expect_false(audit$rate_nonconvergent_path_classified)
  expect_false(audit$general_curved_path_classified)
  expect_false(audit$curved_path_absence_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_match(audit$limitations, "nonvanishing residuals", fixed = TRUE)
})

test_that("negative constant-rate completion makes no curved-path claim", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_affine_transport_fixture("negative")
  audit <-
    fit$config$boundary_audit$gpcm_asymptotically_affine_transport

  expect_identical(audit$state, "no_positive_path_to_transport")
  expect_true(audit$evaluated)
  expect_identical(audit$source_path_count, 0L)
  expect_identical(audit$transported_path_count, 0L)
  expect_false(audit$same_boundary_limit_certified)
  expect_false(audit$positive_transport_theorem_applied)
  expect_false(audit$curved_path_absence_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_match(
    audit$reason_codes, "curved_path_absence_not_certified", fixed = TRUE
  )
})

test_that("slope-only certificates use the same positive transport theorem", {
  fit <- fit_jml_gpcm_affine_transport_fixture("positive")
  config <- fit$config
  classification <-
    config$boundary_audit$gpcm_fixed_objective_classification
  classification$state <- "certified_slope_only_recession"
  classification$positive_boundary_certificate_present <- TRUE
  slope <- list(
    state = "certified_monotone_boundary_path",
    scope_complete = TRUE,
    certificates = data.frame(
      PairId = "SP00001",
      StrictRows = 1L,
      BoundaryLogLikelihood = -2,
      BoundaryImprovement = 0.5,
      Certified = TRUE,
      stringsAsFactors = FALSE
    )
  )
  joint <- list(
    state = "none_certified",
    scope_complete = TRUE,
    rate_scope_complete = TRUE,
    certificates = data.frame(),
    rate_certificates = data.frame()
  )
  audit <- mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
    config, slope, joint, classification
  )

  expect_identical(
    audit$state, "certified_vanishing_residual_curved_neighborhood"
  )
  expect_identical(audit$transported_path_count, 1L)
  expect_identical(
    audit$path_registry$SourceFamily, "slope_only_constant_rate"
  )
  expect_true(audit$path_registry$StrictCertificate)
  expect_true(audit$path_registry$TransportCertified)
})

test_that("production slope-only paths are transportable", {
  fit <- fit_jml_gpcm_affine_slope_transport_fixture()
  audit <-
    fit$config$boundary_audit$gpcm_asymptotically_affine_transport

  expect_identical(
    audit$state, "certified_vanishing_residual_curved_neighborhood"
  )
  expect_identical(audit$transported_path_count, 1L)
  expect_identical(
    audit$path_registry$SourceFamily, "slope_only_constant_rate"
  )
  expect_true(audit$path_registry$StrictCertificate)
  expect_true(audit$path_registry$TransportCertified)
  expect_true(is.finite(audit$path_registry$BoundaryLogLikelihood))
})

test_that("vanishing curved residuals preserve the general-rate boundary", {
  rates <- c(3, -1, -2)
  slope_index <- c(1L, 2L, 3L, 3L)
  score <- rep(2L, 4L)
  direction <- rbind(
    c(0, 1, 2),
    c(0, -1, -2),
    c(0, 1, 2),
    c(0, -1, -2)
  )
  residual_direction <- rbind(
    c(0, -0.2, 0.1),
    c(0, 0.1, -0.1),
    c(0, -0.1, 0.2),
    c(0, 0.2, -0.2)
  )
  distances <- c(4, 8, 16, 32)
  curved <- vapply(distances, function(distance) {
    residual <- sin(distance) / (1 + distance)
    log_slope_residual <- residual * c(1, -2, 1)
    slopes <- exp(rates * distance + log_slope_residual)
    utilities <- distance * direction + residual * residual_direction
    sum(vapply(seq_len(4L), function(observation) {
      gpcm_affine_transport_log_probability(
        slopes[slope_index[observation]] * utilities[observation, ],
        score[observation]
      )
    }, numeric(1)))
  }, numeric(1))
  boundary <- -3 * log(3)

  expect_equal(sum(rates), 0)
  expect_gt(curved[length(curved)], curved[1])
  expect_lt(abs(curved[length(curved)] - boundary), 1e-10)
  expect_lt(
    abs(curved[length(curved)] - boundary),
    abs(curved[1] - boundary)
  )
})

test_that("nonvanishing zero-rate residuals can destroy a common limit", {
  rates <- c(1, 0, -1)
  slope_index <- 1:3
  score <- rep(2L, 3L)
  base <- rbind(c(0, 0, 0), c(0, 1, 2), c(0, 0, 0))
  direction <- rbind(c(0, 1, 2), c(0, 0, 0), c(0, 0, 0))
  n <- 8
  distances <- c(2 * pi * n + pi / 2, 2 * pi * n + 3 * pi / 2)
  observed <- vapply(distances, function(distance) {
    oscillation <- sin(distance)
    residual <- oscillation * c(-0.25, 0.5, -0.25)
    slopes <- exp(rates * distance + residual)
    utilities <- base + distance * direction
    sum(vapply(seq_len(3L), function(observation) {
      gpcm_affine_transport_log_probability(
        slopes[slope_index[observation]] * utilities[observation, ],
        score[observation]
      )
    }, numeric(1)))
  }, numeric(1))
  expected <- vapply(c(1, -1), function(oscillation) {
    zero_slope <- exp(0.5 * oscillation)
    zero_log_probability <- gpcm_affine_transport_log_probability(
      zero_slope * c(0, 1, 2), 2L
    )
    zero_log_probability - log(3)
  }, numeric(1))

  expect_equal(sum(rates), 0)
  expect_equal(observed, expected, tolerance = 1e-10)
  expect_gt(abs(diff(observed)), 0.1)
})

test_that("source incoherence and estimator scope fail closed", {
  fit <- fit_jml_gpcm_affine_transport_fixture("positive")
  config <- fit$config
  slope <- config$boundary_audit$gpcm_slope_boundary
  joint <- config$boundary_audit$gpcm_joint_boundary
  classification <-
    config$boundary_audit$gpcm_fixed_objective_classification

  incoherent_joint <- joint
  incoherent_joint$certificates$BoundaryLogLikelihood[
    incoherent_joint$certificates$Certified
  ] <- NA_real_
  incoherent <-
    mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
      config, slope, incoherent_joint, classification
    )
  expect_identical(incoherent$state, "indeterminate_source_certificate")
  expect_false(incoherent$evaluated)
  expect_false(incoherent$same_boundary_limit_certified)

  mismatch <- classification
  mismatch$contract_version <- "mismatched"
  mismatched <- mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
    config, slope, joint, mismatch
  )
  expect_identical(mismatched$state, "not_evaluated_objective_identity")
  expect_false(mismatched$evaluated)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
    mml_config, slope, joint, classification
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
    pcm_config, slope, joint, classification
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$n_params <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_asymptotically_affine_transport(
    unit_config, slope, joint, classification
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$general_curved_path_classified)
})
