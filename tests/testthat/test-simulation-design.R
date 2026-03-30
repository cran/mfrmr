test_that("simulate_mfrm_data returns long-format data with truth attributes", {
  sim <- simulate_mfrm_data(
    n_person = 30,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 101
  )

  expect_true(is.data.frame(sim))
  expect_named(sim, c("Study", "Person", "Rater", "Criterion", "Score"))
  expect_equal(length(unique(sim$Person)), 30)
  expect_equal(length(unique(sim$Rater)), 4)
  expect_equal(length(unique(sim$Criterion)), 3)
  expect_true(all(sim$Score %in% 1:4))

  truth <- attr(sim, "mfrm_truth")
  expect_true(is.list(truth))
  expect_true(all(c("person", "facets", "steps") %in% names(truth)))
  expect_equal(length(truth$person), 30)
  expect_equal(length(truth$facets$Rater), 4)
  expect_equal(length(truth$facets$Criterion), 3)
})

test_that("build_mfrm_sim_spec returns reusable simulation metadata", {
  spec <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = data.frame(
      StepFacet = rep(c("C01", "C02", "C03", "C04"), each = 3),
      Step = rep(paste0("Step_", 1:3), times = 4),
      Estimate = c(-1.2, 0, 1.2, -1.0, 0.1, 1.0, -0.8, 0.2, 0.9, -1.1, 0.0, 1.1)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$assignment, "rotating")
  expect_equal(spec$model, "PCM")
  expect_true(is.data.frame(spec$threshold_table))
  expect_equal(length(unique(spec$threshold_table$StepFacet)), 4)
})

test_that("simulate_mfrm_data accepts mfrm_sim_spec with step-facet-specific thresholds", {
  spec <- build_mfrm_sim_spec(
    n_person = 16,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = data.frame(
      StepFacet = rep(c("C01", "C02", "C03", "C04"), each = 3),
      StepIndex = rep(1:3, times = 4),
      Estimate = c(-1.1, 0, 1.1, -0.9, 0.1, 1.0, -0.8, 0.2, 0.9, -1.0, 0.0, 1.2)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 515)
  truth <- attr(sim, "mfrm_truth")
  sim_spec <- attr(sim, "mfrm_simulation_spec")

  expect_true(is.data.frame(sim))
  expect_true(is.data.frame(truth$step_table))
  expect_equal(sort(unique(truth$step_table$StepFacet)), c("C01", "C02", "C03", "C04"))
  expect_equal(sim_spec$model, "PCM")
  expect_equal(sim_spec$assignment, "rotating")
})

test_that("simulate_mfrm_data uses PCM step-facet thresholds when sampling scores", {
  spec <- build_mfrm_sim_spec(
    n_person = 500,
    n_rater = 2,
    n_criterion = 2,
    raters_per_person = 2,
    score_levels = 4,
    theta_sd = 0,
    rater_sd = 0,
    criterion_sd = 0,
    noise_sd = 0,
    assignment = "crossed",
    thresholds = data.frame(
      StepFacet = rep(c("C01", "C02"), each = 3),
      StepIndex = rep(1:3, times = 2),
      Estimate = c(-1.5, -0.4, 0.4, 0.4, 1.2, 2.0)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 919)
  mean_by_criterion <- tapply(sim$Score, sim$Criterion, mean)

  expect_gt(unname(mean_by_criterion["C01"]), unname(mean_by_criterion["C02"]))
})

test_that("simulate_mfrm_data can include group-linked signals", {
  sim <- simulate_mfrm_data(
    n_person = 32,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 2,
    group_levels = c("A", "B"),
    dif_effects = data.frame(Group = "B", Criterion = "C04", Effect = 1.0),
    interaction_effects = data.frame(Rater = "R04", Criterion = "C04", Effect = -1.0),
    seed = 111
  )

  expect_true("Group" %in% names(sim))
  expect_setequal(unique(sim$Group), c("A", "B"))

  truth <- attr(sim, "mfrm_truth")
  expect_true(is.list(truth$signals))
  expect_true(is.data.frame(truth$signals$dif_effects))
  expect_true(is.data.frame(truth$signals$interaction_effects))
})

test_that("extract_mfrm_sim_spec captures fitted threshold and assignment metadata", {
  toy <- load_mfrmr_data("example_core")
  toy_people <- unique(toy$Person)[1:18]
  toy <- toy[match(toy$Person, toy_people, nomatch = 0L) > 0L, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(fit)

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$n_person, length(unique(toy$Person)))
  expect_equal(spec$n_rater, length(unique(toy$Rater)))
  expect_equal(spec$n_criterion, length(unique(toy$Criterion)))
  expect_true(spec$assignment %in% c("crossed", "rotating"))
  expect_true(is.data.frame(spec$source_summary$observed_raters_per_person))
  expect_true(is.data.frame(spec$threshold_table))
})

test_that("extract_mfrm_sim_spec can activate empirical latent support and resampled assignment", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "resampled",
    latent_distribution = "empirical"
  )

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$assignment, "resampled")
  expect_equal(spec$latent_distribution, "empirical")
  expect_true(is.list(spec$empirical_support))
  expect_true(all(c("person", "rater", "criterion") %in% names(spec$empirical_support)))
  expect_true(is.data.frame(spec$assignment_profiles))
  expect_true(all(c("TemplatePerson", "Rater") %in% names(spec$assignment_profiles)))
})

test_that("simulate_mfrm_data supports empirical latent draws and resampled assignment profiles", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "resampled",
    latent_distribution = "empirical"
  )
  spec_n30 <- simulation_override_spec_design(
    spec,
    n_person = 30,
    n_rater = spec$n_rater,
    n_criterion = spec$n_criterion,
    raters_per_person = spec$raters_per_person
  )

  sim <- simulate_mfrm_data(sim_spec = spec_n30, seed = 902)
  truth <- attr(sim, "mfrm_truth")
  sim_meta <- attr(sim, "mfrm_simulation_spec")

  expect_equal(length(unique(sim$Person)), 30)
  expect_true(is.list(truth))
  expect_equal(sim_meta$assignment, "resampled")
  expect_equal(sim_meta$latent_distribution, "empirical")
})

test_that("resampled assignment specs reject unsupported design changes", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "resampled",
    latent_distribution = "empirical"
  )

  expect_error(
    simulation_override_spec_design(
      spec,
      n_person = spec$n_person,
      n_rater = spec$n_rater + 1L,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person
    ),
    "supports changing `n_person` only",
    fixed = TRUE
  )
})

test_that("extract_mfrm_sim_spec can record an observed design skeleton", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "skeleton",
    latent_distribution = "empirical",
    source_data = toy,
    person = "Person",
    group = "Group"
  )

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$assignment, "skeleton")
  expect_true(is.data.frame(spec$design_skeleton))
  expect_true(all(c("TemplatePerson", "Rater", "Criterion") %in% names(spec$design_skeleton)))
})

test_that("simulate_mfrm_data supports observed design skeleton reuse", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "skeleton",
    latent_distribution = "empirical",
    source_data = toy,
    person = "Person",
    group = "Group"
  )
  spec_n30 <- simulation_override_spec_design(
    spec,
    n_person = 30,
    n_rater = spec$n_rater,
    n_criterion = spec$n_criterion,
    raters_per_person = spec$raters_per_person
  )

  sim <- simulate_mfrm_data(sim_spec = spec_n30, seed = 903)
  sim_meta <- attr(sim, "mfrm_simulation_spec")

  expect_equal(length(unique(sim$Person)), 30)
  expect_equal(sim_meta$assignment, "skeleton")
  expect_true(is.data.frame(sim_meta$design_skeleton))
})

test_that("observed response skeleton can carry Group and Weight metadata", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  person_groups <- setNames(rep(c("A", "B"), length.out = length(keep_people)), keep_people)
  toy$Group <- unname(person_groups[toy$Person])
  toy$Weight <- rep(c(1, 2), length.out = nrow(toy))

  fit <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      weight = "Weight",
      method = "MML",
      quad_points = 5,
      maxit = 15
    )
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "skeleton",
    latent_distribution = "empirical",
    source_data = toy,
    person = "Person",
    group = "Group"
  )
  spec_n24 <- simulation_override_spec_design(
    spec,
    n_person = 24,
    n_rater = spec$n_rater,
    n_criterion = spec$n_criterion,
    raters_per_person = spec$raters_per_person
  )

  expect_true(all(c("TemplatePerson", "Rater", "Criterion", "Group", "Weight") %in% names(spec$design_skeleton)))

  sim <- simulate_mfrm_data(sim_spec = spec_n24, seed = 904)
  sim_meta <- attr(sim, "mfrm_simulation_spec")

  expect_true(all(c("Group", "Weight") %in% names(sim)))
  expect_true(all(sim$Group %in% c("A", "B")))
  expect_true(all(sim$Weight > 0))
  expect_true(is.data.frame(sim_meta$design_skeleton))
  expect_true(all(c("Group", "Weight") %in% names(sim_meta$design_skeleton)))

  eval_obj <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = 24,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      sim_spec = spec
    )
  )
  expect_s3_class(eval_obj, "mfrm_design_evaluation")
})

test_that("extract_mfrm_sim_spec checks person-level group mapping when source_data is supplied", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  toy$Weight <- rep(c(1, 2), length.out = nrow(toy))
  toy$Group <- ifelse(seq_len(nrow(toy)) %% 2 == 0, "A", "B")

  fit <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      weight = "Weight",
      method = "MML",
      quad_points = 5,
      maxit = 15
    )
  )

  toy_bad <- toy
  toy_bad$Group[toy_bad$Person == keep_people[1]][1] <- "C"

  expect_error(
    extract_mfrm_sim_spec(
      fit,
      assignment = "skeleton",
      latent_distribution = "empirical",
      source_data = toy_bad,
      person = "Person",
      group = "Group"
    ),
    "at most one `group` label per person",
    fixed = TRUE
  )
})

test_that("skeleton assignment specs reject unsupported design changes", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "MML", quad_points = 5, maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(
    fit,
    assignment = "skeleton",
    latent_distribution = "empirical"
  )

  expect_error(
    simulation_override_spec_design(
      spec,
      n_person = spec$n_person,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion + 1L,
      raters_per_person = spec$raters_per_person
    ),
    "supports changing `n_person` only",
    fixed = TRUE
  )
})

test_that("seeded simulation helpers preserve caller RNG state", {
  set.seed(999)
  sim_seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  sim <- simulate_mfrm_data(
    n_person = 24,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 123
  )
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), sim_seed_before)
  expect_true(is.data.frame(sim))

  set.seed(1001)
  design_seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  design_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = 24,
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      maxit = 10,
      seed = 234
    )
  )
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), design_seed_before)
  expect_s3_class(design_eval, "mfrm_design_evaluation")

  set.seed(1003)
  signal_seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  signal_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      n_person = 24,
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      seed = 345
    )
  )
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), signal_seed_before)
  expect_s3_class(signal_eval, "mfrm_signal_detection")
})

test_that("design recovery metrics align location before RMSE and bias", {
  metrics <- design_eval_recovery_metrics(
    est_levels = c("L1", "L2", "L3"),
    est_values = c(0.2, 1.2, 2.2),
    truth_vec = c(L1 = 0, L2 = 1, L3 = 2)
  )

  expect_equal(metrics$raw_bias, 0.2)
  expect_equal(metrics$raw_rmse, 0.2)
  expect_equal(metrics$aligned_bias, 0)
  expect_equal(metrics$aligned_rmse, 0)
})

test_that("evaluate_mfrm_design returns usable summary and plot data", {
  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = c(30, 40),
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      maxit = 15,
      seed = 202
    )
  )

  expect_s3_class(sim_eval, "mfrm_design_evaluation")
  expect_true(is.data.frame(sim_eval$results))
  expect_true(is.data.frame(sim_eval$rep_overview))
  expect_true(all(c("Person", "Rater", "Criterion") %in% unique(sim_eval$results$Facet)))
  expect_true(all(c("SeverityRMSERaw", "SeverityBiasRaw") %in% names(sim_eval$results)))
  expect_true(all(c("GeneratorModel", "GeneratorStepFacet", "FitModel", "FitStepFacet",
                    "RecoveryComparable", "RecoveryBasis") %in% names(sim_eval$results)))
  expect_true(all(sim_eval$results$SeverityRMSE <= sim_eval$results$SeverityRMSERaw | is.na(sim_eval$results$SeverityRMSERaw)))

  s <- summary(sim_eval)
  expect_s3_class(s, "summary.mfrm_design_evaluation")
  expect_true(is.data.frame(s$overview))
  expect_true(is.data.frame(s$design_summary))
  expect_true(all(c("Facet", "MeanSeparation", "MeanSeverityRMSE", "ConvergenceRate",
                    "McseSeparation", "McseSeverityRMSE", "McseConvergenceRate") %in% names(s$design_summary)))
  expect_true(all(c("MeanSeverityRMSERaw", "MeanSeverityBiasRaw") %in% names(s$design_summary)))
  expect_true(all(c("RecoveryComparableRate", "RecoveryBasis") %in% names(s$design_summary)))
  expect_true(is.list(s$ademp))
  expect_true(all(c("aims", "data_generating_mechanism", "estimands", "methods", "performance_measures") %in% names(s$ademp)))

  p <- plot(sim_eval, facet = "Rater", metric = "separation", x_var = "n_person", draw = FALSE)
  expect_true(is.list(p))
  expect_true(is.data.frame(p$data))
  expect_equal(p$facet, "Rater")
  expect_equal(p$metric_col, "MeanSeparation")
})

test_that("recommend_mfrm_design returns threshold tables", {
  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = c(30, 50),
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      reps = 1,
      maxit = 15,
      seed = 303
    )
  )

  rec <- recommend_mfrm_design(
    sim_eval,
    min_separation = 1.5,
    min_reliability = 0.7,
    max_severity_rmse = 1.5,
    max_misfit_rate = 0.5,
    min_convergence_rate = 0
  )

  expect_true(is.list(rec))
  expect_true(is.data.frame(rec$facet_table))
  expect_true(is.data.frame(rec$design_table))
  expect_true(all(c("Pass", "MinSeparation", "MaxSeverityRMSE") %in% names(rec$design_table)))
  expect_true(all(c("SeparationPass", "ReliabilityPass", "Pass") %in% names(rec$facet_table)))
})

test_that("evaluate_mfrm_design accepts sim_spec and carries ADEMP metadata", {
  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating"
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = c(18, 20),
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 808
    )
  )

  expect_s3_class(sim_eval, "mfrm_design_evaluation")
  expect_true(inherits(sim_eval$settings$sim_spec, "mfrm_sim_spec"))
  expect_true(is.list(sim_eval$ademp))
  expect_equal(sim_eval$ademp$data_generating_mechanism$source, "manual")
  expect_equal(sim_eval$ademp$data_generating_mechanism$assignment, "rotating")
  expect_equal(sim_eval$settings$recovery_comparable, TRUE)
})

test_that("evaluate_mfrm_design carries PCM step_facet into fitted recovery contract", {
  spec <- build_mfrm_sim_spec(
    n_person = 24,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = data.frame(
      StepFacet = rep(c("C01", "C02", "C03", "C04"), each = 3),
      StepIndex = rep(1:3, times = 4),
      Estimate = c(-1.1, 0, 1.1, -0.9, 0.1, 1.0, -0.8, 0.2, 0.9, -1.0, 0.0, 1.2)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = 24,
      n_rater = 3,
      n_criterion = 4,
      raters_per_person = 2,
      reps = 1,
      model = "PCM",
      maxit = 10,
      sim_spec = spec,
      seed = 810
    )
  )

  expect_equal(sim_eval$settings$step_facet, "Criterion")
  expect_equal(unique(sim_eval$results$FitStepFacet), "Criterion")
  expect_true(all(sim_eval$results$RecoveryComparable))
})

test_that("evaluate_mfrm_design suppresses recovery metrics when generator and fit contracts differ", {
  spec <- build_mfrm_sim_spec(
    n_person = 36,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    model = "RSM"
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = 36,
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      reps = 1,
      model = "PCM",
      step_facet = "Criterion",
      maxit = 10,
      sim_spec = spec,
      seed = 811
    )
  )

  expect_true(all(!sim_eval$results$RecoveryComparable))
  expect_true(all(sim_eval$results$RecoveryBasis == "generator_fit_model_mismatch"))
  expect_true(all(is.na(sim_eval$results$SeverityRMSE)))
  expect_true(all(is.na(sim_eval$results$SeverityBias)))
})

test_that("evaluate_mfrm_design rejects incompatible step-facet count changes under sim_spec", {
  spec <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = data.frame(
      StepFacet = rep(c("C01", "C02", "C03", "C04"), each = 3),
      StepIndex = rep(1:3, times = 4),
      Estimate = c(-1.1, 0, 1.1, -0.9, 0.1, 1.0, -0.8, 0.2, 0.9, -1.0, 0.0, 1.2)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  expect_error(
    evaluate_mfrm_design(
      n_person = 18,
      n_rater = 3,
      n_criterion = 5,
      raters_per_person = 2,
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 809
    ),
    "design-specific simulation specification",
    fixed = TRUE
  )
})

test_that("evaluate_mfrm_signal_detection returns usable detection summaries", {
  sig_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      n_person = c(36, 48),
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      reps = 1,
      dif_effect = 1.2,
      bias_effect = -1.2,
      maxit = 15,
      seed = 404
    )
  )

  expect_s3_class(sig_eval, "mfrm_signal_detection")
  expect_true(is.data.frame(sig_eval$results))
  expect_true(is.data.frame(sig_eval$rep_overview))
  expect_true(all(c("DIFDetected", "BiasDetected", "BiasScreenMetricAvailable",
                    "DIFFalsePositiveRate", "BiasScreenFalsePositiveRate") %in%
                    names(sig_eval$results)))

  s_sig <- summary(sig_eval)
  expect_s3_class(s_sig, "summary.mfrm_signal_detection")
  expect_true(is.data.frame(s_sig$overview))
  expect_true(is.data.frame(s_sig$detection_summary))
  expect_true(all(c("DIFPower", "BiasScreenRate",
                    "BiasScreenFalsePositiveRate",
                    "BiasScreenMetricAvailabilityRate",
                    "McseDIFPower", "McseBiasScreenRate",
                    "MeanTargetContrast", "MeanTargetBias") %in%
                    names(s_sig$detection_summary)))
  expect_true(is.list(s_sig$ademp))

  p_sig <- plot(sig_eval, signal = "dif", metric = "power", x_var = "n_person", draw = FALSE)
  expect_true(is.list(p_sig))
  expect_true(is.data.frame(p_sig$data))
  expect_equal(p_sig$metric_col, "DIFPower")
  expect_equal(p_sig$display_metric, "DIF target-flag rate")
  expect_match(p_sig$interpretation_note, "DIF-side rates summarize target/non-target flagging behavior", fixed = TRUE)

  p_sig_bias <- plot(sig_eval, signal = "bias", metric = "power", x_var = "n_person", draw = FALSE)
  expect_equal(p_sig_bias$metric_col, "BiasScreenRate")
  expect_equal(p_sig_bias$display_metric, "Bias screening hit rate")
  expect_match(p_sig_bias$interpretation_note, "not formal inferential power or alpha estimates", fixed = TRUE)

  expect_true(any(sig_eval$results$DIFDetected, na.rm = TRUE))
  expect_true(all(is.finite(s_sig$detection_summary$BiasScreenMetricAvailabilityRate)))
  expect_true(any(grepl("Bias-side rates are screening summaries", s_sig$notes, fixed = TRUE)))
})

test_that("evaluate_mfrm_signal_detection accepts sim_spec and keeps signal injection explicit", {
  spec <- build_mfrm_sim_spec(
    n_person = 24,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    group_levels = c("A", "B")
  )

  sig_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      n_person = c(24, 28),
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      reps = 1,
      dif_effect = 1.0,
      bias_effect = -1.0,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 810
    )
  )

  expect_s3_class(sig_eval, "mfrm_signal_detection")
  expect_true(inherits(sig_eval$settings$sim_spec, "mfrm_sim_spec"))
  expect_true(is.list(sig_eval$ademp))
  expect_equal(sig_eval$ademp$data_generating_mechanism$source, "manual")
  expect_true(all(sig_eval$results$BiasTargetCriterion %in% sprintf("C%02d", 1:4)))
})

test_that("predict_mfrm_population returns scenario-level forecast from sim_spec", {
  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating"
  )

  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      n_person = 24,
      reps = 1,
      maxit = 10,
      seed = 811
    )
  )

  expect_s3_class(pred, "mfrm_population_prediction")
  expect_true(is.data.frame(pred$forecast))
  expect_true(is.data.frame(pred$overview))
  expect_true(inherits(pred$sim_spec, "mfrm_sim_spec"))
  expect_equal(pred$sim_spec$n_person, 24L)
  expect_true(is.list(pred$ademp))
  expect_equal(pred$settings$source, "mfrm_sim_spec")

  s_pred <- summary(pred)
  expect_s3_class(s_pred, "summary.mfrm_population_prediction")
  expect_true(is.data.frame(s_pred$forecast))
  expect_true(all(c("Facet", "MeanSeparation", "McseSeparation") %in% names(s_pred$forecast)))
})

test_that("predict_mfrm_population can derive its specification from a fitted model", {
  toy <- load_mfrmr_data("example_core")
  toy_people <- unique(toy$Person)[1:18]
  toy <- toy[match(toy$Person, toy_people, nomatch = 0L) > 0L, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 15)
  )

  pred <- suppressWarnings(
    predict_mfrm_population(
      fit = fit,
      n_person = 20,
      reps = 1,
      maxit = 10,
      seed = 812
    )
  )

  expect_s3_class(pred, "mfrm_population_prediction")
  expect_equal(pred$settings$source, "fit_mfrm")
  expect_equal(pred$sim_spec$source, "fit_mfrm")
  expect_true(all(pred$forecast$Facet %in% c("Person", "Rater", "Criterion")))
})

test_that("predict_mfrm_population requires exactly one source", {
  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating"
  )
  toy <- load_mfrmr_data("example_core")
  toy_people <- unique(toy$Person)[1:18]
  toy <- toy[match(toy$Person, toy_people, nomatch = 0L) > 0L, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 15)
  )

  expect_error(
    predict_mfrm_population(sim_spec = spec, fit = fit, reps = 1, seed = 813),
    "exactly one",
    fixed = TRUE
  )
})
