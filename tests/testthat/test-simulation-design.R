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

test_that("build_mfrm_sim_spec accepts compact step-facet threshold shortcuts", {
  spec_list <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = list(
      C01 = c(-1.2, 0.0, 1.1),
      C02 = c(-0.9, 0.2, 1.0),
      C03 = c(-0.8, 0.3, 1.2)
    ),
    model = "PCM",
    step_facet = "Criterion"
  )

  expect_equal(sort(unique(spec_list$threshold_table$StepFacet)), c("C01", "C02", "C03"))
  expect_equal(spec_list$threshold_table$Estimate[spec_list$threshold_table$StepFacet == "C02"], c(-0.9, 0.2, 1.0))

  mat <- rbind(
    C01 = c(-1.2, 0.0, 1.1),
    C02 = c(-0.9, 0.2, 1.0),
    C03 = c(-0.8, 0.3, 1.2)
  )
  spec_matrix <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating",
    thresholds = mat,
    model = "PCM",
    step_facet = "Criterion"
  )

  expect_equal(spec_matrix$threshold_table, spec_list$threshold_table)
  expect_silent(simulate_mfrm_data(sim_spec = spec_matrix, seed = 2026))
})

test_that("build_mfrm_sim_spec accepts custom public facet names", {
  spec <- build_mfrm_sim_spec(
    n_person = 16,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  expect_identical(unname(spec$facet_names), c("Judge", "Task"))
  expect_identical(spec$facet_levels$rater, c("J01", "J02", "J03"))
  expect_identical(spec$facet_levels$criterion, c("T01", "T02"))
  expect_true(is.list(spec$planning_scope))
  expect_identical(spec$planning_scope$planner_contract, "role_based_two_non_person_facets")
  expect_identical(spec$planning_scope$role_labels, c("Judge", "Task"))
  expect_false(spec$planning_scope$supports_arbitrary_facet_planning)
  expect_true(spec$planning_scope$supports_arbitrary_facet_estimation)
  expect_identical(spec$planning_scope$future_planner_contract, "arbitrary_facet_planning_scaffold")
  expect_identical(spec$planning_scope$future_planner_stage, "schema_only")
  expect_identical(spec$planning_scope$future_branch_input_contract, "design$facets(named counts)")
  expect_true(is.data.frame(spec$planning_scope$facet_manifest))
  expect_identical(spec$planning_scope$facet_manifest$facet, c("Person", "Judge", "Task"))
  expect_identical(spec$planning_scope$facet_manifest$planning_count_alias, c("n_person", "n_judge", "n_task"))
  expect_equal(spec$planning_scope$facet_manifest$level_count, c(16, 3, 2))
  expect_true(is.list(spec$planning_constraints))
  expect_identical(spec$planning_constraints$mutable_design_variables, c("n_person", "n_rater", "n_criterion", "raters_per_person"))
  expect_identical(spec$planning_constraints$locked_design_variables, character(0))
  expect_true(is.list(spec$planning_schema))
  expect_identical(spec$planning_schema$planner_contract, "role_based_two_non_person_facets")
  expect_identical(spec$planning_schema$future_planner_contract, "arbitrary_facet_planning_scaffold")
  expect_identical(spec$planning_schema$future_branch_input_contract, "design$facets(named counts)")
  expect_true(is.data.frame(spec$planning_schema$facet_manifest))
  expect_identical(spec$planning_schema$facet_manifest$facet, c("Person", "Judge", "Task"))
  expect_true(is.data.frame(spec$planning_schema$future_facet_table))
  expect_identical(spec$planning_schema$future_facet_table$future_facet_key, c("person", "judge", "task"))
  expect_identical(spec$planning_schema$future_facet_table$current_planning_count_alias, c("n_person", "n_judge", "n_task"))
  expect_identical(spec$planning_schema$future_facet_table$future_axis_class, c("person_count", "facet_level_count", "facet_level_count"))
  expect_true(is.list(spec$planning_schema$future_design_template))
  expect_identical(spec$planning_schema$future_design_template$facets, list(person = 16L, judge = 3L, task = 2L))
  expect_identical(spec$planning_schema$future_design_template$assignment, 2L)
  expect_true(is.list(spec$planning_schema$future_branch_schema))
  expect_identical(spec$planning_schema$future_branch_schema$planner_contract, "arbitrary_facet_planning_scaffold")
  expect_identical(spec$planning_schema$future_branch_schema$planner_stage, "schema_only")
  expect_identical(spec$planning_schema$future_branch_schema$input_contract, "design$facets(named counts)")
  expect_true(is.data.frame(spec$planning_schema$future_branch_schema$facet_table))
  expect_identical(spec$planning_schema$future_branch_schema$facet_table$future_facet_key, c("person", "judge", "task"))
  expect_true(is.list(spec$planning_schema$future_branch_schema$design_template))
  expect_identical(spec$planning_schema$future_branch_schema$design_template$facets, list(person = 16L, judge = 3L, task = 2L))
  expect_identical(spec$planning_schema$future_branch_schema$design_template$assignment, 2L)
  expect_true(is.list(spec$planning_schema$future_branch_schema$assignment_axis))
  expect_identical(spec$planning_schema$future_branch_schema$assignment_axis$current_planning_count_alias, "judge_per_person")
  expect_identical(spec$planning_schema$future_branch_schema$assignment_axis$future_input_key, "assignment")
  expect_true(is.list(spec$planning_schema$future_branch_schema$design_schema))
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$schema_contract, "arbitrary_facet_design_schema")
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$schema_stage, "schema_only")
  expect_true(is.data.frame(spec$planning_schema$future_branch_schema$design_schema$facet_axes))
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$facet_axes$input_key, c("person", "judge", "task"))
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$assignment_axis$current_planning_count_alias, "judge_per_person")
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$default_design$facets, list(person = 16L, judge = 3L, task = 2L))
  expect_true(is.list(spec$planning_schema$future_branch_schema$design_schema$grid_semantics))
  expect_identical(
    spec$planning_schema$future_branch_schema$design_schema$grid_semantics$canonical_columns,
    c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person")
  )
  expect_true(is.list(spec$planning_schema$future_branch_schema$preview))
  expect_true(spec$planning_schema$future_branch_schema$preview$preview_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$grid_contract))
  expect_true(is.list(spec$planning_schema$future_branch_schema$grid_bundle))
  expect_true(is.list(spec$planning_schema$future_branch_schema$grid_context))
  expect_identical(
    spec$planning_schema$future_branch_schema$grid_contract$contract,
    "arbitrary_facet_design_grid_contract"
  )
  expect_true(spec$planning_schema$future_branch_schema$grid_contract$preview_available)
  expect_identical(
    spec$planning_schema$future_branch_schema$grid_bundle$bundle_contract,
    "arbitrary_facet_design_grid_bundle"
  )
  expect_true(spec$planning_schema$future_branch_schema$grid_bundle$grid_available)
  expect_true(spec$planning_schema$future_branch_schema$grid_context$grid_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$grid_semantics))
  expect_true(is.list(spec$planning_schema$future_branch_preview))
  expect_true(spec$planning_schema$future_branch_preview$preview_available)
  expect_true(is.list(spec$planning_schema$future_branch_preview$grid_semantics))
  expect_true(is.list(spec$planning_schema$future_branch_grid_semantics))
  expect_true(is.list(spec$planning_schema$future_branch_grid_contract))
  expect_true(is.list(spec$planning_schema$future_branch_grid_bundle))
  expect_true(is.list(spec$planning_schema$future_branch_grid_context))
  expect_identical(
    spec$planning_schema$future_branch_grid_contract$planner_contract,
    "arbitrary_facet_planning_scaffold"
  )
  expect_identical(
    spec$planning_schema$future_branch_grid_bundle$planner_contract,
    "arbitrary_facet_planning_scaffold"
  )
  expect_identical(
    spec$planning_schema$future_branch_grid_semantics$branch_columns,
    c("design_id", "person", "judge", "task", "assignment")
  )
  expect_identical(
    spec$planning_schema$future_branch_grid_context$context_contract,
    "arbitrary_facet_design_grid_context"
  )
  expect_true(spec$planning_schema$future_branch_grid_context$grid_available)
  expect_identical(
    names(spec$planning_schema$future_branch_preview$canonical),
    c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person")
  )
  expect_identical(
    spec$planning_schema$future_branch_preview$canonical[1, c("n_person", "n_rater", "n_criterion", "raters_per_person")],
    tibble::tibble(n_person = 16L, n_rater = 3L, n_criterion = 2L, raters_per_person = 2L)
  )
  expect_true(is.data.frame(spec$planning_schema$role_table))
  expect_identical(spec$planning_schema$role_table$alias, c("n_person", "n_judge", "n_task", "judge_per_person"))
  expect_true(all(spec$planning_schema$role_table$mutable))

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 111)
  expect_named(sim, c("Study", "Person", "Judge", "Task", "Score"))
})

test_that("build_mfrm_sim_spec accepts role-based and alias-based design input", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(person = 14, n_judge = 3, n_task = 2, judge_per_person = 2),
    assignment = "rotating"
  )

  expect_equal(spec$n_person, 14)
  expect_equal(spec$n_rater, 3)
  expect_equal(spec$n_criterion, 2)
  expect_equal(spec$raters_per_person, 2)
  expect_identical(unname(spec$facet_names), c("Judge", "Task"))

  expect_error(
    build_mfrm_sim_spec(
      n_person = 12,
      design = list(person = 14),
      assignment = "rotating"
    ),
    "Do not supply the same design variable through both scalar arguments and `design`",
    fixed = TRUE
  )
})

test_that("design$facets accepts schema-only future facet-count keys", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  expect_equal(spec$n_person, 14)
  expect_equal(spec$n_rater, 3)
  expect_equal(spec$n_criterion, 2)
  expect_equal(spec$raters_per_person, 2)
  expect_identical(spec$planning_scope$future_branch_input_contract, "design$facets(named counts)")
  expect_identical(spec$planning_schema$future_facet_table$future_facet_key, c("person", "judge", "task"))
  expect_identical(spec$planning_schema$future_design_template$facets, list(person = 14L, judge = 3L, task = 2L))
  expect_identical(spec$planning_schema$future_design_template$assignment, 2L)
  expect_identical(spec$planning_schema$future_branch_schema$design_template$facets, list(person = 14L, judge = 3L, task = 2L))
  expect_identical(spec$planning_schema$future_branch_schema$assignment_axis$current_planning_count_alias, "judge_per_person")
  expect_identical(spec$planning_schema$future_branch_schema$design_schema$default_design$facets, list(person = 14L, judge = 3L, task = 2L))
  expect_true(spec$planning_schema$future_branch_preview$preview_available)
  expect_true(all(spec$planning_schema$future_branch_preview$canonical$n_person == 14L))

  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      design = list(facets = c(person = 15, judge = 3, task = 2), assignment = 2),
      reps = 1,
      maxit = 10,
      seed = 115
    )
  )
  expect_true(all(pred$design$n_person == 15))
  expect_true(all(pred$design$n_rater == 3))
  expect_true(all(pred$design$n_criterion == 2))

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      sim_spec = spec,
      design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
      reps = 1,
      maxit = 10,
      seed = 116
    )
  )
  expect_true(all(sim_eval$design_grid$n_person == 14))
  expect_true(all(sim_eval$design_grid$n_rater == 3))
  expect_true(all(sim_eval$design_grid$n_criterion == 2))
  expect_false(sim_eval$settings$progress)

  expect_error(
    build_mfrm_sim_spec(
      facet_names = c("Judge", "Task"),
      design = list(n_judge = 3, facets = c(judge = 3, person = 14, task = 2), assignment = 2),
      assignment = "rotating"
    ),
    "same design variable through both top-level names and `design$facets`",
    fixed = TRUE
  )
})

test_that("future_branch_schema normalizes nested facet-count design stubs", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  parse_branch <- getFromNamespace("simulation_parse_future_branch_design", "mfrmr")
  coerce_branch <- getFromNamespace("simulation_coerce_future_branch_design_schema", "mfrmr")

  schema <- coerce_branch(spec$planning_schema$future_branch_schema)
  expect_identical(schema$schema_contract, "arbitrary_facet_design_schema")
  expect_true(is.data.frame(schema$facet_axes))
  expect_identical(schema$facet_axes$input_key, c("person", "judge", "task"))
  expect_identical(schema$assignment_axis$current_planning_count_alias, "judge_per_person")
  expect_identical(schema$default_design$facets, list(person = 14L, judge = 3L, task = 2L))
  expect_true(is.list(schema$grid_semantics))
  expect_identical(schema$grid_semantics$feasibility_rule, "raters_per_person <= n_rater")

  parsed <- parse_branch(
    design = list(facets = c(person = 15, judge = 4, task = 2), assignment = 3),
    future_branch_schema = spec$planning_schema$future_branch_schema
  )
  expect_identical(
    parsed,
    list(n_person = 15, n_rater = 4, n_criterion = 2, raters_per_person = 3)
  )

  parsed_assignment_only <- parse_branch(
    design = list(assignment = 2),
    future_branch_schema = spec$planning_schema$future_branch_schema
  )
  expect_identical(parsed_assignment_only, list(raters_per_person = 2))

  expect_error(
    parse_branch(
      design = list(facets = c(panel = 3)),
      future_branch_schema = spec$planning_schema$future_branch_schema
    ),
    "must use names such as person, judge, task",
    fixed = TRUE
  )
})

test_that("future branch design schema builds schema-only grids", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  build_future_grid <- getFromNamespace("simulation_build_future_branch_design_grid", "mfrmr")

  grid <- build_future_grid(
    design_schema = spec$planning_schema$future_branch_schema,
    design = list(
      facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
      assignment = c(1, 2)
    ),
    id_prefix = "F"
  )

  expect_true(is.data.frame(grid$canonical))
  expect_true(is.data.frame(grid$public))
  expect_true(is.data.frame(grid$branch))
  expect_true(is.list(grid$design_schema))
  expect_true(is.list(grid$grid_semantics))
  expect_identical(names(grid$canonical), c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person"))
  expect_identical(names(grid$branch), c("design_id", "person", "judge", "task", "assignment"))
  expect_identical(grid$grid_semantics$public_columns, c("design_id", "n_person", "n_judge", "n_task", "judge_per_person"))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(grid$public)))
  expect_equal(sort(unique(grid$canonical$n_person)), c(12, 14))
  expect_equal(sort(unique(grid$canonical$n_criterion)), c(2, 4))
  expect_equal(sort(unique(grid$canonical$raters_per_person)), c(1, 2))
  expect_true(all(grid$canonical$raters_per_person <= grid$canonical$n_rater))

  expect_error(
    build_future_grid(
      design_schema = spec$planning_schema$future_branch_schema,
      design = list(facets = list(person = 10, judge = 2, task = 2), assignment = 3)
    ),
    "No valid future-branch design rows remain",
    fixed = TRUE
  )
})

test_that("future branch grid contract can be coerced and rebuilt independently", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  coerce_grid_contract <- getFromNamespace("simulation_coerce_future_branch_grid_contract", "mfrmr")
  build_from_contract <- getFromNamespace("simulation_build_future_branch_design_grid_from_contract", "mfrmr")

  contract <- coerce_grid_contract(spec$planning_schema)
  expect_identical(contract$contract, "arbitrary_facet_design_grid_contract")
  expect_true(contract$preview_available)
  expect_true(is.list(contract$design_schema))
  expect_true(is.list(contract$grid_semantics))

  rebuilt_default <- build_from_contract(spec$planning_schema$future_branch_grid_contract)
  expect_identical(
    rebuilt_default$canonical,
    spec$planning_schema$future_branch_grid_contract$canonical
  )
  expect_identical(
    rebuilt_default$grid_semantics$branch_columns,
    c("design_id", "person", "judge", "task", "assignment")
  )

  rebuilt_custom <- build_from_contract(
    spec$planning_schema$future_branch_grid_contract,
    design = list(facets = list(person = c(12, 14), judge = 3, task = 2), assignment = c(1, 2))
  )
  expect_true(is.data.frame(rebuilt_custom$canonical))
  expect_equal(sort(unique(rebuilt_custom$canonical$n_person)), c(12, 14))
  expect_equal(sort(unique(rebuilt_custom$canonical$raters_per_person)), c(1, 2))
  expect_true(is.list(rebuilt_custom$grid_contract))
})

test_that("future branch grid bundle bundles live branch grids through one object", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  build_bundle <- getFromNamespace("simulation_future_branch_grid_bundle", "mfrmr")
  coerce_bundle <- getFromNamespace("simulation_coerce_future_branch_grid_bundle", "mfrmr")

  bundle <- build_bundle(spec$planning_schema)
  expect_identical(bundle$bundle_contract, "arbitrary_facet_design_grid_bundle")
  expect_true(bundle$grid_available)
  expect_true(is.data.frame(bundle$canonical))
  expect_true(is.list(bundle$grid_contract))
  expect_identical(
    bundle$canonical,
    spec$planning_schema$future_branch_grid_contract$canonical
  )

  coerced <- coerce_bundle(spec$planning_schema)
  expect_identical(coerced$bundle_contract, "arbitrary_facet_design_grid_bundle")
  expect_true(coerced$grid_available)

  custom_bundle <- build_bundle(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = 2), assignment = c(1, 2))
  )
  expect_equal(sort(unique(custom_bundle$canonical$n_person)), c(12, 14))
  expect_equal(sort(unique(custom_bundle$canonical$raters_per_person)), c(1, 2))
})

test_that("future branch grid bundle can be materialized and viewed from live objects", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  materialize_bundle <- getFromNamespace("simulation_materialize_future_branch_grid_bundle", "mfrmr")
  grid_view <- getFromNamespace("simulation_future_branch_grid_view", "mfrmr")

  from_spec <- materialize_bundle(spec)
  expect_identical(from_spec$bundle_contract, "arbitrary_facet_design_grid_bundle")
  expect_true(from_spec$grid_available)
  expect_identical(
    from_spec$branch,
    spec$planning_schema$future_branch_grid_bundle$branch
  )

  from_schema <- materialize_bundle(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = 2), assignment = c(1, 2))
  )
  expect_true(is.data.frame(from_schema$canonical))
  expect_equal(sort(unique(from_schema$canonical$n_person)), c(12, 14))

  expect_identical(
    names(grid_view(spec$planning_schema, view = "branch")),
    c("design_id", "person", "judge", "task", "assignment")
  )
  expect_true(all(
    c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person",
      "n_judge", "n_task", "judge_per_person") %in%
      names(grid_view(spec$planning_schema, view = "public"))
  ))
})

test_that("future branch grid context summarizes varying and fixed axes", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  grid_context <- getFromNamespace("simulation_future_branch_grid_context", "mfrmr")

  context_default <- grid_context(spec$planning_schema)
  expect_identical(context_default$context_contract, "arbitrary_facet_design_grid_context")
  expect_true(context_default$grid_available)
  expect_true(is.data.frame(context_default$axis_table))
  expect_identical(context_default$varying_canonical, character(0))
  expect_identical(sort(context_default$fixed_canonical), sort(c("n_person", "n_rater", "n_criterion", "raters_per_person")))
  expect_true(all(context_default$axis_table$n_values == 1L))

  context_custom <- grid_context(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2))
  )
  expect_true(context_custom$grid_available)
  expect_equal(sort(context_custom$varying_canonical), c("n_criterion", "n_person", "raters_per_person"))
  expect_identical(context_custom$fixed_canonical, "n_rater")
  expect_equal(
    context_custom$axis_table$fixed_value[match("n_rater", context_custom$axis_table$canonical_design_variable)],
    3L
  )
  expect_true(is.list(context_custom$grid_bundle))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  context_unavailable <- grid_context(branch)
  expect_false(context_unavailable$grid_available)
  expect_true(all(is.na(context_unavailable$axis_table$n_values)))
})

test_that("future branch summary and baseline recommendation use schema-only grid context", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  grid_summary <- getFromNamespace("simulation_future_branch_grid_summary", "mfrmr")
  grid_recommendation <- getFromNamespace("simulation_future_branch_grid_recommendation", "mfrmr")

  summary_default <- grid_summary(spec$planning_schema)
  expect_identical(summary_default$summary_contract, "arbitrary_facet_design_grid_summary")
  expect_true(summary_default$grid_available)
  expect_identical(summary_default$n_designs, 1L)
  expect_identical(summary_default$n_varying_axes, 0L)
  expect_identical(summary_default$n_fixed_axes, 4L)
  expect_identical(summary_default$fixed_values[["n_person"]], 14L)
  expect_true(is.list(summary_default$grid_context))

  summary_custom <- grid_summary(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2))
  )
  expect_true(summary_custom$grid_available)
  expect_identical(summary_custom$n_designs, 8L)
  expect_equal(sort(summary_custom$varying_canonical), c("n_criterion", "n_person", "raters_per_person"))
  expect_identical(summary_custom$fixed_values[["n_rater"]], 3L)

  rec <- grid_recommendation(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    prefer = c("assignment", "task", "person")
  )
  expect_identical(rec$recommendation_contract, "arbitrary_facet_design_grid_recommendation")
  expect_true(rec$recommendation_available)
  expect_identical(rec$prefer, c("raters_per_person", "n_criterion", "n_person"))
  expect_identical(rec$recommended_canonical$n_person[[1]], 12L)
  expect_identical(rec$recommended_canonical$n_rater[[1]], 3L)
  expect_identical(rec$recommended_canonical$n_criterion[[1]], 2L)
  expect_identical(rec$recommended_canonical$raters_per_person[[1]], 1L)
  expect_identical(rec$recommended_branch$assignment[[1]], 1L)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  summary_unavailable <- grid_summary(branch)
  expect_false(summary_unavailable$grid_available)
  expect_identical(summary_unavailable$n_designs, 0L)

  rec_unavailable <- grid_recommendation(branch)
  expect_false(rec_unavailable$recommendation_available)
})

test_that("future branch table and plot payload expose schema-only branch views", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  grid_table <- getFromNamespace("simulation_future_branch_grid_table", "mfrmr")
  grid_plot <- getFromNamespace("simulation_future_branch_grid_plot_payload", "mfrmr")

  table_obj <- grid_table(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    prefer = c("assignment", "task", "person"),
    view = "public"
  )
  expect_identical(table_obj$table_contract, "arbitrary_facet_design_grid_table")
  expect_true(table_obj$grid_available)
  expect_true(is.data.frame(table_obj$table))
  expect_true(all(c("design_id", "n_person", "n_judge", "n_task", "judge_per_person", "recommended") %in% names(table_obj$table)))
  expect_identical(sum(table_obj$table$recommended), 1L)

  plot_obj <- grid_plot(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    x_var = "assignment",
    group_var = "task",
    prefer = c("assignment", "task", "person"),
    view = "branch"
  )
  expect_s3_class(plot_obj, "mfrm_plot_data")
  expect_identical(plot_obj$name, "future_branch_grid_schema")
  expect_true(isTRUE(plot_obj$data$plot_available))
  expect_identical(plot_obj$data$x_var, "raters_per_person")
  expect_identical(plot_obj$data$group_var, "n_criterion")
  expect_identical(plot_obj$data$x_label, "assignment")
  expect_identical(plot_obj$data$group_label, "task")
  expect_true(all(c("design_id", "x_value", "group_value", "recommended") %in% names(plot_obj$data$data)))
  expect_identical(sum(plot_obj$data$data$recommended), 1L)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_plot <- grid_plot(branch)
  expect_s3_class(unavailable_plot, "mfrm_plot_data")
  expect_false(unavailable_plot$data$plot_available)
  expect_equal(nrow(unavailable_plot$data$data), 0)
})

test_that("future branch report bundle combines schema-only consumers", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  grid_report <- getFromNamespace("simulation_future_branch_report_bundle", "mfrmr")

  report <- grid_report(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(report$report_contract, "arbitrary_facet_design_grid_report_bundle")
  expect_true(report$report_available)
  expect_identical(report$recommended_design_id, "F01")
  expect_true(is.data.frame(report$overview_table))
  expect_true(all(c("canonical", "public", "branch") %in% names(report$tables)))
  expect_true(all(c("canonical", "public", "branch") %in% names(report$plots)))
  expect_identical(report$tables$public$recommended_design_id, "F01")
  expect_s3_class(report$plots$branch, "mfrm_plot_data")
  expect_true(report$plots$branch$data$plot_available)
  expect_identical(report$plots$branch$data$recommended_design_id, "F01")
  expect_true(is.list(spec$planning_schema$future_branch_report_bundle))
  expect_true(spec$planning_schema$future_branch_report_bundle$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_bundle))
  expect_true(spec$planning_schema$future_branch_schema$report_bundle$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  report_unavailable <- grid_report(branch)
  expect_false(report_unavailable$report_available)
  expect_identical(nrow(report_unavailable$tables$public$table), 0L)
  expect_false(report_unavailable$plots$public$data$plot_available)
})

test_that("future branch report summary exposes a compact component index", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_summary <- getFromNamespace("simulation_future_branch_report_summary", "mfrmr")

  summary_obj <- report_summary(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(summary_obj$report_summary_contract, "arbitrary_facet_design_grid_report_summary")
  expect_true(summary_obj$report_available)
  expect_identical(summary_obj$n_designs, 8L)
  expect_identical(summary_obj$recommended_design_id, "F01")
  expect_equal(summary_obj$available_table_views, c("canonical", "public", "branch"))
  expect_equal(summary_obj$available_plot_views, c("canonical", "public", "branch"))
  expect_true(is.data.frame(summary_obj$component_index))
  expect_identical(nrow(summary_obj$component_index), 8L)
  expect_true(all(c("component_type", "view", "available", "n_rows", "x_var", "group_var", "recommended_design_id") %in% names(summary_obj$component_index)))
  expect_true(is.list(spec$planning_schema$future_branch_report_summary))
  expect_true(spec$planning_schema$future_branch_report_summary$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_summary))
  expect_true(spec$planning_schema$future_branch_schema$report_summary$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  summary_unavailable <- report_summary(branch)
  expect_false(summary_unavailable$report_available)
  expect_identical(summary_unavailable$n_designs, 0L)
  expect_true(is.data.frame(summary_unavailable$component_index))
})

test_that("future branch report overview exposes compact metrics and axis state", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_overview <- getFromNamespace("simulation_future_branch_report_overview_table", "mfrmr")

  overview_obj <- report_overview(
    spec$planning_schema,
    design = list(facets = list(person = c(12, 14), judge = 3, task = c(2, 4)), assignment = c(1, 2)),
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(overview_obj$overview_contract, "arbitrary_facet_design_report_overview_table")
  expect_true(overview_obj$report_available)
  expect_true(is.data.frame(overview_obj$metrics_table))
  expect_true(is.data.frame(overview_obj$axis_overview_table))
  expect_true(is.data.frame(overview_obj$component_index))
  expect_identical(overview_obj$metrics_table$n_designs[[1]], 8L)
  expect_identical(overview_obj$metrics_table$recommended_design_id[[1]], "F01")
  expect_identical(overview_obj$metrics_table$n_table_views[[1]], 3L)
  expect_identical(overview_obj$metrics_table$n_plot_views[[1]], 3L)
  expect_true(all(c(
    "axis_source", "input_key", "canonical_design_variable", "public_design_alias",
    "axis_class", "facet", "axis_state", "n_values", "value_summary", "fixed_value"
  ) %in% names(overview_obj$axis_overview_table)))
  expect_true(all(c("varying", "fixed") %in% overview_obj$axis_overview_table$axis_state))
  expect_true(is.list(spec$planning_schema$future_branch_report_overview_table))
  expect_true(spec$planning_schema$future_branch_report_overview_table$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_overview_table))
  expect_true(spec$planning_schema$future_branch_schema$report_overview_table$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  overview_unavailable <- report_overview(branch)
  expect_false(overview_unavailable$report_available)
  expect_identical(overview_unavailable$metrics_table$n_designs[[1]], 0L)
  expect_true(all(overview_unavailable$axis_overview_table$axis_state == "unavailable"))
  expect_true(is.data.frame(overview_unavailable$component_index))
})

test_that("future branch report overview view exposes selected component tables", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  overview_view <- getFromNamespace("simulation_future_branch_report_overview_view", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  metrics_view <- overview_view(
    spec$planning_schema,
    component = "metrics",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(metrics_view$overview_view_contract, "arbitrary_facet_design_report_overview_view")
  expect_true(metrics_view$report_available)
  expect_identical(metrics_view$component, "metrics")
  expect_identical(metrics_view$component_label, "report metrics")
  expect_identical(metrics_view$table$recommended_design_id[[1]], "F01")

  axes_view <- overview_view(spec$planning_schema, component = "axes", design = design_override)
  expect_identical(axes_view$component_label, "axis overview")
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(axes_view$table)))

  components_view <- overview_view(spec$planning_schema, component = "components", design = design_override)
  expect_identical(components_view$component_label, "component index")
  expect_true(all(c("component_type", "available") %in% names(components_view$table)))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_view <- overview_view(branch, component = "axes")
  expect_false(unavailable_view$report_available)
  expect_true(all(unavailable_view$table$axis_state == "unavailable"))
})

test_that("future branch report catalog enumerates overview surfaces", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_catalog <- getFromNamespace("simulation_future_branch_report_catalog", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  catalog <- report_catalog(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(catalog$catalog_contract, "arbitrary_facet_design_report_catalog")
  expect_true(catalog$report_available)
  expect_identical(catalog$recommended_design_id, "F01")
  expect_true(is.data.frame(catalog$surface_index))
  expect_equal(catalog$surface_index$component, c("metrics", "axes", "components"))
  expect_true(all(catalog$surface_index$available))
  expect_true(is.list(catalog$views))
  expect_true(all(c("metrics", "axes", "components") %in% names(catalog$views)))
  expect_true(is.list(spec$planning_schema$future_branch_report_catalog))
  expect_true(spec$planning_schema$future_branch_report_catalog$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_catalog))
  expect_true(spec$planning_schema$future_branch_schema$report_catalog$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_catalog <- report_catalog(branch)
  expect_false(unavailable_catalog$report_available)
  expect_true(all(!unavailable_catalog$surface_index$available))
})

test_that("future branch report digest exposes headline branch-side report state", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_digest <- getFromNamespace("simulation_future_branch_report_digest", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  digest <- report_digest(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(digest$digest_contract, "arbitrary_facet_design_report_digest")
  expect_true(digest$report_available)
  expect_identical(digest$recommended_design_id, "F01")
  expect_true(is.data.frame(digest$digest_table))
  expect_identical(digest$digest_table$n_designs[[1]], 8L)
  expect_identical(digest$digest_table$n_available_surfaces[[1]], 3L)
  expect_equal(digest$available_surfaces, c("metrics", "axes", "components"))
  expect_equal(sort(digest$varying_axes), sort(c("n_person", "n_criterion", "raters_per_person")))
  expect_equal(digest$fixed_axes, "n_rater")
  expect_true(is.list(spec$planning_schema$future_branch_report_digest))
  expect_true(spec$planning_schema$future_branch_report_digest$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_digest))
  expect_true(spec$planning_schema$future_branch_schema$report_digest$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_digest <- report_digest(branch)
  expect_false(unavailable_digest$report_available)
  expect_identical(unavailable_digest$digest_table$n_designs[[1]], 0L)
  expect_identical(unavailable_digest$digest_table$n_available_surfaces[[1]], 0L)
  expect_length(unavailable_digest$available_surfaces, 0L)
})

test_that("future branch report surface exposes selected compact report layers", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_surface <- getFromNamespace("simulation_future_branch_report_surface", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  digest_surface <- report_surface(
    spec$planning_schema,
    surface = "digest",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(digest_surface$surface_contract, "arbitrary_facet_design_report_surface")
  expect_identical(digest_surface$surface, "digest")
  expect_true(digest_surface$report_available)
  expect_identical(digest_surface$recommended_design_id, "F01")
  expect_true(all(c("n_designs", "available_surfaces", "varying_axes", "fixed_axes") %in% names(digest_surface$table)))

  catalog_surface <- report_surface(spec$planning_schema, surface = "catalog", design = design_override)
  expect_identical(catalog_surface$surface_label, "report catalog")
  expect_true(all(c("surface", "available", "n_rows") %in% names(catalog_surface$table)))

  axes_surface <- report_surface(spec$planning_schema, surface = "axes", design = design_override)
  expect_identical(axes_surface$surface_label, "axis overview")
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(axes_surface$table)))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_surface <- report_surface(branch, surface = "digest")
  expect_false(unavailable_surface$report_available)
  expect_identical(unavailable_surface$table$n_designs[[1]], 0L)
})

test_that("future branch report surface registry enumerates compact report surfaces", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_registry <- getFromNamespace("simulation_future_branch_report_surface_registry", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  registry <- report_registry(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(registry$registry_contract, "arbitrary_facet_design_report_surface_registry")
  expect_true(registry$report_available)
  expect_identical(registry$recommended_design_id, "F01")
  expect_true(is.data.frame(registry$surface_index))
  expect_equal(registry$surface_index$surface, c("digest", "catalog", "metrics", "axes", "components"))
  expect_true(all(registry$surface_index$available))
  expect_true(is.list(registry$surfaces))
  expect_true(all(c("digest", "catalog", "metrics", "axes", "components") %in% names(registry$surfaces)))
  expect_true(is.list(spec$planning_schema$future_branch_report_surface_registry))
  expect_true(spec$planning_schema$future_branch_report_surface_registry$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_surface_registry))
  expect_true(spec$planning_schema$future_branch_schema$report_surface_registry$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_registry <- report_registry(branch)
  expect_false(unavailable_registry$report_available)
  expect_true(all(!unavailable_registry$surface_index$available))
})

test_that("future branch report panel combines one selected surface with headline metadata", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_panel <- getFromNamespace("simulation_future_branch_report_panel", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  default_panel <- report_panel(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(default_panel$panel_contract, "arbitrary_facet_design_report_panel")
  expect_identical(default_panel$surface, "digest")
  expect_true(default_panel$report_available)
  expect_identical(default_panel$recommended_design_id, "F01")
  expect_true(all(c("n_designs", "available_surfaces") %in% names(default_panel$digest_table)))
  expect_true(all(c("surface", "available", "n_rows") %in% names(default_panel$surface_index)))
  expect_true(all(c("n_designs", "available_surfaces", "varying_axes", "fixed_axes") %in% names(default_panel$selected_table)))

  axes_panel <- report_panel(spec$planning_schema, surface = "axes", design = design_override)
  expect_identical(axes_panel$surface_label, "axis overview")
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(axes_panel$selected_table)))

  expect_true(is.list(spec$planning_schema$future_branch_report_panel))
  expect_true(spec$planning_schema$future_branch_report_panel$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_panel))
  expect_true(spec$planning_schema$future_branch_schema$report_panel$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_panel <- report_panel(branch)
  expect_false(unavailable_panel$report_available)
  expect_identical(unavailable_panel$selected_table$n_designs[[1]], 0L)
})

test_that("future branch report operation combines digest, registry, and selected surface", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_operation <- getFromNamespace("simulation_future_branch_report_operation", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  operation <- report_operation(
    spec$planning_schema,
    surface = "axes",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(operation$operation_contract, "arbitrary_facet_design_report_operation")
  expect_identical(operation$surface, "axes")
  expect_true(operation$report_available)
  expect_identical(operation$recommended_design_id, "F01")
  expect_true(all(c("n_designs", "available_surfaces") %in% names(operation$digest_table)))
  expect_true(all(c("surface", "available", "n_rows") %in% names(operation$surface_index)))
  expect_true(all(c("recommended_design_id", "n_table_views") %in% names(operation$metrics_table)))
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(operation$axis_overview_table)))
  expect_true(all(c("component_type", "available") %in% names(operation$component_index)))
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(operation$selected_table)))
  expect_true(is.list(spec$planning_schema$future_branch_report_operation))
  expect_true(spec$planning_schema$future_branch_report_operation$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_operation))
  expect_true(spec$planning_schema$future_branch_schema$report_operation$report_available)
  expect_true(all(c("digest_table", "surface_index", "metrics_table") %in% names(
    spec$planning_schema$future_branch_schema$report_operation
  )))
  expect_false("report_panel" %in% names(spec$planning_schema$future_branch_schema$report_operation))
  expect_false("report_surface_registry" %in% names(spec$planning_schema$future_branch_schema$report_operation))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_operation <- report_operation(branch)
  expect_false(unavailable_operation$report_available)
  expect_identical(unavailable_operation$digest_table$n_designs[[1]], 0L)
  expect_true(all(c("surface", "available", "n_rows") %in% names(unavailable_operation$surface_index)))
})

test_that("future branch report snapshot exposes compact headline tables without nested operation state", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_snapshot <- getFromNamespace("simulation_future_branch_report_snapshot", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  snapshot <- report_snapshot(
    spec$planning_schema,
    surface = "axes",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(snapshot$snapshot_contract, "arbitrary_facet_design_report_snapshot")
  expect_identical(snapshot$surface, "axes")
  expect_true(snapshot$report_available)
  expect_identical(snapshot$recommended_design_id, "F01")
  expect_identical(snapshot$n_designs, 8L)
  expect_equal(snapshot$available_surfaces, c("digest", "catalog", "metrics", "axes", "components"))
  expect_equal(sort(snapshot$varying_axes), sort(c("n_person", "n_criterion", "raters_per_person")))
  expect_equal(snapshot$fixed_axes, "n_rater")
  expect_true(all(c("n_designs", "available_surfaces") %in% names(snapshot$digest_table)))
  expect_true(all(c("surface", "available", "n_rows") %in% names(snapshot$surface_index)))
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(snapshot$selected_table)))
  expect_false("report_panel" %in% names(snapshot))
  expect_false("report_surface_registry" %in% names(snapshot))
  expect_true(is.list(spec$planning_schema$future_branch_report_snapshot))
  expect_true(spec$planning_schema$future_branch_report_snapshot$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_snapshot))
  expect_true(spec$planning_schema$future_branch_schema$report_snapshot$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_snapshot <- report_snapshot(branch)
  expect_false(unavailable_snapshot$report_available)
  expect_identical(unavailable_snapshot$n_designs, 0L)
  expect_length(unavailable_snapshot$available_surfaces, 0L)
})

test_that("future branch report brief exposes one selected surface without snapshot payload", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_brief <- getFromNamespace("simulation_future_branch_report_brief", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  brief <- report_brief(
    spec$planning_schema,
    surface = "components",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(brief$brief_contract, "arbitrary_facet_design_report_brief")
  expect_identical(brief$surface, "components")
  expect_true(brief$report_available)
  expect_true(is.data.frame(brief$headline_table))
  expect_identical(brief$headline_table$n_designs[[1]], 8L)
  expect_identical(brief$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(
    brief$headline_table$available_surfaces[[1]],
    "digest, catalog, metrics, axes, components"
  )
  expect_identical(
    sort(strsplit(brief$headline_table$varying_axes[[1]], ", ", fixed = TRUE)[[1]]),
    sort(c("n_person", "n_criterion", "raters_per_person"))
  )
  expect_identical(brief$headline_table$fixed_axes[[1]], "n_rater")
  expect_true(all(c("component_type", "available") %in% names(brief$selected_table)))
  expect_true(all(c("surface", "available", "n_rows") %in% names(brief$surface_index)))
  expect_false("digest_table" %in% names(brief))
  expect_false("metrics_table" %in% names(brief))
  expect_false("axis_overview_table" %in% names(brief))
  expect_true(is.list(spec$planning_schema$future_branch_report_brief))
  expect_true(spec$planning_schema$future_branch_report_brief$report_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_brief))
  expect_true(spec$planning_schema$future_branch_schema$report_brief$report_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_brief <- report_brief(branch)
  expect_false(unavailable_brief$report_available)
  expect_identical(unavailable_brief$headline_table$n_designs[[1]], 0L)
  expect_identical(unavailable_brief$headline_table$available_surfaces[[1]], "")
})

test_that("future branch report consumer dispatches across brief snapshot and operation modes", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  report_consume <- getFromNamespace("simulation_future_branch_report_consume", "mfrmr")
  mode_registry <- getFromNamespace("simulation_future_branch_report_mode_registry", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  registry <- mode_registry(spec$planning_schema)
  expect_true(is.data.frame(registry))
  expect_identical(registry$mode, c("brief", "snapshot", "operation"))
  expect_identical(registry$default_surface, rep("digest", 3))
  expect_identical(registry$carries_nested_graph, c(FALSE, FALSE, TRUE))

  brief_consumer <- report_consume(
    spec$planning_schema,
    mode = "brief",
    surface = "axes",
    design = design_override,
    prefer = c("assignment", "task", "person"),
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(brief_consumer$consumer_contract, "arbitrary_facet_design_report_consumer")
  expect_identical(brief_consumer$mode, "brief")
  expect_identical(brief_consumer$payload_contract, "arbitrary_facet_design_report_brief")
  expect_true(brief_consumer$report_available)
  expect_identical(brief_consumer$recommended_design_id, "F01")
  expect_identical(brief_consumer$n_designs, 8L)
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(brief_consumer$selected_table)))
  expect_true(all(c("surface", "available", "n_rows") %in% names(brief_consumer$surface_index)))
  expect_true(is.list(brief_consumer$payload))
  expect_identical(brief_consumer$payload$brief_contract, "arbitrary_facet_design_report_brief")

  snapshot_consumer <- report_consume(spec$planning_schema, mode = "snapshot", surface = "metrics", design = design_override)
  expect_identical(snapshot_consumer$payload_contract, "arbitrary_facet_design_report_snapshot")
  expect_true(is.list(snapshot_consumer$payload))
  expect_identical(snapshot_consumer$payload$snapshot_contract, "arbitrary_facet_design_report_snapshot")

  operation_consumer <- report_consume(spec$planning_schema, mode = "operation", surface = "components", design = design_override)
  expect_identical(operation_consumer$payload_contract, "arbitrary_facet_design_report_operation")
  expect_true(is.list(operation_consumer$payload))
  expect_identical(operation_consumer$payload$operation_contract, "arbitrary_facet_design_report_operation")

  expect_true(is.data.frame(spec$planning_schema$future_branch_report_mode_registry))
  expect_true(is.list(spec$planning_schema$future_branch_report_consumer))
  expect_identical(
    spec$planning_schema$future_branch_report_consumer$payload_contract,
    "arbitrary_facet_design_report_brief"
  )
  expect_true(is.data.frame(spec$planning_schema$future_branch_schema$report_mode_registry))
  expect_true(is.list(spec$planning_schema$future_branch_schema$report_consumer))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_consumer <- report_consume(branch, mode = "brief")
  expect_false(unavailable_consumer$report_available)
  expect_identical(unavailable_consumer$n_designs, 0L)
})

test_that("future branch pilot bundles one grid view with one report consumer", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  future_pilot <- getFromNamespace("simulation_future_branch_pilot", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  pilot <- future_pilot(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "branch",
    mode = "brief",
    surface = "axes",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(pilot$pilot_contract, "arbitrary_facet_planning_pilot")
  expect_identical(pilot$pilot_stage, "pilot_active")
  expect_true(pilot$pilot_available)
  expect_identical(pilot$view, "branch")
  expect_identical(pilot$mode, "brief")
  expect_identical(pilot$surface, "axes")
  expect_identical(pilot$recommended_design_id, "F01")
  expect_identical(pilot$n_designs, 8L)
  expect_true(is.data.frame(pilot$grid_table))
  expect_true("recommended" %in% names(pilot$grid_table))
  expect_s3_class(pilot$plot_payload, "mfrm_plot_data")
  expect_true(is.list(pilot$report_consumer))
  expect_identical(
    pilot$report_consumer$payload_contract,
    "arbitrary_facet_design_report_brief"
  )
  expect_true(is.list(spec$planning_schema$future_branch_pilot))
  expect_true(spec$planning_schema$future_branch_pilot$pilot_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$pilot))
  expect_true(spec$planning_schema$future_branch_schema$pilot$pilot_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_pilot <- future_pilot(branch)
  expect_false(unavailable_pilot$pilot_available)
  expect_identical(unavailable_pilot$pilot_stage, "schema_only")
  expect_identical(unavailable_pilot$n_designs, 0L)
})

test_that("future branch pilot summary table and plot expose compact active consumers", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  pilot_summary <- getFromNamespace("simulation_future_branch_pilot_summary", "mfrmr")
  pilot_table <- getFromNamespace("simulation_future_branch_pilot_table", "mfrmr")
  pilot_plot <- getFromNamespace("simulation_future_branch_pilot_plot", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  summary_obj <- pilot_summary(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(summary_obj$pilot_summary_contract, "arbitrary_facet_planning_pilot_summary")
  expect_true(summary_obj$pilot_available)
  expect_true(is.data.frame(summary_obj$headline_table))
  expect_identical(summary_obj$headline_table$n_designs[[1]], 8L)
  expect_identical(summary_obj$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(summary_obj$headline_table$report_payload_contract[[1]], "arbitrary_facet_design_report_snapshot")

  grid_table_obj <- pilot_table(
    spec$planning_schema,
    component = "grid",
    design = design_override,
    view = "branch",
    mode = "brief",
    surface = "axes"
  )
  expect_identical(grid_table_obj$pilot_table_contract, "arbitrary_facet_planning_pilot_table")
  expect_true(grid_table_obj$pilot_available)
  expect_true(all(c("design_id", "recommended") %in% names(grid_table_obj$table)))

  report_table_obj <- pilot_table(
    spec$planning_schema,
    component = "report",
    design = design_override,
    view = "branch",
    mode = "brief",
    surface = "axes"
  )
  expect_true(all(c("axis_state", "canonical_design_variable") %in% names(report_table_obj$table)))

  plot_obj <- pilot_plot(
    spec$planning_schema,
    design = design_override,
    view = "branch",
    mode = "brief",
    surface = "axes",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(plot_obj$pilot_plot_contract, "arbitrary_facet_planning_pilot_plot")
  expect_true(plot_obj$pilot_available)
  expect_s3_class(plot_obj$plot, "mfrm_plot_data")

  expect_true(is.list(spec$planning_schema$future_branch_pilot_summary))
  expect_true(is.list(spec$planning_schema$future_branch_pilot_table))
  expect_true(is.list(spec$planning_schema$future_branch_pilot_plot))
  expect_true(is.list(spec$planning_schema$future_branch_schema$pilot_summary))
  expect_true(is.list(spec$planning_schema$future_branch_schema$pilot_table))
  expect_true(is.list(spec$planning_schema$future_branch_schema$pilot_plot))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$pilot_summary))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$pilot_table))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$pilot_plot))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_summary <- pilot_summary(branch)
  expect_false(unavailable_summary$pilot_available)
  expect_identical(unavailable_summary$headline_table$n_designs[[1]], 0L)
})

test_that("future branch active branch bundles pilot summary table and plot", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_branch <- getFromNamespace("simulation_future_branch_active_branch", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  branch_obj <- active_branch(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(branch_obj$branch_contract, "arbitrary_facet_planning_active_branch")
  expect_identical(branch_obj$branch_stage, "pilot_active")
  expect_true(branch_obj$branch_available)
  expect_identical(branch_obj$recommended_design_id, "F01")
  expect_identical(branch_obj$n_designs, 8L)
  expect_true(is.list(branch_obj$summary))
  expect_true(is.list(branch_obj$table))
  expect_true(is.list(branch_obj$plot))
  expect_identical(
    branch_obj$summary$pilot_summary_contract,
    "arbitrary_facet_planning_pilot_summary"
  )
  expect_identical(
    branch_obj$table$pilot_table_contract,
    "arbitrary_facet_planning_pilot_table"
  )
  expect_identical(
    branch_obj$plot$pilot_plot_contract,
    "arbitrary_facet_planning_pilot_plot"
  )
  expect_true(is.list(spec$planning_schema$future_branch_active_branch))
  expect_true(spec$planning_schema$future_branch_active_branch$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch))
  expect_true(spec$planning_schema$future_branch_schema$active_branch$branch_available)
  expect_true(is.data.frame(spec$planning_schema$future_branch_schema$active_branch$canonical_grid))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$active_branch$summary))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$active_branch$table))
  expect_false("pilot" %in% names(spec$planning_schema$future_branch_schema$active_branch$plot))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_branch_obj <- active_branch(branch)
  expect_false(unavailable_branch_obj$branch_available)
  expect_identical(unavailable_branch_obj$branch_stage, "schema_only")
  expect_identical(unavailable_branch_obj$n_designs, 0L)
})

test_that("future branch active profile exposes deterministic design bookkeeping metrics", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_profile <- getFromNamespace("simulation_future_branch_active_branch_profile", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  profile <- active_profile(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(profile$profile_contract, "arbitrary_facet_planning_active_branch_profile")
  expect_true(profile$branch_available)
  expect_true(is.data.frame(profile$profile_table))
  expect_true(is.data.frame(profile$metric_registry))
  expect_true(is.data.frame(profile$profile_summary_table))
  expect_true(all(c(
    "design_id", "total_observations", "observations_per_person",
    "observations_per_criterion", "expected_observations_per_rater",
    "assignment_fraction", "recommended"
  ) %in% names(profile$profile_table)))
  expect_equal(
    profile$metric_registry$metric,
    c(
      "total_observations",
      "observations_per_person",
      "observations_per_criterion",
      "expected_observations_per_rater",
      "assignment_fraction"
    )
  )
  expect_equal(
    profile$metric_registry$basis_class,
    c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "balanced_expectation",
      "density_ratio"
    )
  )
  expect_true(all(profile$metric_registry$psychometric == FALSE))
  row_f01 <- profile$profile_table[profile$profile_table$design_id == "F01", , drop = FALSE]
  expect_identical(nrow(row_f01), 1L)
  expect_equal(row_f01$total_observations[[1]], 24)
  expect_equal(row_f01$observations_per_person[[1]], 2)
  expect_equal(row_f01$observations_per_criterion[[1]], 12)
  expect_equal(row_f01$expected_observations_per_rater[[1]], 8)
  expect_equal(row_f01$assignment_fraction[[1]], 1 / 3)
  expect_true(row_f01$recommended[[1]])
  summary_row <- profile$profile_summary_table[
    profile$profile_summary_table$metric == "expected_observations_per_rater",
    , drop = FALSE
  ]
  expect_identical(nrow(summary_row), 1L)
  expect_equal(summary_row$recommended_value[[1]], 8)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_profile))
  expect_true(spec$planning_schema$future_branch_active_branch_profile$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_profile))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_profile$branch_available)
  expect_true(all(c("metric_registry", "profile_summary_table", "profile_table") %in% names(
    spec$planning_schema$future_branch_schema$active_branch_profile
  )))
  expect_false("active_branch" %in% names(spec$planning_schema$future_branch_schema$active_branch_profile))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_profile <- active_profile(branch)
  expect_false(unavailable_profile$branch_available)
  expect_true(is.data.frame(unavailable_profile$metric_registry))
  expect_true(is.data.frame(unavailable_profile$profile_summary_table))
  expect_identical(nrow(unavailable_profile$profile_table), 0L)
})

test_that("future branch active overview combines branch headline and metric summaries", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_overview <- getFromNamespace("simulation_future_branch_active_branch_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- active_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(overview$overview_contract, "arbitrary_facet_planning_active_branch_overview")
  expect_true(overview$branch_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$metric_registry))
  expect_true(is.data.frame(overview$metric_summary_table))
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$n_metrics[[1]], 5L)
  expect_true(all(c("metric", "basis_class", "formula", "interpretation") %in% names(overview$metric_registry)))
  expect_true(all(c("metric", "min", "max", "mean", "recommended_value") %in% names(overview$metric_summary_table)))
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_overview$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_overview$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- active_overview(branch)
  expect_false(unavailable_overview$branch_available)
  expect_identical(unavailable_overview$headline_table$n_designs[[1]], 0L)
  expect_identical(nrow(unavailable_overview$metric_summary_table), 0L)
})

test_that("future branch active load-balance diagnostics separate expectations from integer balance", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_load_balance <- getFromNamespace("simulation_future_branch_active_branch_load_balance", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  diagnostics <- active_load_balance(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(diagnostics$diagnostics_contract, "arbitrary_facet_planning_active_branch_load_balance")
  expect_true(diagnostics$branch_available)
  expect_true(is.data.frame(diagnostics$metric_registry))
  expect_true(is.data.frame(diagnostics$diagnostic_summary_table))
  expect_true(is.data.frame(diagnostics$diagnostic_table))
  expect_true(all(c(
    "design_id", "expected_observations_per_rater",
    "expected_person_assignments_per_rater", "observation_load_floor",
    "observation_load_ceiling", "observation_load_remainder",
    "perfect_integer_observation_balance", "recommended"
  ) %in% names(diagnostics$diagnostic_table)))
  expect_equal(
    diagnostics$metric_registry$metric,
    c(
      "expected_observations_per_rater",
      "expected_person_assignments_per_rater",
      "observation_load_floor",
      "observation_load_ceiling",
      "observation_load_remainder",
      "perfect_integer_observation_balance"
    )
  )
  expect_equal(
    diagnostics$metric_registry$basis_class,
    c(
      "balanced_expectation",
      "balanced_expectation",
      "integer_balance_bound",
      "integer_balance_bound",
      "exact_identity",
      "exact_identity"
    )
  )
  row_f01 <- diagnostics$diagnostic_table[
    diagnostics$diagnostic_table$design_id == "F01", , drop = FALSE
  ]
  expect_identical(nrow(row_f01), 1L)
  expect_equal(row_f01$expected_observations_per_rater[[1]], 8)
  expect_equal(row_f01$expected_person_assignments_per_rater[[1]], 4)
  expect_equal(row_f01$observation_load_floor[[1]], 8)
  expect_equal(row_f01$observation_load_ceiling[[1]], 8)
  expect_equal(row_f01$observation_load_remainder[[1]], 0)
  expect_equal(row_f01$perfect_integer_observation_balance[[1]], 1)
  expect_true(row_f01$recommended[[1]])
  summary_row <- diagnostics$diagnostic_summary_table[
    diagnostics$diagnostic_summary_table$metric == "expected_person_assignments_per_rater",
    , drop = FALSE
  ]
  expect_identical(nrow(summary_row), 1L)
  expect_equal(summary_row$recommended_value[[1]], 4)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_load_balance))
  expect_true(spec$planning_schema$future_branch_active_branch_load_balance$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_load_balance))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_load_balance$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_diagnostics <- active_load_balance(branch)
  expect_false(unavailable_diagnostics$branch_available)
  expect_true(is.data.frame(unavailable_diagnostics$metric_registry))
  expect_identical(nrow(unavailable_diagnostics$diagnostic_table), 0L)
})

test_that("future branch active load-balance overview summarizes divisibility patterns", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  load_balance_overview <- getFromNamespace("simulation_future_branch_active_branch_load_balance_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- load_balance_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    overview$overview_contract,
    "arbitrary_facet_planning_active_branch_load_balance_overview"
  )
  expect_true(overview$branch_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$metric_registry))
  expect_true(is.data.frame(overview$diagnostic_summary_table))
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$n_metrics[[1]], 6L)
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(overview$headline_table$n_perfect_integer_balance[[1]], 4L)
  expect_identical(overview$headline_table$n_nondivisible_designs[[1]], 4L)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_load_balance_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_load_balance_overview$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_load_balance_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_load_balance_overview$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- load_balance_overview(branch)
  expect_false(unavailable_overview$branch_available)
  expect_identical(unavailable_overview$headline_table$n_designs[[1]], 0L)
  expect_identical(nrow(unavailable_overview$diagnostic_summary_table), 0L)
})

test_that("future branch active coverage diagnostics summarize scored-cell overlap structure", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_coverage <- getFromNamespace("simulation_future_branch_active_branch_coverage", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  diagnostics <- active_coverage(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(diagnostics$diagnostics_contract, "arbitrary_facet_planning_active_branch_coverage")
  expect_true(diagnostics$branch_available)
  expect_true(is.data.frame(diagnostics$metric_registry))
  expect_true(is.data.frame(diagnostics$diagnostic_summary_table))
  expect_true(is.data.frame(diagnostics$diagnostic_table))
  expect_true(all(c(
    "design_id", "person_criterion_cells", "criterion_replications_per_person",
    "rater_pair_overlap_per_cell", "total_rater_pair_overlaps",
    "pair_coverage_fraction_per_cell", "redundant_scoring", "recommended"
  ) %in% names(diagnostics$diagnostic_table)))
  expect_equal(
    diagnostics$metric_registry$metric,
    c(
      "person_criterion_cells",
      "criterion_replications_per_person",
      "rater_pair_overlap_per_cell",
      "total_rater_pair_overlaps",
      "pair_coverage_fraction_per_cell",
      "redundant_scoring"
    )
  )
  expect_equal(
    diagnostics$metric_registry$basis_class,
    c(
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_identity",
      "exact_ratio",
      "exact_identity"
    )
  )
  row_f01 <- diagnostics$diagnostic_table[
    diagnostics$diagnostic_table$design_id == "F01", , drop = FALSE
  ]
  expect_identical(nrow(row_f01), 1L)
  expect_equal(row_f01$person_criterion_cells[[1]], 24)
  expect_equal(row_f01$criterion_replications_per_person[[1]], 1)
  expect_equal(row_f01$rater_pair_overlap_per_cell[[1]], 0)
  expect_equal(row_f01$total_rater_pair_overlaps[[1]], 0)
  expect_equal(row_f01$pair_coverage_fraction_per_cell[[1]], 0)
  expect_equal(row_f01$redundant_scoring[[1]], 0)
  expect_true(row_f01$recommended[[1]])
  summary_row <- diagnostics$diagnostic_summary_table[
    diagnostics$diagnostic_summary_table$metric == "person_criterion_cells",
    , drop = FALSE
  ]
  expect_identical(nrow(summary_row), 1L)
  expect_equal(summary_row$recommended_value[[1]], 24)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_coverage))
  expect_true(spec$planning_schema$future_branch_active_branch_coverage$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_coverage))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_coverage$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_diagnostics <- active_coverage(branch)
  expect_false(unavailable_diagnostics$branch_available)
  expect_true(is.data.frame(unavailable_diagnostics$metric_registry))
  expect_identical(nrow(unavailable_diagnostics$diagnostic_table), 0L)
})

test_that("future branch active coverage overview summarizes overlap patterns", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  coverage_overview <- getFromNamespace("simulation_future_branch_active_branch_coverage_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- coverage_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    overview$overview_contract,
    "arbitrary_facet_planning_active_branch_coverage_overview"
  )
  expect_true(overview$branch_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$metric_registry))
  expect_true(is.data.frame(overview$diagnostic_summary_table))
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$n_metrics[[1]], 6L)
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(overview$headline_table$n_redundant_designs[[1]], 4L)
  expect_identical(overview$headline_table$n_single_rater_designs[[1]], 4L)
  expect_identical(overview$headline_table$n_pair_connected_designs[[1]], 4L)
  expect_identical(overview$headline_table$n_zero_pair_overlap_designs[[1]], 4L)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_coverage_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_coverage_overview$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_coverage_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_coverage_overview$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- coverage_overview(branch)
  expect_false(unavailable_overview$branch_available)
  expect_identical(unavailable_overview$headline_table$n_designs[[1]], 0L)
  expect_identical(nrow(unavailable_overview$diagnostic_summary_table), 0L)
})

test_that("future branch active guardrails classify exact design regimes", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_guardrails <- getFromNamespace("simulation_future_branch_active_branch_guardrails", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  guardrails <- active_guardrails(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(guardrails$guardrail_contract, "arbitrary_facet_planning_active_branch_guardrails")
  expect_true(guardrails$branch_available)
  expect_true(is.data.frame(guardrails$guardrail_registry))
  expect_true(is.data.frame(guardrails$guardrail_summary_table))
  expect_true(is.data.frame(guardrails$guardrail_table))
  expect_equal(
    guardrails$guardrail_registry$guardrail,
    c(
      "rater_linking_regime",
      "pair_coverage_regime",
      "integer_balance_regime",
      "redundancy_regime"
    )
  )
  expect_true(all(guardrails$guardrail_registry$basis_class == "exact_classification"))
  row_f01 <- guardrails$guardrail_table[
    guardrails$guardrail_table$design_id == "F01", , drop = FALSE
  ]
  expect_identical(nrow(row_f01), 1L)
  expect_identical(row_f01$rater_linking_regime[[1]], "single_rater")
  expect_identical(row_f01$pair_coverage_regime[[1]], "no_pair_coverage")
  expect_identical(row_f01$integer_balance_regime[[1]], "integer_balanced")
  expect_identical(row_f01$redundancy_regime[[1]], "single_rater_only")
  expect_true(row_f01$recommended[[1]])
  summary_row <- guardrails$guardrail_summary_table[
    guardrails$guardrail_summary_table$guardrail == "rater_linking_regime", , drop = FALSE
  ]
  expect_identical(nrow(summary_row), 1L)
  expect_identical(summary_row$recommended_level[[1]], "single_rater")
  expect_true(grepl("partial_overlap", summary_row$observed_levels[[1]], fixed = TRUE))
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_guardrails))
  expect_true(spec$planning_schema$future_branch_active_branch_guardrails$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_guardrails))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_guardrails$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_guardrails <- active_guardrails(branch)
  expect_false(unavailable_guardrails$branch_available)
  expect_true(is.data.frame(unavailable_guardrails$guardrail_registry))
  expect_identical(nrow(unavailable_guardrails$guardrail_table), 0L)
})

test_that("future branch active guardrail overview summarizes regime counts", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  guardrail_overview <- getFromNamespace("simulation_future_branch_active_branch_guardrail_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- guardrail_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    overview$overview_contract,
    "arbitrary_facet_planning_active_branch_guardrail_overview"
  )
  expect_true(overview$branch_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$guardrail_registry))
  expect_true(is.data.frame(overview$guardrail_summary_table))
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$n_guardrails[[1]], 4L)
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(overview$headline_table$n_single_rater_designs[[1]], 4L)
  expect_identical(overview$headline_table$n_partial_overlap_designs[[1]], 4L)
  expect_identical(overview$headline_table$n_fully_crossed_designs[[1]], 0L)
  expect_identical(overview$headline_table$n_integer_balanced_designs[[1]], 4L)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_guardrail_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_guardrail_overview$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_guardrail_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_guardrail_overview$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- guardrail_overview(branch)
  expect_false(unavailable_overview$branch_available)
  expect_identical(unavailable_overview$headline_table$n_designs[[1]], 0L)
  expect_identical(nrow(unavailable_overview$guardrail_summary_table), 0L)
})

test_that("future branch active readiness summarizes exact overlap preconditions", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_readiness <- getFromNamespace("simulation_future_branch_active_branch_readiness", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  readiness <- active_readiness(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(readiness$readiness_contract, "arbitrary_facet_planning_active_branch_readiness")
  expect_true(readiness$branch_available)
  expect_true(is.data.frame(readiness$indicator_registry))
  expect_true(is.data.frame(readiness$readiness_summary_table))
  expect_true(is.data.frame(readiness$readiness_table))
  expect_equal(
    readiness$indicator_registry$indicator,
    c(
      "supports_multi_rater_cells",
      "supports_pair_overlap",
      "supports_integer_balanced_load",
      "supports_full_pair_coverage"
    )
  )
  expect_true(all(readiness$indicator_registry$basis_class == "exact_indicator"))
  row_f01 <- readiness$readiness_table[
    readiness$readiness_table$design_id == "F01", , drop = FALSE
  ]
  expect_identical(nrow(row_f01), 1L)
  expect_equal(row_f01$supports_multi_rater_cells[[1]], 0)
  expect_equal(row_f01$supports_pair_overlap[[1]], 0)
  expect_equal(row_f01$supports_integer_balanced_load[[1]], 1)
  expect_equal(row_f01$supports_full_pair_coverage[[1]], 0)
  expect_identical(row_f01$structural_tier[[1]], "single_rater_only")
  expect_true(row_f01$recommended[[1]])
  summary_row <- readiness$readiness_summary_table[
    readiness$readiness_summary_table$indicator == "supports_integer_balanced_load",
    , drop = FALSE
  ]
  expect_identical(nrow(summary_row), 1L)
  expect_equal(summary_row$recommended_value[[1]], 1)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_readiness))
  expect_true(spec$planning_schema$future_branch_active_branch_readiness$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_readiness))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_readiness$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_readiness <- active_readiness(branch)
  expect_false(unavailable_readiness$branch_available)
  expect_true(is.data.frame(unavailable_readiness$indicator_registry))
  expect_identical(nrow(unavailable_readiness$readiness_table), 0L)
})

test_that("future branch active readiness overview summarizes exact structural tiers", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  readiness_overview <- getFromNamespace("simulation_future_branch_active_branch_readiness_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- readiness_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    overview$overview_contract,
    "arbitrary_facet_planning_active_branch_readiness_overview"
  )
  expect_true(overview$branch_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$indicator_registry))
  expect_true(is.data.frame(overview$readiness_summary_table))
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$n_indicators[[1]], 4L)
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F01")
  expect_identical(overview$headline_table$n_single_rater_only_tiers[[1]], 4L)
  expect_identical(overview$headline_table$n_partial_overlap_balanced_tiers[[1]], 2L)
  expect_identical(overview$headline_table$n_partial_overlap_unbalanced_tiers[[1]], 2L)
  expect_identical(overview$headline_table$n_full_overlap_tiers[[1]], 0L)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_readiness_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_readiness_overview$branch_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_readiness_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_readiness_overview$branch_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- readiness_overview(branch)
  expect_false(unavailable_overview$branch_available)
  expect_identical(unavailable_overview$headline_table$n_designs[[1]], 0L)
  expect_identical(nrow(unavailable_overview$readiness_summary_table), 0L)
})

test_that("future branch active recommendation ranks designs by exact structural tier then counts", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  active_recommendation <- getFromNamespace("simulation_future_branch_active_branch_recommendation", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  recommendation <- active_recommendation(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    recommendation$recommendation_contract,
    "arbitrary_facet_planning_active_branch_recommendation"
  )
  expect_true(recommendation$recommendation_available)
  expect_identical(recommendation$recommended_design_id, "F05")
  expect_identical(recommendation$recommended_tier, "partial_overlap_balanced")
  expect_true(is.data.frame(recommendation$recommendation_table))
  row_f05 <- recommendation$recommendation_table[
    recommendation$recommendation_table$design_id == "F05", , drop = FALSE
  ]
  expect_identical(nrow(row_f05), 1L)
  expect_true(row_f05$recommended[[1]])
  expect_identical(row_f05$structural_priority[[1]], 3L)
  expect_equal(row_f05$total_observations[[1]], 48)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_recommendation))
  expect_true(spec$planning_schema$future_branch_active_branch_recommendation$recommendation_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_recommendation))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_recommendation$recommendation_available)
  expect_true("recommendation_table" %in% names(spec$planning_schema$future_branch_schema$active_branch_recommendation))
  expect_false("active_branch_readiness" %in% names(spec$planning_schema$future_branch_schema$active_branch_recommendation))
  expect_false("active_branch_profile" %in% names(spec$planning_schema$future_branch_schema$active_branch_recommendation))

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_recommendation <- active_recommendation(branch)
  expect_false(unavailable_recommendation$recommendation_available)
  expect_identical(nrow(unavailable_recommendation$recommendation_table), 0L)
})

test_that("future branch active recommendation overview summarizes conservative selection", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  recommendation_overview <- getFromNamespace("simulation_future_branch_active_branch_recommendation_overview", "mfrmr")
  design_override <- list(
    facets = list(person = c(12, 14), judge = 3, task = c(2, 4)),
    assignment = c(1, 2)
  )

  overview <- recommendation_overview(
    spec$planning_schema,
    design = design_override,
    prefer = c("assignment", "task", "person"),
    view = "public",
    mode = "snapshot",
    surface = "metrics",
    table_component = "report",
    x_var = "assignment",
    group_var = "task"
  )
  expect_identical(
    overview$overview_contract,
    "arbitrary_facet_planning_active_branch_recommendation_overview"
  )
  expect_true(overview$recommendation_available)
  expect_true(is.data.frame(overview$headline_table))
  expect_true(is.data.frame(overview$recommendation_table))
  expect_identical(overview$headline_table$n_designs[[1]], 8L)
  expect_identical(overview$headline_table$recommended_design_id[[1]], "F05")
  expect_identical(overview$headline_table$recommended_tier[[1]], "partial_overlap_balanced")
  expect_identical(overview$headline_table$n_top_tier_candidates[[1]], 2L)
  expect_true(is.list(spec$planning_schema$future_branch_active_branch_recommendation_overview))
  expect_true(spec$planning_schema$future_branch_active_branch_recommendation_overview$recommendation_available)
  expect_true(is.list(spec$planning_schema$future_branch_schema$active_branch_recommendation_overview))
  expect_true(spec$planning_schema$future_branch_schema$active_branch_recommendation_overview$recommendation_available)

  branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")(facet_names = c("Judge", "Task"))
  unavailable_overview <- recommendation_overview(branch)
  expect_false(unavailable_overview$recommendation_available)
  expect_identical(nrow(unavailable_overview$recommendation_table), 0L)
})

test_that("future branch grid can be materialized from live planning objects", {
  spec <- build_mfrm_sim_spec(
    facet_names = c("Judge", "Task"),
    design = list(facets = c(person = 14, judge = 3, task = 2), assignment = 2),
    assignment = "rotating"
  )

  materialize_branch_grid <- getFromNamespace("simulation_materialize_future_branch_grid", "mfrmr")

  from_spec <- materialize_branch_grid(spec)
  expect_identical(
    from_spec$canonical,
    spec$planning_schema$future_branch_grid_contract$canonical
  )

  from_schema <- materialize_branch_grid(
    spec$planning_schema,
    design = list(facets = list(person = 16, judge = 4, task = 2), assignment = 3)
  )
  expect_true(is.data.frame(from_schema$canonical))
  expect_equal(sort(unique(from_schema$canonical$n_person)), 16)
  expect_equal(sort(unique(from_schema$canonical$n_rater)), 4)
  expect_equal(sort(unique(from_schema$canonical$raters_per_person)), 3)
  expect_true(is.list(from_schema$grid_contract))

  wrapped <- list(settings = list(planning_schema = spec$planning_schema))
  from_wrapped <- materialize_branch_grid(wrapped)
  expect_identical(
    from_wrapped$branch,
    spec$planning_schema$future_branch_grid_contract$branch
  )
})

test_that("future branch preview reports unavailable defaults for facet-name-only schemas", {
  build_branch <- getFromNamespace("simulation_future_branch_schema", "mfrmr")

  branch <- build_branch(facet_names = c("Judge", "Task"))
  expect_true(is.list(branch$preview))
  expect_false(branch$preview$preview_available)
  expect_true(grepl("does not yet have complete default facet counts", branch$preview$reason, fixed = TRUE))
  expect_true(is.list(branch$grid_contract))
  expect_false(branch$grid_contract$preview_available)
  expect_identical(branch$grid_contract$contract, "arbitrary_facet_design_grid_contract")
  expect_true(is.list(branch$grid_bundle))
  expect_false(branch$grid_bundle$grid_available)
  expect_identical(branch$grid_bundle$bundle_contract, "arbitrary_facet_design_grid_bundle")
  expect_true(is.list(branch$grid_context))
  expect_false(branch$grid_context$grid_available)
  expect_identical(branch$grid_context$context_contract, "arbitrary_facet_design_grid_context")
})

test_that("build_mfrm_sim_spec can store a latent-regression population generator", {
    population_covariates <- data.frame(
    TemplatePerson = sprintf("TP%02d", 1:20),
    X = seq(-1, 1, length.out = 20),
    G = rep(c("A", "B"), each = 10),
    stringsAsFactors = FALSE
  )

  spec <- build_mfrm_sim_spec(
    n_person = 24,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    population_formula = ~ X + G,
    population_coefficients = c(`(Intercept)` = 0.2, X = 0.8, GB = 0.5),
    population_sigma2 = 0.15,
    population_covariates = population_covariates
  )

  expect_true(is.list(spec$population))
  expect_true(isTRUE(spec$population$active))
  expect_identical(names(spec$population$coefficients), c("(Intercept)", "X", "GB"))
  expect_equal(spec$population$sigma2, 0.15)
  expect_true(all(c("TemplatePerson", "X", "G") %in% names(spec$population$covariate_template)))
  expect_identical(spec$population$xlevels$G, c("A", "B"))
  expect_true("G" %in% names(spec$population$contrasts))
})

test_that("design-grid helper returns canonical and public design metadata", {
  spec <- build_mfrm_sim_spec(
    n_person = 16,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )
  build_grid <- getFromNamespace("simulation_build_design_grid", "mfrmr")

  grid_meta <- build_grid(
    n_person = c(12, 14),
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    sim_spec = spec,
    id_prefix = "X"
  )

  expect_named(grid_meta, c("canonical", "public", "aliases", "descriptor"))
  expect_identical(names(grid_meta$canonical), c("design_id", "n_person", "n_rater", "n_criterion", "raters_per_person"))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(grid_meta$public)))
  expect_identical(grid_meta$descriptor$role, c("person", "rater", "criterion", "assignment"))
  expect_identical(grid_meta$descriptor$alias, c("n_person", "n_judge", "n_task", "judge_per_person"))
  expect_identical(grid_meta$canonical$design_id, c("X01", "X02"))

  expect_error(
    build_grid(
      n_person = 12,
      n_rater = 3,
      n_criterion = 2,
      raters_per_person = 4,
      sim_spec = spec,
      id_prefix = "X"
    ),
    "judge_per_person <= n_judge",
    fixed = TRUE
  )
})

test_that("manual custom-name specs normalize assignment profiles and design skeleton inputs", {
  assignment_profiles <- data.frame(
    TemplatePerson = c("TP1", "TP1", "TP2", "TP2", "TP3", "TP3"),
    Judge = c("J01", "J02", "J02", "J03", "J01", "J03"),
    stringsAsFactors = FALSE
  )
  resampled_spec <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    facet_names = c("Judge", "Task"),
    assignment = "resampled",
    assignment_profiles = assignment_profiles
  )

  expect_true(all(c("TemplatePerson", "Rater") %in% names(resampled_spec$assignment_profiles)))
  expect_setequal(unique(resampled_spec$assignment_profiles$Rater), c("J01", "J02", "J03"))

  design_skeleton <- expand.grid(
    TemplatePerson = c("TP1", "TP2"),
    Judge = c("J01", "J02", "J03"),
    Task = c("T01", "T02"),
    stringsAsFactors = FALSE
  )
  skeleton_spec <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 3,
    facet_names = c("Judge", "Task"),
    assignment = "skeleton",
    design_skeleton = design_skeleton
  )

  expect_true(all(c("TemplatePerson", "Rater", "Criterion") %in% names(skeleton_spec$design_skeleton)))
  expect_setequal(unique(skeleton_spec$design_skeleton$Rater), c("J01", "J02", "J03"))
  expect_setequal(unique(skeleton_spec$design_skeleton$Criterion), c("T01", "T02"))
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

test_that("RSM simulation reduces exactly to PCM simulation under common thresholds", {
  base_args <- list(
    n_person = 24,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 4,
    theta_sd = 0.8,
    rater_sd = 0.25,
    criterion_sd = 0.2,
    noise_sd = 0,
    assignment = "rotating",
    thresholds = c(-1.1, 0.1, 1.0),
    step_facet = "Criterion"
  )
  spec_rsm <- do.call(build_mfrm_sim_spec, c(base_args, list(model = "RSM")))
  spec_pcm <- do.call(build_mfrm_sim_spec, c(base_args, list(model = "PCM")))

  sim_rsm <- simulate_mfrm_data(sim_spec = spec_rsm, seed = 3141)
  sim_pcm <- simulate_mfrm_data(sim_spec = spec_pcm, seed = 3141)

  visible_cols <- c("Study", "Person", "Rater", "Criterion", "Score")
  expect_identical(as.list(sim_rsm[visible_cols]), as.list(sim_pcm[visible_cols]))

  truth_rsm <- attr(sim_rsm, "mfrm_truth")
  truth_pcm <- attr(sim_pcm, "mfrm_truth")
  expect_equal(truth_rsm$person, truth_pcm$person, tolerance = 1e-12)
  expect_equal(truth_rsm$facets, truth_pcm$facets, tolerance = 1e-12)
  expect_equal(truth_rsm$steps, truth_pcm$steps, tolerance = 1e-12)
  expect_equal(truth_rsm$step_table, truth_pcm$step_table)

  step_cum <- c(0, cumsum(truth_rsm$steps))
  step_cum_mat <- matrix(
    rep(step_cum, times = length(unique(sim_rsm$Criterion))),
    nrow = length(unique(sim_rsm$Criterion)),
    byrow = TRUE
  )
  step_levels <- sort(unique(as.character(sim_rsm$Criterion)))
  criterion_idx <- match(as.character(sim_rsm$Criterion), step_levels)
  eta <- unname(
    truth_rsm$person[sim_rsm$Person] -
      truth_rsm$facets$Rater[sim_rsm$Rater] -
      truth_rsm$facets$Criterion[sim_rsm$Criterion]
  )

  probs_rsm <- mfrmr:::category_prob_rsm(eta, step_cum)
  probs_pcm <- mfrmr:::category_prob_pcm(
    eta = eta,
    step_cum_mat = step_cum_mat,
    criterion_idx = criterion_idx
  )
  expect_equal(unname(probs_rsm), unname(probs_pcm), tolerance = 1e-12)
})

test_that("GPCM direct and sim-spec generators carry slope-aware truth metadata", {
  sim_direct <- simulate_mfrm_data(
    n_person = 18,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 4,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 0.8, C02 = 1.0, C03 = 1.25),
    seed = 717
  )

  truth_direct <- attr(sim_direct, "mfrm_truth")
  spec_direct <- attr(sim_direct, "mfrm_simulation_spec")

  expect_true(all(sim_direct$Score %in% seq_len(4)))
  expect_true(is.data.frame(truth_direct$slope_table))
  expect_equal(spec_direct$model, "GPCM")
  expect_equal(spec_direct$slope_facet, "Criterion")
  expect_setequal(truth_direct$slope_table$SlopeFacet, c("C01", "C02", "C03"))

  spec <- build_mfrm_sim_spec(
    n_person = 18,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 4,
    assignment = "rotating",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 0.8, C02 = 1.0, C03 = 1.25)
  )

  expect_equal(spec$model, "GPCM")
  expect_equal(spec$step_facet, "Criterion")
  expect_equal(spec$slope_facet, "Criterion")
  expect_true(is.data.frame(spec$slope_table))
  expect_true(all(spec$slope_table$Estimate > 0))
  expect_equal(exp(mean(log(spec$slope_table$Estimate))), 1, tolerance = 1e-12)

  sim_spec <- simulate_mfrm_data(sim_spec = spec, seed = 718)
  truth_spec <- attr(sim_spec, "mfrm_truth")

  expect_true(all(sim_spec$Score %in% seq_len(spec$score_levels)))
  expect_true(is.data.frame(truth_spec$slope_table))
  expect_equal(attr(sim_spec, "mfrm_simulation_spec")$model, "GPCM")
  expect_equal(attr(sim_spec, "mfrm_simulation_spec")$slope_facet, "Criterion")
  expect_equal(truth_spec$slope_table$Estimate, spec$slope_table$Estimate, tolerance = 1e-12)
})

test_that("unit-slope GPCM simulation reduces exactly to PCM simulation", {
  thresholds <- data.frame(
    StepFacet = rep(c("C01", "C02", "C03"), each = 3),
    StepIndex = rep(1:3, times = 3),
    Estimate = c(
      -1.2, -0.1, 1.1,
      -0.8,  0.3, 1.3,
      -1.5,  0.0, 0.9
    )
  )

  base_args <- list(
    n_person = 20,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 4,
    theta_sd = 0.7,
    rater_sd = 0.2,
    criterion_sd = 0.15,
    noise_sd = 0,
    assignment = "rotating",
    thresholds = thresholds,
    step_facet = "Criterion"
  )
  spec_pcm <- do.call(build_mfrm_sim_spec, c(base_args, list(model = "PCM")))
  spec_gpcm <- do.call(
    build_mfrm_sim_spec,
    c(
      base_args,
      list(
        model = "GPCM",
        slope_facet = "Criterion",
        slopes = c(C01 = 1, C02 = 1, C03 = 1)
      )
    )
  )

  sim_pcm <- simulate_mfrm_data(sim_spec = spec_pcm, seed = 2718)
  sim_gpcm <- simulate_mfrm_data(sim_spec = spec_gpcm, seed = 2718)

  visible_cols <- c("Study", "Person", "Rater", "Criterion", "Score")
  expect_identical(as.list(sim_gpcm[visible_cols]), as.list(sim_pcm[visible_cols]))

  truth_pcm <- attr(sim_pcm, "mfrm_truth")
  truth_gpcm <- attr(sim_gpcm, "mfrm_truth")
  expect_equal(truth_gpcm$person, truth_pcm$person, tolerance = 1e-12)
  expect_equal(truth_gpcm$facets, truth_pcm$facets, tolerance = 1e-12)
  expect_equal(truth_gpcm$step_table, truth_pcm$step_table)
  expect_equal(truth_gpcm$slope_table$Estimate, rep(1, 3), tolerance = 1e-12)

  step_table <- truth_pcm$step_table
  step_lookup <- split(step_table$Estimate, step_table$StepFacet)
  step_levels <- unique(as.character(step_table$StepFacet))
  step_cum_mat <- t(vapply(
    step_levels,
    function(level) c(0, cumsum(step_lookup[[level]])),
    numeric(spec_pcm$score_levels)
  ))
  criterion_idx <- match(as.character(sim_pcm$Criterion), step_levels)
  eta <- unname(
    truth_pcm$person[sim_pcm$Person] -
      truth_pcm$facets$Rater[sim_pcm$Rater] -
      truth_pcm$facets$Criterion[sim_pcm$Criterion]
  )

  probs_pcm <- mfrmr:::category_prob_pcm(
    eta = eta,
    step_cum_mat = step_cum_mat,
    criterion_idx = criterion_idx
  )
  probs_gpcm <- mfrmr:::category_prob_gpcm(
    eta = eta,
    step_cum_mat = step_cum_mat,
    criterion_idx = criterion_idx,
    slopes = rep(1, length(step_levels)),
    slope_idx = criterion_idx
  )

  expect_equal(unname(probs_gpcm), unname(probs_pcm), tolerance = 1e-12)
})

test_that("GPCM simulation slopes are normalized to the identified scale", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 2,
    n_criterion = 3,
    raters_per_person = 2,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 0.6, C02 = 1.1, C03 = 1.8)
  )

  expect_equal(exp(mean(log(spec$slope_table$Estimate))), 1, tolerance = 1e-12)
  expect_equal(
    spec$slope_table$Estimate,
    exp(log(c(0.6, 1.1, 1.8)) - mean(log(c(0.6, 1.1, 1.8)))),
    tolerance = 1e-12
  )

  spec_df <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 2,
    n_criterion = 3,
    raters_per_person = 2,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = data.frame(
      SlopeFacet = c("C03", "C01", "C02"),
      Estimate = c(1.8, 0.6, 1.1)
    )
  )
  expect_equal(spec_df$slope_table$Estimate, spec$slope_table$Estimate, tolerance = 1e-12)
})

test_that("simulate_mfrm_data accepts role-based design input in the direct path", {
  sim <- simulate_mfrm_data(
    design = list(person = 12, rater = 3, criterion = 2, assignment = 2),
    seed = 1818
  )

  expect_true(is.data.frame(sim))
  expect_equal(length(unique(sim$Person)), 12)
  expect_equal(length(unique(sim$Rater)), 3)
  expect_equal(length(unique(sim$Criterion)), 2)

  expect_error(
    simulate_mfrm_data(
      n_person = 10,
      design = list(person = 12),
      seed = 1818
    ),
    "Do not supply the same design variable through both scalar arguments and `design`",
    fixed = TRUE
  )
})

test_that("simulate_mfrm_data generates person covariates under latent regression", {
  population_covariates <- data.frame(
    TemplatePerson = sprintf("TP%02d", 1:30),
    X = seq(-1.5, 1.5, length.out = 30),
    stringsAsFactors = FALSE
  )
  spec <- build_mfrm_sim_spec(
    n_person = 36,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating",
    population_formula = ~ X,
    population_coefficients = c(`(Intercept)` = 0.1, X = 1.1),
    population_sigma2 = 0.1,
    population_covariates = population_covariates
  )

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 5151)
  truth <- attr(sim, "mfrm_truth")
  pop_data <- attr(sim, "mfrm_population_data")

  expect_true(is.list(pop_data))
  expect_true(isTRUE(pop_data$active))
  expect_true(is.data.frame(pop_data$person_data))
  expect_true(all(c("Person", "X") %in% names(pop_data$person_data)))
  expect_equal(nrow(pop_data$person_data), length(unique(sim$Person)))
  expect_gt(stats::cor(unname(truth$person[pop_data$person_data$Person]), pop_data$person_data$X), 0.6)
  expect_equal(as.character(pop_data$population_formula), c("~", "X"))
})

test_that("simulate_mfrm_data preserves categorical population coding", {
  population_covariates <- data.frame(
    TemplatePerson = c("TP01", "TP02"),
    Group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  spec <- build_mfrm_sim_spec(
    n_person = 2,
    n_rater = 2,
    n_criterion = 2,
    raters_per_person = 2,
    score_levels = 2,
    assignment = "crossed",
    population_formula = ~ Group,
    population_coefficients = c(`(Intercept)` = 0, GroupB = 1),
    population_sigma2 = 0,
    population_covariates = population_covariates
  )

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 5152)
  truth <- attr(sim, "mfrm_truth")
  pop_data <- attr(sim, "mfrm_population_data")

  expect_true(isTRUE(spec$population$active))
  expect_identical(spec$population$xlevels$Group, c("A", "B"))
  expect_true("Group" %in% names(spec$population$contrasts))
  one_person <- mfrmr:::simulation_generate_population_person_data(spec$population, "P001")
  expect_true(is.factor(one_person$person_data$Group))
  expect_identical(levels(one_person$person_data$Group), c("A", "B"))
  expect_true(is.factor(pop_data$person_data$Group))
  expect_identical(levels(pop_data$person_data$Group), c("A", "B"))
  expect_identical(pop_data$xlevels$Group, c("A", "B"))
  expect_true("Group" %in% names(pop_data$contrasts))
  expect_identical(truth$population$xlevels$Group, c("A", "B"))
  expect_true("Group" %in% names(truth$population$contrasts))
  expect_identical(as.character(pop_data$population_formula), c("~", "Group"))
  expect_identical(pop_data$design_columns, c("(Intercept)", "GroupB"))
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

test_that("extract_mfrm_sim_spec captures bounded GPCM slope metadata", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:14]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit_gpcm <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = 5,
      maxit = 20
    )
  )

  spec <- extract_mfrm_sim_spec(fit_gpcm)

  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$model, "GPCM")
  expect_equal(spec$step_facet, "Criterion")
  expect_equal(spec$slope_facet, "Criterion")
  expect_true(is.data.frame(spec$slope_table))
  expect_true(all(c("SlopeFacet", "Estimate") %in% names(spec$slope_table)))
  expect_true(all(is.finite(spec$slope_table$Estimate)))
  expect_true(all(spec$slope_table$Estimate > 0))
  expect_setequal(spec$slope_table$SlopeFacet, unique(toy$Criterion))
})

test_that("fit-derived GPCM specs generate data but remain blocked in planning helpers", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:14]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit_gpcm <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = 5,
      maxit = 20
    )
  )
  spec <- extract_mfrm_sim_spec(fit_gpcm)

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 719)
  truth <- attr(sim, "mfrm_truth")
  expect_true(all(sim$Score %in% seq_len(spec$score_levels)))
  expect_true(is.data.frame(truth$slope_table))
  expect_equal(attr(sim, "mfrm_simulation_spec")$model, "GPCM")

  expect_error(
    evaluate_mfrm_design(
      n_person = spec$n_person,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 720
    ),
    "`evaluate_mfrm_design()` does not yet support bounded `GPCM` simulation specifications.",
    fixed = TRUE
  )
  expect_error(
    evaluate_mfrm_design(
      n_person = spec$n_person,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 720
    ),
    "role-based person x rater-like x criterion-like",
    fixed = TRUE
  )
  expect_error(
    evaluate_mfrm_signal_detection(
      n_person = spec$n_person,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 721
    ),
    "`evaluate_mfrm_signal_detection()` does not yet support bounded `GPCM` simulation specifications.",
    fixed = TRUE
  )
  expect_error(
    evaluate_mfrm_signal_detection(
      n_person = spec$n_person,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 721
    ),
    "role-based person x rater-like x criterion-like",
    fixed = TRUE
  )
  expect_error(
    predict_mfrm_population(sim_spec = spec, n_person = spec$n_person, reps = 1, seed = 722),
    "`predict_mfrm_population()` is not yet validated for `GPCM` simulation specifications.",
    fixed = TRUE
  )
  expect_error(
    predict_mfrm_population(sim_spec = spec, n_person = spec$n_person, reps = 1, seed = 722),
    "role-based person x rater-like x criterion-like",
    fixed = TRUE
  )
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
  spec_n30 <- mfrmr:::simulation_override_spec_design(
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
  expect_identical(spec$planning_constraints$mutable_design_variables, c("n_person", "n_criterion"))
  expect_identical(spec$planning_constraints$locked_design_variables, c("n_rater", "raters_per_person"))
  expect_true(all(c("n_rater", "raters_per_person") %in% names(spec$planning_constraints$lock_reasons)))

  expect_error(
    mfrmr:::simulation_override_spec_design(
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
  spec_n30 <- mfrmr:::simulation_override_spec_design(
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
  spec_n24 <- mfrmr:::simulation_override_spec_design(
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
  expect_identical(spec$planning_constraints$mutable_design_variables, "n_person")
  expect_identical(spec$planning_constraints$locked_design_variables, c("n_rater", "n_criterion", "raters_per_person"))

  expect_error(
    mfrmr:::simulation_override_spec_design(
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

test_that("fit-derived specs preserve custom facet names across planning helpers", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  names(toy)[names(toy) == "Rater"] <- "Judge"
  names(toy)[names(toy) == "Criterion"] <- "Task"

  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Judge", "Task"), "Score", method = "JML", maxit = 15)
  )

  spec <- extract_mfrm_sim_spec(fit)
  expect_identical(unname(spec$facet_names), c("Judge", "Task"))

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 905)
  expect_named(sim, c("Study", "Person", "Judge", "Task", "Score"))
  expect_setequal(unique(sim$Task), unique(toy$Task))

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = c(12, 14),
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 906
    )
  )
  expect_true(all(c("Person", "Judge", "Task") %in% unique(sim_eval$results$Facet)))

  p_custom <- plot(sim_eval, facet = "Judge", metric = "separation", x_var = "n_person", draw = FALSE)
  expect_equal(p_custom$facet, "Judge")

  rec <- recommend_mfrm_design(sim_eval)
  expect_true(all(c("Judge", "Task") %in% unique(rec$facet_table$Facet)))

  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      n_person = 14,
      reps = 1,
      maxit = 10,
      seed = 907
    )
  )
  expect_true(all(pred$forecast$Facet %in% c("Person", "Judge", "Task")))

  sig_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      n_person = 12,
      n_rater = spec$n_rater,
      n_criterion = spec$n_criterion,
      raters_per_person = spec$raters_per_person,
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 908
    )
  )
  expect_s3_class(sig_eval, "mfrm_signal_detection")
  expect_equal(nrow(sig_eval$results), 1)
  expect_true(sig_eval$results$DIFTargetLevel[1] %in% unique(toy$Task))
})

test_that("manual specs preserve custom facet names across planning helpers", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = c(12, 14),
      n_rater = c(spec$n_rater, spec$n_rater + 1L),
      n_criterion = c(spec$n_criterion, spec$n_criterion + 1L),
      raters_per_person = c(spec$raters_per_person, spec$raters_per_person + 1L),
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 910
    )
  )
  expect_true(all(c("Person", "Judge", "Task") %in% unique(sim_eval$results$Facet)))
  expect_identical(
    unname(sim_eval$settings$design_variable_aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")]),
    c("n_person", "n_judge", "n_task", "judge_per_person")
  )
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sim_eval$design_grid)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sim_eval$results)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sim_eval$rep_overview)))
  expect_true(is.data.frame(sim_eval$design_descriptor))
  expect_identical(sim_eval$design_descriptor$role, c("person", "rater", "criterion", "assignment"))
  expect_identical(sim_eval$design_descriptor$alias, c("n_person", "n_judge", "n_task", "judge_per_person"))
  expect_true(is.list(sim_eval$planning_scope))
  expect_identical(sim_eval$planning_scope$planner_contract, "role_based_two_non_person_facets")
  expect_identical(sim_eval$planning_scope$role_labels, c("Judge", "Task"))
  expect_false(sim_eval$planning_scope$supports_arbitrary_facet_planning)
  expect_true(sim_eval$planning_scope$supports_arbitrary_facet_estimation)
  expect_true(grepl("arbitrary facet counts", sim_eval$planning_scope$note, fixed = TRUE))
  expect_true(is.list(sim_eval$planning_constraints))
  expect_identical(sim_eval$planning_constraints$mutable_design_variables, c("n_person", "n_rater", "n_criterion", "raters_per_person"))
  expect_true(is.list(sim_eval$planning_schema))
  expect_true(is.data.frame(sim_eval$planning_schema$role_table))
  expect_identical(sim_eval$planning_schema$role_table$alias, c("n_person", "n_judge", "n_task", "judge_per_person"))
  expect_true(is.data.frame(sim_eval$planning_schema$facet_manifest))
  expect_identical(sim_eval$planning_schema$facet_manifest$facet, c("Person", "Judge", "Task"))
  expect_true(is.data.frame(sim_eval$planning_schema$future_facet_table))
  expect_identical(sim_eval$planning_schema$future_facet_table$future_facet_key, c("person", "judge", "task"))
  expect_true(is.list(sim_eval$planning_schema$future_branch_schema))
  expect_identical(sim_eval$planning_schema$future_branch_schema$input_contract, "design$facets(named counts)")
  expect_identical(sim_eval$planning_schema$future_branch_schema$design_template$assignment, 2L)
  expect_true(is.list(sim_eval$planning_schema$future_branch_preview))
  expect_true(sim_eval$planning_schema$future_branch_preview$preview_available)
  expect_true(is.list(sim_eval$planning_schema$future_branch_grid_semantics))
  expect_true(is.list(sim_eval$planning_schema$future_branch_grid_contract))

  s_custom <- summary(sim_eval)
  expect_equal(s_custom$facet_names[["rater"]], "Judge")
  expect_equal(s_custom$design_variable_aliases[["n_rater"]], "n_judge")
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(s_custom$design_summary)))
  expect_true(is.data.frame(s_custom$design_descriptor))
  expect_identical(s_custom$design_descriptor$canonical, c("n_person", "n_rater", "n_criterion", "raters_per_person"))
  expect_true(is.list(s_custom$planning_scope))
  expect_identical(s_custom$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(any(grepl("arbitrary facet counts", s_custom$notes, fixed = TRUE)))
  expect_true(is.list(s_custom$planning_constraints))
  expect_true(any(grepl("All current design variables remain mutable", s_custom$notes, fixed = TRUE)))
  expect_true(is.list(s_custom$planning_schema))
  expect_true(is.data.frame(s_custom$planning_schema$role_table))
  expect_true(is.list(s_custom$planning_schema$future_branch_preview))
  expect_true(is.list(s_custom$planning_schema$future_branch_grid_semantics))
  expect_true(is.list(s_custom$planning_schema$future_branch_grid_contract))

  p_custom <- plot(
    sim_eval,
    facet = "Judge",
    metric = "separation",
    x_var = "rater",
    group_var = "criterion",
    draw = FALSE
  )
  expect_equal(p_custom$facet, "Judge")
  expect_equal(p_custom$x_var, "n_rater")
  expect_equal(p_custom$x_label, "n_judge")
  expect_equal(p_custom$group_var, "n_criterion")
  expect_equal(p_custom$group_label, "n_task")
  expect_true(is.data.frame(p_custom$design_descriptor))
  expect_true(is.list(p_custom$planning_scope))
  expect_identical(p_custom$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(is.list(p_custom$planning_constraints))
  expect_true(is.list(p_custom$planning_schema))

  rec <- recommend_mfrm_design(
    sim_eval,
    prefer = c("rater", "assignment"),
    min_convergence_rate = 0
  )
  expect_identical(rec$thresholds$prefer, c("n_rater", "raters_per_person"))
  expect_equal(rec$design_variable_aliases[["n_rater"]], "n_judge")
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(rec$facet_table)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(rec$design_table)))
  expect_true(is.data.frame(rec$design_descriptor))
  expect_identical(rec$design_descriptor$role, c("person", "rater", "criterion", "assignment"))
  expect_true(is.list(rec$planning_scope))
  expect_identical(rec$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(is.list(rec$planning_constraints))
  expect_true(is.list(rec$planning_schema))
  expect_true(is.data.frame(rec$planning_schema$role_table))

  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      n_person = 14,
      reps = 1,
      maxit = 10,
      seed = 911
    )
  )
  expect_true(all(pred$forecast$Facet %in% c("Person", "Judge", "Task")))
  expect_true(is.list(pred$planning_scope))
  expect_identical(pred$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(is.list(pred$planning_constraints))
  expect_true(is.list(pred$planning_schema))
  expect_true(is.data.frame(pred$planning_schema$facet_manifest))
  expect_identical(pred$planning_schema$future_planner_contract, "arbitrary_facet_planning_scaffold")
  expect_true(is.data.frame(pred$planning_schema$future_facet_table))
  expect_identical(pred$planning_schema$future_facet_table$future_facet_key, c("person", "judge", "task"))
  expect_true(is.list(pred$planning_schema$future_branch_schema))
  expect_identical(pred$planning_schema$future_branch_schema$planner_stage, "schema_only")

  sig_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      n_person = c(12, 14),
      n_rater = c(spec$n_rater, spec$n_rater + 1L),
      n_criterion = c(spec$n_criterion, spec$n_criterion + 1L),
      raters_per_person = c(spec$raters_per_person, spec$raters_per_person + 1L),
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 912
    )
  )
  s_sig <- summary(sig_eval)
  expect_equal(s_sig$design_variable_aliases[["n_rater"]], "n_judge")
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sig_eval$design_grid)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sig_eval$results)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sig_eval$rep_overview)))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(s_sig$detection_summary)))
  expect_true(is.data.frame(sig_eval$design_descriptor))
  expect_true(is.data.frame(s_sig$design_descriptor))
  expect_true(is.list(sig_eval$planning_scope))
  expect_true(is.list(s_sig$planning_scope))
  expect_identical(sig_eval$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(any(grepl("arbitrary facet counts", s_sig$notes, fixed = TRUE)))
  expect_true(is.list(sig_eval$planning_constraints))
  expect_true(is.list(s_sig$planning_constraints))
  expect_true(is.list(sig_eval$planning_schema))
  expect_true(is.list(s_sig$planning_schema))
  expect_true(is.data.frame(sig_eval$planning_schema$facet_manifest))
  expect_identical(sig_eval$planning_schema$facet_manifest$facet, c("Person", "Judge", "Task"))
  expect_true(is.data.frame(sig_eval$planning_schema$future_facet_table))
  expect_identical(sig_eval$planning_schema$future_facet_table$future_facet_key, c("person", "judge", "task"))
  expect_true(is.list(sig_eval$planning_schema$future_branch_schema))
  p_sig <- plot(sig_eval, signal = "dif", metric = "power", x_var = "rater", group_var = "criterion", draw = FALSE)
  expect_equal(p_sig$x_var, "n_rater")
  expect_equal(p_sig$x_label, "n_judge")
  expect_equal(p_sig$group_var, "n_criterion")
  expect_equal(p_sig$group_label, "n_task")
  expect_true(is.data.frame(p_sig$design_descriptor))
  expect_true(is.list(p_sig$planning_scope))
  expect_identical(p_sig$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(is.list(p_sig$planning_constraints))
  expect_true(is.list(p_sig$planning_schema))
})

test_that("role-based design input is accepted in spec override and forecast helpers", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  spec_over <- mfrmr:::simulation_override_spec_design(
    spec,
    design = list(person = 14, n_judge = spec$n_rater, n_task = spec$n_criterion, judge_per_person = spec$raters_per_person)
  )
  expect_equal(spec_over$n_person, 14L)
  expect_equal(spec_over$n_rater, spec$n_rater)
  expect_equal(spec_over$n_criterion, spec$n_criterion)
  expect_equal(spec_over$raters_per_person, spec$raters_per_person)

  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      design = list(person = 14),
      reps = 1,
      maxit = 10,
      seed = 913
    )
  )
  expect_equal(pred$design$n_person[1], 14)
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(pred$design)))

  expect_error(
    predict_mfrm_population(
      sim_spec = spec,
      n_person = 14,
      design = list(person = 15),
      reps = 1,
      maxit = 10,
      seed = 914
    ),
    "Do not supply the same design variable",
    fixed = TRUE
  )
})

test_that("role-based design input is accepted in design and signal evaluation helpers", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      design = list(person = c(12, 14), n_judge = spec$n_rater, n_task = spec$n_criterion, judge_per_person = spec$raters_per_person),
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 915
    )
  )
  expect_equal(sort(unique(sim_eval$design_grid$n_person)), c(12, 14))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sim_eval$design_grid)))

  expect_error(
    evaluate_mfrm_design(
      n_person = c(12, 14),
      design = list(person = 15),
      reps = 1,
      maxit = 10,
      sim_spec = spec,
      seed = 916
    ),
    "Do not supply the same design variable",
    fixed = TRUE
  )

  sig_eval <- suppressWarnings(
    evaluate_mfrm_signal_detection(
      design = list(person = 12, n_judge = spec$n_rater, n_task = spec$n_criterion, judge_per_person = spec$raters_per_person),
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 917
    )
  )
  expect_equal(unique(sig_eval$design_grid$n_person), 12)
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(sig_eval$design_grid)))

  expect_error(
    evaluate_mfrm_signal_detection(
      n_person = 12,
      design = list(person = 14),
      reps = 1,
      maxit = 10,
      bias_max_iter = 1,
      sim_spec = spec,
      seed = 918
    ),
    "Do not supply the same design variable",
    fixed = TRUE
  )
})

test_that("fit-derived PCM specs honor custom step facet names and stored levels", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  names(toy)[names(toy) == "Rater"] <- "Judge"
  names(toy)[names(toy) == "Criterion"] <- "Task"

  fit <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Judge", "Task"), "Score",
      method = "JML",
      model = "PCM",
      step_facet = "Task",
      maxit = 15
    )
  )

  spec <- extract_mfrm_sim_spec(fit)
  expect_equal(spec$model, "PCM")
  expect_equal(spec$step_facet, "Task")
  expect_setequal(unique(spec$threshold_table$StepFacet), unique(toy$Task))

  sim <- simulate_mfrm_data(sim_spec = spec, seed = 909)
  expect_named(sim, c("Study", "Person", "Judge", "Task", "Score"))
  expect_setequal(unique(sim$Task), unique(toy$Task))
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
  metrics <- mfrmr:::design_eval_recovery_metrics(
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
  expect_s3_class(s$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")
  expect_true(is.data.frame(s$future_branch_active_summary$overview))
  expect_true(is.data.frame(s$future_branch_active_summary$readiness_summary))
  expect_true(is.data.frame(s$future_branch_active_summary$recommendation_table))
  printed <- capture.output(print(summary(sim_eval)))
  expect_true(any(grepl("mfrmr Design Evaluation Summary", printed, fixed = TRUE)))
  expect_true(any(grepl("Future arbitrary-facet planning scaffold", printed, fixed = TRUE)))

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
  expect_identical(
    unname(rec$design_variable_aliases[c("n_person", "n_rater", "n_criterion", "raters_per_person")]),
    c("n_person", "n_rater", "n_criterion", "raters_per_person")
  )
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

test_that("evaluate_mfrm_design and predict_mfrm_population accept latent-regression sim_spec in MML", {
  population_covariates <- data.frame(
    TemplatePerson = sprintf("TP%02d", 1:24),
    X = seq(-1.2, 1.2, length.out = 24),
    stringsAsFactors = FALSE
  )
  spec <- build_mfrm_sim_spec(
    n_person = 24,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating",
    population_formula = ~ X,
    population_coefficients = c(`(Intercept)` = 0.0, X = 0.9),
    population_sigma2 = 0.2,
    population_covariates = population_covariates
  )

  sim_eval <- suppressWarnings(
    evaluate_mfrm_design(
      n_person = 24,
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      fit_method = "MML",
      maxit = 10,
      quad_points = 7,
      sim_spec = spec,
      seed = 1808
    )
  )
  pred <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      n_person = 26,
      reps = 1,
      fit_method = "MML",
      maxit = 10,
      quad_points = 7,
      seed = 1811
    )
  )

  expect_s3_class(sim_eval, "mfrm_design_evaluation")
  expect_true(all(sim_eval$rep_overview$RunOK))
  expect_true(isTRUE(sim_eval$settings$sim_spec$population$active))
  expect_s3_class(pred, "mfrm_population_prediction")
  expect_true(isTRUE(pred$sim_spec$population$active))
  expect_true(all(pred$forecast$Facet %in% c("Person", "Rater", "Criterion")))

  expect_error(
    evaluate_mfrm_design(
      n_person = 24,
      n_rater = 4,
      n_criterion = 3,
      raters_per_person = 2,
      reps = 1,
      fit_method = "JML",
      maxit = 10,
      sim_spec = spec,
      seed = 1809
    ),
    "require `fit_method = \"MML\"`",
    fixed = TRUE
  )
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
  expect_s3_class(s_sig$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")
  expect_true(is.data.frame(s_sig$future_branch_active_summary$overview))
  expect_true(is.data.frame(s_sig$future_branch_active_summary$recommendation_table))
  printed_sig <- capture.output(print(summary(sig_eval)))
  expect_true(any(grepl("mfrmr Signal Detection Summary", printed_sig, fixed = TRUE)))
  expect_true(any(grepl("Future arbitrary-facet planning scaffold", printed_sig, fixed = TRUE)))

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
  p_sig_bias_screen <- plot(sig_eval, signal = "bias", metric = "screen_rate", x_var = "n_person", draw = FALSE)
  expect_equal(p_sig_bias_screen$metric_col, "BiasScreenRate")
  expect_equal(p_sig_bias_screen$display_metric, "Bias screening hit rate")

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
    assignment = "rotating",
    facet_names = c("Judge", "Task")
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
  expect_identical(pred$facet_names, c("Judge", "Task"))
  expect_equal(pred$design_variable_aliases[["n_rater"]], "n_judge")
  expect_equal(pred$design_variable_aliases[["n_criterion"]], "n_task")
  expect_equal(pred$design_variable_aliases[["raters_per_person"]], "judge_per_person")
  expect_true(is.data.frame(pred$design_descriptor))
  expect_true(is.list(pred$planning_scope))
  expect_identical(pred$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(is.list(pred$planning_constraints))
  expect_true(is.list(pred$planning_schema))
  expect_true(all(c("n_judge", "n_task", "judge_per_person") %in% names(pred$design)))

  s_pred <- summary(pred)
  expect_s3_class(s_pred, "summary.mfrm_population_prediction")
  expect_true(is.data.frame(s_pred$forecast))
  expect_identical(s_pred$facet_names, c("Judge", "Task"))
  expect_equal(s_pred$design_variable_aliases[["n_rater"]], "n_judge")
  expect_true(is.data.frame(s_pred$design_descriptor))
  expect_true(is.list(s_pred$planning_scope))
  expect_identical(s_pred$planning_scope$role_labels, c("Judge", "Task"))
  expect_true(any(grepl("arbitrary facet counts", s_pred$notes, fixed = TRUE)))
  expect_true(is.list(s_pred$planning_constraints))
  expect_true(is.list(s_pred$planning_schema))
  expect_s3_class(s_pred$future_branch_active_summary, "summary.mfrm_future_branch_active_branch")
  expect_true(is.data.frame(s_pred$future_branch_active_summary$overview))
  expect_true(is.data.frame(s_pred$future_branch_active_summary$recommendation_table))
  expect_true(any(grepl("All current design variables remain mutable", s_pred$notes, fixed = TRUE)))
  expect_true(all(c("Facet", "MeanSeparation", "McseSeparation") %in% names(s_pred$forecast)))
  printed_pred <- capture.output(print(summary(pred)))
  expect_true(any(grepl("mfrmr Population Prediction Summary", printed_pred, fixed = TRUE)))
  expect_true(any(grepl("Future arbitrary-facet planning scaffold", printed_pred, fixed = TRUE)))
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

test_that("extract_mfrm_sim_spec preserves latent-regression population metadata", {
  population_covariates <- data.frame(
    Person = sprintf("P%02d", 1:24),
    X = seq(-1, 1, length.out = 24),
    Group = factor(rep(c("A", "B"), length.out = 24), levels = c("A", "B")),
    stringsAsFactors = FALSE
  )
  sim_spec <- build_mfrm_sim_spec(
    n_person = 24,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    assignment = "rotating",
    population_formula = ~ X + Group,
    population_coefficients = c(`(Intercept)` = 0.0, X = 0.8, GroupB = 0.3),
    population_sigma2 = 0.2,
    population_covariates = transform(population_covariates, TemplatePerson = Person)[, c("TemplatePerson", "X", "Group")]
  )
  sim <- simulate_mfrm_data(sim_spec = sim_spec, seed = 1911)
  fit <- suppressWarnings(
    fit_mfrm(
      sim,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      maxit = 20,
      quad_points = 7,
      population_formula = ~ X + Group,
      person_data = attr(sim, "mfrm_population_data")$person_data,
      person_id = "Person"
    )
  )

  spec_from_fit <- extract_mfrm_sim_spec(fit)

  expect_true(isTRUE(spec_from_fit$population$active))
  expect_identical(names(spec_from_fit$population$coefficients), c("(Intercept)", "X", "GroupB"))
  expect_true(all(c("TemplatePerson", "X", "Group") %in% names(spec_from_fit$population$covariate_template)))
  expect_identical(spec_from_fit$population$xlevels$Group, c("A", "B"))
  expect_true("Group" %in% names(spec_from_fit$population$contrasts))
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

test_that("future arbitrary-facet active branch exposes summary and plot contracts", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  expect_s3_class(active, "mfrm_future_branch_active_branch")

  s_active <- summary(active)
  expect_s3_class(s_active, "summary.mfrm_future_branch_active_branch")
  expect_true(is.data.frame(s_active$overview))
  expect_true(is.data.frame(s_active$table_index))
  expect_true(is.data.frame(s_active$profile_summary))
  expect_true(is.data.frame(s_active$role_summary))
  expect_true(is.data.frame(s_active$table_profile))
  expect_true(is.data.frame(s_active$load_balance_summary))
  expect_true(is.data.frame(s_active$coverage_summary))
  expect_true(is.data.frame(s_active$guardrail_summary))
  expect_true(is.data.frame(s_active$readiness_summary))
  expect_true(is.data.frame(s_active$recommendation_table))
  expect_true(is.data.frame(s_active$plot_index))
  expect_true(is.data.frame(s_active$table_catalog))
  expect_true(is.data.frame(s_active$appendix_presets))
  expect_true(is.data.frame(s_active$appendix_role_summary))
  expect_true(is.data.frame(s_active$appendix_section_summary))
  expect_true(is.data.frame(s_active$selection_summary))
  expect_true(is.data.frame(s_active$selection_table_summary))
  expect_true(is.data.frame(s_active$selection_table_preset_summary))
  expect_true(is.data.frame(s_active$selection_handoff_table_summary))
  expect_true(is.data.frame(s_active$selection_handoff_preset_summary))
  expect_true(is.data.frame(s_active$selection_handoff_summary))
  expect_true(is.data.frame(s_active$selection_handoff_bundle_summary))
  expect_true(is.data.frame(s_active$selection_handoff_role_summary))
  expect_true(is.data.frame(s_active$selection_handoff_role_section_summary))
  expect_true(is.data.frame(s_active$selection_role_summary))
  expect_true(is.data.frame(s_active$selection_section_summary))
  expect_true(is.data.frame(s_active$reporting_map))
  expect_true(all(c("Table", "Rows", "Role", "Description") %in% names(s_active$table_index)))
  expect_true(all(c("Table", "PlotReady", "NumericColumns", "DefaultPlotTypes") %in% names(s_active$plot_index)))
  expect_true(all(c("Role", "Tables", "TotalRows", "TotalCols") %in% names(s_active$role_summary)))
  expect_true(all(c("Table", "Rows", "Cols", "NumericColumns", "MissingValues") %in% names(s_active$table_profile)))
  expect_true("recommended_design_id" %in% names(s_active$overview))
  expect_true(all(c("RecommendedAppendixTables", "CompactAppendixTables", "NumericTables", "AnyNumericTable") %in% names(s_active$overview)))
  expect_true(all(c(
    "future_branch_overview",
    "future_branch_profile",
    "future_branch_selection_table_presets",
    "future_branch_selection_handoff_tables",
    "future_branch_selection_handoff_presets",
    "future_branch_selection_handoff",
    "future_branch_selection_handoff_bundles",
    "future_branch_selection_handoff_roles",
    "future_branch_selection_handoff_role_sections",
    "future_branch_selection_tables",
    "future_branch_selection_summary",
    "future_branch_reporting_map"
  ) %in% s_active$table_index$Table))
  expect_true(any(grepl("selection_tables", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_handoff_presets", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_handoff", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_handoff_bundles", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_handoff_roles", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_handoff_role_sections", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("selection_roles", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(any(grepl("profile_metrics", s_active$plot_index$DefaultPlotTypes, fixed = TRUE)))
  expect_true(all(c("Preset", "Bundle", "TablesAvailable", "TablesSelected", "SelectionFraction", "PlotReadyFraction", "NumericFraction") %in% names(s_active$selection_summary)))
  expect_true(all(c("Table", "PresetsSelected", "Rows") %in% names(s_active$selection_table_summary)))
  expect_true(all(c("Preset", "Table", "Rows") %in% names(s_active$selection_table_preset_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Role", "Bundle", "Table", "Rows", "NumericColumns", "PlotReady", "ExportReady", "ApaTableReady") %in%
                    names(s_active$selection_handoff_table_summary)))
  expect_true(all(c("recommended", "compact") %in% s_active$selection_table_preset_summary$Preset))
  expect_true(all(c("Preset", "SectionsCovered", "PlotReadyTables", "PlotReadyFraction", "NumericFraction", "KeySections") %in% names(s_active$selection_handoff_preset_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Tables", "PlotReadyFraction", "NumericFraction", "KeyTables") %in% names(s_active$selection_handoff_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Bundle", "PlotReadyTables", "PlotReadyFraction", "NumericFraction", "RolesCovered", "KeyTables") %in% names(s_active$selection_handoff_bundle_summary)))
  expect_true(all(c("Preset", "Role", "PlotReadyTables", "PlotReadyFraction", "NumericFraction", "KeyTables") %in% names(s_active$selection_handoff_role_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Role", "PlotReadyTables", "PlotReadyFraction", "NumericFraction", "KeyTables") %in% names(s_active$selection_handoff_role_section_summary)))
  expect_true(all(c("Preset", "Role", "Tables", "PlotReadyFraction", "NumericFraction") %in% names(s_active$selection_role_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Tables", "PlotReadyFraction", "NumericFraction") %in% names(s_active$selection_section_summary)))
  expect_true(all(is.na(s_active$selection_summary$SelectionFraction) | (s_active$selection_summary$SelectionFraction >= 0 & s_active$selection_summary$SelectionFraction <= 1)))
  expect_true(all(is.na(s_active$selection_handoff_summary$PlotReadyFraction) | (s_active$selection_handoff_summary$PlotReadyFraction >= 0 & s_active$selection_handoff_summary$PlotReadyFraction <= 1)))
  printed_active <- capture.output(print(s_active))
  expect_true(any(grepl("mfrmr Future Arbitrary-Facet Planning Summary", printed_active, fixed = TRUE)))

  p_profile <- plot(active, type = "profile_metrics", draw = FALSE)
  expect_s3_class(p_profile, "mfrm_plot_data")
  expect_identical(p_profile$name, "future_branch_active_branch")
  expect_identical(p_profile$data$plot, "profile_metrics")

  p_tiers <- plot(active, type = "readiness_tiers", draw = FALSE)
  expect_s3_class(p_tiers, "mfrm_plot_data")
  expect_identical(p_tiers$data$plot, "readiness_tiers")

  p_sections <- plot(active, type = "appendix_sections", draw = FALSE)
  expect_s3_class(p_sections, "mfrm_plot_data")
  expect_identical(p_sections$name, "summary_table_bundle")
  expect_identical(p_sections$data$plot, "appendix_sections")

  p_roles <- plot(active, type = "appendix_roles", draw = FALSE)
  expect_s3_class(p_roles, "mfrm_plot_data")
  expect_identical(p_roles$name, "summary_table_bundle")
  expect_identical(p_roles$data$plot, "appendix_roles")

  p_presets <- plot(active, type = "appendix_presets", draw = FALSE)
  expect_s3_class(p_presets, "mfrm_plot_data")
  expect_identical(p_presets$name, "summary_table_bundle")
  expect_identical(p_presets$data$plot, "appendix_presets")

  p_sel_tables <- plot(active, type = "selection_tables", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_tables, "mfrm_plot_data")
  expect_identical(p_sel_tables$name, "future_branch_active_branch")
  expect_identical(p_sel_tables$data$plot, "selection_tables")
  expect_identical(p_sel_tables$data$appendix_preset, "recommended")

  p_sel_handoff_presets <- plot(active, type = "selection_handoff_presets", appendix_preset = "all", draw = FALSE)
  expect_s3_class(p_sel_handoff_presets, "mfrm_plot_data")
  expect_identical(p_sel_handoff_presets$name, "future_branch_active_branch")
  expect_identical(p_sel_handoff_presets$data$plot, "selection_handoff_presets")
  expect_identical(p_sel_handoff_presets$data$appendix_preset, "all")

  p_sel_handoff <- plot(active, type = "selection_handoff", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_handoff, "mfrm_plot_data")
  expect_identical(p_sel_handoff$name, "future_branch_active_branch")
  expect_identical(p_sel_handoff$data$plot, "selection_handoff")
  expect_identical(p_sel_handoff$data$appendix_preset, "recommended")

  p_sel_handoff_fraction <- plot(active, type = "selection_handoff", appendix_preset = "recommended", selection_value = "fraction", draw = FALSE)
  expect_s3_class(p_sel_handoff_fraction, "mfrm_plot_data")
  expect_identical(p_sel_handoff_fraction$data$plot, "selection_handoff")
  expect_identical(p_sel_handoff_fraction$data$selection_value, "fraction")
  expect_true(all(is.finite(as.numeric(p_sel_handoff_fraction$data$table$PlotReadyFraction))))

  p_sel_handoff_bundles <- plot(active, type = "selection_handoff_bundles", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_handoff_bundles, "mfrm_plot_data")
  expect_identical(p_sel_handoff_bundles$name, "future_branch_active_branch")
  expect_identical(p_sel_handoff_bundles$data$plot, "selection_handoff_bundles")
  expect_identical(p_sel_handoff_bundles$data$appendix_preset, "recommended")

  p_sel_handoff_roles <- plot(active, type = "selection_handoff_roles", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_handoff_roles, "mfrm_plot_data")
  expect_identical(p_sel_handoff_roles$name, "future_branch_active_branch")
  expect_identical(p_sel_handoff_roles$data$plot, "selection_handoff_roles")
  expect_identical(p_sel_handoff_roles$data$appendix_preset, "recommended")

  p_sel_handoff_role_sections <- plot(active, type = "selection_handoff_role_sections", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_handoff_role_sections, "mfrm_plot_data")
  expect_identical(p_sel_handoff_role_sections$name, "future_branch_active_branch")
  expect_identical(p_sel_handoff_role_sections$data$plot, "selection_handoff_role_sections")
  expect_identical(p_sel_handoff_role_sections$data$appendix_preset, "recommended")

  p_sel_roles <- plot(active, type = "selection_roles", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(p_sel_roles, "mfrm_plot_data")
  expect_identical(p_sel_roles$name, "future_branch_active_branch")
  expect_identical(p_sel_roles$data$plot, "selection_roles")
  expect_identical(p_sel_roles$data$appendix_preset, "recommended")

  p_sel_sections <- plot(active, type = "selection_sections", appendix_preset = "compact", draw = FALSE)
  expect_s3_class(p_sel_sections, "mfrm_plot_data")
  expect_identical(p_sel_sections$name, "future_branch_active_branch")
  expect_identical(p_sel_sections$data$plot, "selection_sections")
  expect_identical(p_sel_sections$data$appendix_preset, "compact")

  expect_error(
    plot(active, type = "selection_tables", appendix_preset = "recommended", selection_value = "fraction", draw = FALSE),
    "not available for `type = \"selection_tables\"`",
    fixed = TRUE
  )
})
