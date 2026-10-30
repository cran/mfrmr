test_that("population simulation retains the template transformation basis", {
  template <- data.frame(TemplatePerson = paste0("T", 1:12),
    X = seq(-2, 3, length.out = 12), G = rep(c("a", "b"), 6))
  for (formula in list(~ scale(X), ~ poly(X, 2), ~ splines::ns(X, df = 3),
                       ~ scale(X) * G)) {
    mf <- model.frame(formula, template)
    tt <- attr(mf, "terms")
    mm <- model.matrix(tt, mf)
    beta <- setNames(seq_len(ncol(mm)) / 5, colnames(mm))
    spec <- build_mfrm_sim_spec(n_person = 5, n_rater = 2, n_criterion = 2,
      population_formula = formula, population_coefficients = beta,
      population_sigma2 = 0, population_covariates = template)
    sim <- simulate_mfrm_data(sim_spec = spec, seed = 1121)
    pop <- attr(sim, "mfrm_population_data")
    expected <- model.matrix(tt, model.frame(tt, pop$person_data,
      xlev = spec$population$xlevels), contrasts.arg = spec$population$contrasts)
    mu <- drop(expected %*% beta)
    expect_equal(unname(attr(sim, "mfrm_truth")$person), unname(mu), tolerance = 1e-13)
    refit <- prepare_mfrm_population_scaffold(sim, "Person", pop$population_formula,
      pop$person_data, "Person", require_full_rank = FALSE)
    expect_equal(unname(refit$design_matrix), unname(expected), tolerance = 1e-13,
      ignore_attr = TRUE)
    replay <- unserialize(serialize(spec, NULL))
    replay_sim <- simulate_mfrm_data(sim_spec = replay, seed = 1121)
    expect_identical(lapply(replay_sim, identity), lapply(sim, identity))
    expect_identical(attr(replay_sim, "mfrm_truth")$person, attr(sim, "mfrm_truth")$person)
  }
})

test_that("population simulation and refitting retain categorical contrasts", {
  withr::local_options(contrasts = c("contr.sum", "contr.poly"))
  template <- data.frame(TemplatePerson = paste0("T", 1:12),
    G = rep(c("a", "b", "c"), 4))
  spec <- build_mfrm_sim_spec(n_person = 20, n_rater = 2, n_criterion = 2,
    population_formula = ~ G, population_coefficients = c(0, .5, 1),
    population_sigma2 = 0, population_covariates = template)
  options(contrasts = c("contr.treatment", "contr.poly"))
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 34)
  pop <- attr(sim, "mfrm_population_data")
  refit <- prepare_mfrm_population_scaffold(sim, "Person", pop$population_formula,
    pop$person_data, "Person")
  expected <- model.matrix(~ G, pop$person_data, contrasts.arg = list(G = "contr.sum"))
  expect_identical(colnames(refit$design_matrix), colnames(expected))
  expect_equal(unname(refit$design_matrix), unname(expected), ignore_attr = TRUE)
  expect_equal(unname(attr(sim, "mfrm_truth")$person),
    unname(drop(refit$design_matrix %*% spec$population$coefficients)))
})

test_that("factor calls retain their template levels and contrasts during refitting", {
  withr::local_options(contrasts = c("contr.sum", "contr.poly"))
  template <- data.frame(TemplatePerson = paste0("T", 1:3), G = c("a", "b", "c"))
  spec <- build_mfrm_sim_spec(n_person = 2, n_rater = 2, n_criterion = 2,
    population_formula = ~ factor(G), population_coefficients = c(0, .5, 1),
    population_sigma2 = 0, population_covariates = template)
  options(contrasts = c("contr.treatment", "contr.poly"))
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 34)
  pop <- attr(sim, "mfrm_population_data")
  refit <- prepare_mfrm_population_scaffold(sim, "Person", pop$population_formula,
    pop$person_data, "Person", require_full_rank = FALSE)
  expected <- c(a = .5, b = 1, c = -1.5)[pop$person_data$G]
  expect_identical(refit$design_columns, spec$population$design_columns)
  expect_identical(refit$xlevels, spec$population$xlevels)
  expect_equal(unname(attr(sim, "mfrm_truth")$person), unname(expected))
  expect_equal(unname(drop(refit$design_matrix %*% spec$population$coefficients)), unname(expected))
})

simulation_template_fit <- function(owner = "Rater", custom = FALSE) {
  facets <- if (custom) c("Judge", "Task") else c("Rater", "Criterion")
  levels <- list(Person = paste0("P", 1:12), r = paste0("J", 1:3), c = paste0("T", 1:2))
  names(levels) <- c("Person", facets)
  data <- expand.grid(levels, stringsAsFactors = FALSE)
  data$Score <- rep(1:3, length.out = nrow(data))
  owner <- if (owner == "Rater") facets[1] else facets[2]
  structure(list(prep = list(data = data, levels = levels),
    summary = data.frame(Model = "GPCM"),
    config = list(model = "GPCM", facet_names = facets, step_facet = owner, slope_facet = owner),
    facets = list(person = data.frame(Estimate = seq(-1, 1, length.out = 12)),
      others = data.frame(Facet = rep(facets, c(3, 2)), Estimate = c(-.5, 0, .5, -.2, .2))),
    steps = data.frame(StepFacet = rep(levels[[owner]], each = 2),
      Step = rep(1:2, length(levels[[owner]])), Estimate = rep(c(-1, 1), length(levels[[owner]]))),
    slopes = data.frame(SlopeFacet = levels[[owner]],
      Estimate = exp(seq(-.4, .4, length.out = length(levels[[owner]]))))), class = "mfrm_fit")
}

test_that("fit-derived slopes follow their owner and refresh planning metadata", {
  for (owner in c("Rater", "Criterion")) for (custom in c(FALSE, TRUE)) {
    fit <- simulation_template_fit(owner, custom)
    spec <- extract_mfrm_sim_spec(fit)
    expect_equal(spec$slope_table, tibble::as_tibble(fit$slopes))
    expect_identical(spec$planning_scope$role_labels, fit$config$facet_names)
    expect_equal(spec$planning_schema, simulation_planning_schema(spec))
    sim <- simulate_mfrm_data(sim_spec = spec, seed = 67)
    expect_setequal(attr(sim, "mfrm_truth")$slope_table$SlopeFacet,
      fit$prep$levels[[fit$config$slope_facet]])
  }
})

test_that("fit-derived population basis survives complete-case selection", {
  fit <- simulation_template_fit()
  fit$summary$Model <- fit$config$model <- "RSM"
  fit$config$slope_facet <- NULL
  fit$steps <- data.frame(Step = 1:2, Estimate = c(-1, 1))
  people <- data.frame(Person = paste0("P", 1:12), X = 1:12, Z = c(NA, rep(0:1, 5), 0))
  pop <- prepare_mfrm_population_scaffold(people, "Person", ~ scale(X) + Z,
    people, "Person", population_policy = "omit")
  pop$coefficients <- setNames(c(.2, .5, .8), pop$design_columns)
  pop$sigma2 <- 0
  fit$population <- pop
  fit$prep$data <- fit$prep$data[fit$prep$data$Person != "P1", ]
  fit$prep$levels$Person <- pop$included_persons
  spec <- extract_mfrm_sim_spec(fit)
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 80)
  p <- attr(sim, "mfrm_population_data")$person_data
  expected <- .2 + .5 * ((p$X - mean(people$X)) / sd(people$X)) + .8 * p$Z
  expect_equal(unname(attr(sim, "mfrm_truth")$person), expected)
  fit$population$terms <- NULL
  expect_equal(extract_mfrm_sim_spec(fit)$population$terms, spec$population$terms)
})

test_that("legacy transformed templates fail explicitly while plain ones still work", {
  people <- data.frame(TemplatePerson = paste0("P", 1:6), X = 1:6)
  for (formula in list(~ scale(X), ~ X)) {
    spec <- build_mfrm_sim_spec(population_formula = formula,
      population_coefficients = c(0, 1), population_sigma2 = .4,
      population_covariates = people)
    spec$population$terms <- NULL
    if (identical(as.character(formula)[2], "X")) {
      expect_s3_class(simulate_mfrm_data(sim_spec = spec, seed = 3), "data.frame")
    } else {
      expect_error(simulate_mfrm_data(sim_spec = spec, seed = 3), "Rebuild.*specification")
    }
  }
})
