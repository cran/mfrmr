assignment_template_spec <- function(assignment = "skeleton", reuse = FALSE,
                                     population_assignment = "matched", group = TRUE) {
  people <- data.frame(TemplatePerson = paste0("T", 1:4), X = 1:4,
    Group = rep(c("A", "B"), each = 2))
  roster <- merge(people[c("TemplatePerson", "Group")],
    data.frame(Criterion = c("C01", "C02")))
  roster$Rater <- ifelse(roster$Group == "A", "R01", "R02")
  if (reuse) roster$TemplatePersonReuse <- TRUE
  if (!group) people$Group <- NULL
  build_mfrm_sim_spec(n_person = if (reuse) 4 else 16, n_rater = 2, n_criterion = 2,
    raters_per_person = 1, assignment = assignment,
    assignment_profiles = unique(roster[c("TemplatePerson", "Rater", "Group")]),
    design_skeleton = roster, population_formula = ~ X,
    population_coefficients = c(0, 1), population_sigma2 = 0,
    population_covariates = people, population_assignment = population_assignment)
}

test_that("population covariates and response Group refer to the same Persons", {
  spec <- build_mfrm_sim_spec(n_person = 20, n_rater = 2, n_criterion = 2,
    group_levels = c("A", "B"), population_formula = ~ Group,
    population_coefficients = c(0, 2), population_sigma2 = 0,
    population_covariates = data.frame(TemplatePerson = c("a", "b"), Group = c("A", "B")))
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 71)
  people <- attr(sim, "mfrm_population_data")$person_data
  expect_equal(sim$Group, as.character(people$Group[match(sim$Person, people$Person)]))
  expect_equal(unname(attr(sim, "mfrm_truth")$person), 2 * (people$Group == "B"))
})

test_that("matched templates preserve covariate, group and assignment associations", {
  for (assignment in c("skeleton", "resampled")) {
    spec <- assignment_template_spec(assignment)
    sim <- simulate_mfrm_data(sim_spec = spec, seed = 19)
    pop <- attr(sim, "mfrm_population_data")$person_data
    expect_equal(sim$Group, as.character(pop$Group[match(sim$Person, pop$Person)]))
    expect_equal(sim$Rater, ifelse(sim$Group == "A", "R01", "R02"))
    expected_x <- pop$X[match(sim$Person, pop$Person)]
    expect_equal(sim$Group, ifelse(expected_x <= 2, "A", "B"))
    truth <- attr(sim, "mfrm_truth")
    expect_identical(truth$population$template_person, truth$design$template_person)
    expect_equal(unname(truth$person), pop$X)
    expect_identical(attr(sim, "mfrm_simulation_spec")$population_assignment, "matched")
  }
})

test_that("fixed skeleton order preserves repeated rows and matched covariates", {
  spec <- assignment_template_spec(reuse = TRUE)
  spec$design_skeleton <- rbind(spec$design_skeleton, spec$design_skeleton[1:2, ])
  spec$population$covariate_template <- spec$population$covariate_template[4:1, ]
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 84)
  expect_equal(nrow(sim), nrow(spec$design_skeleton))
  truth <- attr(sim, "mfrm_truth")
  expect_identical(unname(truth$design$template_person), unique(spec$design_skeleton$TemplatePerson))
  expect_equal(unname(truth$person), as.numeric(sub("T", "", truth$design$template_person)))
  expect_equal(table(sim$Person), setNames(table(spec$design_skeleton$TemplatePerson), unique(sim$Person)),
    ignore_attr = TRUE)
})

test_that("independent sampling is explicit and conflicting Group definitions fail", {
  spec <- assignment_template_spec(population_assignment = "independent", group = FALSE)
  spec$population$covariate_template$TemplatePerson <- paste0("U", 1:4)
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 19)
  truth <- attr(sim, "mfrm_truth")
  expect_true(all(grepl("U", truth$population$template_person)))
  expect_true(all(grepl("T", truth$design$template_person)))
  expect_identical(attr(sim, "mfrm_simulation_spec")$population_assignment, "independent")
  spec$population_assignment <- "matched"
  expect_error(simulate_mfrm_data(sim_spec = spec, seed = 19), "same.*TemplatePerson")
  expect_error(assignment_template_spec(population_assignment = "independent"), "[Ii]ndependent.*Group")
  spec <- assignment_template_spec()
  spec$population$covariate_template$Group[1] <- "B"
  expect_error(simulate_mfrm_data(sim_spec = spec, seed = 19), "Group.*agree")
})

test_that("saved independent specs keep their sampling law and invalid groups fail early", {
  spec <- assignment_template_spec(population_assignment = "independent", group = FALSE)
  expected <- simulate_mfrm_data(sim_spec = spec, seed = 23)
  spec$population_assignment <- NULL
  legacy <- simulate_mfrm_data(sim_spec = spec, seed = 23)
  expect_identical(lapply(legacy, identity), lapply(expected, identity))
  expect_identical(attr(legacy, "mfrm_truth"), attr(expected, "mfrm_truth"))
  spec$group_levels <- "A"
  expect_error(simulate_mfrm_data(sim_spec = spec), "Template Group.*group_levels")
  spec <- assignment_template_spec()
  spec$design_skeleton$Group[1] <- NA
  expect_error(simulate_mfrm_data(sim_spec = spec), "non-missing")
  spec <- assignment_template_spec(reuse = TRUE)
  spec$design_skeleton$TemplatePersonReuse[1] <- FALSE
  expect_error(simulate_mfrm_data(sim_spec = spec), "uniformly")
})

test_that("skeleton normalization and fit extraction retain observation multiplicity", {
  skeleton <- data.frame(TemplatePerson = rep(c("a", "b"), each = 3),
    Rater = rep(c("R01", "R02"), each = 3), Criterion = rep(c("C01", "C01", "C02"), 2),
    Weight = c(1, 1, 2, 1, 1, 2), TemplatePersonReuse = TRUE)
  spec <- build_mfrm_sim_spec(n_person = 2, n_rater = 2, n_criterion = 2,
    raters_per_person = 1, assignment = "skeleton", design_skeleton = skeleton)
  expect_equal(nrow(spec$design_skeleton), 6)
  sim <- simulate_mfrm_data(sim_spec = spec, seed = 6)
  expect_equal(nrow(sim), 6)
  expect_equal(sim$Weight, skeleton$Weight)
  # Extraction contract only: supplied facet/step tables are controlled values.
  fit <- structure(list(summary = data.frame(Model = "RSM"),
    config = list(model = "RSM", facet_names = c("Rater", "Criterion")),
    prep = list(data = sim, levels = list(Person = unique(sim$Person), Rater = c("R01", "R02"),
      Criterion = c("C01", "C02"))), facets = list(person = data.frame(Estimate = c(-1, 1)),
      others = data.frame(Facet = rep(c("Rater", "Criterion"), each = 2), Estimate = rep(c(-.2, .2), 2))),
    steps = data.frame(Step = 1:3, Estimate = c(-1, 0, 1))), class = "mfrm_fit")
  sim$Score <- fit$prep$data$Score <- rep(1:3, 2)
  fit$steps <- data.frame(Step = 1:2, Estimate = c(-1, 1))
  extracted <- extract_mfrm_sim_spec(fit, assignment = "skeleton")
  expect_equal(nrow(extracted$design_skeleton), 6)
})

test_that("matched groups and repeated rows follow the declared adjacent logits", {
  captured <- NULL
  implementations <- lapply(c("category_prob_rsm", "category_prob_pcm", "category_prob_gpcm"),
    function(n) get(n, envir = asNamespace("mfrmr")))
  local_mocked_bindings(
    category_prob_rsm = function(...) { captured <<- implementations[[1]](...); captured },
    category_prob_pcm = function(...) { captured <<- implementations[[2]](...); captured },
    category_prob_gpcm = function(...) { captured <<- implementations[[3]](...); captured },
    .package = "mfrmr")
  for (model in c("RSM", "PCM", "GPCM")) {
    spec <- assignment_template_spec(reuse = TRUE)
    spec$model <- model
    spec$design_skeleton <- rbind(spec$design_skeleton, spec$design_skeleton[1:2, ])
    spec$dif_effects <- data.frame(Group = "B", Criterion = "C02", Effect = .7)
    spec$interaction_effects <- data.frame(Rater = "R01", Criterion = "C01", Effect = -.3)
    if (model == "GPCM") {
      spec$slope_facet <- "Criterion"
      spec$slope_table <- data.frame(SlopeFacet = c("C01", "C02"), Estimate = c(.5, 2))
    }
    sim <- simulate_mfrm_data(sim_spec = spec, seed = 88)
    truth <- attr(sim, "mfrm_truth")
    expected <- t(vapply(seq_len(nrow(sim)), function(i) {
      row <- sim[i, ]
      theta <- as.numeric(sub("T", "", truth$design$template_person[row$Person]))
      eta <- theta - truth$facets$Rater[row$Rater] - truth$facets$Criterion[row$Criterion] +
        .7 * (row$Group == "B" && row$Criterion == "C02") -
        .3 * (row$Rater == "R01" && row$Criterion == "C01")
      slope <- if (model == "GPCM") c(C01 = .5, C02 = 2)[row$Criterion] else 1
      weights <- c(1, cumprod(exp(slope * (eta - spec$threshold_table$Estimate))))
      weights / sum(weights)
    }, numeric(spec$score_levels)))
    expect_equal(unname(captured), unname(expected), tolerance = 1e-13)
    expect_equal(nrow(sim), nrow(spec$design_skeleton))
    fit_args <- simulation_add_fit_score_support(list(), sim)
    expect_identical(fit_args$repeated_ratings, "independent")
    expect_equal(c(fit_args$rating_min, fit_args$rating_max), c(1, spec$score_levels))
  }
})
