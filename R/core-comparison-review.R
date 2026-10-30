# Prospective interaction comparisons, separate from residual/refit DFF effects.
mfrmr_comparison_rows <- function(data, facets) {
  out <- as.data.frame(data[c("Person", facets, "Score", "Weight")])
  for (f in c("Person", facets)) out[[f]] <- as.character(out[[f]])
  out$Score <- as.numeric(out$Score); out$Weight <- as.numeric(out$Weight)
  out <- out[do.call(order, c(unname(out), list(method = "radix"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Interpret existing map checks without promoting conditional rank to marginal identification.
mfrmr_comparison_basis <- function(table) {
  identified <- table$ConditionalState %in% c("all_parameters_identified", "target_identified_only")
  free <- table$FreePersonState %in% c("all_parameters_identified", "target_identified_only")
  basis <- rep("not_evaluated", nrow(table))
  basis[table$State == "unsupported"] <- "unsupported"
  basis[table$ConditionalState == "target_unidentified"] <- "unidentified"
  basis[identified & table$Method == "CML"] <- "conditional_support"
  basis[identified & (table$Method == "JML" | (table$Method == "MML" & free))] <- "free_person_response_model"
  basis[identified & table$Method == "MML" & table$FreePersonState == "target_unidentified"] <- "population_restriction"
  detail <- c(
    not_evaluated = "The required map checks are unavailable or numerically sensitive; no basis is assigned.",
    unsupported = "The requested model is unsupported; no substitute model was reviewed.",
    unidentified = "The reviewed conditional law permits a change in this target without changing response probabilities.",
    conditional_support = "The target is identified by conditional response patterns at the observed Person totals and rosters, without an ability distribution.",
    free_person_response_model = "The declared response map identifies this target with unrestricted Person locations; this is conditional on its response-model restrictions.",
    population_restriction = "The target has an alias with unrestricted Persons that is removed by the declared population restrictions in the conditional map.")
  table$IdentificationBasis <- basis
  table$IdentificationDetail <- unname(detail[basis])
  mml <- table$Method == "MML" & identified
  table$IdentificationDetail[mml] <- paste(table$IdentificationDetail[mml],
    "MML estimation still assumes its declared population; global marginal identification is not established.")
  table
}

mfrmr_comparison_rowspace <- function(design, targets) {
  nr <- nrow(targets); p <- ncol(design)
  empty <- list(state = rep("not_evaluated", nr), rank = NA_integer_, dimension = p,
    null_fraction = rep(NA_real_, nr), reason = "Numerical map exceeds the bounded dense row-space review.")
  if (as.double(nrow(design)) * p > 2e6 || as.double(p) * p > 2e6) return(empty)
  X <- as.matrix(design)
  norms <- sqrt(colSums(X^2)); norms[norms == 0] <- 1
  X <- sweep(X, 2L, norms, "/")
  C <- sweep(as.matrix(targets), 2L, norms, "/")
  if (any(!is.finite(C))) stop("Comparison coefficients overflow the model map; rescale the contrast weights.", call. = FALSE)
  # Row rescaling must not turn large or tiny nonzero targets into zero rows.
  row_scale <- apply(abs(C), 1L, max); row_scale[row_scale == 0] <- 1
  C <- C / row_scale
  length <- sqrt(rowSums(C^2)); length[length == 0] <- 1
  C <- C / length
  s <- svd(X, nu = 0L, nv = p)
  tolerances <- c(1e-12, 1e-10, 1e-8)
  ranks <- vapply(tolerances, function(tol) sum(s$d > max(s$d) * tol), integer(1))
  states <- fractions <- matrix(NA_real_, nr, length(tolerances))
  for (i in seq_along(tolerances)) {
    rank <- ranks[i]
    null <- if (rank < p) s$v[, seq.int(rank + 1L, p), drop = FALSE] else matrix(0, p, 0)
    fractions[, i] <- sqrt(rowSums((C %*% null)^2))
    states[, i] <- fractions[, i] < 1e-9
  }
  # A separate augmented QR check protects against a misleading SVD threshold.
  qr_rank <- function(A) qr(A, tol = 1e-10, LAPACK = FALSE)$rank
  base_rank <- qr_rank(X)
  qr_identified <- apply(C, 1L, function(c) qr_rank(rbind(X, c)) == base_rank)
  stable <- length(unique(ranks)) == 1L & base_rank == ranks[2L] &
    apply(states, 1L, function(v) length(unique(v)) == 1L) & qr_identified == states[, 2L]
  state <- ifelse(states[, 2L], if (ranks[2L] == p) "all_parameters_identified" else
    "target_identified_only", "target_unidentified")
  state[!stable] <- "not_evaluated"
  list(state = state, rank = ranks[2L], dimension = p,
    null_fraction = fractions[, 2L], reason = ifelse(stable,
      "Constrained linear response map; complete null space and augmented QR agree.",
      "Rank or row-space membership is numerically sensitive; no conclusion assigned."))
}

mfrmr_comparison_model <- function(prep, s) {
  population <- prepare_mfrm_population_scaffold(prep$data, "Person",
    population_formula = s$population_formula,
    person_data = unique(prep$data[c("Person", s$group)]))
  interactions <- build_facet_interaction_specs(prep,
    paste(s$facet, s$group, sep = ":"), interaction_policy = "silent")
  build_estimation_config(prep, s$model, s$method, step_facet = s$step_facet,
    slope_facet = NULL, interaction_specs = interactions, weight_col = NULL,
    facet_signs = build_facet_signs(prep$facet_names, s$positive_facets)$signs,
    positive_facets = s$positive_facets, noncenter_facet = "Person", dummy_facets = s$group,
    anchor_df = NULL, group_anchor_df = NULL, population = population, centering_scope = "facet")
}

mfrmr_build_comparison_review <- function(prep, specification) {
  if (is.null(specification)) return(NULL)
  if (!isTRUE(prep$keep_original)) stop("Comparison review requires category_policy = 'preserve' (or keep_original = TRUE).", call. = FALSE)
  allowed <- c("facet", "group", "groups", "contrasts", "model", "method", "step_facet",
    "population_formula", "positive_facets")
  if (!is.list(specification) || is.null(names(specification)) ||
      anyDuplicated(names(specification)) || any(!names(specification) %in% allowed)) {
    stop("`comparisons` must be a named specification; see describe_mfrm_data().", call. = FALSE)
  }
  s <- specification
  for (field in c("facet", "group")) {
    value <- s[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% prep$facet_names) stop("Comparison ", field, " must name a reviewed facet.", call. = FALSE)
  }
  if (s$facet == s$group) stop("Comparison facet and group must differ.", call. = FALSE)
  levels <- as.character(prep$levels[[s$facet]])
  groups <- as.character(prep$levels[[s$group]])
  if (!is.character(s$groups) || length(s$groups) != 2L || anyNA(s$groups) ||
      anyDuplicated(s$groups) || !all(s$groups %in% groups)) {
    stop("Comparison `groups` must name two distinct observed origin groups, in subtraction order.", call. = FALSE)
  }
  W <- s$contrasts
  if (!is.matrix(W) || !is.numeric(W) || is.complex(W) || !nrow(W) ||
      any(!is.finite(W)) || is.null(rownames(W)) || anyNA(rownames(W)) || any(!nzchar(rownames(W))) ||
      anyDuplicated(rownames(W)) || is.null(colnames(W)) || anyDuplicated(colnames(W)) ||
      !setequal(colnames(W), levels) || any(rowSums(abs(W)) == 0)) {
    stop("Comparison contrasts need unique row labels, every facet level as a named column, and finite nonzero rows summing to zero.", call. = FALSE)
  }
  scaled_weights <- W / apply(abs(W), 1L, max)
  if (any(abs(rowSums(scaled_weights)) > 1e-12 * rowSums(abs(scaled_weights))))
    stop("Comparison contrasts need finite nonzero rows summing to zero.", call. = FALSE)
  W <- W[, levels, drop = FALSE]
  origin <- unique(prep$data[c("Person", s$group)])
  if (anyDuplicated(origin$Person)) stop("A Person/recording must retain one origin group across ratings.", call. = FALSE)
  if (any(prep$data$Weight != 1)) stop("Comparison review currently requires unit-weight ratings.", call. = FALSE)
  s$model <- s$model %||% "RSM"; s$method <- s$method %||% "JML"
  if (!is.character(s$model) || length(s$model) != 1L || !s$model %in% c("RSM", "PCM", "GPCM") ||
      !is.character(s$method) || length(s$method) != 1L || !s$method %in% c("JML", "MML", "CML")) {
    stop("Comparison model must be RSM, PCM or GPCM and method must be JML, MML or CML.", call. = FALSE)
  }
  if (s$model == "PCM" && (!is.character(s$step_facet) || length(s$step_facet) != 1L ||
      is.na(s$step_facet) || !s$step_facet %in% prep$facet_names)) {
    stop("PCM comparison review requires an explicit reviewed `step_facet`.", call. = FALSE)
  }
  if (s$model == "RSM" && !is.null(s$step_facet)) stop("RSM uses shared steps; omit comparison step_facet.", call. = FALSE)
  if (s$method == "JML" && !is.null(s$population_formula)) stop("JML comparison review has free Person locations, not a population formula.", call. = FALSE)
  if (s$method == "CML" && !is.null(s$population_formula)) stop("CML conditions out Person abilities and estimates no population model.", call. = FALSE)
  if (!is.null(s$population_formula) && (!inherits(s$population_formula, "formula") ||
      any(!all.vars(s$population_formula) %in% s$group))) {
    stop("Comparison population_formula may use only the declared origin group (or an intercept).", call. = FALSE)
  }
  s$positive_facets <- s$positive_facets %||% character()
  if (!is.character(s$positive_facets) || anyNA(s$positive_facets) ||
      any(!s$positive_facets %in% prep$facet_names)) stop("Invalid comparison positive_facets.", call. = FALSE)
  cells <- expand.grid(Level = levels, Group = groups, stringsAsFactors = FALSE)
  cells$N <- cells$NPersons <- integer(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    take <- prep$data[[s$facet]] == cells$Level[i] & prep$data[[s$group]] == cells$Group[i]
    cells$N[i] <- sum(take); cells$NPersons[i] <- length(unique(prep$data$Person[take]))
  }
  fit_arguments <- list(person = prep$source_columns$person %||% "Person",
    facets = prep$facet_names, score = prep$source_columns$score %||% "Score",
    rating_min = prep$rating_min, rating_max = prep$rating_max, category_policy = "preserve",
    model = s$model, method = s$method, step_facet = s$step_facet,
    dummy_facets = s$group, noncenter_facet = "Person", centering_scope = "facet",
    facet_interactions = paste(s$facet, s$group, sep = ":"), positive_facets = s$positive_facets,
    population_formula = s$population_formula)
  if (s$method == "CML") fit_arguments$population_formula <- NULL
  if (!is.null(s$population_formula)) {
    fit_arguments$person_data <- origin
    names(fit_arguments$person_data)[1L] <- fit_arguments$person
  }
  table <- data.frame(Comparison = rownames(W), Group1 = s$groups[1], Group2 = s$groups[2],
    Facet = s$facet, Model = s$model, Method = s$method,
    StepFacet = if (s$model == "RSM") "shared" else s$step_facet %||% NA_character_,
    SlopeFacet = NA_character_,
    State = "not_evaluated", ConditionalState = "not_evaluated", FreePersonState = "not_evaluated",
    Basis = "not_evaluated", EvidenceScope = "not_evaluated", Rank = NA_integer_,
    Dimension = NA_integer_, NullFraction = NA_real_, Estimate = NA_real_,
    FormalInferenceEligible = FALSE, stringsAsFactors = FALSE)
  for (j in 1:2) {
    table[[paste0("N_Group", j)]] <- vapply(seq_len(nrow(W)), function(i) {
      sum(prep$data[[s$facet]] %in% levels[W[i, ] != 0] & prep$data[[s$group]] == s$groups[j])
    }, integer(1))
    table[[paste0("NPersons_Group", j)]] <- vapply(seq_len(nrow(W)), function(i) {
      length(unique(prep$data$Person[prep$data[[s$facet]] %in% levels[W[i, ] != 0] &
        prep$data[[s$group]] == s$groups[j]]))
    }, integer(1))
  }
  notes <- c("Prospective facet-by-origin interaction comparisons; these are not residual or subgroup-refit contrasts.",
    "Each target is sum of the declared facet weights times interaction(Group1) minus interaction(Group2).",
    "All conclusions condition on the declared response model and constraints. Identifiability does not guarantee finite estimates, convergence or precision.",
    "Group is a zero main facet; other facet locations and interaction margins are centered. Anchors and additional interactions are outside this review.",
    "Comparison counts cover nonzero-weight facet levels; cells retains all observed and empty facet-by-origin cells. Repeated ratings are not independent Persons.")
  output <- list(table = table, weights = W, cells = cells, specification = s,
    fit_arguments = fit_arguments, notes = notes,
    source = list(facets = prep$facet_names, score_map = prep$score_map,
      rows = mfrmr_comparison_rows(prep$data, prep$facet_names)),
    contract_version = 1L)
  if (s$model == "GPCM") {
    output$table$State <- "unsupported"
    output$table$Detail <- "GPCM does not support the requested facet-by-group interaction model; no additive-rank substitute was used."
    output$table <- mfrmr_comparison_basis(output$table)
    return(output)
  }
  built <- mfrmr_comparison_model(prep, s)
  cfg <- built$config; sizes <- built$sizes
  idx <- build_indices(prep, step_facet = s$step_facet, interaction_specs = cfg$interaction_specs)
  adjacent <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes)
  # Targets are formed in expanded interaction cells, then mapped to actual
  # optimizer coordinates; do not use truncated diagnostic null directions.
  spec <- cfg$interaction_specs[[1L]]
  weights <- matrix(0, nrow(W), spec$n_a * spec$n_b)
  for (a in seq_len(spec$n_a)) for (b in seq_len(spec$n_b)) {
    lev <- if (spec$facet_a == s$facet) spec$levels_a[a] else spec$levels_b[b]
    grp <- if (spec$facet_a == s$group) spec$levels_a[a] else spec$levels_b[b]
    weights[, a + (b - 1L) * spec$n_a] <- W[, lev] *
      (as.numeric(grp == s$groups[1]) - as.numeric(grp == s$groups[2]))
  }
  block <- mfrmr_interaction_jacobian_sparse(spec)
  target_map <- function(map) {
    C <- matrix(0, nrow(W), nrow(map))
    C[, match(block$map$Coordinate, map$Coordinate)] <- as.matrix(weights %*% block$jacobian)
    C
  }
  if (s$method == "CML") {
    problem <- mfrmr_cml_problem(prep, idx, cfg, sizes)
    support <- mfrmr_cml_support_design(problem)
    adjacent <- list(design = support, map = problem$map)
    if (is.null(support)) {
      output$table$Detail <- "Conditional support exceeds the bounded dense row-space review."
      output$table <- mfrmr_comparison_basis(output$table)
      return(output)
    }
  }
  audit <- mfrmr_comparison_rowspace(adjacent$design, target_map(adjacent$map))
  output$table$ConditionalState <- audit$state
  output$table$State <- audit$state
  output$table$Rank <- audit$rank; output$table$Dimension <- audit$dimension
  output$table$NullFraction <- audit$null_fraction
  output$table$Detail <- audit$reason
  output$table$EvidenceScope <- "conditional_linear_response_map"
  output$table$Basis <- "observed_design_under_declared_response_model"
  output$table$FreePersonState <- audit$state
  if (s$method == "CML") {
    output$table$EvidenceScope <- "conditional_support_at_observed_totals"
    output$table$Basis <- "person_totals_and_observed_rosters"
    output$table$Detail <- sub("Constrained linear response map", "Conditional support differences", audit$reason)
    counter <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes,
      include_person = TRUE, include_population_beta = FALSE)
    output$table$FreePersonState <- mfrmr_comparison_rowspace(counter$design, target_map(counter$map))$state
    output$notes <- c(output$notes,
      "CML identification uses all possible response patterns at the observed Person totals and rosters; extreme totals carry no conditional information.",
      "Person-constant main effects disappear on conditioning. The zero origin main facet is coding, not an assumption of equal group abilities; no population SD is estimated.")
  }
  if (s$method == "MML") {
    counter <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes,
      include_person = TRUE, include_population_beta = FALSE)
    counter_audit <- mfrmr_comparison_rowspace(counter$design, target_map(counter$map))
    output$table$FreePersonState <- counter_audit$state
    # A linear same-law alias survives mixing; a positive conditional result
    # alone does not certify identification after marginalization.
    negative <- audit$state == "target_unidentified"
    output$table$State[!negative] <- "not_evaluated"
    output$table$EvidenceScope <- ifelse(negative, "same_conditional_law_alias_survives_mixing",
      "conditional_map_only_marginal_identification_not_evaluated")
    output$table$Basis <- if (isTRUE(cfg$population_spec$active)) "declared_normal_population_regression" else "fixed_common_standard_normal_population"
    output$notes <- c(output$notes, "Positive conditional MML ranks are not global marginal identification certificates; no fitted-point local or uncertainty calculation is performed here.")
    output$notes <- c(output$notes, if (!isTRUE(cfg$population_spec$active))
      "The model fixes one common N(0,1) ability population across origins; equal origin distributions are a substantive restriction, not scale coding." else if (
        identical(cfg$population_spec$design_columns, "(Intercept)"))
      "The model estimates one common normal ability mean and SD across origins; equal origin distributions are a substantive restriction, not scale coding." else
      "The declared normal regression models ability means by origin with a common residual SD; identification depends on the specified mean design, not just the MML label.")
    output$notes <- c(output$notes,
      "Interpret the zero origin main effect jointly with the population mean restrictions. A free-Person response-map basis does not make MML estimation distribution-free.")
  }
  output$table <- mfrmr_comparison_basis(output$table)
  output$parameter_map <- adjacent$map
  output$target_map <- target_map(adjacent$map)
  output
}

mfrmr_print_comparison_review <- function(review) {
  if (is.null(review)) return(invisible(NULL))
  cat("\nProspective interaction-comparison review\n")
  table <- mfrmr_comparison_basis(review$table)
  print(table[c("Comparison", "Group1", "Group2", "State", "ConditionalState", "IdentificationBasis")], row.names = FALSE)
  if (any(table$IdentificationBasis == "population_restriction"))
    cat("Population-restriction rows lose target identification with unrestricted Persons; global marginal identification remains unverified.\n")
  cat("Separate from displayed residual/refit effects; no estimates or formal inference supplied.\n")
  invisible(review)
}

mfrmr_bind_comparison_review <- function(review, fit, data, facet, group, focal = NULL) {
  if (is.null(review)) return(NULL)
  if (inherits(review, "mfrm_data_description")) review <- review$comparison_review
  if (!is.list(review) || !identical(review$contract_version, 1L) ||
      !is.data.frame(review$table) || !is.list(review$source)) {
    stop("Supply the comparison_review from describe_mfrm_data(comparisons = ...).", call. = FALSE)
  }
  s <- review$specification
  if (!identical(s$facet, facet) || !identical(s$group, group) ||
      (!is.null(focal) && !identical(s$groups[1], focal)) ||
      !setequal(review$source$facets, union(fit$prep$facet_names, group)) ||
      !identical(review$source$score_map, fit$prep$score_map)) {
    stop("Comparison review facet, group direction, score coding or response facets do not match this DFF request.", call. = FALSE)
  }
  matched <- .dff_match_observations(fit, data, group)
  rows <- matched$data
  rows[[group]] <- rows$.group_var
  rows <- mfrmr_comparison_rows(rows, review$source$facets)
  if (!identical(rows, review$source$rows)) stop(
    "Comparison review does not describe the same retained responses and origin groups; rerun describe_mfrm_data() on the requested ratings.", call. = FALSE)
  review
}

mfrmr_match_comparison_model <- function(prep, cfg, review) {
  if (inherits(review, "mfrm_data_description")) review <- review$comparison_review
  if (!is.list(review) || !identical(review$contract_version, 1L) ||
      !is.list(review$specification) || !is.list(review$source)) {
    stop("Supply the comparison_review from describe_mfrm_data(comparisons = ...).", call. = FALSE)
  }
  if (!setequal(review$source$facets, prep$facet_names) ||
      !identical(review$source$score_map, prep$score_map) ||
      !identical(review$source$rows, mfrmr_comparison_rows(prep$data, review$source$facets))) {
    stop("Comparison review does not match the joint fit's retained responses, origin groups or score coding.", call. = FALSE)
  }
  # Rebuild in the fit's level order; a saved target map uses its own coordinates.
  current <- mfrmr_build_comparison_review(prep, review$specification)
  expected <- mfrmr_comparison_model(prep, current$specification)$config
  keys <- c("model", "method", "n_cat", "n_person", "person_levels", "facet_names", "facet_levels",
    "step_facet", "slope_facet", "facet_signs",
    "theta_spec", "facet_specs", "step_specs", "facet_interactions")
  same <- function(x, y) isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
  interaction_keys <- c("name", "facet_a", "facet_b", "levels_a", "levels_b", "n_a", "n_b",
    "n_params", "identification")
  interaction_contract <- function(x) lapply(x$interaction_specs, `[`, interaction_keys)
  population_design <- function(x) {
    p <- x$population_spec
    if (!isTRUE(p$active)) return(NULL)
    p$design_matrix[p$person_lookup, , drop = FALSE]
  }
  if (!same(cfg[keys], expected[keys]) ||
      !same(interaction_contract(cfg), interaction_contract(expected)) ||
      !same(population_design(cfg), population_design(expected))) {
    stop("The joint fit must use the review's model, constraints, interaction, step owner, orientation and population design.", call. = FALSE)
  }
  current
}

mfrmr_fitted_comparison_table <- function(fit, review) {
  cfg <- fit$config; prep <- fit$prep
  if (inherits(fit, "mfrm_imported_fit") || !cfg$model %in% c("RSM", "PCM") ||
      !cfg$method %in% c("JML", "MML") || mfrm_has_jml_adjustment(fit) ||
      !identical(cfg$facet_shrinkage, "none")) {
    stop("Fitted comparisons require a native, unadjusted RSM/PCM joint fit without shrinkage.", call. = FALSE)
  }
  current <- mfrmr_match_comparison_model(prep, cfg, review)
  same <- function(x, y) isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
  sizes <- build_param_sizes(cfg)
  par <- fit$opt$par
  if (length(par) != sum(unlist(sizes)) || any(!is.finite(par))) {
    stop("The joint fit must retain a finite parameter vector with matching coordinates.", call. = FALSE)
  }
  params <- expand_params(par, sizes, cfg)
  cells <- build_interaction_effect_table(cfg, prep, params)
  saved <- fit$interactions$effects
  cell_keys <- c("Interaction", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level", "Estimate")
  if (!all(cell_keys %in% names(saved)) || !same(saved[cell_keys], cells[cell_keys])) {
    stop("Saved interaction effects do not match the joint fit's retained parameters.", call. = FALSE)
  }
  map <- current$parameter_map
  if (is.null(map) || anyNA(map$OptimizerIndex) || any(map$OptimizerIndex > length(par))) {
    stop("The comparison review has no matching optimizer-coordinate map.", call. = FALSE)
  }
  table <- current$table
  table$OptimizerEstimate <- drop(current$target_map %*% par[map$OptimizerIndex])
  ready <- mfrmr_get_readiness_record(fit)
  for (key in c("FitReadiness", "InputState", "EstimabilityState", "CategoryState", "BoundaryState", "NumericalState")) {
    table[[key]] <- ready$fit[[key]][1L]
  }
  table$EstimationConverged <- identical(as.integer(fit$opt$convergence), 0L)
  table$EstimateStatus <- "unavailable"
  table$EstimateReason <- "Full-model identification or retained fit readiness does not permit extraction."
  table$SE <- table$Lower <- table$Upper <- NA_real_

  idx <- build_indices(prep, cfg$step_facet, interaction_specs = cfg$interaction_specs)
  mml <- identical(cfg$method, "MML")
  q <- if (mml) cfg$estimation_control$quad_points else 1L
  if (length(q) != 1L || !is.finite(q) || q < 1 || q != floor(q)) {
    stop("The joint fit must retain its integration order.", call. = FALSE)
  }
  evaluate <- function(order) {
    cache <- make_param_cache(sizes, cfg, idx, is_mml = mml)
    evaluator <- make_mfrm_direct_evaluator(cfg$method, cache, idx, cfg, sizes,
      if (mml) gauss_hermite_normal(order) else NULL)
    c(Objective = evaluator$value(par), GradientMaxAbs = max(abs(evaluator$gradient(par))))
  }
  numerical <- tryCatch({
    base <- evaluate(q)
    fine <- if (mml) evaluate(2L * q + 1L) else base
    c(base, ObjectiveDifference = unname(abs(base[1L] - fit$opt$value)),
      FinerGradientMaxAbs = unname(fine[2L]),
      QuadratureChangePerPerson = unname(abs(fine[1L] - base[1L]) / cfg$n_person))
  }, error = identity)
  table$GradientMaxAbs <- table$FinerGradientMaxAbs <- table$QuadratureChangePerPerson <- NA_real_
  table$QuadratureState <- if (mml) "not_evaluated" else "not_applicable"
  if (!inherits(numerical, "error")) {
    table$GradientMaxAbs <- unname(numerical["GradientMaxAbs"])
    table$FinerGradientMaxAbs <- if (mml) unname(numerical["FinerGradientMaxAbs"]) else NA_real_
    table$QuadratureChangePerPerson <- if (mml) unname(numerical["QuadratureChangePerPerson"]) else NA_real_
    gradient_tolerance <- max(1e-4, 10 * (cfg$estimation_control$reltol %||% 1e-9))
    numerical_ok <- all(is.finite(numerical)) && numerical[3L] <= 1e-6 &&
      numerical[2L] <= gradient_tolerance
    quadrature_ok <- !mml || (numerical[4L] <= gradient_tolerance && numerical[5L] <= 1e-6)
    if (mml) table$QuadratureState <- if (quadrature_ok) "stable_at_retained_point" else "review"
    eligible <- table$ConditionalState == "all_parameters_identified" &
      table$EstimationConverged & table$NumericalState == "ready" &
      table$FitReadiness %in% c("ready", "review") &
      table$InputState %in% c("pass", "review") &
      table$CategoryState %in% c("adequate", "weak_information", "not_applicable") &
      table$BoundaryState %in% c("finite", "not_applicable", "not_evaluated") &
      numerical_ok & quadrature_ok
    table$Estimate[eligible] <- table$OptimizerEstimate[eligible]
    table$EstimateStatus[eligible] <- "descriptive_only"
    table$EstimateReason[eligible] <- if (mml) {
      "Fitted under the declared population model; global marginal identification and formal uncertainty are not certified."
    } else {
      "Identified conditional-model contrast at the converged fit; formal uncertainty and finite-MLE existence are not certified by this table."
    }
    if (!numerical_ok || !quadrature_ok) table$EstimateReason <-
      "Retained objective/gradient or finer-grid integration check failed; optimizer contrasts are diagnostic traces only."
  } else {
    table$EstimateReason <- paste("Retained-point numerical check failed:", conditionMessage(numerical))
  }
  notes <- c(
    "Estimate is a descriptive joint-fit contrast after source, full conditional rank and numerical checks; OptimizerEstimate is a trace, including when Estimate is unavailable.",
    "No partial model is fitted. SE and interval endpoints remain missing; FormalInferenceEligible remains false.",
    "Finer MML integration is checked at the retained point without refitting; stability is not a global identification or precision certificate.")
  checks <- list(readiness = ready, numerical = if (inherits(numerical, "error"))
    conditionMessage(numerical) else numerical, quadrature = if (mml) c(q, 2L*q+1L) else NULL,
    source = mfrm_checkpoint_fingerprint(mfrm_gpcm_inference_source(fit), canonical = TRUE))
  out <- tibble::as_tibble(table)
  attr(out, "comparison_review") <- current
  attr(out, "fit_checks") <- checks
  attr(out, "notes") <- notes
  out
}
