# Shared by both exported runners; no mfrmr installation is needed for replay.
mfrm_population_label <- function(population = NULL) {
  if (is.null(population)) return("Fixed N(0,1)")
  if (!is.list(population) || !identical(population$mode, "estimated_normal") ||
      !identical(population$formula, "~1") ||
      !identical(population$parameterization, "noncentered"))
    stop("Unsupported exported population specification.")
  "Normal population with estimated mean and SD"
}

mfrm_parameter_map <- function(bundle) {
  d <- bundle$stan_data
  population_label <- mfrm_population_label(bundle$population)
  free_population <- !is.null(bundle$population)
  rows <- list()
  add <- function(variable, quantity, facet, level, fixed, value, anchored, reference, direction, transition = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(variable, quantity, facet, level, transition,
      value_role = if (!fixed) "estimated" else if (anchored) "fixed_anchor" else "fixed_by_constraint",
      fixed_value = if (fixed) value else NA_real_, reference, direction,
      stringsAsFactors = FALSE)
  }
  for (i in seq_len(nrow(bundle$levels))) {
    r <- bundle$levels[i, ]
    person <- startsWith(r$Parameter, "theta[")
    fixed <- !person && all(d$location_basis[r$Index, ] == 0)
    value <- if (fixed && !is.null(d$location_offset)) d$location_offset[r$Index] else 0
    spec <- bundle$constraints$facets[[r$Facet]]
    anchored <- !is.null(spec) && !is.na(spec$anchors[match(r$Level, spec$levels)])
    reference <- if (person && free_population) "Estimated normal population; native relative-slope metric; observed Persons are not centered." else
      if (person) "Population N(0,1); observed Persons are not forced to average zero." else
      if (r$Centered) "Saved facet centering constraints; zero is not a performance standard." else
        if (free_population) "Uncentered location on the native relative-slope metric." else
          "Uncentered location on the fixed N(0,1) ability scale."
    add(r$Parameter, if (person) "ability" else "location", r$Facet, r$Level,
      fixed, value, anchored, reference,
      if (r$Sign == 1) "Higher values favor higher coded scores." else "Higher values favor lower coded scores.")
  }
  if (free_population) {
    add("population_mean", "population_mean", "Population", "Mean", FALSE, NA_real_, FALSE,
      population_label, "Higher values favor higher coded scores; not necessarily better performance.")
    add("population_sd", "population_sd", "Population", "SD", FALSE, NA_real_, FALSE,
      "Native metric with geometric-mean-one relative slopes; not an individual Person SE.",
      "Larger values indicate more between-Person variation in the assumed population.")
  }
  for (g in seq_len(d$G)) for (k in seq_len(d$C - 1L)) {
    index <- (g - 1L) * (d$C - 1L) + k
    basis_row <- if (isTRUE(bundle$affine)) index else k
    fixed <- all(d$step_basis[basis_row, ] == 0)
    value <- if (fixed && !is.null(d$step_offset)) d$step_offset[index] else 0
    spec <- bundle$constraints$steps[[g]]
    anchored <- !is.null(spec) && !is.na(spec$anchors[k])
    add(paste0("steps[", g, ",", k, "]"), "step_offset", "Step owner", bundle$step_owners[g],
      fixed, value, anchored, paste0("Transition ", bundle$categories$Score[k], " to ",
        bundle$categories$Score[k + 1L], "; saved within-owner sum constraint."),
      "Higher offsets reduce the adjacent-category odds of the higher coded score.",
      paste(bundle$categories$Score[k], "->", bundle$categories$Score[k + 1L]))
  }
  if (!is.null(bundle$slope_levels)) for (i in seq_len(nrow(bundle$slope_levels))) {
    r <- bundle$slope_levels[i, ]
    fixed <- all(d$slope_basis[r$Index, ] == 0)
    reference <- if (r$ScaleReference == "geometric_mean_one")
      "Owner geometric mean slope = 1; one is a reference, not a quality standard." else
        "Component scale identified by the fixed N(0,1) ability population."
    add(r$Parameter, "slope", r$Owner, r$Level, fixed, 1, FALSE, reference,
      "Steepness of the response relationship, not rater severity.")
    add(r$LogParameter, "log_slope", r$Owner, r$Level, fixed, 0, FALSE, reference,
      "Log component slope, not rater severity.")
  }
  do.call(rbind, rows)
}

mfrm_sampling_review <- function(parameters, sampler, chains) {
  checks <- list()
  add <- function(check, value, criterion, ok, action) {
    checks[[length(checks) + 1L]] <<- data.frame(Check = check, Value = value,
      Criterion = criterion, Status = if (is.na(ok)) "not_available" else if (ok) "no_flag" else "review",
      NextAction = if (isTRUE(ok)) "No action from this check." else action, stringsAsFactors = FALSE)
  }
  p <- parameters[!parameters$value_role %in% c("fixed_anchor", "fixed_by_constraint"), ]
  add("Chains", chains, ">= 4", chains >= 4, "Use at least four chains for the final sampling review.")
  for (name in c("rhat", "ess_bulk", "ess_tail")) {
    v <- p[[name]]
    known <- length(v) > 0L && all(is.finite(v))
    value <- if (!known) NA_real_ else if (name == "rhat") max(v) else min(v)
    add(name, value, if (name == "rhat") "< 1.01 for all varying coordinates" else paste0(">= ", 100 * chains, " (100 per chain)"),
      if (!known) NA else if (name == "rhat") value < 1.01 else value >= 100 * chains,
      "Inspect all varying coordinates and chain traces; diagnose poor mixing before adding draws.")
  }
  complete <- nrow(sampler) == chains && "chain" %in% names(sampler) &&
    !anyNA(sampler$chain) && !anyDuplicated(sampler$chain)
  for (name in c("num_divergent", "num_max_treedepth", "ebfmi")) {
    v <- sampler[[name]]
    known <- complete && length(v) == chains && all(is.finite(v)) && all(v >= 0)
    value <- if (!known) NA_real_ else if (name == "ebfmi") min(v) else sum(v)
    add(name, value, if (name == "ebfmi") ">= 0.3 in every chain" else "0 across chains",
      if (!known) NA else if (name == "ebfmi") value >= 0.3 else value == 0,
      if (name == "num_max_treedepth") "Review sampling efficiency and geometry; tree-depth hits alone do not establish bias." else
        "Inspect sampler diagnostics and model geometry before interpreting posterior intervals.")
  }
  do.call(rbind, checks)
}

summarize_mfrm_draws <- function(draws, bundle, sampler, backend) {
  draws <- posterior::as_draws_array(draws)
  map <- mfrm_parameter_map(bundle)
  if (anyDuplicated(map$variable) || !all(map$variable %in% posterior::variables(draws)))
    stop("Draws do not contain the complete exported parameter map. Keep fit.rds and check the matching bundle.")
  # Fixed status comes from the model map, never from a zero sample variance.
  for (i in which(!is.na(map$fixed_value))) {
    v <- draws[, , map$variable[i]]
    if (any(!is.finite(v)) || any(abs(v - map$fixed_value[i]) > 1e-7 * max(1, abs(map$fixed_value[i]))))
      stop("Draws disagree with fixed constraints. Check that fit.rds and model-data.rds belong together.")
  }
  p <- as.data.frame(posterior::summarise_draws(draws,
    mean = base::mean, median = stats::median, sd = stats::sd,
    lower = function(x) unname(stats::quantile(x, .025)),
    upper = function(x) unname(stats::quantile(x, .975)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail,
    mcse_lower = function(x) unname(posterior::mcse_quantile(x, .025)),
    mcse_upper = function(x) unname(posterior::mcse_quantile(x, .975))))
  m <- map[match(p$variable, map$variable), setdiff(names(map), "variable"), drop = FALSE]
  unknown <- is.na(m$value_role)
  m$value_role[unknown] <- "sampling_coordinate"
  m$quantity[unknown] <- "internal_coordinate"
  m$reference[unknown] <- "Internal sampling coordinate; use mapped parameters for substantive reporting."
  p <- cbind(p, m); rownames(p) <- NULL
  p$role <- ifelse(unknown, "free_coordinate", "reported_parameter")
  fixed <- p$value_role %in% c("fixed_anchor", "fixed_by_constraint")
  p$interval_computed <- is.finite(p$lower) & is.finite(p$upper) & !fixed
  # A fixed value has no estimated uncertainty; raw constant draws remain in fit.rds.
  p[fixed, c("sd", "lower", "upper", "rhat", "ess_bulk", "ess_tail", "mcse_lower", "mcse_upper")] <- NA_real_
  width <- p$upper - p$lower
  p$endpoint_mcse_fraction <- ifelse(is.finite(width) & width > 0,
    pmax(p$mcse_lower, p$mcse_upper) / width, NA_real_)
  checks <- mfrm_sampling_review(p, sampler, dim(draws)[2L])
  flagged <- checks$Status != "no_flag"
  status <- if (any(flagged)) "sampling_review_required" else "no_sampling_flags_detected"
  p$sampling_review <- status
  p$interval_status <- ifelse(fixed, "fixed_value_no_interval", ifelse(!p$interval_computed,
    "unavailable", if (any(flagged)) "computed_sampling_review_required" else "computed_model_review_required"))
  first <- c("variable", "quantity", "facet", "level", "transition", "value_role", "mean", "median", "sd",
    "lower", "upper", "interval_status", "sampling_review", "fixed_value", "reference", "direction")
  p <- p[, c(first, setdiff(names(p), first))]
  reasons <- c(Chains = "Fewer than four chains", rhat = "Chain mixing (R-hat) needs review",
    ess_bulk = "Effective sample size for central summaries needs review",
    ess_tail = "Effective sample size for tail summaries needs review",
    num_divergent = "Divergent transitions", num_max_treedepth = "Tree-depth limit reached",
    ebfmi = "Energy exploration needs review")
  why <- unname(reasons[checks$Check])
  why[checks$Status == "not_available"] <- paste("Diagnostic unavailable:", checks$Check[checks$Status == "not_available"])
  decision <- data.frame(Status = status,
    Interpretation = if (any(flagged)) "Review sampling before interpreting intervals" else "No sampling flags detected; model review remains",
    Why = if (any(flagged)) paste(why[flagged], collapse = "; ") else "Listed sampling checks detected no flags.",
    NextAction = if (any(flagged)) "Read sampling-review.csv before interpreting estimates or intervals." else
      "Review model fit, priors and target-specific Monte Carlo error before reporting.")
  result <- list(schema_version = if (is.null(bundle$population)) 1L else 2L, chains = as.integer(dim(draws)[2L]),
    model_data_md5 = if (file.exists("model-data.rds")) unname(tools::md5sum("model-data.rds")) else NA_character_,
    parameters = p, sampler = sampler, decision = decision, checks = checks,
    backend = backend, interval = "equal-tailed posterior credible interval for varying parameters", level = .95,
    notes = c("Sampling checks are a screen, not a certificate of validity or calibrated coverage.",
      "Fixed anchors and constraint-determined values have no estimated interval; calibration uncertainty is not propagated.",
      "Computed intervals remain visible during review; examine all varying coordinates and chains.",
      "Endpoint MCSE is simulation error, not posterior uncertainty. Choose sufficient precision for the intended claim.",
      "Higher coded scores need not mean better performance; retain the original rubric direction.",
      "Review prior sensitivity and posterior predictive fit; Bayesian intervals do not validate MML coverage.",
      "Screening criteria: https://mc-stan.org/learn-stan/diagnostics-warnings.html"))
  for (name in c("parameters", "sampler", "checks", "decision")) {
    file <- c(parameters = "parameter-summary.csv", sampler = "sampler-diagnostics.csv",
      checks = "sampling-review.csv", decision = "summary-decision.csv")[[name]]
    utils::write.csv(result[[name]], file, row.names = FALSE, fileEncoding = "UTF-8")
  }
  saveRDS(result, "posterior-summary.rds")
  cat("Posterior review (", backend, ")\n", sep = "")
  cat(decision$Interpretation, "\nWhy: ", decision$Why, "\nNext: ", decision$NextAction, "\n", sep = "")
  cat("Fixed values have no estimated interval. Sampling checks do not establish model validity.\n")
  main <- p[p$role == "reported_parameter" & p$quantity != "log_slope", ]
  preview <- which(!duplicated(main[, c("quantity", "value_role")]))
  print(utils::head(main[preview, c("quantity", "facet", "level", "transition", "value_role", "mean", "lower", "upper", "interval_status")], 8L), row.names = FALSE)
  cat(nrow(main), "mapped parameters; all rows and diagnostic coordinates are in parameter-summary.csv.\n")
  invisible(result)
}

mfrm_brms_draws <- function(fit, parameters) {
  # Transform each preserved chain/iteration in the same chain-major order used by extract_mfrm_parameters().
  raw <- posterior::as_draws_array(fit)
  parts <- list()
  for (name in c("theta", "locations", "steps", "component_slopes", "log_component_slopes")) {
    x <- parameters[[name]]
    if (is.null(x) || length(x) == 0L) next
    dims <- dim(x)
    mat <- matrix(x, nrow = dims[1L])
    indices <- if (length(dims) == 3L) apply(expand.grid(seq_len(dims[2L]), seq_len(dims[3L])), 1L, paste, collapse = ",") else seq_len(ncol(mat))
    colnames(mat) <- paste0(name, "[", indices, "]")
    parts[[name]] <- mat
  }
  mapped <- do.call(cbind, parts)
  if (nrow(mapped) != prod(dim(raw)[1:2])) stop("Mapped brms draws lost their chain/iteration alignment.")
  transformed <- posterior::as_draws_array(array(mapped, c(dim(raw)[1:2], ncol(mapped)),
    dimnames = list(dimnames(raw)[[1]], dimnames(raw)[[2]], colnames(mapped))))
  internal <- posterior::subset_draws(raw, variable = grep("^(b_|lp__$)", posterior::variables(raw), value = TRUE))
  posterior::bind_draws(transformed, internal, along = "variable")
}

mfrm_brms_sampler <- function(fit) {
  nuts <- brms::nuts_params(fit)
  depth <- brms::control_params(fit)$max_treedepth
  chains <- posterior::nchains(fit)
  do.call(rbind, lapply(seq_len(chains), function(chain) {
    get <- function(name) {
      x <- nuts[nuts$Chain == chain & nuts$Parameter == name, ]
      x$Value[order(x$Iteration)]
    }
    divergent <- get("divergent__"); tree <- get("treedepth__"); energy <- get("energy__")
    expected <- posterior::niterations(fit)
    complete <- function(x) length(x) == expected && all(is.finite(x))
    data.frame(chain = chain,
      num_divergent = if (complete(divergent)) sum(divergent) else NA_real_,
      num_max_treedepth = if (complete(tree) && length(depth) == 1L && is.finite(depth)) sum(tree >= depth) else NA_real_,
      ebfmi = if (complete(energy) && length(energy) > 1L && stats::var(energy) > 0)
        mean(diff(energy)^2) / stats::var(energy) else NA_real_)
  }))
}
