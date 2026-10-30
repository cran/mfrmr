# Source from the bundle directory; mfrmr is not needed for fitting or replay.
if (!requireNamespace("brms", quietly = TRUE)) stop("Install brms first.")
bundle <- readRDS("model-data.rds")
d <- bundle$stan_data
model_data <- data.frame(y = d$y, row_id = seq_len(d$N), constant = 1,
  known = as.vector(d$eta_offset),
  person = factor(paste0("p", d$person), levels = paste0("p", seq_len(d$J))))
z_names <- if (d$P) paste0("z", seq_len(d$P)) else character()
for (p in seq_len(d$P)) model_data[[z_names[p]]] <- d$X[, p]
step_names <- if (d$T) paste0("step", seq_len(d$T), "offset") else character()
slope_names <- if (d$S) paste0("slope", seq_len(d$S), "offset") else character()
coordinate_names <- c(step_names, slope_names)

# Closures retain only the model maps required by prediction after reloading.
family <- local({
  basis <- d$step_basis
  offset <- as.vector(d$step_offset)
  owner <- as.vector(d$owner)
  slope_X <- d$slope_X
  category_count <- d$C
  step_dpars <- step_names
  slope_dpars <- slope_names
  log_probabilities <- function(i, prep) {
    row <- prep$data$vint1[i]
    if (length(row) != 1L || is.na(row) || row < 1L || row > length(owner))
      stop("Prediction requires a retained model row identifier.")
    mu <- brms::get_dpar(prep, "mu", i = i)
    rows <- (owner[row] - 1L) * (category_count - 1L) + seq_len(category_count - 1L)
    steps <- matrix(rep(offset[rows], each = length(mu)), length(mu))
    for (q in seq_along(step_dpars))
      steps <- steps + outer(brms::get_dpar(prep, step_dpars[q], i = i), basis[rows, q])
    log_slope <- rep(0, length(mu))
    for (q in seq_along(slope_dpars))
      log_slope <- log_slope + brms::get_dpar(prep, slope_dpars[q], i = i) * slope_X[row, q]
    logits <- matrix(0, length(mu), category_count)
    for (k in 2:category_count)
      logits[, k] <- logits[, k - 1L] + mu - steps[, k - 1L]
    logits <- logits * exp(log_slope)
    logits <- logits - apply(logits, 1L, max)
    logits - log(rowSums(exp(logits)))
  }
  probabilities <- function(i, prep) exp(log_probabilities(i, prep))
  brms::custom_family("mfrmraffine", dpars = c("mu", coordinate_names),
    links = rep("identity", 1L + length(coordinate_names)), type = "int",
    vars = c("vint1[n]", "mfrmr_basis", "mfrmr_offset", "mfrmr_owner", "mfrmr_slopeX", "mfrmr_C"),
    log_lik = function(i, prep) log_probabilities(i, prep)[, prep$data$Y[i]],
    posterior_predict = function(i, prep, ...) apply(probabilities(i, prep), 1L,
      function(p) sample.int(length(p), 1L, prob = p)),
    posterior_epred = function(prep) {
      n <- length(prep$data$Y)
      first <- probabilities(1L, prep)
      out <- array(NA_real_, c(nrow(first), n, category_count))
      out[, 1L, ] <- first
      if (n > 1L) for (i in 2:n) out[, i, ] <- probabilities(i, prep)
      out
    })
})
person_term <- "person"
if (d$J == 1L) { model_data$personp1 <- 1; person_term <- "personp1" }
main_formula <- stats::reformulate(c(person_term, z_names, "offset(known)"),
  response = "y | vint(row_id)", intercept = FALSE)
coordinate_formulas <- lapply(coordinate_names, function(s) stats::as.formula(paste(s, "~ 0 + constant")))
formula <- do.call(brms::bf, c(list(main_formula), coordinate_formulas, list(family = family)))
prior <- brms::set_prior("normal(0, 1)", class = "b")
for (p in seq_len(d$P)) prior <- c(prior,
  brms::set_prior(sprintf("normal(0, %.17g)", d$location_sd[p]), class = "b", coef = z_names[p]))
coordinate_sd <- c(d$step_sd, d$log_slope_sd)
for (q in seq_along(coordinate_names)) prior <- c(prior,
  brms::set_prior(sprintf("normal(0, %.17g)", coordinate_sd[q]),
    class = "b", coef = "constant", dpar = coordinate_names[q]))
arguments <- paste(c("int y", "real mu", paste("real", coordinate_names), "int obs",
  "matrix basis", "vector step_offset", "array[] int owner", "matrix slopeX", "int C"), collapse = ", ")
# paste() on a zero-length vector would otherwise introduce an empty real argument.
if (!length(coordinate_names)) arguments <- paste("int y, real mu, int obs, matrix basis,",
  "vector step_offset, array[] int owner, matrix slopeX, int C")
stan_functions <- paste(c(
  paste0("real mfrmraffine_lpmf(", arguments, ") {"),
  "  vector[C] logits;", "  real log_slope = 0;", "  logits[1] = 0;",
  if (d$S) paste0("  log_slope += ", slope_names, " * slopeX[obs, ", seq_len(d$S), "];") else character(),
  "  for (k in 2:C) {", "    int row = (owner[obs] - 1) * (C - 1) + k - 1;",
  "    real step = step_offset[row];",
  if (d$T) paste0("    step += ", step_names, " * basis[row, ", seq_len(d$T), "];") else character(),
  "    logits[k] = logits[k - 1] + mu - step;", "  }",
  "  return categorical_logit_lpmf(y | exp(log_slope) * logits);", "}"), collapse = "\n")
stanvars <- brms::stanvar(scode = stan_functions, block = "functions") +
  brms::stanvar(d$C, name = "mfrmr_C", scode = "int<lower=2> mfrmr_C;") +
  brms::stanvar(d$step_basis, name = "mfrmr_basis") +
  brms::stanvar(as.vector(d$step_offset), name = "mfrmr_offset") +
  brms::stanvar(as.integer(d$owner), name = "mfrmr_owner") +
  brms::stanvar(d$slope_X, name = "mfrmr_slopeX")

extract_mfrm_parameters <- function(fit) {
  draws <- unclass(posterior::as_draws_matrix(fit))
  theta <- draws[, paste0("b_personp", seq_len(d$J)), drop = FALSE]
  free <- if (d$P) draws[, paste0("b_", z_names), drop = FALSE] else matrix(numeric(), nrow(draws), 0L)
  locations <- sweep(free %*% t(d$location_basis), 2L, d$location_offset, `+`)
  step <- if (d$T) draws[, paste0("b_", step_names, "_constant"), drop = FALSE] else matrix(numeric(), nrow(draws), 0L)
  step <- sweep(step %*% t(d$step_basis), 2L, d$step_offset, `+`)
  steps <- array(0, c(nrow(draws), d$G, d$C - 1L))
  for (g in seq_len(d$G)) steps[, g, ] <- step[, (g - 1L) * (d$C - 1L) + seq_len(d$C - 1L), drop = FALSE]
  slope <- if (d$S) draws[, paste0("b_", slope_names, "_constant"), drop = FALSE] else matrix(numeric(), nrow(draws), 0L)
  log_slopes <- slope %*% t(d$slope_basis)
  list(theta = theta, locations = locations, steps = steps,
    log_component_slopes = log_slopes, component_slopes = exp(log_slopes),
    levels = bundle$levels, step_owners = bundle$step_owners, slope_levels = bundle$slope_levels)
}
