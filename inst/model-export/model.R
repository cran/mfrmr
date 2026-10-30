# Source from the bundle directory. No mfrmr installation is needed.
if (!requireNamespace("brms", quietly = TRUE)) stop("Install brms first.")
bundle <- readRDS("model-data.rds")
d <- bundle$stan_data
model_data <- data.frame(y = d$y,
                        person = factor(paste0("p", d$person), levels = paste0("p", seq_len(d$J))),
                        owner = factor(paste0("g", d$owner), levels = paste0("g", seq_len(d$G))))
z_names <- paste0("z", seq_len(d$P))
if (!d$P) z_names <- character()
for (p in seq_len(d$P)) model_data[[z_names[p]]] <- d$X[, p]
step_names <- if (d$Q) paste0("step", seq_len(d$Q), "offset") else character()

# Closure state is retained when the fitted brms object is saved and reloaded.
family <- local({
  basis <- d$step_basis
  category_count <- d$C
  step_dpars <- step_names
  log_probabilities <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    offsets <- matrix(0, length(mu), category_count - 1L)
    for (q in seq_along(step_dpars)) {
      offsets <- offsets + outer(brms::get_dpar(prep, step_dpars[q], i = i), basis[, q])
    }
    logits <- matrix(0, length(mu), category_count)
    for (k in 2:category_count) logits[, k] <- logits[, k - 1L] + mu - offsets[, k - 1L]
    logits <- logits - apply(logits, 1L, max)
    logits - log(rowSums(exp(logits)))
  }
  probabilities <- function(i, prep) exp(log_probabilities(i, prep))
  brms::custom_family(
    "mfrmrpcm", dpars = c("mu", step_dpars),
    links = rep("identity", 1L + length(step_dpars)), type = "int",
    vars = c("mfrmr_basis", "mfrmr_C"),
    log_lik = function(i, prep) log_probabilities(i, prep)[, prep$data$Y[i]],
    posterior_predict = function(i, prep, ...) {
      prob <- probabilities(i, prep)
      apply(prob, 1L, function(p) sample.int(length(p), 1L, prob = p))
    },
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
if (d$J == 1L) {
  model_data$personp1 <- 1
  person_term <- "personp1"
}
main_formula <- stats::reformulate(c(person_term, z_names), response = "y", intercept = FALSE)
step_formulas <- lapply(step_names, function(s) {
  stats::as.formula(paste(s, if (d$G == 1L) "~ 1" else "~ 0 + owner"))
})
formula <- do.call(brms::bf, c(list(main_formula), step_formulas, list(family = family)))
prior <- brms::set_prior("normal(0, 1)", class = "b")
for (p in seq_len(d$P)) {
  prior <- c(prior, brms::set_prior(sprintf("normal(0, %.17g)", d$location_sd[p]),
                                   class = "b", coef = z_names[p]))
}
for (s in step_names) {
  prior <- c(prior, brms::set_prior(sprintf("normal(0, %.17g)", d$step_sd),
                                   class = if (d$G == 1L) "Intercept" else "b", dpar = s))
}
arguments <- paste(c("int y", "real mu", paste("real", step_names), "matrix basis", "int C"), collapse = ", ")
if (!length(step_names)) arguments <- "int y, real mu, matrix basis, int C"
coordinates <- if (d$Q) paste0("  coordinates[", seq_len(d$Q), "] = ", step_names, ";") else character()
stan_functions <- paste(c(
  paste0("real mfrmrpcm_lpmf(", arguments, ") {"),
  paste0("  vector[", d$Q, "] coordinates;"),
  "  vector[C - 1] steps;", "  vector[C] logits;", coordinates,
  "  steps = basis * coordinates;", "  logits[1] = 0;",
  "  for (k in 2:C) logits[k] = logits[k - 1] + mu - steps[k - 1];",
  "  return categorical_logit_lpmf(y | logits);", "}"), collapse = "\n")
stanvars <- brms::stanvar(scode = stan_functions, block = "functions") +
  brms::stanvar(d$C, name = "mfrmr_C", scode = "int<lower=2> mfrmr_C;") +
  brms::stanvar(d$step_basis, name = "mfrmr_basis",
                scode = paste0("matrix[mfrmr_C - 1, ", d$Q, "] mfrmr_basis;"))

# Convert posterior coefficients to the original facet/step parameterization.
extract_mfrm_parameters <- function(fit) {
  draws <- unclass(posterior::as_draws_matrix(fit))
  theta <- draws[, paste0("b_personp", seq_len(d$J)), drop = FALSE]
  free <- if (d$P) draws[, paste0("b_", z_names), drop = FALSE] else
    matrix(numeric(), nrow(draws), 0L)
  locations <- free %*% t(d$location_basis)
  steps <- array(0, c(nrow(draws), d$G, d$C - 1L))
  for (g in seq_len(d$G)) {
    cols <- if (d$G == 1L) paste0("b_", step_names, "_Intercept") else
      paste0("b_", step_names, "_ownerg", g)
    if (d$Q) steps[, g, ] <- draws[, cols, drop = FALSE] %*% t(d$step_basis)
  }
  list(theta = theta, locations = locations, steps = steps,
       levels = bundle$levels, step_owners = bundle$step_owners)
}
