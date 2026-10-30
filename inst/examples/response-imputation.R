# Regenerate the joint RSM completions used by the assigned-score MI tutorial.
# Requires cmdstanr, an already installed CmdStan toolchain, and posterior.
# No software is installed by this script. Run from an installed mfrmr session:
# source(system.file("examples", "response-imputation.R", package = "mfrmr"))
# example <- response_imputation_example("mi-example-output")
response_imputation_example <- function(directory,
  stan_file = system.file("examples", "response-imputation.stan", package = "mfrmr")) {
  stopifnot(requireNamespace("cmdstanr", quietly = TRUE),
    requireNamespace("posterior", quietly = TRUE), file.exists(stan_file))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(directory, all.files = TRUE, no.. = TRUE)))
    stop("Use an empty output directory to preserve earlier runs.")
  ratings <- mfrmr::load_mfrmr_data("example_core")
  ratings$Event <- paste0("E", seq_len(nrow(ratings)))
  ratings$Assigned <- TRUE
  unassigned <- which(ratings$Rater == "R04" & ratings$Criterion == "Content")[1:6]
  ratings$Assigned[unassigned] <- FALSE
  r01 <- tapply(ratings$Score[ratings$Rater == "R01"],
    ratings$Person[ratings$Rater == "R01"], mean)
  eligible <- which(ratings$Assigned & ratings$Rater == "R04" &
    ratings$Criterion %in% c("Content", "Language"))
  probability <- plogis(-.5 + (2.5 - r01[ratings$Person[eligible]]))
  set.seed(26230901)
  missing <- eligible[runif(length(eligible)) < probability]
  held_out <- ratings$Score[missing]
  ratings$Score[c(unassigned, missing)] <- NA_integer_
  observed <- which(ratings$Assigned & !is.na(ratings$Score))
  ids <- lapply(ratings[c("Person", "Rater", "Criterion")], function(x) sort(unique(x)))
  index <- lapply(names(ids), function(n) match(ratings[[n]], ids[[n]]))
  names(index) <- names(ids)
  orthonormal <- function(n) {
    q <- stats::contr.helmert(n)
    sweep(q, 2, sqrt(colSums(q^2)), "/")
  }
  data <- list(P = length(ids$Person), R = length(ids$Rater), C = length(ids$Criterion),
    K = 4L, N_obs = length(observed), N_mis = length(missing),
    person_obs = index$Person[observed], rater_obs = index$Rater[observed],
    criterion_obs = index$Criterion[observed], y_obs = ratings$Score[observed],
    person_mis = index$Person[missing], rater_mis = index$Rater[missing],
    criterion_mis = index$Criterion[missing], Q_rater = orthonormal(length(ids$Rater)),
    Q_criterion = orthonormal(length(ids$Criterion)), Q_step = orthonormal(3L), prior_sd = 2.5)
  settings <- list(chains = 4L, parallel_chains = 4L, iter_warmup = 1000L,
    iter_sampling = 2000L, seed = 26230902L, adapt_delta = .9, max_treedepth = 10L,
    refresh = 500L, sig_figs = 18L, output_dir = normalizePath(directory),
    output_basename = "joint-rsm")
  specification <- list(method = "Joint RSM posterior predictive draws",
    data = data, ids = ids, code = readLines(stan_file),
    settings = settings[setdiff(names(settings), "output_dir")],
    software = c(cmdstanr = as.character(utils::packageVersion("cmdstanr")),
      cmdstan = as.character(cmdstanr::cmdstan_version()),
      posterior = as.character(utils::packageVersion("posterior"))),
    selection = list(eligible = ratings$Event[eligible], probability = unname(probability),
      seed = 26230901L, assumption = "MAR conditional on observed R01 scores"))
  saveRDS(list(ratings = ratings, missing = missing, unassigned = unassigned,
    held_out = held_out, specification = specification), file.path(directory, "input.rds"))
  model_path <- file.path(directory, "response-imputation.stan")
  writeLines(specification$code, model_path)
  model <- cmdstanr::cmdstan_model(model_path, cpp_options = list(PRECOMPILED_HEADERS = "false"))
  fit <- do.call(model$sample, c(list(data = data), settings))
  fit$save_object(file.path(directory, "posterior.rds"))
  variables <- c("theta", "rater_free", "criterion_free", "step_free",
    "severity", "difficulty", "step")
  diagnostics <- fit$summary(variables)
  sampler <- fit$diagnostic_summary()
  saveRDS(list(diagnostics = diagnostics, sampler = sampler), file.path(directory, "diagnostics.rds"))
  good <- all(fit$return_codes() == 0) &&
    all(is.finite(diagnostics$rhat) & diagnostics$rhat < 1.01) &&
    all(diagnostics$ess_bulk >= 400 & diagnostics$ess_tail >= 400) &&
    all(sampler$num_divergent == 0 & sampler$num_max_treedepth == 0 & sampler$ebfmi >= .3)
  if (!isTRUE(good)) stop("Review the retained sampling diagnostics before using any completions.")
  draws <- posterior::as_draws_df(fit$draws(c("y_mis", "probability_mis")))
  set.seed(26230903)
  selected <- sample.int(nrow(draws), 40L)
  scores <- t(as.matrix(as.data.frame(draws)[selected, paste0("y_mis[", seq_along(missing), "]")]))
  storage.mode(scores) <- "integer"
  specification$diagnostics <- diagnostics
  specification$sampler <- sampler
  specification$draw_ids <- as.data.frame(draws)[selected, c(".chain", ".iteration", ".draw")]
  specification$completion_seed <- 26230903L
  # Use all posterior probability draws for the held-out check, not just 40 categories.
  probability_mean <- vapply(1:4, function(k) vapply(seq_along(missing), function(i)
    mean(draws[[sprintf("probability_mis[%d,%d]", i, k)]]), numeric(1)), numeric(length(missing)))
  severity <- as.matrix(fit$draws("severity", format = "draws_matrix"))
  difference <- as.numeric(severity[, sprintf("severity[%d]", match("R04", ids$Rater))] -
    severity[, sprintf("severity[%d]", match("R01", ids$Rater))])
  reference <- data.frame(Method = "Observed-score Bayes", Estimate = mean(difference),
    SE = stats::sd(difference), Lower = unname(stats::quantile(difference, .025)),
    Upper = unname(stats::quantile(difference, .975)))
  example <- list(ratings = ratings, missing = missing, unassigned = unassigned,
    held_out = held_out, scores = scores, probability_mean = probability_mean,
    imputation_model = specification, reference = reference)
  saveRDS(example, file.path(directory, "completions.rds"), compress = "xz")
  example
}
