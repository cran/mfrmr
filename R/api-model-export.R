#' Export a measurement model to Stan or brms
#'
#' Write a self-contained Bayesian model bundle without running an estimator.
#' Supports ordinary RSM/PCM and fixed-population one-family GPCM MML fits in
#' both formats. Estimated-population one-family GPCM and two-family GPCM use
#' Stan. Fixed facet/group and step anchors are retained.
#' Sources must be unit-weight, with a fixed standard-normal Person population
#' or (Stan only, one-family GPCM) an estimated intercept-only normal population,
#' and no facet interactions or
#' post-fit shrinkage. Unsupported specifications fail before files are written.
#'
#' @param fit A native RSM/PCM or GPCM MML `mfrm_fit` returned by
#'   [fit_mfrm()].
#' @param output_dir Destination directory. The bundle includes observed ratings
#'   and original level labels; review it before sharing.
#' @param format Either `"stan"` (CmdStanR) or `"brms"` (custom family).
#' @param location_prior_sd Positive marginal prior SD for each varying facet
#'   location; fixed or uniquely determined locations have zero prior SD.
#' @param step_prior_sd Positive marginal prior SD for each varying step offset.
#' @param overwrite Replace bundle files that already exist? Defaults to `FALSE`.
#' @param acknowledge_sensitive Suppress the reminder that ratings and original
#'   identifiers are included. Defaults to `FALSE`.
#' @param log_slope_prior_sd Positive marginal prior SD of each component log
#'   slope for GPCM export. Defaults to 1; inspect prior
#'   sensitivity. Not used by RSM/PCM models.
#' @param population_prior Normal prior parameters for an estimated Person
#'   population: a named list with `mean = c(location, scale)` for its mean
#'   and `log_sd = c(location, scale)` for the log of its standard deviation.
#'   Each scale must be positive. Defaults give `mu ~ normal(0, 1)` and
#'   `log(sigma) ~ normal(0, 0.5)`. Unused for fixed-population exports.
#'
#' @details
#' The default one-family GPCM estimates its population mean and variance.
#' Its intercept-only population can be exported to Stan, retaining the native
#' metric, relative-slope constraint and supported fixed anchors. The generated
#' model uses `theta = mu + sigma * z`, with independent standard-normal `z`.
#' The population mean and SD are sampled, not fixed at their MML estimates.
#' `population_prior` adds explicit Bayesian priors: by default sigma has median
#' 1 and a central 95% prior interval of approximately 0.38 to 2.66. These are
#' starting assumptions, not a population-size recommendation; inspect prior
#' predictions and sensitivity for your rubric and anchors. The displayed
#' Person prior SD includes uncertainty in both mu and sigma, not only sigma.
#' Population regression, estimated-mean models with a fixed SD, and estimated-population
#' brms or two-family exports are not supported. Choosing
#' `gpcm_mml_identification = "fixed_standard_normal"` changes the model;
#' it is not a file-format conversion or a workaround for these exclusions.
#'
#' Both formats retain the category ladder, facet signs, facet centering and
#' RSM/PCM step ownership of the source specification. RSM/PCM adjacent-category log odds
#' are `theta + signed facet locations - step`. Steps follow their declared
#' within-owner sum targets (zero by default) and are **not ordered**. A binary
#' ladder has its single offset fixed at its target. Fixed-population ability
#' parameters have independent standard-normal priors; realized abilities are
#' not centered across the observed Persons in either population mode.
#'
#' All supported exports preserve the native logistic metric (L-metric),
#' with scaling constant \eqn{D = 1}. No extra 1.7 or 1.702 factor is used.
#' Standard-normal ability priors do not imply a normal-ogive response link.
#' The "Logistic metric (D = 1)" section of [fit_mfrm()] explains the scale
#' and the adjustments needed before comparing other parameterizations.
#'
#' Centered locations and steps use an orthonormal contrast basis. Independent
#' normal coordinates are scaled so each constrained component has the requested
#' marginal SD. A centered singleton is fixed at zero. Uncentered locations have
#' independent normal priors. These proper structural priors are additions to
#' MML: the Bayesian posterior is not a reproduction of MML point estimates or
#' frequentist confidence intervals. Review prior sensitivity for the application.
#'
#' Two-family GPCM is available with `format = "stan"` only. Its effective
#' slope is the product of the two component slopes and multiplies the entire
#' adjacent-category predictor, including locations and steps. The first
#' family's log slopes sum to zero; the second family's log slopes are free
#' on the fixed N(0,1) ability scale. Normal priors on orthonormal/free log-slope
#' coordinates give each component the specified marginal log-slope SD.
#' `slope-levels.csv` maps the generated `component_slopes` and
#' `log_component_slopes` to their owners and levels. This preserves the source
#' likelihood, not an MML interval or the exact priors/program of Uto and Ueno
#' (2020). In particular, the extra 1.7 multiplier in their Appendix is absent.
#' One-family GPCM supports both formats. Its single owner's log slopes sum
#' to zero (geometric mean one); its slope also multiplies the entire predictor.
#' A singleton slope owner is fixed at one.
#'
#' Fixed non-Person facet anchors, group means and step anchors are exact affine
#' constraints. The export preserves whole-facet versus remainder centering and
#' each step owner's target. Priors on the remaining free contrasts are centered
#' so each free level shares its block's remaining total equally. These priors
#' are defined on the constrained subspace, not obtained by conditioning an
#' earlier unanchored posterior. `levels.csv` and `steps.csv` in these bundles
#' record the affine prior means, marginal SDs and fixed-coordinate indicators;
#' `model-data.rds$constraints` retains the source constraints. One-family bundles
#' use the same affine representation even without anchors. Anchors are fixed
#' inputs: zero posterior spread does not establish certainty about their earlier
#' calibration. Uncertain external calibration is not propagated. Retained step
#' constraints are exported; this does not add a step-anchor fitting argument.
#' Native two-family fitting retains its anchor exclusion. Person anchors
#' require JML in native fitting and remain outside this MML exporter.
#'
#' `model-data.rds` records the design matrices, priors and mappings used by both
#' formats. `levels.csv` maps original labels to generated parameters;
#' `categories.csv` maps Stan/brms categories to the source score scale.
#' `run.R` fits the exported model and saves draws, diagnostics and replicated
#' ratings. Run it from the bundle directory, after installing CmdStanR/CmdStan
#' and, for the brms format, brms. Neither is required to export a bundle.
#' Both formats support `Rscript run.R --summarize` to regenerate summaries
#' from `fit.rds` without compiling or sampling. Keep `model-data.rds`,
#' `summarize.R` and the brms `model.R` with the saved fit.
#' Read `summary-decision.csv` first, then `sampling-review.csv` for the checks,
#' criteria and next actions. `parameter-summary.csv` joins original labels to
#' estimates, references and score-direction notes. `value_role` distinguishes
#' estimated values, fixed anchors and constraint-determined values. Fixed
#' values have no estimated SD or interval (`NA`); they do not establish zero
#' external calibration uncertainty. Varying parameters retain equal-tailed
#' 95% credible bounds, endpoint MCSE, R-hat and bulk/tail ESS, even during
#' review, with `interval_status` beside the values. `interval_computed` means
#' finite bounds for a varying parameter, not reliable sampling or coverage.
#' `sampler-diagnostics.csv` contains per-chain diagnostics;
#' `posterior-summary.rds` retains all four tables and interpretation notes.
#'
#' The screen checks at least four chains, R-hat below 1.01, bulk/tail ESS at
#' least 100 per chain, zero divergences/tree-depth hits and E-BFMI at least
#' 0.3. Missing diagnostics require review; fixed coordinates are excluded
#' using model constraints, not sample variance. Tree-depth hits flag efficiency
#' separately from divergences. These checks follow the
#' [Stan diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html)
#' and are not a certificate of validity. Review all varying coordinates,
#' model fit and priors. The larger endpoint MCSE divided by interval width
#' measures simulation precision; no universal acceptable fraction is imposed.
#' Read the matching summary and model-data files with [read_mfrm_posterior()]
#' to use `summary()`, `plot()`, stored credible intervals and [mfrm_report()]
#' without either sampling backend. Keep `fit.rds` for joint-draw calculations.
#' Generated brms callbacks support `log_lik()`, `posterior_predict()` and
#' `posterior_epred()`; the last returns category probabilities.
#'
#' Log likelihood and replicated ratings condition on the sampled abilities of
#' existing Persons. They do not provide leave-one-Person-out validation or
#' predictions for new Persons. Exporting does not qualify the original fit's
#' convergence, inference, or an external sampler's diagnostics. Repeated ratings
#' retain the source conditional-independence assumption.
#'
#' @return Invisibly, a list with `format`, `model`, `output_dir` and named
#'   `files`. No model is fitted and no software is installed.
#' @seealso [export_mfrm_results()], [fit_mfrm()]
#' @references
#' Uto, M., & Ueno, M. (2020). A generalized many-facet Rasch model and its
#' Bayesian estimation using Hamiltonian Monte Carlo. *Behaviormetrika, 47*,
#' 469-496. \doi{10.1007/s41237-020-00115-7}. Equation 9, Sections 6-7 and Appendix.
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("export_mfrm_model", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' # Illustrative fixed reference; applications need justified calibration values.
#' anchor <- data.frame(Facet = "Rater",
#'                      Level = as.character(ratings$Rater[1]), Anchor = 0)
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                  model = "GPCM", method = "MML", step_facet = "Criterion",
#'                  slope_facet = "Criterion", anchors = anchor,
#'                  anchor_policy = "silent")
#' # Keep the default estimated population. The example data are synthetic.
#' # Export needs neither Stan nor brms installed and runs no sampler.
#' stan <- export_mfrm_model(fit, tempfile("mfrmr-stan-"), "stan",
#'                            acknowledge_sensitive = TRUE)
#' names(stan$files)
#' bundle <- readRDS(stan$files[["model-data.rds"]])
#' bundle$categories # Generated categories mapped to original scores
#' bundle$population # Estimated mean/SD and their priors
#' bundle$priors
#' # Run this bundle's run.R separately if external sampling is wanted.
#' unlink(stan$output_dir, recursive = TRUE)
#' }
#' @export
export_mfrm_model <- function(fit, output_dir, format = c("stan", "brms"),
                              location_prior_sd = 1, step_prior_sd = 1,
                              overwrite = FALSE, acknowledge_sensitive = FALSE,
                              log_slope_prior_sd = 1,
                              population_prior = list(mean = c(0, 1), log_sd = c(0, 0.5))) {
  format <- match.arg(format)
  for (x in c("overwrite", "acknowledge_sensitive")) {
    value <- get(x)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", x, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be one nonempty path.", call. = FALSE)
  }
  product <- mfrm_has_product_slopes(fit) && identical(fit$config$model, "GPCM")
  if (product && format != "stan") stop(
    "Two-family GPCM model export supports format = 'stan' only.", call. = FALSE)
  bundle <- mfrm_model_export_data(fit, location_prior_sd, step_prior_sd,
                                  log_slope_prior_sd, population_prior)
  free_population <- !is.null(bundle$population)
  if (free_population && format != "stan") stop(paste(
    "Estimated-population one-family GPCM export supports format = 'stan' only.",
    "Keep the original fit; choosing fixed_standard_normal changes the model."
  ), call. = FALSE)
  template_dir <- system.file("model-export", package = "mfrmr")
  model_file <- if (format == "stan") "model.stan" else "model.R"
  templates <- c(if (free_population) "model-population.stan" else
                   if (isTRUE(bundle$affine)) paste0("model-affine.", if (format == "stan") "stan" else "R") else
                   if (product) "model-gmfrm.stan" else model_file,
                 paste0("run-", format, ".R"), if (isTRUE(bundle$affine)) "README-affine.txt" else
                   if (product) "README-gmfrm.txt" else "README.txt", "summarize.R")
  source_files <- file.path(template_dir, templates)
  if (!all(file.exists(source_files))) {
    stop("Model export templates are missing from this installation.", call. = FALSE)
  }
  filenames <- c("model-data.rds", "levels.csv", "categories.csv",
                 model_file, "run.R", "README.txt", if (!is.null(bundle$slope_levels)) "slope-levels.csv",
                 if (isTRUE(bundle$affine)) "steps.csv", "summarize.R")
  files <- stats::setNames(file.path(output_dir, filenames), filenames)
  if (any(dir.exists(files)) || (!overwrite && any(file.exists(files)))) {
    stop("Bundle files already exist; use a new directory or `overwrite = TRUE`.",
         call. = FALSE)
  }
  if (!acknowledge_sensitive) {
    warning("The model bundle includes observed ratings and original level labels. Review it before sharing.",
            call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create `output_dir`.", call. = FALSE)
  }
  saveRDS(bundle, files[["model-data.rds"]])
  utils::write.csv(bundle$levels, files[["levels.csv"]], row.names = FALSE,
                   fileEncoding = "UTF-8")
  utils::write.csv(bundle$categories, files[["categories.csv"]], row.names = FALSE)
  if (!is.null(bundle$slope_levels)) utils::write.csv(bundle$slope_levels, files[["slope-levels.csv"]],
    row.names = FALSE, fileEncoding = "UTF-8")
  if (isTRUE(bundle$affine)) utils::write.csv(bundle$step_levels, files[["steps.csv"]],
    row.names = FALSE, fileEncoding = "UTF-8")
  if (!all(file.copy(source_files, files[c(model_file, "run.R", "README.txt", "summarize.R")], overwrite = overwrite))) {
    stop("Could not copy all model export templates.", call. = FALSE)
  }
  invisible(list(format = format, model = bundle$model,
                 output_dir = normalizePath(output_dir), files = files))
}

mfrm_export_basis <- function(n, centered = TRUE) {
  if (!centered) return(diag(n))
  if (n == 1L) return(matrix(numeric(), 1L, 0L))
  h <- stats::contr.helmert(n)
  sweep(h, 2L, sqrt(colSums(h^2)), "/")
}

mfrm_model_export_data <- function(fit, location_prior_sd, step_prior_sd,
                                   log_slope_prior_sd = 1,
                                   population_prior = list(mean = c(0, 1), log_sd = c(0, 0.5))) {
  if (!inherits(fit, "mfrm_fit") ||
      inherits(fit, c("mfrm_testlet", "mfrm_random_rater"))) {
    stop("`fit` must be an ordinary mfrm_fit.", call. = FALSE)
  }
  for (x in c("location_prior_sd", "step_prior_sd", "log_slope_prior_sd")) {
    value <- get(x)
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("`", x, "` must be a finite positive number.", call. = FALSE)
    }
  }
  cfg <- fit$config
  product <- identical(cfg$model, "GPCM") && mfrm_has_product_slopes(fit)
  if (!cfg$model %in% c("RSM", "PCM", "GPCM") || !identical(cfg$method, "MML") ||
      inherits(fit, "mfrm_imported_fit")) {
    stop("Model export supports native ordinary RSM/PCM/GPCM MML specifications only.",
         call. = FALSE)
  }
  single <- identical(cfg$model, "GPCM") && !product
  if (single && (length(cfg$slope_facet) != 1L ||
      !identical(cfg$gpcm_spec$levels, cfg$facet_levels[[cfg$slope_facet]]) ||
      !identical(cfg$gpcm_spec$identification, "sum_to_zero_log_slopes") ||
      !identical(cfg$gpcm_spec$scale_reference, "geometric_mean_one") ||
      cfg$gpcm_spec$n_params != max(length(cfg$gpcm_spec$levels) - 1L, 0L))) {
    stop("One-family slope levels or constraints disagree with the fitted specification.", call. = FALSE)
  }
  if (product) {
    mfrm_gpcm_slope_roles(cfg)
    if (any(lengths(cfg$gpcm_spec$component_levels) < 2L) ||
        !identical(cfg$gpcm_spec$component_levels, cfg$facet_levels[cfg$slope_facet]) ||
        !identical(cfg$step_facet, cfg$slope_facet[2L])) stop(
      "Two-family slope levels or step ownership disagree with the fitted specification.", call. = FALSE)
  }
  pop <- cfg$population_spec
  free_population <- isTRUE(pop$active)
  if (free_population) {
    Xp <- pop$design_matrix
    if (!single || !is.null(pop$fixed_sd) || !is.matrix(Xp) ||
        !is.numeric(Xp) || ncol(Xp) != 1L || !nrow(Xp) ||
        anyNA(Xp) || any(Xp != 1) ||
        !identical(pop$design_columns, "(Intercept)")) {
      stop("Estimated-population export requires one-family GPCM with an intercept-only normal population and free SD; regression and two-family populations are unsupported.", call. = FALSE)
    }
    ids <- as.character(fit$population$person_table[[pop$person_id]])
    lookup <- pop$person_lookup
    if (length(ids) != nrow(Xp) || anyNA(ids) || anyDuplicated(ids) ||
        !is.numeric(lookup) || length(lookup) != length(cfg$person_levels) ||
        anyNA(lookup) || any(!is.finite(lookup)) || any(lookup != as.integer(lookup)) ||
        any(lookup < 1 | lookup > length(ids)) ||
        !identical(ids[lookup], as.character(cfg$person_levels))) {
      stop("Population rows do not match the fitted Person identities.", call. = FALSE)
    }
    if (!is.list(population_prior) || length(population_prior) != 2L ||
        !setequal(names(population_prior), c("mean", "log_sd")) ||
        !all(vapply(population_prior, function(x) is.numeric(x) && length(x) == 2L &&
          all(is.finite(x)) && x[2L] > 0, logical(1)))) {
      stop("`population_prior` must contain mean and log_sd, each c(location, positive scale).", call. = FALSE)
    }
  }
  if (length(cfg$interaction_specs) || length(cfg$dummy_facets) ||
      !identical(cfg$facet_shrinkage %||% "none", "none")) {
    stop("Model export does not support interactions, dummy facets or shrinkage.",
         call. = FALSE)
  }
  specs <- c(list(cfg$theta_spec), cfg$facet_specs)
  anchored <- vapply(specs, function(s) {
    any(!is.na(s$anchors)) || any(!is.na(s$groups) & nzchar(s$groups))
  }, logical(1))
  step_specs <- mfrmr_step_specs(cfg)
  if (anchored[1L]) {
    stop("MML model export does not support Person anchors.", call. = FALSE)
  }
  affine <- single || any(anchored) || any(vapply(step_specs, function(s) {
    any(!is.na(s$anchors)) || !isTRUE(s$target == 0)
  }, logical(1)))
  d <- fit$prep$data
  if (!is.data.frame(d) || !nrow(d) || anyNA(d$Weight) ||
      any(d$Weight != 1)) {
    stop("Model export requires retained observations with unit weights.", call. = FALSE)
  }
  facets <- cfg$facet_names
  maps <- lapply(facets, function(f) {
    s <- cfg$facet_specs[[f]]
    if (!identical(as.character(s$levels), as.character(cfg$facet_levels[[f]]))) {
      stop("Facet levels and constraints disagree.", call. = FALSE)
    }
    mfrm_export_affine_map(s, location_prior_sd)
  })
  basis <- lapply(maps, `[[`, "basis")
  nr <- vapply(basis, nrow, integer(1))
  nc <- vapply(basis, ncol, integer(1))
  B <- matrix(0, sum(nr), sum(nc))
  X <- matrix(0, nrow(d), sum(nc))
  location_sd <- numeric(sum(nc))
  location_offset <- numeric(sum(nr))
  eta_offset <- numeric(nrow(d))
  levels <- list(data.frame(Facet = "Person", Level = cfg$person_levels,
                            Index = seq_along(cfg$person_levels),
                            Parameter = paste0("theta[", seq_along(cfg$person_levels), "]"),
                            Sign = 1, Centered = FALSE))
  row_start <- col_start <- 0L
  for (i in seq_along(facets)) {
    f <- facets[i]
    rows <- row_start + seq_len(nr[i])
    cols <- col_start + seq_len(nc[i])
    index <- match(as.character(d[[f]]), cfg$facet_levels[[f]])
    sign <- cfg$facet_signs[[f]]
    if (anyNA(index) || length(sign) != 1L || !sign %in% c(-1, 1)) {
      stop("Invalid retained facet indices or signs.", call. = FALSE)
    }
    B[rows, cols] <- basis[[i]]
    X[, cols] <- sign * basis[[i]][index, , drop = FALSE]
    centered <- isTRUE(cfg$facet_specs[[f]]$centered)
    location_sd[cols] <- maps[[i]]$sd
    location_offset[rows] <- maps[[i]]$offset
    eta_offset <- eta_offset + sign * maps[[i]]$offset[index]
    levels[[length(levels) + 1L]] <- data.frame(
      Facet = f, Level = cfg$facet_levels[[f]], Index = rows,
      Parameter = paste0("locations[", rows, "]"), Sign = sign, Centered = centered)
    row_start <- row_start + nr[i]
    col_start <- col_start + nc[i]
  }
  C <- as.integer(cfg$n_cat)
  H <- mfrm_export_basis(C - 1L)
  owner <- if (cfg$model == "RSM") rep(1L, nrow(d)) else
    match(as.character(d[[cfg$step_facet]]), cfg$facet_levels[[cfg$step_facet]])
  person <- match(as.character(d$Person), cfg$person_levels)
  y <- as.integer(d$score_k) + 1L
  if (anyNA(person) || anyNA(owner) || anyNA(y) || any(y < 1L | y > C) ||
      any(d$score_k != y - 1L)) {
    stop("Invalid retained Person, step-owner or category indices.", call. = FALSE)
  }
  out <- list(model = cfg$model, source_method = "MML", package_version =
         as.character(utils::packageVersion("mfrmr")),
       stan_data = list(N = nrow(d), J = length(cfg$person_levels), C = C,
                        P = ncol(X), L = nrow(B), G = length(step_specs),
                        Q = ncol(H), person = as.array(person),
                        owner = as.array(owner), y = as.array(y),
                        X = X, location_basis = B, step_basis = H,
                        location_sd = as.array(location_sd),
                        step_sd = if (C > 2L) step_prior_sd * sqrt((C - 1) / (C - 2)) else step_prior_sd),
       levels = do.call(rbind, levels), step_owners = names(step_specs),
       categories = data.frame(Category = seq_len(C),
                               InternalScore = cfg$rating_min + seq_len(C) - 1L,
                               Score = cfg$score_map$OriginalScore[match(
                                 cfg$rating_min + seq_len(C) - 1L,
                                 cfg$score_map$InternalScore)]),
       priors = list(ability = "independent normal(0, 1)",
                     location_marginal_sd = location_prior_sd,
                     step_marginal_sd = step_prior_sd),
       source_columns = cfg$source_columns,
       repeated_ratings = cfg$repeated_ratings)
  if (product) {
    owners <- cfg$slope_facet
    n <- lengths(cfg$gpcm_spec$component_levels)
    first <- mfrm_export_basis(n[1L])
    slope_basis <- matrix(0, sum(n), sum(n) - 1L)
    first_cols <- seq_len(n[1L] - 1L)
    second_cols <- n[1L] - 1L + seq_len(n[2L])
    slope_basis[seq_len(n[1L]), first_cols] <- first
    slope_basis[n[1L] + seq_len(n[2L]), second_cols] <- diag(n[2L])
    ids <- lapply(owners, function(f) match(as.character(d[[f]]), cfg$facet_levels[[f]]))
    if (anyNA(unlist(ids))) stop("Invalid retained slope-owner indices.", call. = FALSE)
    slope_X <- cbind(first[ids[[1L]], , drop = FALSE], diag(n[2L])[ids[[2L]], , drop = FALSE])
    out$stan_data <- c(out$stan_data, list(S = ncol(slope_basis), A = nrow(slope_basis),
      slope_basis = slope_basis, slope_X = slope_X,
      log_slope_sd = as.array(c(rep(log_slope_prior_sd * sqrt(n[1L] / (n[1L] - 1L)), n[1L] - 1L),
        rep(log_slope_prior_sd, n[2L])))))
    out$slope_levels <- data.frame(Owner = rep(owners, n),
      Level = unlist(cfg$gpcm_spec$component_levels, use.names = FALSE),
      Index = seq_len(sum(n)), Parameter = paste0("component_slopes[", seq_len(sum(n)), "]"),
      LogParameter = paste0("log_component_slopes[", seq_len(sum(n)), "]"),
      ScaleReference = rep(c("geometric_mean_one", "fixed_standard_normal"), n))
    out$priors$log_slope_marginal_sd <- log_slope_prior_sd
  }
  if (single) {
    f <- cfg$slope_facet
    labels <- cfg$gpcm_spec$levels
    n <- length(labels)
    sb <- mfrm_export_basis(n)
    index <- match(as.character(d[[f]]), labels)
    if (anyNA(index)) stop("Invalid retained slope-owner indices.", call. = FALSE)
    out$stan_data <- c(out$stan_data, list(S = ncol(sb), A = n,
      slope_basis = sb, slope_X = sb[index, , drop = FALSE],
      log_slope_sd = as.array(rep(log_slope_prior_sd * if (n > 1L) sqrt(n / (n - 1L)) else 1, ncol(sb)))))
    out$slope_levels <- data.frame(Owner = f, Level = labels, Index = seq_len(n),
      Parameter = paste0("component_slopes[", seq_len(n), "]"),
      LogParameter = paste0("log_component_slopes[", seq_len(n), "]"),
      ScaleReference = "geometric_mean_one")
    out$priors$log_slope_marginal_sd <- log_slope_prior_sd
  }
  if (affine) {
    sm <- lapply(step_specs, mfrm_export_affine_map, prior_sd = step_prior_sd, step = TRUE)
    counts <- vapply(sm, function(x) ncol(x$basis), integer(1))
    K <- C - 1L
    step_basis <- matrix(0, length(sm) * K, sum(counts))
    start <- 0L
    for (g in seq_along(sm)) {
      step_basis[(g - 1L) * K + seq_len(K), start + seq_len(counts[g])] <- sm[[g]]$basis
      start <- start + counts[g]
    }
    out$affine <- TRUE
    out$constraints <- list(facets = cfg$facet_specs, steps = step_specs)
    out$stan_data$Q <- NULL
    out$stan_data$T <- sum(counts)
    out$stan_data$step_basis <- step_basis
    out$stan_data$step_sd <- as.array(unlist(lapply(sm, `[[`, "sd"), use.names = FALSE))
    out$stan_data$step_offset <- as.array(unlist(lapply(sm, `[[`, "offset"), use.names = FALSE))
    out$stan_data$location_offset <- as.array(location_offset)
    out$stan_data$eta_offset <- as.array(eta_offset)
    if (is.null(out$stan_data$S)) out$stan_data <- c(out$stan_data, list(S = 0L, A = 0L,
      slope_basis = matrix(0, 0, 0), slope_X = matrix(0, nrow(d), 0), log_slope_sd = array(numeric())))
    out$levels$PriorMean <- c(rep(0, out$stan_data$J), location_offset)
    out$levels$PriorSD <- c(rep(1, out$stan_data$J), sqrt(rowSums(sweep(B, 2L, location_sd, `*`)^2)))
    out$levels$Fixed <- out$levels$PriorSD == 0
    out$step_levels <- data.frame(Owner = rep(names(step_specs), each = K), Transition = rep(seq_len(K), length(sm)),
      Parameter = unlist(lapply(seq_along(sm), function(g) paste0("steps[", g, ",", seq_len(K), "]"))),
      PriorMean = as.vector(out$stan_data$step_offset),
      PriorSD = sqrt(rowSums(sweep(step_basis, 2L, out$stan_data$step_sd, `*`)^2)))
    out$step_levels$Fixed <- out$step_levels$PriorSD == 0
    out$priors$anchor_basis <- "Conditional on fixed anchors and group means; normal free contrasts about the affine prior means."
  }
  if (free_population) {
    prior <- lapply(population_prior[c("mean", "log_sd")], as.numeric)
    person_sd <- sqrt(prior$mean[2L]^2 + exp(2 * prior$log_sd[1L] + 2 * prior$log_sd[2L]^2))
    if (!is.finite(person_sd)) stop("The marginal Person prior SD exceeds the numeric range; review population_prior.", call. = FALSE)
    out$schema_version <- 2L
    out$population <- list(mode = "estimated_normal", formula = "~1",
      parameterization = "noncentered", prior = prior,
      person_levels = as.character(cfg$person_levels),
      source_identification = pop$identification_role)
    out$stan_data$population_mean_prior <- as.array(prior$mean)
    out$stan_data$population_log_sd_prior <- as.array(prior$log_sd)
    out$priors$ability <- "normal(population_mean, population_sd); conditionally independent"
    out$priors$population <- prior
    is_person <- out$levels$Facet == "Person"
    out$levels$PriorMean[is_person] <- prior$mean[1L]
    out$levels$PriorSD[is_person] <- person_sd
  }
  out
}

# Use the native affine constraints, but symmetric priors on free contrasts.
mfrm_export_affine_map <- function(spec, prior_sd, step = FALSE) {
  if (step) {
    checked <- build_step_constraint(spec$n_steps, spec$anchors[!is.na(spec$anchors)], spec$target, spec$scope)
    if (!isTRUE(all.equal(checked, spec))) stop("Inconsistent step constraints for export.", call. = FALSE)
    offset <- expand_step_constraint(numeric(spec$n_params), spec)
    blocks <- list(spec$free_index)
    centered <- TRUE
  } else {
    checked <- build_facet_constraint(spec$levels, spec$anchors[!is.na(spec$anchors)],
      spec$groups[!is.na(spec$groups)], spec$group_values, spec$centered, spec$centering_scope %||% "remainder")
    if (!isTRUE(all.equal(checked, spec))) stop("Inconsistent facet constraints for export.", call. = FALSE)
    offset <- as.numeric(expand_facet_with_constraints(numeric(spec$n_params), spec))
    free <- which(is.na(spec$anchors))
    ids <- unique(stats::na.omit(spec$groups[free]))
    ids <- ids[nzchar(ids)]
    blocks <- lapply(ids, function(g) which(is.na(spec$anchors) & !is.na(spec$groups) & spec$groups == g))
    blocks <- c(blocks, list(free[is.na(spec$groups[free]) | !nzchar(spec$groups[free])]))
    centered <- c(rep(TRUE, length(ids)), isTRUE(spec$centered))
  }
  B <- matrix(0, length(offset), spec$n_params)
  sd <- numeric(spec$n_params)
  used <- 0L
  for (i in seq_along(blocks)) {
    rows <- blocks[[i]]
    if (!length(rows)) next
    h <- mfrm_export_basis(length(rows), centered[i])
    cols <- used + seq_len(ncol(h))
    if (centered[i]) offset[rows] <- mean(offset[rows])
    B[rows, cols] <- h
    sd[cols] <- prior_sd * if (centered[i] && length(rows) > 1L) sqrt(length(rows) / (length(rows) - 1L)) else 1
    used <- used + ncol(h)
  }
  stopifnot(used == spec$n_params, all(is.finite(offset)))
  list(offset = offset, basis = B, sd = sd)
}
