# Explicit-order root solving for the internal shared-owner JML equation.
# This is not an mfrm_fit: ordinary likelihood audits/information cannot be
# applied to an adjusted-equation root. Public output integration is separate.
mfrm_jml_adjustment_fit <- function(problem, order, starts = NULL, maxit = 150L,
    sampling = c("fixed_rosters", "random_rosters")) {
  sampling <- match.arg(sampling)
  if (!is.numeric(order) || is.complex(order) || length(order) != 1L || !is.finite(order) ||
      order < 1 || order > .Machine$integer.max || order != floor(order))
    stop("Supply an explicit positive integer correction order; automatic selection is not available.")
  if (!is.numeric(maxit) || is.complex(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
      maxit < 1 || maxit > .Machine$integer.max || maxit != floor(maxit))
    stop("maxit must be a positive integer.")
  if (!requireNamespace("nleqslv", quietly = TRUE))
    stop("Install the optional package 'nleqslv' to solve the adjusted JML equations.")
  metadata <- problem$parameters; p <- nrow(metadata)
  if (is.null(p) || p < 1L || !is.function(problem$mean_score) ||
      !is.function(problem$evaluate)) stop("Supply an observed-data JML adjustment problem.")
  if (is.null(starts)) {
    neutral <- numeric(p)
    step <- which(metadata$Type == "step")
    if (length(step)) {
      centered_steps <- seq(-.5, .5, length.out = problem$specification$rating_max)
      neutral[step] <- centered_steps[metadata$Step[step]]
    }
    starts <- list(neutral = neutral,
      perturbed = neutral + rep(c(-.15, .15), length.out = p))
  }
  if (!is.list(starts) || length(starts) < 2L ||
      any(!vapply(starts, function(x) is.numeric(x) && !is.complex(x) &&
        length(x) == p && all(is.finite(x)), logical(1))) ||
      anyDuplicated(lapply(starts, as.numeric)))
    stop("Supply at least two distinct, finite starting vectors in the declared structural coordinates.")
  if (is.null(names(starts))) names(starts) <- paste0("start_", seq_along(starts))
  if (anyNA(names(starts)) || any(!nzchar(names(starts))) || anyDuplicated(names(starts)))
    stop("Starting vectors must have distinct nonempty names.")
  policy <- list(version = "explicit_order_v1", solver = "nleqslv",
    solver_version = as.character(utils::packageVersion("nleqslv")),
    maxit = as.integer(maxit), ftol = 1e-10, xtol = 1e-11,
    equation_tolerance = 1e-7, root_step_tolerance = 1e-5,
    root_agreement_tolerance = 1e-6,
    jacobian_step = 5e-5, starts = starts,
    stages = data.frame(method = c("Broyden", "Newton", "Newton"),
      global = c("cline", "dbldog", "dbldog"), stepmax = c(1, 1, .25)))
  fn <- function(b) problem$mean_score(b, order)
  attempts <- list(); roots <- integer(0)
  for (s in seq_along(starts)) for (stage in seq_len(nrow(policy$stages))) {
    setting <- policy$stages[stage, ]
    warnings <- character(0)
    attempt <- list(start = names(starts)[s], stage = stage, settings = setting,
      solver = NULL, beta = NULL, residual = NA_real_, jacobian = NULL,
      singular_values = numeric(0), newton_step = NA_real_,
      available = FALSE, reason = "solver_failed")
    attempt <- tryCatch(withCallingHandlers({
      z <- nleqslv::nleqslv(starts[[s]], fn, method = setting$method,
        global = setting$global, control = list(maxit = policy$maxit,
          ftol = policy$ftol, xtol = policy$xtol, stepmax = setting$stepmax,
          allowSingular = FALSE))
      attempt$solver <- z; attempt$beta <- z$x
      attempt$residual <- max(abs(fn(z$x)))
      attempt$reason <- "equation_not_solved"
      if (is.finite(attempt$residual) && attempt$residual <= policy$equation_tolerance) {
        derivative <- mfrmr_numeric_transformation_jacobian(fn, z$x,
          relative_step = policy$jacobian_step)
        attempt$jacobian <- derivative$jacobian
        attempt$reason <- "equation_jacobian_unresolved"
        if (derivative$valid) {
          d <- svd(derivative$jacobian, nu = 0L, nv = 0L)$d
          attempt$singular_values <- d
          # Numerical rank is a point-root check. Covariance has its separate
          # derivative-stability, conditioning and empirical-score checks.
          attempt$available <- length(d) == p &&
            min(d) > max(d) * p * .Machine$double.eps
          attempt$reason <- if (attempt$available) "local_root" else "equation_rank_deficient"
          if (attempt$available) {
            attempt$newton_step <- max(abs(solve(derivative$jacobian, fn(z$x))))
            attempt$available <- is.finite(attempt$newton_step) &&
              attempt$newton_step <= policy$root_step_tolerance
            if (!attempt$available) attempt$reason <- "root_step_unresolved"
          }
        }
      }
      attempt
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      attempt$available <- FALSE
      attempt$reason <- if (is.null(attempt$solver)) "solver_failed" else "equation_review_failed"
      attempt$error <- conditionMessage(e); attempt
    })
    attempt$warnings <- warnings
    attempts[[length(attempts) + 1L]] <- attempt
    # Covariance cannot trigger a refit or change which root is retained.
    if (isTRUE(attempt$available)) {
      roots <- c(roots, length(attempts)); break
    }
  }
  point <- list(available = FALSE, status = "no_local_root", beta = NULL,
    theta = NULL, raw_profile_criterion = NA_real_, equation_residual = NA_real_,
    selected_attempt = NA_integer_, successful_starts = length(roots),
    requested_starts = length(starts), root_difference = NA_real_)
  covariance <- list(available = FALSE, status = "point_unavailable", reason = "",
    result = NULL)
  if (length(roots)) {
    b <- attempts[[roots[1L]]]$beta
    spread <- max(vapply(roots, function(i) max(abs(attempts[[i]]$beta - b)), numeric(1)))
    point$root_difference <- if (length(roots) > 1L) spread else NA_real_
    if (spread > policy$root_agreement_tolerance) {
      # Neither a smaller SE nor a larger raw likelihood chooses among roots.
      point$status <- "multiple_roots"
    } else {
      point$available <- TRUE
      point$status <- if (length(roots) == length(starts)) "consistent_roots" else "root_with_unresolved_starts"
      point$beta <- b; point$selected_attempt <- roots[1L]
      z <- problem$evaluate(b, order)
      point$theta <- z$theta; point$raw_profile_criterion <- sum(z$q)
      point$equation_residual <- max(abs(colMeans(z$value)))
      covariance <- tryCatch(list(available = TRUE, status = "local_root_covariance",
        reason = "", result = mfrm_jml_adjustment_covariance(problem, b, order, sampling)),
        error = function(e) list(available = FALSE, status = "covariance_unavailable",
          reason = conditionMessage(e), result = NULL))
    }
  }
  list(estimator = list(name = "finite_MLE_plugin_profile_score_adjustment",
      order = as.integer(order), order_selection = "explicit", sampling = sampling,
      formal_structural_intervals = FALSE,
      target = "Adjusted-equation root; residual structural bias may remain."),
    specification = problem$specification, parameters = metadata,
    persons = problem$persons, levels = problem$levels,
    policy = policy, point = point, covariance = covariance, attempts = attempts)
}

# Reuse the ordinary parameter-table builders, with explicit root uncertainty.
# No mfrm_fit class, likelihood criteria, or Person SE is manufactured here.
mfrm_jml_adjustment_tables <- function(result) {
  if (!isTRUE(result$point$available))
    stop("Parameter tables require an available, unambiguous adjusted-equation root.")
  spec <- result$specification; levels <- result$levels[spec$facets]
  ng <- length(levels[[spec$owner]]); nstep <- spec$rating_max
  sizes <- c(lengths(levels), rep(nstep, ng), ng)
  p <- length(result$point$beta)
  if (8 * sum(sizes) * p > getOption("mfrmr.max_information_bytes", 256 * 1024^2))
    stop("Parameter transformation exceeds mfrmr.max_information_bytes.")
  H <- as.matrix(Matrix::bdiag(lapply(sizes, sum_zero_jacobian)))
  expanded <- drop(H %*% result$point$beta)
  root_se <- if (isTRUE(result$covariance$available)) {
    sqrt(pmax(0, rowSums((H %*% result$covariance$result$vcov) * H)))
  } else rep(NA_real_, length(expanded))
  ends <- cumsum(sizes); begin <- c(1L, head(ends, -1L) + 1L)
  values <- lapply(seq_along(sizes), function(j) expanded[begin[j]:ends[j]])
  nf <- length(levels)
  params <- list(facets = setNames(values[seq_len(nf)], names(levels)),
    steps_mat = do.call(rbind, values[nf + seq_len(ng)]),
    log_slopes = tail(values, 1L)[[1L]], slopes = exp(tail(values, 1L)[[1L]]))
  config <- list(model = "GPCM", facet_names = spec$facets, n_cat = nstep + 1L,
    step_facet = spec$owner, slope_facet = spec$owner)
  prep <- list(levels = levels)
  tables <- list(locations = build_other_facet_table(config, prep, params),
    steps = build_step_table(config, prep, params),
    slopes = build_slope_table(config, prep, params))
  offset <- 0L
  for (name in names(tables)) {
    rows <- offset + seq_len(nrow(tables[[name]])); offset <- max(rows)
    if (name == "slopes") {
      tables[[name]]$LogRootSE <- root_se[rows]
      tables[[name]]$RootSE <- tables[[name]]$Estimate * root_se[rows]
    } else tables[[name]]$RootSE <- root_se[rows]
    tables[[name]]$Estimator <- result$estimator$name
    tables[[name]]$CorrectionOrder <- result$estimator$order
    tables[[name]]$PointStatus <- result$point$status
    tables[[name]]$CovarianceStatus <- result$covariance$status
    tables[[name]]$CovarianceReason <- result$covariance$reason
    tables[[name]]$UncertaintyTarget <- result$estimator$target
    tables[[name]]$CIEligible <- FALSE
  }
  person <- tibble::tibble(Person = result$persons, Estimate = result$point$theta,
    SE = NA_real_, Extreme = ifelse(is.infinite(result$point$theta),
      ifelse(result$point$theta < 0, "low", "high"), "none"),
    Estimator = result$estimator$name, CorrectionOrder = result$estimator$order,
    PointStatus = result$point$status,
    EstimateTarget = "Person profile at adjusted structural estimates")
  list(facets = list(person = person, others = tables$locations),
    steps = tables$steps, slopes = tables$slopes, estimator = result$estimator)
}
