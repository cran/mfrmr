# Local derivative of the actual two-stage estimator, in its affine free
# coordinates. Category k has unnormalized log probability X_k beta + O_k a.
# With Hs = sum Cov(Xs), B = sum Cov(Xs, Os), anchor map A = da / d beta_c,
# and target derivatives Gs, Ga, the event loadings are
#   Xs Hs^{-1} Gs' + Xc Hc^{-1} A' (Ga - Gs Hs^{-1} B)'.
# Map both terms onto the SAME union event before taking categorical covariance.
# The union probabilities can differ from either observed stage's fitted model:
# center the event loadings under those probabilities, not under the stage model.

mfrm_fair_delta_design <- function(fit, idx = NULL) {
  cfg <- fit$config
  if (is.null(idx)) idx <- build_indices(fit$prep, cfg$step_facet,
    cfg$slope_facet, cfg$interaction_specs)
  sizes <- build_param_sizes(cfg)
  adjacent <- mfrmr_estimability_adjacent_design(fit$prep, idx, cfg, sizes)
  p <- length(fit$opt$par)
  if (!identical(sort(adjacent$map$OptimizerIndex), seq_len(p)))
    stop("The delta design does not cover the optimizer coordinates.", call. = FALSE)
  order <- order(adjacent$map$OptimizerIndex)
  n <- length(idx$score_k)
  categories <- vector("list", cfg$n_cat)
  categories[[1]] <- mfrmr_empty_sparse_matrix(n, p)
  for (k in seq_len(cfg$n_cat - 1L)) categories[[k + 1L]] <- categories[[k]] +
    adjacent$design[seq_len(n) + (k - 1L) * n, order, drop = FALSE]
  list(category = categories, idx = idx, parameter_map = adjacent$map[order, , drop = FALSE])
}

mfrm_fair_delta_anchor_map <- function(source, anchors) {
  sizes <- build_param_sizes(source$config); slices <- build_param_slices(sizes)
  A <- mfrmr_empty_sparse_matrix(nrow(anchors), length(source$opt$par))
  for (facet in unique(anchors$Facet)) {
    rows <- which(anchors$Facet == facet)
    spec <- if (facet == "Person") source$config$theta_spec else source$config$facet_specs[[facet]]
    cols <- slices[[if (facet == "Person") "theta" else facet]]
    jac <- mfrmr_constraint_jacobian_sparse(spec, facet)$jacobian
    if (length(cols)) A[rows, cols] <- jac[match(anchors$Level[rows], spec$levels), , drop = FALSE]
  }
  A
}

mfrm_fair_delta_anchor_design <- function(fit, idx, anchors, anchor_map) {
  n <- length(idx$score_k); cfg <- fit$config
  eta <- mfrmr_empty_sparse_matrix(n, nrow(anchors))
  for (facet in unique(anchors$Facet)) {
    rows <- which(anchors$Facet == facet)
    spec <- if (facet == "Person") cfg$theta_spec else cfg$facet_specs[[facet]]
    selected <- match(anchors$Level[rows], spec$levels)
    if (anyNA(selected) || any(!is.finite(spec$anchors[selected])))
      stop("An anchor is absent from the scoring constraint map.", call. = FALSE)
    if (any(!is.na(spec$groups) & nzchar(spec$groups)))
      stop("Delta anchor propagation does not support group constraints.", call. = FALSE)
    C <- Matrix::sparseMatrix(i = selected, j = rows, x = 1,
      dims = c(length(spec$levels), nrow(anchors)))
    free <- which(is.na(spec$anchors))
    if (isTRUE(spec$centered) && identical(spec$centering_scope, "facet")) {
      if (length(free)) {
        # The final free level absorbs the changing sum of all direct anchors.
        C[tail(free, 1L), rows] <- -1
      } else if (length(anchor_map) &&
                 any(abs(Matrix::colSums(anchor_map[rows, , drop = FALSE])) > 1e-12)) {
        stop("Estimated anchors do not preserve the fully anchored facet's centering constraint.", call. = FALSE)
      }
    }
    lookup <- if (facet == "Person") idx$person else idx$facets[[facet]]
    sign <- if (facet == "Person") 1 else cfg$facet_signs[[facet]]
    eta <- eta + sign * C[lookup, , drop = FALSE]
  }
  lapply(seq_len(cfg$n_cat) - 1L, function(k) k * eta)
}

mfrm_fair_delta_center <- function(category, probs) {
  mean <- category[[1]] * 0
  for (k in seq_along(category)) mean <- mean + Matrix::Diagonal(x = probs[, k]) %*% category[[k]]
  lapply(category, function(x) x - mean)
}

mfrm_fair_delta_information <- function(left, probs, right = left) {
  out <- mfrmr_empty_sparse_matrix(ncol(left[[1]]), ncol(right[[1]]))
  for (k in seq_along(left)) out <- out + Matrix::crossprod(left[[k]],
    Matrix::Diagonal(x = probs[, k]) %*% right[[k]])
  out
}

# Solve for selected target loadings; do not form the full parameter inverse.
# Diagonal scaling improves numerical units without a ridge or rank truncation.
mfrm_fair_delta_solve <- function(information, rhs, stage) {
  p <- nrow(information)
  if (!p) return(list(solution = matrix(0, 0L, ncol(rhs)), residual = 0))
  diagonal <- as.numeric(Matrix::diag(information))
  if (any(!is.finite(diagonal)) || any(diagonal <= 0))
    stop(stage, " information has a nonpositive/nonfinite diagonal.", call. = FALSE)
  D <- Matrix::Diagonal(x = 1 / sqrt(diagonal))
  H <- Matrix::forceSymmetric(D %*% information %*% D)
  factor <- tryCatch(Matrix::Cholesky(H, LDL = FALSE, perm = TRUE),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(factor)) stop(stage, " information is not numerically positive definite; no regularization was applied.", call. = FALSE)
  scaled_rhs <- D %*% rhs
  z <- Matrix::solve(factor, scaled_rhs, system = "A")
  residual <- max(abs(as.matrix(H %*% z - scaled_rhs))) / max(1, max(abs(as.matrix(scaled_rhs))))
  if (!is.finite(residual) || residual > 1e-8 || any(!is.finite(z)))
    stop(stage, " information solve failed its residual check.", call. = FALSE)
  list(solution = as.matrix(D %*% z), residual = residual)
}

mfrm_fair_delta_target <- function(fit, plan, parm, anchors, anchor_map) {
  cfg <- fit$config; m <- nrow(plan)
  rows <- rep(seq_len(m), times = length(parm))
  idx <- list(score_k = rep(0L, length(rows)),
    person = rep(match(parm, cfg$person_levels), each = m),
    facets = lapply(cfg$facet_names, function(f) match(plan[[f]][rows], cfg$facet_levels[[f]])),
    interactions = list())
  names(idx$facets) <- cfg$facet_names
  idx$step_idx <- if (cfg$model == "RSM") rep(1L, length(rows)) else idx$facets[[cfg$step_facet]]
  design <- mfrm_fair_delta_design(fit, idx)$category
  kernel <- mfrm_fair_reference_kernel(fit, plan)
  theta <- expand_params(fit$opt$par, build_param_sizes(cfg), cfg)$theta[idx$person]
  eta <- theta + kernel$offsets[rows]
  probs <- category_prob_pcm(eta, kernel$cumulative, rows)
  weighted_derivative <- function(category) {
    center <- mfrm_fair_delta_center(category, probs)
    derivative <- center[[1]] * 0
    for (k in seq_along(center)) derivative <- derivative +
      Matrix::Diagonal(x = (k - 1L) * probs[, k]) %*% center[[k]]
    aggregate <- Matrix::sparseMatrix(i = rep(seq_along(parm), each = m),
      j = seq_along(rows), x = rep(plan$NormalizedWeight, times = length(parm)),
      dims = c(length(parm), length(rows)))
    as.matrix(aggregate %*% derivative)
  }
  list(free = weighted_derivative(design), anchors = weighted_derivative(
    mfrm_fair_delta_anchor_design(fit, idx, anchors, anchor_map)))
}

mfrm_fair_delta_covariance <- function(fit, calibration, plan, parm, L, anchors, maps, probs) {
  two_stage <- !is.null(calibration)
  for (f in if (two_stage) list(fit, calibration) else list(fit)) {
    theta <- f$facets$person$PrimaryEstimate %||% f$facets$person$Estimate
    if (any(!is.finite(theta)))
      stop("Delta propagation requires finite primary estimates for every Person in both stages; boundary Persons need the refit route.", call. = FALSE)
  }
  propagated <- if (two_stage) anchors else anchors[FALSE, , drop = FALSE]
  A <- if (two_stage) mfrm_fair_delta_anchor_map(calibration, anchors) else
    mfrmr_empty_sparse_matrix(0L, 0L)
  scoring <- mfrm_fair_delta_design(fit)
  ps <- mfrm_fair_generator_probs(fit)
  Ds <- mfrm_fair_delta_center(scoring$category, ps)
  Hs <- mfrm_fair_delta_information(Ds, ps)
  penalized <- mfrm_is_jeffreys_fit(fit)
  penalized_s <- penalized_c <- NULL
  if (penalized) {
    penalized_s <- mfrm_fair_jeffreys_jacobian(fit, scoring, ps, propagated, A)
    Hs <- penalized_s$information
  }
  target <- mfrm_fair_delta_target(fit, plan, parm, propagated, A)
  Gs <- L %*% target$free; Ga <- L %*% target$anchors
  solved_s <- mfrm_fair_delta_solve(Hs, t(Gs), "Scoring")
  loading_s <- solved_s$solution
  loading_c <- matrix(0, 0L, nrow(L)); solved_c <- list(residual = 0)
  calibration_design <- NULL; coupling <- matrix(0, nrow(L), 0L)
  if (two_stage) {
    O <- mfrm_fair_delta_anchor_design(fit, scoring$idx, anchors, A)
    B <- mfrm_fair_delta_information(Ds, ps, mfrm_fair_delta_center(O, ps))
    if (penalized) B <- penalized_s$anchors
    coupling <- as.matrix(Ga - t(loading_s) %*% B)
    calibration_design <- mfrm_fair_delta_design(calibration)
    pc <- mfrm_fair_generator_probs(calibration)
    Dc <- mfrm_fair_delta_center(calibration_design$category, pc)
    Hc <- mfrm_fair_delta_information(Dc, pc)
    if (penalized) {
      penalized_c <- mfrm_fair_jeffreys_jacobian(calibration, calibration_design, pc,
        anchors[FALSE, , drop = FALSE], mfrmr_empty_sparse_matrix(0L, 0L))
      Hc <- penalized_c$information
    }
    solved_c <- mfrm_fair_delta_solve(Hc, Matrix::t(A) %*% t(coupling), "Calibration")
    loading_c <- solved_c$solution
  }
  n <- nrow(probs); q <- nrow(L)
  scoring_events <- calibration_events <- vector("list", ncol(probs))
  for (k in seq_len(ncol(probs))) {
    s <- c <- matrix(0, n, q)
    s[maps$scoring, ] <- as.matrix(scoring$category[[k]] %*% loading_s)
    if (two_stage) c[maps$calibration, ] <- as.matrix(calibration_design$category[[k]] %*% loading_c)
    scoring_events[[k]] <- s; calibration_events[[k]] <- c
  }
  center <- function(events) {
    mean <- matrix(0, n, q)
    for (k in seq_len(ncol(probs))) mean <- mean + probs[, k] * events[[k]]
    lapply(events, function(x) x - mean)
  }
  scoring_events <- center(scoring_events); calibration_events <- center(calibration_events)
  total <- direct <- through <- shared <- matrix(0, q, q, dimnames = list(rownames(L), rownames(L)))
  for (k in seq_len(ncol(probs))) {
    s <- sqrt(probs[, k]) * scoring_events[[k]]
    c <- sqrt(probs[, k]) * calibration_events[[k]]
    total <- total + crossprod(s + c)
    direct <- direct + crossprod(s); through <- through + crossprod(c)
    cross <- crossprod(s, c); shared <- shared + cross + t(cross)
  }
  if (any(!is.finite(total))) stop("Delta target covariance is nonfinite.", call. = FALSE)
  list(covariance = total,
    components = list(direct_scoring = direct, through_calibration = through, shared_events = shared),
    linearization = list(target_free = Gs, target_anchors = Ga, anchor_map = A,
      anchor_coupling = as.matrix(coupling), scoring_loading = loading_s,
      calibration_loading = loading_c, scoring_coordinates = scoring$parameter_map,
      calibration_coordinates = calibration_design$parameter_map,
      score_jacobian = if (penalized) list(scoring = penalized_s$diagnostics,
        calibration = penalized_c$diagnostics) else list(equation = "ordinary_jml_score"),
      solve_residual = c(scoring = solved_s$residual, calibration = solved_c$residual)))
}

mfrm_fair_delta_result <- function(fit, calibration, generator, plan, parm, L,
    estimate, truth, anchors, maps, probs, level, gradient_tol, anchor_uncertainty,
    sources, source_checks, has_contrasts) {
  calculated <- tryCatch({
    if (is.null(probs))
      stop("Delta propagation requires finite primary estimates for every Person in the generating fit.", call. = FALSE)
    mfrm_fair_delta_covariance(fit, calibration, plan, parm, L, anchors, maps, probs)
  }, error = function(e) list(reason = conditionMessage(e)))
  ok <- is.null(calculated$reason)
  covariance <- if (ok) calculated$covariance else matrix(NA_real_, nrow(L), nrow(L),
    dimnames = list(rownames(L), rownames(L)))
  se <- sqrt(diag(covariance))
  z <- stats::qnorm((1 + level) / 2)
  two_stage <- !is.null(calibration)
  out <- data.frame(Target = rownames(L), Estimate = unname(estimate), SE = unname(se),
    CI_Lower = unname(estimate - z * se), CI_Upper = unname(estimate + z * se),
    Level = level, Available = ok, Status = if (ok) "experimental" else "unavailable_information",
    Planned = NA_integer_, Resolved = NA_integer_, Unresolved = NA_integer_,
    ExpectedTailDraws = NA_real_, CIEligible = FALSE,
    Uncertainty = if (two_stage) "calibration_and_scoring_delta" else if (nrow(anchors))
      "scoring_delta_conditional_on_fixed_anchors" else "joint_calibration_and_person_delta",
    InferenceReview = "Experimental pointwise local model-based normal approximation; finite-sample coverage and JML bias remain unqualified.",
    Reason = if (ok) "available" else calculated$reason, row.names = NULL)
  attr(out, "covariance") <- covariance
  attr(out, "covariance_components") <- calculated$components
  attr(out, "linearization") <- calculated$linearization
  attr(out, "generating_target") <- setNames(truth, rownames(L))
  attr(out, "reference_plan") <- plan; attr(out, "response_map") <- maps
  attr(out, "source_checks") <- source_checks
  attr(out, "source_identity") <- lapply(sources, mfrm_fair_fit_identity)
  attr(out, "settings") <- list(parm = parm, contrasts = L, level = level,
    gradient_tol = gradient_tol, anchor_uncertainty = anchor_uncertainty, anchors = anchors,
    method = "local_fixed_person_delta", linearization_at = "observed_stage_fits",
    estimator = if (mfrm_is_jeffreys_fit(fit)) "jeffreys_canonical_joint" else "ordinary_jml",
    event_covariance = "generating_fit_category_probabilities",
    target = if (has_contrasts) "common_plan_expected_average_contrast" else "common_plan_expected_average",
    shared_responses = if (two_stage) length(intersect(maps$calibration, maps$scoring)) else 0L)
  if (mfrm_is_jeffreys_fit(fit)) out$InferenceReview <- paste(out$InferenceReview,
    "Analytic Jeffreys score and anchor derivatives; categorical event covariance, not inverse penalized curvature. Normal-reference limits can substantially under-cover in small samples; see help. Score bias and coverage are unqualified.")
  class(out) <- c("mfrm_fair_intervals", "data.frame")
  out
}
