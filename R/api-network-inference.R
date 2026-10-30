# Pairwise alpha retains negative values and distinguishes undefined agreement.
.mfrm_network_alpha <- function(wide, pairs, level, categories) {
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- wide[[pairs$Rater1[i]]]; b <- wide[[pairs$Rater2[i]]]
    keep <- is.finite(a) & is.finite(b)
    if (!any(keep)) return(data.frame(KrippendorffAlpha = NA_real_,
      AlphaStatus = "no_pairable_units"))
    dat <- data.frame(Rater = rep(c("A", "B"), each = sum(keep)),
      Unit = rep(seq_len(sum(keep)), 2), Score = c(a[keep], b[keep]))
    z <- krippendorff_alpha_table(dat, "Rater", "Unit", "Score",
      level = level, categories = categories)$table
    data.frame(KrippendorffAlpha = z$Alpha, AlphaStatus = z$Status)
  })
  do.call(rbind, rows)
}

# Equations 1--3 use UNFILTERED eligible magnitude comparisons and c = 0.
.mfrm_network_fit_indices <- function(pairs, raters) {
  n <- length(raters)
  w <- t <- matrix(0, n, n, dimnames = list(raters, raters))
  for (k in which(pairs$EligiblePair)) {
    i <- pairs$Rater1[k]; j <- pairs$Rater2[k]
    t[i, j] <- t[j, i] <- pairs$DirectionN[k]
    w[i, j] <- pairs$Rater1HigherSum[k]
    w[j, i] <- pairs$Rater2HigherSum[k]
  }
  prop <- w / t; prop[t == 0] <- 0
  outgoing <- rowSums(prop); incoming <- colSums(prop)
  ratio <- outgoing / incoming
  ratio[outgoing + incoming == 0] <- NA_real_
  rows <- lapply(seq_len(n), function(i) {
    js <- which(t[i, ] > 0)
    valid <- length(js) > 0 && is.finite(ratio[i]) && ratio[i] > 0 &&
      all(is.finite(ratio[js]) & ratio[js] > 0)
    sws <- if (valid) sum(w[i, js] * ratio[js] / ratio[i] +
      w[js, i] * ratio[i] / ratio[js]) / sum(t[i, js]) - 1 else NA_real_
    data.frame(Rater = raters[i], Partners = length(js),
      Comparisons = sum(t[i, js]), OutStrength = outgoing[i],
      InStrength = incoming[i], SeverityIndex = -log(ratio[i]),
      SWS = sws,
      PUS = if (length(js)) mean(w[i, js] * w[js, i] / t[i, js]^2) else NA_real_,
      SWSStatus = if (!length(js)) "no_comparisons" else
        if (valid) "available" else "undefined_strength_ratio",
      PUSStatus = if (length(js)) "available" else "no_comparisons")
  })
  do.call(rbind, rows)
}

.mfrm_halo_contrast <- function(wide_info, pairs, method, min_pair_n,
                                bootstrap, level, seed, reference) {
  use <- is.finite(pairs$Estimate) & pairs$N >= min_pair_n
  if (reference == "different_criteria") use <- use & pairs$Criterion1 != pairs$Criterion2
  fixed <- pairs[use, , drop = FALSE]
  halo <- fixed$EdgeType == "halo"
  target_ok <- any(halo) && any(!halo)
  estimate <- if (target_ok) mean(fixed$Estimate[halo]) - mean(fixed$Estimate[!halo]) else NA_real_
  persons <- wide_info$persons
  person_known <- !is.null(persons) && length(persons) == nrow(wide_info$wide) &&
    !anyNA(persons) && all(nzchar(as.character(persons)))
  person_ok <- person_known && length(unique(persons)) >= 2L
  draws <- data.frame(Replicate = integer(), Estimate = numeric(), Status = character())
  status <- if (!target_ok) "no_comparable_edges" else if (!bootstrap) "not_requested" else
    if (!person_ok) "person_clusters_unavailable" else "available"
  if (status == "available") {
    blocks <- split(seq_len(nrow(wide_info$wide)), persons)
    draws <- with_preserved_rng_seed(seed, {
      do.call(rbind, lapply(seq_len(bootstrap), function(b) {
        # Repeated blocks remain repeated rows; never collapse resampled persons.
        idx <- unlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)], use.names = FALSE)
        z <- vapply(seq_len(nrow(fixed)), function(k) {
          a <- wide_info$wide[[fixed$From[k]]][idx]
          c <- wide_info$wide[[fixed$To[k]]][idx]
          ok <- is.finite(a) & is.finite(c)
          if (sum(ok) < min_pair_n || length(unique(a[ok])) < 2L || length(unique(c[ok])) < 2L)
            return(NA_real_)
          stats::cor(a[ok], c[ok], method = method)
        }, numeric(1))
        available <- all(is.finite(z))
        data.frame(Replicate = b,
          Estimate = if (available) mean(z[halo]) - mean(z[!halo]) else NA_real_,
          Status = if (available) "available" else "undefined_fixed_pair")
      }))
    })
    if (any(!is.finite(draws$Estimate))) status <- "incomplete_bootstrap"
  }
  ci <- if (status == "available") stats::quantile(draws$Estimate,
    c((1 - level) / 2, (1 + level) / 2), names = FALSE) else c(NA_real_, NA_real_)
  list(table = data.frame(Estimate = estimate, Lower = ci[1], Upper = ci[2],
    Level = level, Status = status, Method = "Person-cluster percentile bootstrap",
    Reference = reference, HaloPairs = sum(halo), ReferencePairs = sum(!halo),
    Persons = if (person_known) length(unique(persons)) else NA_integer_,
    Requested = bootstrap, Attempted = nrow(draws), Valid = sum(is.finite(draws$Estimate)),
    Failed = sum(!is.finite(draws$Estimate))), draws = draws, pairs = fixed)
}

#' Fit an ERGM to a declared binary rater network
#'
#' Model network configurations and rater covariates after examining shared
#' rating coverage with [rater_network_analysis()] or
#' [rater_halo_network_analysis()]. Estimation uses the optional `ergm` package.
#' @param x For `fit_rater_network()`: A rater or halo network result. All eligible pair metrics are used,
#'   before graph/plot filtering; missing comparisons never become zero ties.
#'
#' For `plot()`: Output from [fit_rater_network()].
#' @param terms One-sided ERGM formula, for example `~ edges + mutual` for a
#'   directed network, or `~ edges + nodecov("Experience")`. See
#'   `ergm::ergmTerm` for terms appropriate to the declared graph.
#' @param threshold Required finite cutoff for the chosen pair metric.
#' @param relation `"greater"` makes a tie when the metric is strictly greater
#'   than the threshold; `"less"` makes one when strictly smaller. Equal values
#'   are observed non-ties. Negative alpha is not automatically a positive tie.
#' @param node_data Optional data frame with `Node` IDs and complete node
#'   covariates. IDs must match all nodes exactly; order is matched by ID.
#' @param missing `"error"` requires all dyads to be observed. `"model"`
#'   explicitly permits unobserved dyads, represented by `NA`, using ergm's
#'   observed-data likelihood under an ignorable missingness assumption.
#' @param control Optional `ergm::control.ergm()` object.
#' @param gof Optional one-sided formula for simulation-based model checks,
#'   e.g. `~ degree + distance` for an undirected graph. See `ergm::gof.ergm`.
#' @param gof_control Optional `ergm::control.gof.ergm()` object.
#' @param seed Optional integer seed; the caller's RNG state is restored.
#' @details This is a binary ERGM, not a valued/signed ERGM or a rating-response model.
#' A coefficient is a conditional network log-odds coefficient for its change
#' statistic; it is not an MFRM severity, a causal rater effect or an independent
#' rating-level regression. Thresholds, sampling of raters, overlap and node
#' covariates define the question. Prespecify them and assess sensitivity by
#' refitting at substantively plausible thresholds; do not select a threshold
#' for significance. Different thresholds define different response networks:
#' their AIC values do not select a best threshold.
#'
#' A missing dyad includes insufficient overlap and an undefined association.
#' Such missingness need not be ignorable (e.g. constant scores or assignment
#' related to rater characteristics). `missing = "model"` asserts an assumption,
#' not a correction for informative missingness. Self-ties are excluded.
#' Unobserved dyads are not structural zeroes. For structural constraints or
#' valued ERGMs, use the returned `network` with `ergm::ergm()` explicitly.
#'
#' Estimation requests MLE. For dyad-independent terms MPLE is also the MLE;
#' dependent terms require MCMC. A returned fit is not proof of convergence or
#' model adequacy. Inspect `ergm::mcmc.diagnostics(result$fit)` and the saved
#' `gof` object; review warnings, degeneracy, effective sample size and mixing.
#' GOF tail areas are simulation diagnostics, not tests of halo. Intervals in
#' the coefficient table are approximate pointwise normal intervals, conditional
#' on the chosen network model and adequate estimation. Small rater networks
#' may not support complex models. Neither mfrmr nor this wrapper validates
#' arbitrary ERGM specifications or corrects estimated-edge uncertainty.
#' @return An `mfrm_network_ergm` object with the backend `fit`, binary `network`,
#'   complete `dyads`, `nodes`, coefficient/diagnostic/settings tables, backend
#'   warnings and optional `gof` with long-form simulation tables.
#'   `summary()` retains model settings and diagnostic cautions. `coef()` and
#'   `vcov()` extract backend coefficient coordinates; `confint()` returns
#'   approximate pointwise normal intervals without refitting. `plot()` uses
#'   the originally stored 95% intervals or saved simulation GOF.
#'   Use [mfrm_report()] and [export_summary_appendix()] for stored reporting.
#' @references Lamprianou, I. (2018). Investigation of rater effects using social network
#' analysis and exponential random graph models. Educational and Psychological
#' Measurement, 78, 430--459. doi:10.1177/0013164416689696.
#' @seealso [rater_network_analysis()], [rater_halo_network_analysis()],
#'   [plot.mfrm_network_ergm()]
#' @examples # Optional integration example (loads additional modelling or rendering packages).
#' # example("fit_rater_network", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("ergm", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   ratings <- load_mfrmr_data("example_core")
#'   fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'     method = "JML", maxit = 300)
#'   rn <- rater_network_analysis(fit, weight_metric = "KrippendorffAlpha")
#'   # Four synthetic raters: illustrate the handoff, not a complex ERGM.
#'   # Tie = pairwise interval alpha strictly above the declared 0.35 cutoff.
#'   model <- fit_rater_network(rn, ~ edges, threshold = 0.35, seed = 41)
#'   model$dyads
#'   summary(model)$diagnostics
#'   confint(model, level = 0.90) # Approximate, conditional on this binary network
#'   mfrm_report(model)$tables$coefficients
#'   plot(model, draw = FALSE)
#' }
#' }
#' @export
#' @order 1
fit_rater_network <- function(x, terms, threshold, relation = c("greater", "less"),
                              node_data = NULL, missing = c("error", "model"),
                              control = NULL, gof = NULL, gof_control = NULL, seed = NULL) {
  if (!inherits(x, c("mfrm_rater_network", "mfrm_halo_network")))
    stop("`x` must be a rater or halo network result.", call. = FALSE)
  if (!requireNamespace("ergm", quietly = TRUE) || !requireNamespace("network", quietly = TRUE))
    stop("Install the optional ergm and network packages.", call. = FALSE)
  if (!inherits(terms, "formula") || length(terms) != 2L)
    stop("`terms` must be a one-sided ERGM formula.", call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("Supply one finite `threshold`.", call. = FALSE)
  relation <- match.arg(relation); missing <- match.arg(missing)
  halo <- inherits(x, "mfrm_halo_network")
  directed <- !halo && isTRUE(x$settings$directed)
  ids <- as.character(x$node_metrics[[if (halo) "Node" else "Rater"]])
  if (length(ids) < 3L || anyNA(ids) || anyDuplicated(ids))
    stop("At least three unique network node IDs are required.", call. = FALSE)
  p <- x$pair_metrics
  if (halo) {
    dyads <- data.frame(From = p$From, To = p$To, Metric = p$Estimate,
      N = p$N, Eligible = p$N >= x$settings$min_pair_n)
  } else if (directed) {
    suffix <- if (identical(x$settings$weight_metric, "DirectionalMagnitude")) "Mean" else "Prop"
    dyads <- rbind(data.frame(From = p$Rater1, To = p$Rater2,
      Metric = p[[paste0("Rater1Higher", suffix)]], N = p$DirectionN, Eligible = p$EligiblePair),
      data.frame(From = p$Rater2, To = p$Rater1,
      Metric = p[[paste0("Rater2Higher", suffix)]], N = p$DirectionN, Eligible = p$EligiblePair))
  } else {
    dyads <- data.frame(From = p$Rater1, To = p$Rater2,
      Metric = p[[x$settings$weight_metric]], N = p$N, Eligible = p$EligiblePair)
  }
  adj <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  for (k in seq_len(nrow(dyads))) {
    value <- if (isTRUE(dyads$Eligible[k]) && is.finite(dyads$Metric[k]))
      as.numeric(if (relation == "greater") dyads$Metric[k] > threshold else dyads$Metric[k] < threshold) else NA_real_
    adj[dyads$From[k], dyads$To[k]] <- value
    if (!directed) adj[dyads$To[k], dyads$From[k]] <- value
  }
  diag(adj) <- 0
  idx <- which(if (directed) row(adj) != col(adj) else upper.tri(adj), arr.ind = TRUE)
  all_pairs <- data.frame(From = ids[idx[, 1]], To = ids[idx[, 2]], Tie = adj[idx])
  # Keep every possible dyad, including pairs absent from older sparse tables.
  if (!directed) {
    swap <- match(dyads$From, ids) > match(dyads$To, ids)
    tmp <- dyads$From[swap]; dyads$From[swap] <- dyads$To[swap]; dyads$To[swap] <- tmp
  }
  all_pairs <- merge(all_pairs, dyads, by = c("From", "To"), all.x = TRUE, sort = FALSE)
  unknown <- sum(is.na(all_pairs$Tie))
  if (unknown && missing == "error") stop("Unobserved dyads are present; review coverage or explicitly set missing = 'model'.", call. = FALSE)
  observed <- all_pairs$Tie[!is.na(all_pairs$Tie)]
  if (length(unique(observed)) < 2L) stop("Both observed ties and observed non-ties are required for estimation.", call. = FALSE)
  net <- network::network(adj, directed = directed, loops = FALSE, matrix.type = "adjacency")
  net <- network::set.vertex.attribute(net, "vertex.names", ids)
  nodes <- data.frame(Node = ids)
  if (!is.null(node_data)) {
    if (!is.data.frame(node_data) || !"Node" %in% names(node_data) ||
        anyNA(node_data) || anyDuplicated(node_data$Node) || !setequal(node_data$Node, ids))
      stop("`node_data` must contain complete covariates and exactly one row for every Node.", call. = FALSE)
    nodes <- node_data[match(ids, node_data$Node), , drop = FALSE]
    if (any(c("vertex.names", "na") %in% setdiff(names(nodes), "Node")))
      stop("Reserved network attribute name in node_data.", call. = FALSE)
    for (nm in setdiff(names(nodes), "Node")) net <- network::set.vertex.attribute(net, nm, nodes[[nm]])
  }
  env <- new.env(parent = environment(terms)); env$net <- net
  formula <- stats::as.formula(call("~", quote(net), terms[[2L]]), env = env)
  warnings <- character()
  fitted <- with_preserved_rng_seed(seed, withCallingHandlers(
    ergm::ergm(formula, estimate = "MLE", control = control %||% ergm::control.ergm()),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))))
  coef <- stats::coef(fitted); vc <- stats::vcov(fitted)
  se <- sqrt(diag(vc))
  coefficients <- data.frame(Term = names(coef), Estimate = unname(coef),
    SE = unname(se), Lower = unname(coef - stats::qnorm(.975) * se),
    Upper = unname(coef + stats::qnorm(.975) * se), Level = .95)
  bad <- !is.finite(coef) | !is.finite(se) | se <= 0 | isTRUE(fitted$failure)
  coefficients$Status <- ifelse(bad, "unavailable", "approximate_review_required")
  coefficients$Lower[bad] <- coefficients$Upper[bad] <- NA_real_
  simulated <- if (is.null(gof)) NULL else with_preserved_rng_seed(seed,
    ergm::gof(fitted, GOF = gof, control = gof_control %||% ergm::control.gof.ergm()))
  summaries <- if (is.null(simulated)) list() else simulated[grepl("^summary[.]", names(simulated))]
  gof_table <- if (!length(summaries)) data.frame() else do.call(rbind,
    lapply(names(summaries), function(nm) data.frame(Check = sub("^summary[.]", "", nm),
      Statistic = rownames(summaries[[nm]]), as.data.frame(summaries[[nm]]), row.names = NULL)))
  simulations <- if (is.null(simulated)) list() else simulated[grepl("^sim[.]", names(simulated))]
  gof_draws <- if (!length(simulations)) data.frame() else do.call(rbind,
    lapply(names(simulations), function(nm) {
      m <- as.matrix(simulations[[nm]])
      data.frame(Check = sub("^sim[.]", "", nm), Replicate = rep(seq_len(nrow(m)), ncol(m)),
        Statistic = rep(colnames(m) %||% as.character(seq_len(ncol(m))), each = nrow(m)),
        Value = as.vector(m))
    }))
  structure(list(fit = fitted, network = net, nodes = nodes, dyads = all_pairs,
    coefficients = coefficients, gof = simulated, gof_table = gof_table, gof_draws = gof_draws,
    diagnostics = data.frame(BackendFailure = isTRUE(fitted$failure),
      MPLEisMLE = isTRUE(fitted$MPLE_is_MLE), FiniteCoefficients = all(is.finite(coef)),
      MCMCSampleStored = !is.null(fitted$sample),
      Status = if (isTRUE(fitted$failure)) "backend_failure" else "diagnostic_review_required"),
    warnings = data.frame(Message = unique(warnings)),
    source_settings = x$settings,
    settings = data.frame(Directed = directed, Threshold = threshold, Relation = relation,
      Metric = if (halo) x$settings$method else x$settings$weight_metric,
      MinPairN = x$settings$min_pair_n, Missing = missing, MissingDyads = unknown,
      ObservedDyads = length(observed), Ties = sum(observed),
      Seed = if (is.null(seed)) NA_integer_ else as.integer(seed[1L]),
      Terms = paste(deparse(terms), collapse = " "), Estimation = fitted$estimate,
      BackendVersion = as.character(utils::packageVersion("ergm"))),
    notes = c("Binary network likelihood conditional on the declared threshold and observed node set; not an MFRM rating likelihood.",
      "Review backend warnings, MCMC mixing/degeneracy and simulation GOF before interpreting approximate coefficient intervals.",
      "Missing dyads remain unobserved, not zero. Observed-data likelihood assumes ignorable missingness.",
      "Edge-estimation uncertainty and data-driven threshold/model selection are not included.")),
    class = "mfrm_network_ergm")
}

#' @export
print.mfrm_network_ergm <- function(x, ...) {
  print(x$settings, row.names = FALSE)
  print(x$coefficients, row.names = FALSE)
  print(x$diagnostics, row.names = FALSE)
  invisible(x)
}

#' @rdname fit_rater_network
#' @param object An `mfrm_network_ergm` result.
#' @param ... For `summary.mfrm_network_ergm()`: Unused; unexpected arguments are rejected by extractors.
#'
#' For `plot()`: For GOF, arguments passed to ergm's `plot.gof`; for coefficients,
#'   arguments passed to base [graphics::plot()].
#' @export
summary.mfrm_network_ergm <- function(object, ...) {
  object[c("settings", "coefficients", "diagnostics", "gof_table", "warnings", "notes")]
}

#' @rdname fit_rater_network
#' @export
coef.mfrm_network_ergm <- function(object, ...) {
  rlang::check_dots_empty()
  stats::setNames(object$coefficients$Estimate, object$coefficients$Term)
}

#' @rdname fit_rater_network
#' @export
vcov.mfrm_network_ergm <- function(object, ...) {
  rlang::check_dots_empty()
  if (!requireNamespace("ergm", quietly = TRUE))
    stop("Install the optional ergm package to extract the backend covariance.", call. = FALSE)
  stats::vcov(object$fit)
}

#' @rdname fit_rater_network
#' @param parm Coefficient names or positive integer positions; NULL selects all.
#' @param level Confidence level strictly between zero and one.
#' @details `confint()` preserves unavailable coefficients as missing bounds.
#'   Its matrix attributes retain `level`, `method`, `diagnostics` and `notes`.
#'   These intervals remain conditional on the selected binary network and
#'   model; they do not propagate estimated-edge or model-selection uncertainty.
#'   Changing `level` does not change the fitted object or its stored plot.
#' @export
confint.mfrm_network_ergm <- function(object, parm = NULL, level = .95, ...) {
  rlang::check_dots_empty()
  if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) ||
      length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be one finite number strictly between zero and one.", call. = FALSE)
  tab <- object$coefficients
  index <- if (is.null(parm)) seq_len(nrow(tab)) else if (is.character(parm)) {
    match(parm, tab$Term)
  } else if (is.numeric(parm) && !is.complex(parm) && !is.object(parm) &&
             is.null(dim(parm)) && all(is.finite(parm)) && all(parm == floor(parm))) {
    parm
  } else NA_integer_
  if (anyNA(index) || any(index < 1 | index > nrow(tab)))
    stop("`parm` must select existing coefficient names or positions.", call. = FALSE)
  tab <- tab[index, , drop = FALSE]
  available <- tab$Status %in% "approximate_review_required" &
    is.finite(tab$Estimate) & is.finite(tab$SE) & tab$SE > 0 &
    !isTRUE(object$diagnostics$BackendFailure)
  z <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  out <- cbind(tab$Estimate - z * tab$SE, tab$Estimate + z * tab$SE)
  out[!available, ] <- NA_real_
  rownames(out) <- tab$Term
  colnames(out) <- paste0(format(100 * c((1-level)/2, (1+level)/2), trim = TRUE), "%")
  tab$Status[!available] <- "unavailable"
  attr(out, "level") <- unname(level)
  attr(out, "method") <- "Approximate pointwise normal ERGM coefficient intervals"
  attr(out, "diagnostics") <- tab[c("Term", "Estimate", "SE", "Status")]
  attr(out, "notes") <- object$notes
  out
}

#' @rdname fit_rater_network
#' @details `plot()` method: Plot ERGM model checks or coefficient intervals
#' @param type `"coefficients"` or `"gof"` (requires saved GOF).
#' @param draw Whether to draw; `FALSE` returns plot data only.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object containing coefficient rows or
#'   GOF summaries and long-form simulation draws. Use [plot_data()] to export.
#' @details For `plot()`: Coefficient intervals are approximate pointwise normal intervals.
#'   They are not evidence of MCMC convergence. The GOF plot delegates to ergm;
#'   its graphical options remain available through `...`.
#' @export
plot.mfrm_network_ergm <- function(x, type = c("coefficients", "gof"), draw = TRUE, ...) {
  type <- match.arg(type)
  if (type == "gof" && is.null(x$gof)) stop("No saved GOF: supply gof when fitting.", call. = FALSE)
  if (isTRUE(draw)) {
    if (type == "gof") {
      args <- list(x = x$gof, main = "")
      dots <- list(...); args[names(dots)] <- dots
      do.call(graphics::plot, args)
    } else {
      z <- x$coefficients; keep <- is.finite(z$Estimate) & is.finite(z$Lower) & is.finite(z$Upper)
      z <- z[keep, , drop = FALSE]
      if (!nrow(z)) stop("No finite coefficient intervals to plot.", call. = FALSE)
      args <- list(x = z$Estimate, y = seq_len(nrow(z)),
        xlim = range(c(z$Lower, z$Upper)), yaxt = "n", xlab = "Conditional log-odds coefficient", ylab = "", main = "")
      dots <- list(...); args[names(dots)] <- dots
      do.call(graphics::plot, args)
      graphics::segments(z$Lower, seq_len(nrow(z)), z$Upper, seq_len(nrow(z)))
      graphics::axis(2, seq_len(nrow(z)), z$Term, las = 1)
      graphics::abline(v = 0, lty = 2)
    }
  }
  invisible(new_mfrm_plot_data("network_ergm", list(plot = type,
    table = if (type == "gof") x$gof_table else x$coefficients,
    draws = if (type == "gof") x$gof_draws else NULL, notes = x$notes)))
}
