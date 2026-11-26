
mfrm_fair_reference_score_fit <- function(fit) .mfrm_validate_fit(fit)


mfrm_fair_reference_plan <- function(fit, plan) {
    mfrm_fair_reference_score_fit(fit)
    facets <- fit$encoding$base_facets
    levels <- stats::setNames(fit$encoding$facet_levels, fit$encoding$facets)
    if ("Weight" %in% facets)
        stop("The reference plan reserves the column name 'Weight'.", call. = FALSE)
    if (!is.data.frame(plan) || !nrow(plan) || anyDuplicated(names(plan)) || !setequal(setdiff(names(plan),
        "Weight"), facets))
        stop("`reference_plan` must contain every fitted facet and optional Weight only.", call. = FALSE)
    out <- as.data.frame(plan[facets], stringsAsFactors = FALSE)
    for (f in facets) {
        x <- plan[[f]]
        if (!(is.character(x) || is.factor(x) || is.numeric(x)) || !is.null(dim(x)) || anyNA(x) || any(!as.character(x) %in%
            levels[[f]]))
            stop("Unknown or invalid reference-plan level in facet `", f, "`.", call. = FALSE)
        out[[f]] <- as.character(x)
    }
    if (anyDuplicated(out))
        stop("Use weights for duplicate reference-plan cells.", call. = FALSE)
    w <- plan$Weight %||% rep(1, nrow(out))
    if (!is.numeric(w) || !is.null(dim(w)) || any(!is.finite(w)) || any(w <= 0))
        stop("Reference-plan weights must be finite and positive.", call. = FALSE)
    out$Weight <- w
    out$NormalizedWeight <- (w/max(w))/sum(w/max(w))
    out
}

mfrm_fair_reference_scores <- function(fit, plan, level = NULL) {
    .mfrm_validate_fit(fit)
    ready <- .mfrm_readiness(fit)$ready
    nd <- plan[fit$encoding$base_facets]
    eval_plan <- function(theta) {
        out <- rep(NA_real_, length(theta))
        ok <- which(is.finite(theta))
        batch <- max(1L, as.integer(50000/nrow(plan)))
        starts <- if (length(ok))
            seq.int(1L, length(ok), by = batch)
        else integer()
        for (start in starts) {
            ids <- ok[seq.int(start, min(length(ok), start + batch - 1L))]
            block <- nd[rep(seq_len(nrow(plan)), length(ids)), , drop = FALSE]
            block$theta <- rep(theta[ids], each = nrow(plan))
            score <- matrix(predict(fit, block, type = "response"), nrow = nrow(plan))
            out[ids] <- drop(crossprod(plan$NormalizedWeight, score))
        }
        out
    }
    score <- if (ready)
        eval_plan(fit$person$Estimate)
    else rep(NA_real_, nrow(fit$person))
    out <- data.frame(Person = fit$person$Person, Estimate = score, Status = ifelse(is.finite(score), paste0(tolower(fit$score_method),
        "_plug_in_estimate"), "unavailable_calibration"), Target = "common_plan_expected_average", Units = "original_score",
        Uncertainty = "not_computed", stringsAsFactors = FALSE)
    if (!is.null(level)) {
        out$Lower <- out$Upper <- NA_real_
        if (ready) {
            ci <- confint(fit, component = "person", level = level)
            out$Lower <- eval_plan(ci[, 1L])
            out$Upper <- eval_plan(ci[, 2L])
        }
        out$Level <- level
        basis <- if (identical(fit$score_method, "EAP"))
            "conditional_normal_posterior_approximation"
        else "conditional_normal_ability_interval"
        out$Uncertainty <- ifelse(is.finite(out$Lower) & is.finite(out$Upper), basis, "unavailable")
    }
    return(out)
}
