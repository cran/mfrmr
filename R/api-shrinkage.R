

apply_empirical_bayes_shrinkage <- function(fit, facet_prior_sd = NULL, shrink_person = FALSE) {
    .mfrm_validate_fit(fit)
    shrink_vec <- function(est, se, tau2_override = NULL) {
        ok <- is.finite(est) & is.finite(se) & (se > 0)
        if (sum(ok) < 2L)
            return(list(shrunk = est, shrunk_se = se, tau2 = 0))
        mu <- mean(est[ok])
        tau2 <- if (is.null(tau2_override))
            max(0, stats::var(est[ok]) - mean(se[ok]^2))
        else as.numeric(tau2_override)^2
        b_wt <- ifelse(ok & (tau2 + se^2 > 0), tau2/(tau2 + se^2), 0)
        est_sh <- ifelse(ok, mu + b_wt * (est - mu), est)
        se_sh <- ifelse(ok, sqrt(b_wt * se^2), se)
        list(shrunk = est_sh, shrunk_se = se_sh, tau2 = tau2)
    }
    out <- fit
    out$facets$ShrunkEstimate <- rep(NA_real_, nrow(out$facets))
    out$facets$ShrunkSE <- rep(NA_real_, nrow(out$facets))
    for (fnm in unique(out$facets$Facet)) {
        idx <- which(out$facets$Facet == fnm)
        res <- shrink_vec(out$facets$Estimate[idx], out$facets$ModelSE[idx], facet_prior_sd)
        out$facets$ShrunkEstimate[idx] <- res$shrunk
        out$facets$ShrunkSE[idx] <- res$shrunk_se
    }
    if (isTRUE(shrink_person)) {
        pr <- shrink_vec(out$person$Estimate, out$person$ModelSE, facet_prior_sd)
        out$person$ShrunkEstimate <- pr$shrunk
        out$person$ShrunkSE <- pr$shrunk_se
    }
    return(out)
}
