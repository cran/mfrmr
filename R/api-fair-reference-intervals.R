confint.mfrm_fair_average <- function(object, parm = NULL, level = 0.95, ...) {
  rlang::check_dots_empty()
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be in (0, 1).", call. = FALSE)
  fit <- object$reference_fit
  .mfrm_validate_fit(fit)
  plan <- object$reference_plan
  if (is.null(plan)) stop("Supply an explicit reference_plan to fair_average_table().", call. = FALSE)
  plan <- mfrm_fair_reference_plan(fit, plan[setdiff(names(plan), "NormalizedWeight")])
  current <- mfrm_fair_reference_scores(fit, plan)
  saved <- object$reference_scores
  if (!identical(current$Person, saved$Person) ||
      !isTRUE(all.equal(current$Estimate, saved$Estimate, tolerance = 1e-10)))
    stop("The retained calibration or reference plan changed; recreate fair_average_table().", call. = FALSE)
  out <- mfrm_fair_reference_scores(fit, plan, level)
  if (!is.null(parm)) {
    if (!is.character(parm) || anyNA(parm) || anyDuplicated(parm) || any(!parm %in% out$Person))
      stop("`parm` must select distinct known Person identifiers.", call. = FALSE)
    out <- out[match(parm, out$Person), , drop = FALSE]
  }
  out
}
