#' Estimate abilities for people already included in a fitted RSM
#'
#' Use the existing ratings and fitted model to estimate each person's ability.
#' This is the common scoring entry for supported ordinary, shared-rater and
#' testlet rating-scale models (RSMs). It reuses the fitted calibration, such
#' as rater severity and category thresholds; it does not fit the model again.
#'
#' The estimate is the mean of the conditional ability distribution (EAP).
#' Posterior SDs and continuous equal-tail intervals describe its uncertainty
#' with calibration held fixed. Selecting persons changes the returned rows,
#' not the ratings used to condition shared effects.
#' @param fit A numerically ready ordinary RSM MML, testlet or shared-rater fit.
#' @param persons Distinct source Person IDs; `NULL` returns all source Persons,
#'   including Persons whose assigned scores are all missing. For a shared-rater
#'   fit, start with a few actual IDs: scoring everyone can be slow. Selecting
#'   IDs changes output rows, not the data used to account for shared raters.
#' @param level Conditional equal-tail interval probability, between zero and one.
#' @param quad_points Extension quadrature order; `NULL` uses its fit's order.
#'   Ordinary ability integration is continuous and does not use quadrature;
#'   leave this argument `NULL` for ordinary fits.
#' @param cut_scores Optional finite, strictly increasing numeric ability cuts,
#'   with optional distinct names. Ordinary RSM scoring only. Returns a
#'   `decisions` table of probabilities above each cut and in successive bands
#'   (`Lower < ability <= Upper`). Calibration, population parameters and cuts
#'   are held fixed. Labels do not establish CEFR alignment or a pass/fail rule.
#' @details Ordinary fits use the bounded model/population specification of
#'   [mfrm_response_diagnostics()]: additive unanchored severity facets and
#'   known N(0,1) or an estimated intercept-only normal population. Both the
#'   fitted population mean and variance are retained. EAP, posterior SD and
#'   interval endpoints integrate the continuous conditional density; endpoints
#'   are not quadrature-grid quantiles. Ordinary plug-in scoring is unchanged.
#'
#'   For extensions this calls [predict.mfrm_testlet()] or
#'   [score_mfrm_random_rater()] on the complete source roster. Their numerical
#'   checks, approximation limits and uncertainty definitions apply unchanged.
#'   Use those functions directly for a different complete scoring roster.
#'
#'   Intervals condition on fitted calibration and the assumed normal
#'   population. They exclude calibration-estimation uncertainty, do not test
#'   Person differences and do not guarantee frequentist coverage for each
#'   fixed ability. Missing-only Persons remain `prior_only`; zero ability
#'   variance and failed integration return `unavailable`, not zero-width
#'   intervals. No missing scores are imputed.
#' @inheritSection predict_mfrm_units Conditional ability decisions
#' @return Ordinary fits return `mfrm_person_scores` with `table`,
#'   `scoring_data`, data usage, settings and source metadata. Extensions retain
#'   their existing scoring classes. If `cut_scores` is supplied, `decisions`
#'   retains probabilities, boundaries, labels, priors and availability;
#'   `summary(x)$decisions` and `apa_table(x, which = "decisions")` reuse it.
#'   Supply saved scores to [compare_mfrm()]
#'   as `person_scores`; extension scores also support model-aware maps through
#'   [mfrm_results()]. Saved plots/reports do not recompute scores.
#'   In `table`, `Status` and `Reason` distinguish `available_conditional`
#'   (a computed conditional posterior), `prior_only` (no observed responses;
#'   the estimate and interval describe the assumed population alone), and
#'   `unavailable` (the estimate, SD and bounds are withheld, with the cause).
#'   Keep these columns with exported scores. More integration points cannot
#'   supply missing responses or resolve a zero ability-variance restriction.
#' @seealso [score_mfrm_calibration()] for new people under an eligible saved
#'   ordinary RSM/PCM calibration; [predict_mfrm_units()] for supplied response
#'   rows under a fitted calibration; [plot.mfrm_testlet_scores()],
#'   [mfrm_response_diagnostics()]
#' @examples
#' ratings <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score")
#' scores <- score_mfrm_persons(fit, persons = unique(ratings$Person)[1:4],
#'   cut_scores = c(Basic = -0.5, Advanced = 0.8))
#' scores$table
#' scores$decisions
#' apa_table(scores, which = "decisions")
#' plot(scores, draw = FALSE)
#' @export
score_mfrm_persons <- function(fit, persons = NULL, level = .95, quad_points = NULL,
                               cut_scores = NULL) {
  cut_scores <- mfrm_validate_cut_scores(cut_scores)
  if (!is.null(cut_scores) && mfrm_extended_fit(fit))
    stop("Cut-score probabilities currently require ordinary RSM/PCM scoring.", call. = FALSE)
  if (inherits(fit, "mfrm_testlet")) return(predict(fit, persons = persons,
    level = level, quad_points = quad_points %||% fit$settings$quad_points))
  if (inherits(fit, "mfrm_random_rater")) return(score_mfrm_random_rater(fit,
    persons = persons, level = level, quad_points = quad_points %||% fit$settings$quad_points))
  input <- mfrm_ordinary_response_input(fit)
  if (!is.null(quad_points)) stop("Ordinary source scoring uses continuous integration; leave quad_points = NULL.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("Supply 0 < level < 1.", call. = FALSE)
  ids <- unique(input$assigned_data[[input$columns$person]])
  persons <- persons %||% ids
  if (!is.character(persons) || !length(persons) || anyNA(persons) || anyDuplicated(persons) || !all(persons %in% ids)) stop("Choose distinct Person IDs in the source roster.", call. = FALSE)
  probabilities <- matrix(NA_real_, length(persons), 2L*length(cut_scores)+1L)
  table <- lapply(persons, function(id) {
    rows <- which(input$person == id); reason <- ""
    value <- setNames(rep(NA_real_, 4), c("Estimate", "ConditionalSD", "Lower", "Upper"))
    if (input$sd == 0) {
      status <- "unavailable"; reason <- "The fitted ability variance is zero; individual ability scores are withheld."
    } else if (!length(rows)) {
      value <- c(Estimate = input$mean, ConditionalSD = input$sd,
        Lower = stats::qnorm((1-level)/2, input$mean, input$sd),
        Upper = stats::qnorm((1-level)/2, input$mean, input$sd, lower.tail = FALSE))
      status <- "prior_only"
      reason <- "No observed responses for this Person; the estimate, SD and interval describe only the assumed population, not measured individual ability."
      if (!is.null(cut_scores)) probabilities[match(id,persons),] <<- mfrm_cut_probabilities(
        function(a,b) stats::integrate(stats::dnorm,a,b,rel.tol=1e-10,abs.tol=0)$value,
        (unname(cut_scores)-input$mean)/input$sd, 1)
    } else {
      answer <- tryCatch(withCallingHandlers({
        loglikelihood <- function(z) {
          theta <- input$mean + input$sd * z
          loglik <- vapply(rows, function(i) {
            lw <- outer(theta - input$offset[i], 0:length(input$steps)) -
              matrix(c(0,cumsum(input$steps)),length(z),length(input$steps)+1L,byrow=TRUE)
            lw <- lw - row_max_fast(lw)
            lw[, input$y[i]+1L] - mfrm_testlet_logsum_rows(lw)
          }, numeric(length(z)))
          rowSums(matrix(loglik, nrow=length(z)))
        }
        ans <- mfrm_person_posterior_interval(loglikelihood, input$sd, level,
          if (is.null(cut_scores)) NULL else (unname(cut_scores)-input$mean)/input$sd)
        if (!is.null(cut_scores)) probabilities[match(id,persons),] <<- ans[-seq_len(5L)]
        ans <- ans[1:4]
        ans[c("Estimate", "Lower", "Upper")] <- ans[c("Estimate", "Lower", "Upper")] + input$mean
        ans
      }, warning=function(w) stop(conditionMessage(w),call.=FALSE)),
      error=function(e) {reason <<- conditionMessage(e); NULL})
      status <- if (is.null(answer)) "unavailable" else "available_conditional"
      if (!is.null(answer)) value <- answer
    }
    data.frame(Person=id, Observed=length(rows), as.list(value), Status=status, Reason=reason, row.names=NULL)
  })
  out <- structure(list(table=do.call(rbind,table), scoring_data=input$assigned_data,
    data_usage=c(Input=input$input_rows,Observed=nrow(input$data),Omitted=length(input$omitted_rows)),
    omitted_rows=input$omitted_rows, source=mfrm_response_source(fit),
    settings=list(level=level, calibration_uncertainty=FALSE, person_mean=input$mean,
      person_variance=input$sd^2, integration="Continuous normal ability posterior",
      target="Person marginal posterior conditional on calibration and the complete scoring roster")),
    class="mfrm_person_scores")
  if (!is.null(cut_scores)) {
    out$settings$cut_scores <- cut_scores
    out$decisions <- mfrm_decision_table(probabilities, persons, cut_scores,
      input$mean,input$sd,out$table$Status,out$table$Reason)
    out$notes <- mfrm_decision_note()
  }
  out
}

#' @rdname score_mfrm_persons
#' @param object,x An ordinary `mfrm_person_scores` result.
#' @param ... Display controls passed to [plot.mfrm_testlet_scores()]; unused
#'   by print and summary.
#' @export
summary.mfrm_person_scores <- function(object, ...) {
  mfrm_validate_decisions(object)
  out <- summary.mfrm_testlet_scores(object, ...)
  if (!is.null(object$decisions)) {
    out$decisions <- object$decisions
    out$notes <- object$notes
  }
  out
}

#' @rdname score_mfrm_persons
#' @export
print.mfrm_person_scores <- function(x, ...) {
  mfrm_validate_decisions(x)
  cat("Conditional Person scores from ordinary RSM calibration\n")
  print(x$table, row.names=FALSE)
  cat("Intervals exclude calibration-estimation uncertainty.\n")
  if (!is.null(x$decisions)) {
    cat("\nConditional decision probabilities (first 10; full table in $decisions)\n")
    print(utils::head(x$decisions, 10L), row.names=FALSE)
    print_wrapped_line(mfrm_decision_note())
  }
  invisible(x)
}

#' @rdname score_mfrm_persons
#' @export
plot.mfrm_person_scores <- function(x, ...) {
  out <- plot.mfrm_testlet_scores(x, ...)
  out$name <- "person_scores"
  invisible(out)
}

mfrm_validate_person_scores <- function(fit, scores) {
  expected <- if (inherits(fit,"mfrm_testlet")) "mfrm_testlet_scores" else
    if (inherits(fit,"mfrm_random_rater")) "mfrm_random_rater_scores" else "mfrm_person_scores"
  input <- if (mfrm_extended_fit(fit)) fit$input else mfrm_ordinary_response_input(fit)
  columns <- names(input$assigned_data %||% input$data)
  roster <- input$assigned_data %||% input$data
  if (!inherits(scores,expected) || !identical(scores$source,mfrm_response_source(fit)) ||
      is.null(scores$scoring_data) || !identical(mfrm_compare_events(scores$scoring_data,columns),
        mfrm_compare_events(input$assigned_data %||% input$data,columns))) stop(
    "Supply saved Person scores with matching calibration and the complete source roster; regenerate older scores lacking scoring_data.",call.=FALSE)
  if (!identical(scores$settings$calibration_uncertainty,FALSE) ||
      !is.data.frame(scores$table) || !all(c("Person","Observed","Estimate","ConditionalSD","Lower","Upper","Status","Reason") %in% names(scores$table)) ||
      anyDuplicated(scores$table$Person) || !all(scores$table$Person %in% roster[[input$columns$person]])) stop("Saved Person-score identities or conditional uncertainty fields are invalid.",call.=FALSE)
  invisible(TRUE)
}

mfrm_compare_persons <- function(fits, scores) {
  if (!is.list(scores) || length(scores)!=2L) stop("person_scores must contain two saved scoring results in fit order.",call.=FALSE)
  for (i in 1:2) mfrm_validate_person_scores(fits[[i]],scores[[i]])
  a <- scores[[1]]; b <- scores[[2]]
  if (!setequal(a$table$Person,b$table$Person) || !identical(a$settings$level,b$settings$level)) stop("Compare the same requested Persons and conditional interval level.",call.=FALSE)
  origin <- vapply(fits,function(f) if(mfrm_extended_fit(f)) 0 else mfrm_ordinary_response_input(f)$mean,numeric(1))
  ids <- sort(a$table$Person)
  aa <- a$table[match(ids,a$table$Person),]; bb <- b$table[match(ids,b$table$Person),]
  if (!identical(aa$Observed,bb$Observed)) stop("Person observed-row counts differ.",call.=FALSE)
  available <- aa$Status=="available_conditional" & bb$Status=="available_conditional"
  status <- ifelse(available,"available_descriptive",ifelse(aa$Status=="prior_only" | bb$Status=="prior_only","prior_only","unavailable"))
  table <- data.frame(Person=ids,Observed=aa$Observed,
    SourceReference=aa$Estimate,SourceComparison=bb$Estimate,
    OriginReference=origin[1],OriginComparison=origin[2],
    Reference=aa$Estimate-origin[1],Comparison=bb$Estimate-origin[2],
    ConditionalSDReference=aa$ConditionalSD,ConditionalSDComparison=bb$ConditionalSD,
    LowerReference=aa$Lower-origin[1],UpperReference=aa$Upper-origin[1],
    LowerComparison=bb$Lower-origin[2],UpperComparison=bb$Upper-origin[2],
    ReferenceStatus=aa$Status,ComparisonStatus=bb$Status,Status=status,
    Reason=ifelse(available,"",paste0("Reference: ",aa$Status," ",aa$Reason," | Comparison: ",bb$Status," ",bb$Reason)))
  table$Difference <- ifelse(available,table$Comparison-table$Reference,NA_real_)
  list(table=table,level=a$settings$level,
    note="EAPs and conditional endpoints are centered at each fitted population mean, on the unit-slope logit scale. Calibration is fixed; conditional intervals are not intervals for model differences. Prior-only and unavailable differences are withheld; no ranking or Person-difference test is supplied.")
}
