#' Automated-scoring compatibility interface for score comparisons
#'
#' Prefer [score_comparison_table()] for new analyses. This earlier interface
#' retains its argument names, automated/human output columns and default
#' Williamson (2012) numerical screen. The generic interface uses score/reference
#' names and applies no thresholds unless requested.
#'
#' @param automated,human,human_reference Earlier names for `scores`,
#'   `reference`, and `reference2` in [score_comparison_table()].
#' @param categories Declared numeric score categories in increasing order.
#' @param automated_categories Earlier name for `score_categories`.
#' @param smd `"pooled"` or `"human"`; the latter maps to `"reference"`.
#' @param qwk_min,correlation_min,max_degradation,max_abs_smd Numerical screen
#'   limits; see [score_comparison_table()]. These defaults retain the earlier
#'   automated-scoring examples and are not universal acceptance criteria.
#' @param fit,diagnostics,rater_facet Optional saved MFRM source and rater facet;
#'   see [score_comparison_table()].
#' @param automated_rater Earlier name for `rater`.
#' @return The same list structure as [score_comparison_table()], retaining the
#'   earlier automated/human column and category-table names. Saved results
#'   created through the earlier interface can still be read directly.
#' @seealso [score_comparison_table()]
#' @export
automated_scoring_table <- function(automated, human, human_reference = NULL,
                                   categories, automated_categories = automated,
                                   smd = c("pooled", "human"), qwk_min = .70,
                                   correlation_min = .70, max_degradation = .10,
                                   max_abs_smd = .15, fit = NULL,
                                   diagnostics = NULL, rater_facet = NULL,
                                   automated_rater = NULL) {
  smd <- match.arg(smd)
  x <- score_comparison_table(scores = automated, reference = human,
    reference2 = human_reference, categories = categories,
    score_categories = automated_categories,
    smd = if (smd == "human") "reference" else smd, screen = "williamson2012",
    qwk_min = qwk_min, correlation_min = correlation_min,
    max_degradation = max_degradation, max_abs_smd = max_abs_smd,
    fit = fit, diagnostics = diagnostics, rater_facet = rater_facet,
    rater = automated_rater)
  rename <- c(Reference2Supplied = "HumanReferenceSupplied",
    MeanScore = "MeanAutomated", MeanReference = "MeanHuman",
    SDScore = "SDAutomated", SDReference = "SDHuman", ReferenceQWK = "HumanQWK",
    ReferenceCorrelation = "HumanCorrelation", SMDReference = "SMDHuman")
  names(x$table)[match(names(rename), names(x$table))] <- unname(rename)
  x$table$SMDBasis <- x$settings$smd <- smd
  x$criteria$Basis <- sub("_smd_reference$", "_smd_human", x$criteria$Basis)
  names(x$cases) <- sub("Reference2", "HumanReference", names(x$cases))
  names(x$cases) <- sub("ScoreCategory", "AutomatedCategory", names(x$cases))
  names(x$cases) <- sub("MissingScore$", "MissingAutomated", names(x$cases))
  names(x$cases) <- sub("MissingReference$", "MissingHuman", names(x$cases))
  names(x$counts) <- c("automated_human", "human_human")
  for (nm in names(x$counts)) {
    if (!is.null(x$counts[[nm]])) names(dimnames(x$counts[[nm]])) <- c("Scorer", "Human")
  }
  x
}
