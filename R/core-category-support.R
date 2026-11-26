mfrm_category_settings <- function(prep, config = list()) {
    keep <- prep$keep_original %||% config$keep_original %||% config$replay_inputs$keep_original
    policy <- if (identical(keep, TRUE))
        "preserve"
    else if (identical(keep, FALSE))
        "collapse"
    else "not_recorded"
    map <- prep$score_map %||% config$score_map
    recoded <- NA
    if (is.data.frame(map) && nrow(map) > 0L && all(c("OriginalScore", "InternalScore") %in% names(map)) &&
        is.numeric(map$OriginalScore) && is.numeric(map$InternalScore) && all(is.finite(map$OriginalScore)) &&
        all(is.finite(map$InternalScore))) {
        recoded <- any(map$OriginalScore != map$InternalScore)
    }
    list(CategoryPolicy = policy, ScoreRecoded = recoded)
}

mfrm_original_score_labels <- function(table, map) {
    if (!"OriginalScore" %in% names(table) && "Score" %in% names(table) && is.data.frame(map) && all(c("OriginalScore",
        "InternalScore") %in% names(map)) && isTRUE(any(map$OriginalScore != map$InternalScore))) {
        table$OriginalScore <- map$OriginalScore[match(table$Score, map$InternalScore)]
    }
    table
}
