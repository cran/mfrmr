# Display state is scoped to the call and never stored in a fitted model.
mfrmr_start_fit_progress <- function(enabled) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled))
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  previous <- getOption("mfrmr._fit_progress")
  id <- if (enabled) cli_progress_bar(name = "mfrmr",
    status = "Fitting model", total = NA,
    format = "{cli::pb_spin} {cli::pb_name}: {cli::pb_status} [{cli::pb_elapsed}]",
    clear = TRUE, auto_terminate = FALSE, .auto_close = FALSE) else NULL
  options(mfrmr._fit_progress = id)
  function() {
    options(mfrmr._fit_progress = previous)
    if (!is.null(id)) cli_progress_done(id = id, result = "clear")
    invisible(NULL)
  }
}

mfrmr_fit_progress <- function(status) {
  id <- getOption("mfrmr._fit_progress")
  if (!is.null(id)) cli_progress_update(status = status, id = id)
  invisible(NULL)
}
