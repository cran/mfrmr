# Observe the arguments mfrmr sends to graphics, including implicit headings
# drawn inside plot()/barplot(); never intercept graphics' own internal calls.
local_graphics_headings <- function(record, .env = parent.frame()) {
  functions <- c("plot", "barplot", "image", "hist", "dotchart", "plot.new")
  replacements <- lapply(functions, function(name) {
    original <- get(name, envir = asNamespace("mfrmr"), inherits = TRUE)
    function(...) {
      args <- list(...)
      record(args$main)
      do.call(original, args)
    }
  })
  do.call(testthat::local_mocked_bindings,
    c(stats::setNames(replacements, functions), list(.package = "mfrmr", .env = .env)))
}
