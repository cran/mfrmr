# Recreate fitted inputs for examples about tables, plots and reports.
# Sourcing this file defines a function; it does not run estimation.
# Each saved object contains the fit and its matching default diagnostics.
# The estimation and diagnostic entry-point examples still run those methods.
#
# library(mfrmr)
# source(system.file("examples", "ordinary-fits.R", package = "mfrmr"))
# saved <- ordinary_fit_examples()
# for (name in names(saved)) saveRDS(saved[[name]], paste0(name, ".rds"), compress = "xz")
ordinary_fit_examples <- function() {
  make_example <- function(dataset, method) {
    ratings <- mfrmr::load_mfrmr_data(dataset)
    args <- list(data = ratings, person = "Person",
      facets = c("Rater", "Criterion"), score = "Score", method = method,
      model = "RSM")
    if (method == "JML") args$maxit <- 300
    fit <- do.call(mfrmr::fit_mfrm, args)
    # Constant population formulas must not retain the generation session.
    for (field in c("formula", "terms")) {
      if (!is.null(fit$population[[field]]))
        environment(fit$population[[field]]) <- baseenv()
      if (!is.null(fit$config$population_spec[[field]]))
        environment(fit$config$population_spec[[field]]) <- baseenv()
    }
    list(dataset = dataset, fit = fit,
      diagnostics = mfrmr::diagnose_mfrm(fit))
  }
  list(
    `core-jml` = make_example("example_core", "JML"),
    `bias-jml` = make_example("example_bias", "JML"),
    `operational-mml` = make_example("example_operational", "MML")
  )
}
