# Regenerate the synthetic results used by simulation summary and plot examples.
# Sourcing this file only defines a function; generation is an explicit action.
# One replication illustrates object structure, not reliable recovery, power,
# coverage or design recommendations. Inspect the saved run/status tables.
# Unresolved runs remain in the results; generation warnings are retained too.
#
# source(system.file("examples", "simulation-results.R", package = "mfrmr"))
# RNGkind("Mersenne-Twister", "Inversion", "Rejection") # recorded generation RNG
# saved <- simulation_result_examples()
# for (name in setdiff(names(saved), "generation_warnings")) {
#   saveRDS(list(result = saved[[name]], generation_warnings = saved$generation_warnings),
#     paste0("simulation-", gsub("_", "-", name), ".rds"), compress = "xz")
# }
simulation_result_examples <- function() {
  warnings <- character()
  results <- withCallingHandlers(list(
    recovery = mfrmr::simulate_mfrm_recovery(
      n_person = 12, n_rater = 2, n_criterion = 2,
      reps = 1, maxit = 30, seed = 123
    ),
    design = mfrmr::evaluate_mfrm_design(
      n_person = c(8, 12), n_rater = 2, n_criterion = 2,
      raters_per_person = 2, reps = 1, maxit = 30, seed = 123
    ),
    diagnostic_screening = mfrmr::evaluate_mfrm_diagnostic_screening(
      design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
      reps = 1, maxit = 30, include_report = TRUE, seed = 123
    ),
    signal_detection = mfrmr::evaluate_mfrm_signal_detection(
      n_person = 8, n_rater = 2, n_criterion = 2, raters_per_person = 2,
      reps = 1, maxit = 30, bias_max_iter = 1, seed = 123
    ),
    population = mfrmr::predict_mfrm_population(
      sim_spec = mfrmr::build_mfrm_sim_spec(
        n_person = 16, n_rater = 3, n_criterion = 2,
        raters_per_person = 2, assignment = "rotating"
      ),
      design = list(person = 18), reps = 1, maxit = 30, seed = 123
    )
  ), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w)))
  })
  c(results, list(generation_warnings = warnings))
}
