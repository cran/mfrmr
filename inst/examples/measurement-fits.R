# Regenerate synthetic fitted inputs for interval and export help examples.
# Sourcing this file defines a function; it does not fit any models.
# Interval calculations, scoring and export still run in the help examples.
#
# source(system.file("examples", "measurement-fits.R", package = "mfrmr"))
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")
# saved <- measurement_fit_examples()
# saveRDS(saved, "measurement-fits.rds", compress = "xz")
measurement_fit_examples <- function() {
  ratings <- mfrmr::load_mfrmr_data("example_core")
  gpcm <- mfrmr::fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
    model = "GPCM", method = "MML", step_facet = "Criterion",
    slope_facet = "Rater", rating_min = 1, rating_max = 4,
    category_policy = "preserve", quad_points = 61)
  # The automatic ~1 formula needs no fitting-session environment when saved.
  for (field in c("formula", "terms")) {
    environment(gpcm$population[[field]]) <- baseenv()
    environment(gpcm$config$population_spec[[field]]) <- baseenv()
  }

  # Same 60-person, six-item binary case as the ConQuest bundle generator.
  latent <- withr::with_seed(20240403, {
    persons <- sprintf("P%03d", 1:60)
    criteria <- paste0("C", 1:6)
    x <- seq(-1.5, 1.5, length.out = 60)
    theta <- .25 + .80 * x + stats::rnorm(60, sd = sqrt(.36))
    difficulty <- seq(-1, 1, length.out = 6)
    data <- expand.grid(Person = persons, Criterion = criteria,
      stringsAsFactors = FALSE)
    eta <- theta[match(data$Person, persons)] - difficulty[match(data$Criterion, criteria)]
    data$Score <- stats::rbinom(nrow(data), 1, stats::plogis(eta))
    list(data = data, persons = data.frame(Person = persons, X = x))
  })
  conquest <- function(points, iterations, tolerance) {
    mfrmr::fit_mfrm(latent$data, "Person", "Criterion", "Score",
      method = "MML", model = "RSM",
      population_formula = stats::as.formula("~ X", env = baseenv()),
      person_data = latent$persons, person_id = "Person", population_policy = "error",
      quad_points = points, maxit = iterations, reltol = tolerance, mml_engine = "direct")
  }
  list(gpcm = gpcm, conquest = conquest(31, 2000, 1e-10),
    conquest_quick = conquest(7, 40, 1e-9))
}
