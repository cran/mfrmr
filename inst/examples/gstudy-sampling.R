# Regenerate the incomplete-design results used in the sampling-covariance help.
# Sourcing this file defines a function; it does not start a calculation.
# The observed-design moment calculation can take tens of seconds.
#
# source(system.file("examples", "gstudy-sampling.R", package = "mfrmr"))
# example <- gstudy_sampling_example()
# saveRDS(example, "gstudy-sampling.rds", compress = "xz")
gstudy_sampling_example <- function() {
  data <- new.env(parent = emptyenv())
  utils::data("mfrmr_example_gtheory_crossed", package = "mfrmr", envir = data)
  observed <- data$mfrmr_example_gtheory_crossed[-c(1, 8, 23), ]
  g <- mfrmr::mfrm_g_study(observed,
    c("Content", "Organization"), method = "minque0")
  list(study = g,
    variance = stats::vcov(g, score = "Content", assumption = "normal"),
    covariance = stats::vcov(g, score = c("Content", "Organization"),
      assumption = "normal"))
}
