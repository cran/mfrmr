.mfrm_dataset_names <- c(
  example_core = "mfrmr_example_core", example_bias = "mfrmr_example_bias",
  example_operational = "mfrmr_example_operational", study1 = "ej2021_study1",
  study2 = "ej2021_study2", combined = "ej2021_combined",
  study1_itercal = "ej2021_study1_itercal", study2_itercal = "ej2021_study2_itercal",
  combined_itercal = "ej2021_combined_itercal", example_agreement = "mfrmr_example_agreement",
  example_gtheory_crossed = "mfrmr_example_gtheory_crossed",
  example_gtheory_nested = "mfrmr_example_gtheory_nested",
  example_operational_design = "mfrmr_example_operational_design")





































































































NULL
















































































NULL





































NULL




















































list_mfrmr_data <- function(details = FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("`details` must be TRUE or FALSE.", call. = FALSE)
  }
  keys <- names(.mfrm_dataset_names)
  if (!isTRUE(details)) return(keys)

  catalog <- data.frame(
    Key = keys[seq_len(9L)],
    Rows = c(768L, 384L, 282L, 1842L, 3287L, 5129L, 1842L, 3341L, 5183L),
    Persons = c(48L, 48L, 48L, 307L, 206L, 307L, 307L, 206L, 307L),
    Raters = c(4L, 4L, 6L, 18L, 12L, 18L, 18L, 12L, 18L),
    Criteria = c(4L, 4L, 3L, 3L, 9L, 12L, 3L, 9L, 12L),
    CountBasis = c(
      rep("unique labels", 5L),
      "raw labels; 513 persons and 30 raters when Study-prefixed",
      rep("unique labels", 2L),
      "raw labels; 513 persons and 30 raters when Study-prefixed"
    ),
    PrimaryUse = c(
      "Idealized fast examples",
      "DFF and bias demonstrations with planted effects",
      "Beginner applied workflow",
      "Unequal-workload sparse-design review",
      "Larger sparse-design review",
      "Identity/linking design review; not direct fit",
      "Legacy synthetic variant review",
      "Legacy synthetic variant review",
      "Identity/linking sensitivity review; not direct fit"
    ),
    Design = c(
      "Complete crossing; no planned omissions",
      "Balanced two-rater assignment; planted non-null effects",
      "Connected two-rater assignment; six planned omissions",
      "Two raters per person; highly unequal rater workloads",
      "Two raters per person; incomplete criterion coverage",
      "Overlapping IDs; requires explicit anchors/linking for a common scale",
      "Legacy Study 1 variant; rows and scores can differ",
      "Legacy Study 2 variant; rows and scores can differ",
      "Overlapping IDs; requires explicit anchors/linking for a common scale"
    ),
    Empirical = FALSE,
    stringsAsFactors = FALSE
  )
  extras <- keys[-seq_len(9L)]
  purposes <- c("Rater agreement with repeated contexts", "Crossed multivariate G-Study",
    "Nested multivariate G-Study", "Planned rating assignments; contains no scores")
  designs <- c("Long ratings; Person x Task defines the rating unit",
    "Wide Content and Organization outcomes; crossed raters and tasks",
    "Wide Content and Organization outcomes; raters nested within tasks",
    "Complete planned roster for example_operational")
  extra_rows <- lapply(seq_along(extras), function(i) {
    d <- getExportedValue("mfrmr", .mfrm_dataset_names[[extras[i]]])
    count <- function(column) if (column %in% names(d)) length(unique(d[[column]])) else NA_integer_
    data.frame(Key = extras[i], Rows = nrow(d), Persons = count("Person"), Raters = count("Rater"),
      Criteria = count("Criterion"), CountBasis = "unique labels; absent columns reported as NA",
      PrimaryUse = purposes[i], Design = designs[i], Empirical = FALSE)
  })
  dplyr::bind_rows(c(list(catalog), extra_rows))
}





























































load_mfrmr_data <- function(name = "example_core", details = FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("details must be TRUE or FALSE.", call. = FALSE)
  if (is.null(name) || identical(tolower(name), "list")) return(list_mfrmr_data(details))
  if (details) stop("details = TRUE is only available for the dataset list.", call. = FALSE)
  key <- match.arg(tolower(name), choices = names(.mfrm_dataset_names))
  getExportedValue("mfrmr", .mfrm_dataset_names[[key]])
}














































NULL



































NULL
