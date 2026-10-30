resampling_unit_fixture <- function() {
  data.frame(Person = rep(c("a", "b", "c"), c(2, 3, 1)),
    Rater = c("R1", "R1", "R2", "R3", "R3", "R1"),
    Item = c("I1", "I2", "I1", "I1", "I2", "I1"), Score = c(0, 1, 1, 2, 1, 0),
    Region = rep(c("A", "B", "C"), c(2, 3, 1)), stringsAsFactors = FALSE)
}

test_that("person and stratum identities cannot collapse in resampling", {
  dat <- resampling_unit_fixture()
  dat$Person <- rep(c(1, 1 + 1e-15, 2), c(2, 3, 1))
  expect_error(build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score"),
    "distinct text labels")
  dat <- resampling_unit_fixture()
  dat$.mfrm_original_person <- "user data"
  expect_error(build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score"), "reserved")
  dat <- data.frame(Person = letters[1:6], Rater = "R", Item = "I", Score = 0,
    S1 = c("A | B", "A", NA, "(missing)", "", "(empty)"),
    S2 = c("C", "B | C", "C", "C", "C", "C"))
  spec <- build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score",
    strata = c("S1", "S2"), reps = 1, sample_fraction = 1, seed = 11)
  expect_equal(nrow(spec$stratum_overview), 6)
  draws <- draw_mfrm_resamples(spec)
  expect_equal(draws$manifest$StrataRepresented, 6)
  expect_equal(sort(draws$samples[[1]]$Person), letters[1:6])
  expect_equal(sort(draws$stratum_manifest$SelectedPersonClusters), rep(1L, 6))
  names(dat)[1:2] <- c("Rows", "n_strata")
  unusual <- build_mfrm_resampling_spec(dat, "Rows", c("n_strata", "Item"), "Score",
    strata = c("S1", "S2"), reps = 1, sample_fraction = 1, seed = 11)
  expect_equal(sum(unusual$stratum_overview$Rows), nrow(dat))
  expect_equal(draw_mfrm_resamples(unusual)$manifest$OriginalPersons, 6)
  names(dat)[names(dat) == "S1"] <- "sep"
  expect_equal(nrow(build_mfrm_resampling_spec(dat, "Rows", c("n_strata", "Item"),
    "Score", strata = c("sep", "S2"))$stratum_overview), 6)
  dat <- data.frame(Person = rep("p", 2), Rater = "R", Score = 1, Group = c(NA_real_, NaN))
  expect_equal(nrow(build_mfrm_resampling_spec(dat, "Person", "Rater", "Score",
    strata = "Group")$person_table), 1)
})

test_that("top-up adds a cluster only while a facet level remains missing", {
  dat <- resampling_unit_fixture()
  spec <- build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score",
    preserve_facets = "Rater", design = "stratified_bootstrap", seed = 1)
  selected <- spec$person_table[spec$person_table$Person == "a", ]
  selected$.mfrm_topup <- FALSE
  selected$.mfrm_topup_reason <- NA_character_
  topped <- mfrmr:::resampling_topup_units(spec, selected)
  # Person b supplies both missing raters; it must be added only once.
  expect_equal(topped$Person, c("a", "b"))
  expect_equal(sum(topped$.mfrm_topup), 1)
})

test_that("stratum counts describe the complete returned draw after top-up", {
  dat <- resampling_unit_fixture()
  spec <- build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score",
    preserve_facets = "Rater", strata = NULL, sample_n = 1, seed = 24, reps = 8)
  draws <- draw_mfrm_resamples(spec)
  expect_true(any(draws$manifest$TopupPersonClusters > 0))
  for (i in seq_along(draws$samples)) {
    sample <- draws$samples[[i]]
    sm <- draws$stratum_manifest[draws$stratum_manifest$Replicate == i, ]
    expect_equal(sum(sm$SelectedPersonClusters), length(unique(sample$.mfrm_draw_unit)))
    expect_equal(sum(sm$TopupPersonClusters), draws$manifest$TopupPersonClusters[i])
    expect_equal(sum(sm$BasePersonClusters), 1)
    expect_equal(sum(sm$Rows), nrow(sample))
  }
  expect_true(all(grepl("coverage_augmented", draws$manifest$SamplingScheme)))
  legacy <- draws
  legacy$stratum_manifest$SelectedPersonClusters <- legacy$stratum_manifest$BasePersonClusters
  legacy$stratum_manifest$BasePersonClusters <- legacy$stratum_manifest$TopupPersonClusters <- NULL
  old_summary <- summary(legacy)
  affected <- old_summary$full_stratum_manifest$Replicate %in%
    legacy$manifest$Replicate[legacy$manifest$TopupPersonClusters > 0]
  expect_true(all(is.na(old_summary$full_stratum_manifest$SelectedPersonClusters[affected])))
  expect_equal(old_summary$full_stratum_manifest$BasePersonClusters, rep(1, 8))
  expect_true(any(grepl("Legacy stratum counts", old_summary$terminology)))
})

test_that("bootstrap draws preserve unequal-sized clusters and disclose the sampling law", {
  dat <- resampling_unit_fixture()
  spec <- build_mfrm_resampling_spec(dat, "Person", c("Rater", "Item"), "Score",
    design = "stratified_bootstrap", sample_fraction = 1, topup_preserve_facets = FALSE,
    reps = 12, seed = 19)
  set.seed(3)
  rng <- .Random.seed
  draws <- draw_mfrm_resamples(spec)
  expect_identical(.Random.seed, rng)
  expect_identical(spec$sampling_scheme, "stratified_person_bootstrap")
  expect_equal(draws$manifest$PersonClusters, rep(3L, 12))
  expect_true(any(draws$manifest$OriginalPersons < 3))
  for (sample in draws$samples) {
    for (block in split(sample, sample$.mfrm_draw_unit)) {
      original <- dat[dat$Person == block$.mfrm_original_person[1], ]
      expect_equal(block[c("Rater", "Item", "Score", "Region")],
        original[c("Rater", "Item", "Score", "Region")], ignore_attr = TRUE)
      expect_equal(length(unique(block$Person)), 1)
    }
  }
  small <- spec
  small$sample_fraction <- .5
  # Old or edited specifications are rebuilt from their declared inputs.
  smaller <- draw_mfrm_resamples(small, keep_data = FALSE)
  expect_true(all(smaller$manifest$PersonClusters == 2))
  expect_true(all(smaller$manifest$SamplingScheme == "stratified_person_resample_with_replacement"))
  expect_identical(summary(draws)$overview$SamplingScheme, "stratified_person_bootstrap")
})
