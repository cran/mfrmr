# FACETS 4.5.0 hidden batch outputs, directly executed 2026-10-24.
# The reference files contain literal output, not mfrmr-generated expectations.
xtreme_fit <- function(case, model = "RSM", anchored_person = FALSE, cell_weight = NULL,
                       weight_scale = 1) {
  rows <- read.csv(test_path("fixtures", "facets-xtreme-input.csv"))
  rows <- rows[rows$Case == case, ]
  if (!is.null(cell_weight)) rows$Weight[rows$Person==4 & rows$Rater==2 & rows$Task==3] <- cell_weight
  rows$Weight <- rows$Weight * weight_scale
  anchors <- data.frame(Facet = c(rep("Rater", 2), rep("Task", 3)),
    Level = as.character(c(1, 2, 1, 2, 3)), Anchor = c(-.4, .4, -.8, 0, .8))
  if (case == "bias_weighted_anchors") anchors <- rbind(anchors,
    data.frame(Facet = "Person", Level = as.character(1:6), Anchor = ((1:6)-3.5)*.3))
  if (anchored_person) anchors <- rbind(anchors,
    data.frame(Facet = "Person", Level = "5", Anchor = -1))
  fit <- suppressWarnings(fit_mfrm(rows, "Person", c("Rater", "Task"), "Score",
    weight = "Weight", method = "JML", model = model,
    step_facet = if (model == "PCM") "Task" else NULL, anchors = anchors, maxit = 200))
  list(fit = fit, dx = diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy"))
}

test_that("public total-score extrapolation matches native binary counts and weights", {
  refs <- read.csv(test_path("fixtures", "facets-xtreme-reference.csv"))
  for (case in c("binary", "binary_half", "weighted")) {
    x <- xtreme_fit(case); before <- serialize(x, NULL)
    ref <- refs[refs$Case == case, ]
    b <- fair_average_table(x$fit, x$dx, xtreme = ref$Correction[1], xtreme_method = "facets")
    p <- b$raw_by_facet$Person
    ref <- ref[match(p$Level, as.character(ref$Person)), ]
    # Native six-decimal measures/SEs have small iterative extrapolation
    # residuals beyond print rounding. Independent score-equation tests below
    # require much tighter accuracy from this implementation.
    expect_equal(p$Measure, ref$Measure, tolerance = 2e-5)
    expect_equal(p$ModelSE, ref$SE, tolerance = 2e-5)
    expect_lt(max(abs(p$FairM - ref$FairM)), .0051)
    expect_identical(serialize(x, NULL), before)
    expect_true(all(is.infinite(p$PrimaryMeasure[p$ExtremeAdjustment > 0])))
    if (case == "weighted") {
      expect_equal(p$ExtremeAdjustment[match(c("1", "2", "3"), p$Level)], c(.6, .15, .15))
      # The .25 observation on another Person must not set this correction.
      expect_false(any(p$ExtremeAdjustment == .3*.25))
    }
    saved <- unserialize(serialize(b, NULL))
    expect_true(any(grepl("total-score extrapolation", summary(saved)$notes)))
    plot <- plot(x$fit, type = "fair_average", diagnostics = x$dx, facet = "Person",
      xtreme = ref$Correction[1], xtreme_method = "facets", show_ci = TRUE, draw = FALSE)
    expect_true(all(is.na(plot$data$data$CI_SE)))
    expect_true(any(grepl("extrapolation", plot_data(plot, "notes")$Text)))
    broken <- saved; broken$settings$xtreme_method <- NULL
    expect_error(summary(broken), "extreme-score convention")
  }
})

test_that("PCM extrapolation uses actual scales and preserves anchors and primary fits", {
  x <- xtreme_fit("pcm_weighted", model = "PCM")
  pars <- mfrmr:::expand_params(x$fit$opt$par, mfrmr:::build_param_sizes(x$fit$config), x$fit$config)
  # Separate kernel check against native fixed thresholds; source fitting
  # remains independent and is not claimed to estimate these fixed values.
  native_pars <- pars
  native_pars$steps_mat <- rbind(c(-.5,.5), c(-1,1), c(-1.5,1.5))
  native <- mfrmr:::.mfrm_extreme_report(x$fit, x$dx$obs, "Person", as.character(1:3),
    rep("", 3), .3, native_pars)
  ref <- read.csv(test_path("fixtures", "facets-xtreme-reference.csv"))
  ref <- ref[ref$Case == "pcm_weighted" & ref$Person %in% 1:3, ]
  expect_equal(native$measure, ref$Measure, tolerance = 2e-5)
  expect_equal(native$se, ref$SE, tolerance = 2e-5)
  b <- fair_average_table(x$fit, x$dx, xtreme = .3, xtreme_method = "facets")
  p <- b$raw_by_facet$Person
  for (i in which(p$ExtremeAdjustment > 0)) {
    z <- x$dx$obs[as.character(x$dx$obs$Person) == p$Level[i], ]
    task <- match(as.character(z$Task), x$fit$prep$levels$Task)
    rater <- match(as.character(z$Rater), x$fit$prep$levels$Rater)
    eta <- p$Measure[i] - pars$facets$Rater[rater] - pars$facets$Task[task]
    # Independent three-category exp-polynomial; no package probability helper.
    logits <- cbind(0, eta - pars$steps_mat[task,1], 2*eta - rowSums(pars$steps_mat[task,,drop=FALSE]))
    prob <- exp(logits - apply(logits, 1, max)); prob <- prob/rowSums(prob)
    expected <- prob[,2] + 2*prob[,3]
    variance <- rowSums(prob * (matrix(0:2, nrow(z), 3, byrow=TRUE) - expected)^2)
    target <- if (all(z$Observed == 0)) p$ExtremeAdjustment[i] else 2*sum(z$Weight)-p$ExtremeAdjustment[i]
    expect_equal(sum(z$Weight*expected), target, tolerance = 1e-9)
    expect_equal(p$ModelSE[i], 1/sqrt(sum(z$Weight*variance)), tolerance = 1e-9)
  }
  scaled <- fair_average_table(x$fit, x$dx, xtreme=.3, xtreme_method="facets", umean=50, uscale=-10)$raw_by_facet$Person
  at <- match(p$Level, scaled$Level)
  expect_equal(scaled$Measure[at], 50-10*p$Measure)
  expect_equal(scaled$ModelSE[at], 10*p$ModelSE)
  expect_equal(scaled$FairM[at], p$FairM)
  anchor <- xtreme_fit("weighted", anchored_person=TRUE)
  b <- fair_average_table(anchor$fit, anchor$dx, xtreme=.3, xtreme_method="facets")
  a <- subset(b$raw_by_facet$Person, Level == "5")
  expect_equal(a$Measure, -1); expect_equal(a$PrimaryMeasure, -1)
  expect_equal(a$ExtremeAdjustment, 0)
})

test_that("bias extrapolation matches native points while fractional tails remain distinct", {
  x <- xtreme_fit("bias_weighted_anchors")
  b <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5)
  ref <- read.csv(test_path("fixtures", "facets-xtreme-bias-reference.csv"))
  z <- b$table; ref <- ref[match(paste(z$FacetA_Level,z$FacetB_Level),paste(ref$Rater,ref$Task)), ]
  expect_lt(max(abs(z$`Bias Size`-ref$Bias)), .0011)
  expect_lt(max(abs(z$`S.E.`-ref$SE)), .0011)
  integer <- z$`Observd Count` == floor(z$`Observd Count`)
  expect_lt(max(abs(z$`Prob.`[integer]-ref$p[integer])), .00011)
  expect_equal(z$DFDisplay, ref$DF)
  fractional <- which(z$DFUnrounded != floor(z$DFUnrounded))
  expect_length(fractional, 1)
  expect_equal(z$DFUnrounded[fractional], 12.25)
  expect_equal(z$`d.f.`[fractional], 12.25)
  expect_equal(z$`Prob.`, 2*pt(-abs(z$t), df=z$`d.f.`), tolerance=1e-14)
  adjusted <- z$BiasStatus == "extreme_adjusted"
  expect_equal(sum(adjusted), 2)
  expect_true(all(is.infinite(z$PrimaryBias[adjusted])))
  expect_true(all(is.finite(z$`S.E.`[adjusted])))
  expect_false(any(z$SupportsFormalInference))
  expect_equal(z$AfterBiasExpectedScore[adjusted] - z$`Observd Score`[adjusted],
    -sign(z$PrimaryBias[adjusted])*z$ExtremeAdjustment[adjusted], tolerance = 1e-9)
  raw <- estimate_bias(x$fit, x$dx, "Rater", "Task")
  expect_true(all(is.na(raw$table$`S.E.`[adjusted])))
  limited <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5, max_abs=.1)
  expect_false(any(limited$table$BiasStatus == "extreme_adjusted"))
  expect_no_error(summary(unserialize(serialize(b, NULL))))
  expect_no_error(plot(b, draw=FALSE))
  expect_no_error(build_fixed_reports(b))
})

test_that("fractional native df display is separate from the Student tail calculation", {
  ref <- read.csv(test_path("fixtures", "facets-xtreme-fractional-reference.csv"))
  for (i in seq_len(nrow(ref))) {
    x <- xtreme_fit("bias_weighted_anchors", cell_weight=ref$Weight[i])
    z <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5)$table
    z <- subset(z, FacetA_Level=="2" & FacetB_Level=="3")
    expect_equal(z$`Observd Count`, ref$Count[i])
    expect_lt(abs(z$`Bias Size`-ref$Bias[i]), .0011)
    expect_lt(abs(z$`S.E.`-ref$SE[i]), .0011)
    expect_equal(z$DFDisplay, ref$DF[i])
    expect_equal(z$`d.f.`, ref$Count[i]-1)
    expect_equal(z$`Prob.`, 2*pt(-abs(z$t), ref$Count[i]-1), tolerance=1e-14)
    if (ref$Count[i] == floor(ref$Count[i])) {
      expect_lt(abs(z$`Prob.`-ref$p[i]), .00011)
    } else {
      # Preserve the known disagreement, rather than fitting a correction to it.
      expect_gt(abs(z$`Prob.`-ref$p[i]), .002)
    }
  }
})

test_that("extrapolation excludes infinite calibration and refuses coupled models", {
  x <- xtreme_fit("weighted")
  b <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5, omit_extreme=FALSE)
  expect_equal(sum(b$table$`Observd Count`),
    sum(x$dx$obs$Weight[as.character(x$dx$obs$Person) %in% c("4","6")]))
  group <- x$fit; group$config$anchor_summary$GroupAnchors[1] <- 1L
  expect_error(fair_average_table(group, xtreme=.3, xtreme_method="facets"), "group-mean")
  expect_error(estimate_bias(group, x$dx, "Rater", "Task", xtreme=.5), "group-mean")
  corrected <- x$fit; corrected$config$method <- "JML_BC"
  expect_error(fair_average_table(corrected, xtreme=.3, xtreme_method="facets"), "ordinary additive")
  expect_error(estimate_bias(corrected, x$dx, "Rater", "Task", xtreme=.5), "ordinary additive")
})

test_that("retaining a fixed extreme does not retain other unanchored extremes", {
  x <- xtreme_fit("weighted", anchored_person = TRUE)
  kept <- x$dx$obs[as.character(x$dx$obs$Person) %in% c("4", "5", "6"), ]
  for (convention in c("display", "facets")) {
    result <- fair_average_table(x$fit, x$dx, totalscore = FALSE,
      xtreme = .3, xtreme_method = convention)
    for (facet in c("Rater", "Task")) {
      tab <- result$raw_by_facet[[facet]]
      expect_equal(tab$TotalCount, as.numeric(table(kept[[facet]])[tab$Level]))
      expect_equal(tab$WeightdCount, as.numeric(tapply(kept$Weight,kept[[facet]],sum)[tab$Level]))
    }
    fixed <- subset(result$raw_by_facet$Person, Level == "5")
    expect_equal(fixed$Measure, -1)
    expect_equal(fixed$ExtremeAdjustment, 0)
  }
})

test_that("default Person bias retains direct extreme anchors without qualifying their cells", {
  for (model in c("RSM", "PCM")) {
    x <- xtreme_fit(if (model == "RSM") "weighted" else "pcm_weighted",
      model = model, anchored_person = TRUE)
    before <- serialize(x, NULL)
    for (facet in c("Rater", "Task")) {
      raw <- suppressWarnings(estimate_bias(x$fit, x$dx, "Person", facet))
      reference <- suppressWarnings(estimate_bias(x$fit, x$dx, "Person", facet,
        omit_extreme = FALSE))
      fixed <- subset(raw$table, FacetA_Level == "5")
      expected <- subset(reference$table, FacetA_Level == "5")
      expect_equal(nrow(fixed), length(unique(x$dx$obs[[facet]][x$dx$obs$Person == "5"])))
      expect_equal(fixed$`Observd Count`, expected$`Observd Count`)
      expect_equal(fixed$`Bias Size`, expected$`Bias Size`)
      expect_true(all(grepl("search_boundary", fixed$BiasStatus)))
      expect_true(all(is.na(fixed$`S.E.`) & is.na(fixed$`Prob.`)))
      expect_false(any(fixed$SupportsFormalInference))
      expect_false(any(raw$table$FacetA_Level %in% c("1", "2", "3")))
      expect_identical(summary(unserialize(serialize(raw, NULL))), summary(raw))
      expect_no_error(plot(raw, draw = FALSE))
    }
    collection <- suppressWarnings(estimate_all_bias(x$fit, x$dx,
      pairs = list(c("Person", "Rater"))))
    expect_equal(sum(collection$by_pair[[1]]$table$FacetA_Level == "5"),
      length(unique(x$dx$obs$Rater[x$dx$obs$Person == "5"])))
    expect_identical(serialize(x, NULL), before)
  }
})

test_that("fixed extreme main facets survive and inward correction can reverse the observed sign", {
  d <- rbind(data.frame(Person = "p", Rater = rep(c("a", "b"), each = 3), Score = 0),
    data.frame(Person = "q", Rater = c("a", "b"), Score = c(0, 1)))
  anchors <- data.frame(Facet = c("Person", "Person", "Rater", "Rater"),
    Level = c("p", "q", "a", "b"), Anchor = c(qlogis(.01), 0, 0, 0))
  for (model in c("RSM", "PCM")) {
    fit_warnings <- character()
    f <- withCallingHandlers(fit_mfrm(d, "Person", "Rater", "Score", method = "JML", model = model,
      step_facet = if (model == "PCM") "Rater" else NULL, anchors = anchors,
      rating_min = 0, rating_max = 1, category_policy = "preserve",
      repeated_ratings = "independent"), warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    expect_true(any(grepl("Category support", fit_warnings)))
    expect_true(any(grepl("Boundary-constant", fit_warnings)))
    dx <- diagnose_mfrm(f, residual_pca = "none")
    raw <- suppressWarnings(estimate_bias(f, dx, "Person", "Rater"))$table
    adjusted <- suppressWarnings(estimate_bias(f, dx, "Person", "Rater", xtreme = .5))$table
    expect_equal(nrow(raw), 4)
    raw <- subset(raw, FacetA_Level == "p")
    adjusted <- subset(adjusted, FacetA_Level == "p")
    expect_equal(nrow(raw), 2)
    expect_equal(raw$`Observd Count`, c(3, 3))
    expect_true(all(raw$`Bias Size` < 0 & is.na(raw$`S.E.`)))
    expect_equal(adjusted$`Bias Size`, rep(qlogis(.5/3) - qlogis(.01), 2), tolerance = 1e-9)
    expect_equal(adjusted$AfterBiasExpectedScore, c(.5, .5), tolerance = 1e-9)
    expect_true(all(adjusted$PrimaryBias == -Inf))
    expect_true(all(adjusted$`Obs-Exp Average` < 0 & adjusted$`Bias Size` > 0))
    expect_true(all(grepl("reverse the sign", adjusted$ScreenReason)))
    expect_false(any(adjusted$SupportsFormalInference))
  }
})

test_that("fractional bias tails use the continuous density and retain their reference", {
  refs <- read.csv(test_path("fixtures", "facets-bias-fractional-tails.csv"))
  # Public-route checks cover the integer boundary, odd fractional branch,
  # low fractional df, and the df=1 Cauchy special case. The remaining native
  # rows are retained evidence from independently solved scalar cell equations.
  refs <- refs[refs$Case %in% c("near_integer", "half", "df_1.25", "df_1"), ]
  for (i in seq_len(nrow(refs))) {
    ref <- refs[i, ]
    x <- xtreme_fit("bias_weighted_anchors", cell_weight=ref$CellWeight,
      weight_scale=ref$WeightScale)
    result <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5)
    row <- subset(result$table, FacetA_Level=="2" & FacetB_Level=="3")
    expect_equal(row$`Observd Count`, ref$Count, tolerance=1e-13)
    expect_lt(abs(row$`Bias Size`-ref$NativeBias), .00011)
    expect_lt(abs(row$`S.E.`-ref$NativeSE), .0011)
    expect_equal(row$DFUnrounded, ref$Count-1, tolerance=1e-13)
    df <- ref$Count-1
    norm <- lgamma((df+1)/2)-lgamma(df/2)-.5*log(pi*df)
    tail <- 2*integrate(function(t) exp(norm-(df+1)/2*log1p(t*t/df)),
      abs(row$t), Inf, rel.tol=1e-11, abs.tol=1e-12)$value
    expect_equal(row$`Prob.`, tail, tolerance=1e-10)
    expect_identical(row$ProbabilityStatus, "available_conditional")
    expect_match(row$ProbabilityReference, "continuous count-minus-one")
    expect_match(row$ProbabilityReason, "screening approximation")
    if (ref$Case=="near_integer") {
      expect_lt(row$`Prob.`, .05); expect_gt(ref$NativeP, .05)
    }
    if (ref$Case=="df_1.25") expect_gt(abs(row$`Prob.`-ref$NativeP), .24)
    if (ref$Case=="df_1") expect_lt(abs(row$`Prob.`-ref$NativeP), .00006)
    saved <- unserialize(serialize(result,NULL))
    expect_identical(saved$table$ProbabilityStatus, result$table$ProbabilityStatus)
    expect_identical(summary(saved), summary(result))
    expect_no_error(build_fixed_reports(saved))
  }
})

test_that("nonpositive bias reference df is explicit without discarding finite estimates", {
  for (count in c(.5,1)) {
    x <- xtreme_fit("bias_weighted_anchors", weight_scale=count/13.25)
    result <- estimate_bias(x$fit, x$dx, "Rater", "Task", xtreme=.5)
    row <- subset(result$table, FacetA_Level=="2" & FacetB_Level=="3")
    expect_equal(row$DFUnrounded, count-1, tolerance=1e-13)
    expect_true(is.finite(row$`Bias Size`)); expect_true(is.finite(row$`S.E.`))
    expect_true(is.finite(row$t))
    expect_true(all(is.na(row[c("d.f.", "DFDisplay", "Prob.")])))
    expect_identical(row$ProbabilityStatus, "unavailable_df")
    expect_match(row$ProbabilityReason, "nonpositive")
    expect_false(row$SupportsFormalInference)
    s <- summary(result)
    expect_true(any(grepl("ProbabilityStatus",s$notes)))
    expect_true("ProbabilityStatus" %in% names(s$top_rows))
    if (count==.5) {
      expect_equal(s$overview$ScreenAvailable,0)
      expect_true(is.na(s$overview$ScreenPositive))
    }
    expect_no_error(build_fixed_reports(result))
    expect_no_error(plot(result, draw=FALSE))
  }
})
