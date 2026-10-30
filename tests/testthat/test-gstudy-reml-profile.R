component_profile_source <- function(incomplete = TRUE) {
  x <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  if (incomplete) x <- x[(x$Person+x$Task) %% 3 != 0, ]
  mfrm_g_study(x, c("V","W"), rater = NULL,
    method = if (incomplete) "minque0" else "anova")
}

test_that("component profiles preserve original estimates and separate SE availability", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- component_profile_source()
  before <- serialize(g, NULL)
  w <- c(V=-1,W=1)
  ci <- confint(g,weights=w,method="reml_profile",assumption="normal")
  original <- vapply(g$components,function(m)as.numeric(crossprod(w,m %*% w)),numeric(1))
  expect_equal(ci$Estimate,unname(original),tolerance=1e-12)
  expect_true(any(ci$Estimate<0))
  expect_true(all(ci$Available))
  expect_true(any(ci$IntervalCenter==0))
  expect_true(all(ci$Lower<=ci$IntervalCenter & ci$Upper>=ci$IntervalCenter))
  expect_true(all(is.na(ci$SE) & is.na(ci$RawLower) & is.na(ci$RawUpper)))
  expect_true(all(ci$SEStatus=="Not supplied by profile inversion"))
  expect_identical(serialize(g,NULL),before)
  expect_null(attr(ci,"reml_source")$objective)
  expect_identical(names(attr(ci,"diagnostics")),ci$Source)
  for(i in seq_len(nrow(ci))) {
    z <- attr(ci,"diagnostics")[[i]]
    scale <- attr(ci,"reml_source")$score_scale
    expect_equal(z$bounds*scale*scale,c(ci$Lower[i],ci$Upper[i]))
    expect_true(all(z$scan$Valid))
    for(p in z$profiles) expect_equal(p$components[i],p$target,tolerance=1e-10,ignore_attr=TRUE)
  }
  selected <- confint(g,parm=c("Residual","Person"),weights=w,
    method="reml_profile",assumption="normal")
  expect_identical(selected$Source,c("Residual","Person"))
  expect_identical(names(attr(selected,"diagnostics")),selected$Source)
  expect_equal(selected$Lower,ci$Lower[c(3,1)])
  expect_equal(selected$Upper,ci$Upper[c(3,1)])
})

test_that("component profiles respect composite units, shifts and unused scores", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- component_profile_source(); w <- c(V=-.3,W=.7)
  ci <- confint(g,weights=w,method="reml_profile",assumption="normal")
  x <- g$data[nrow(g$data):1,]; x$Combined <- drop(as.matrix(x[c("V","W")]) %*% w)
  h <- mfrm_g_study(x,"Combined",rater=NULL,method="minque0")
  direct <- confint(h,method="reml_profile",assumption="normal")
  transformed <- confint(g,weights=-3*w,method="reml_profile",assumption="normal")
  for(column in c("Estimate","IntervalCenter","Lower","Upper")) {
    expect_equal(ci[[column]],direct[[column]],tolerance=1e-6)
    expect_equal(transformed[[column]],9*ci[[column]],tolerance=1e-6)
  }
  base <- confint(g,score="V",method="reml_profile",assumption="normal")
  x <- g$data; x$V <- (100-3*x$V)*1e-80; x$W <- x$W*1e80
  h <- mfrm_g_study(x,c("V","W"),rater=NULL,method="minque0")
  tiny <- confint(h,score="V",method="reml_profile",assumption="normal")
  for(column in c("Estimate","IntervalCenter","Lower","Upper"))
    expect_equal(tiny[[column]]/9e-160,base[[column]],tolerance=1e-6)
})

test_that("profile component plots retain points, intervals and replayable diagnostics", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- component_profile_source(); w <- c(V=-1,W=1)
  p <- plot(g,weights=w,intervals=TRUE,method="reml_profile",assumption="normal",draw=FALSE)
  ci <- plot_data(p,"intervals")$Composite; tab <- plot_data(p,"table")
  expect_equal(tab$Estimate,ci$Estimate)
  expect_identical(tab$IntervalCenter,ci$IntervalCenter)
  expect_identical(tab$IntervalAvailable,ci$Available)
  expect_true(all(is.na(tab$SE)))
  expect_match(paste(plot_data(p,"notes")$Text,collapse=" "),"Points remain original")
  f <- tempfile(fileext=".rds"); on.exit(unlink(f),add=TRUE)
  saveRDS(list(plot=p,intervals=ci),f); restored <- readRDS(f)
  expect_identical(restored$plot,p); expect_identical(restored$intervals,ci)
  local_mocked_bindings(.mfrm_mvgt_profile_source=function(...)stop("Unexpected refit"),
    confint.mfrm_multivariate_gstudy=function(...)stop("Unexpected interval calculation"))
  if(requireNamespace("ggplot2",quietly=TRUE)) {
    gg <- as_ggplot(restored$plot)
    expect_identical(plot_data(gg,"table"),tab)
    expect_identical(plot_data(gg,"intervals")$Composite,ci)
    expect_null(gg$labels$title); expect_null(gg$labels$caption)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
})

test_that("component profiles reject unsupported inputs and retain degenerate outcomes", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- component_profile_source()
  expect_error(confint(g,method="reml_profile"),"assumption")
  expect_error(confint(g,parm="absent",method="reml_profile",assumption="normal"),"parm")
  expect_error(confint(g,weights=c(V=1),method="reml_profile",assumption="normal"),"weights")
  expect_error(confint(g,score=c("V","V"),method="reml_profile",assumption="normal"),"distinct")
  expect_error(confint(g,method="reml_profile",assumption="normal",draws=1000),"without draws")
  expect_error(confint(g,method="auto",assumption="normal"),"balanced ANOVA")
  broken <- g; broken$data <- NULL
  expect_error(confint(broken,method="reml_profile",assumption="normal"),"used ratings")
  broken <- g; names(broken$components)[1] <- "Unknown"
  expect_error(confint(broken,method="reml_profile",assumption="normal"),"identities")
  x <- g$data; x$V <- 1
  zero <- mfrm_g_study(x,"V",rater=NULL,method="minque0")
  z <- confint(zero,method="reml_profile",assumption="normal")
  expect_true(all(z$Estimate==0)); expect_false(any(z$Available))
  expect_true(all(is.na(z$Lower) & is.na(z$Upper) & is.na(z$SE)))
  expect_true(all(grepl("Degenerate",z$Status)))
  p <- plot(zero,intervals=TRUE,method="reml_profile",assumption="normal",draw=FALSE)
  expect_equal(nrow(plot_data(p,"unavailable_intervals")),3)
})
