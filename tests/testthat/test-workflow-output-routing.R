test_that("component selection cannot turn unsupported views into misleading bars", {
  # Installing ggplot2 cannot enable an unsupported view. Check the capability
  # before the optional dependency, including when a component is requested.
  local_mocked_bindings(.require_mfrmr_ggplot2 = function() {
    stop("optional ggplot2 dependency required", call. = FALSE)
  }, .package = "mfrmr")
  # These tables have enough numeric/label columns for the old generic fallback
  # to draw a bar chart while discarding axes, groups, intervals or denominators.
  views <- c("screening_performance",
    "multivariate_d_comparison", "d_study")
  table <- data.frame(ID = c("A", "B"), Estimate = c(-.5, .5),
    Lower = c(-1, 0), Upper = c(0, 1), OtherAxis = c(2, -2))
  for (view in views) {
    x <- new_mfrm_plot_data(view, list(table = table, matrix = matrix(1:4, 2)))
    for (component in list(NULL, "table", "matrix")) {
      expect_error(as_ggplot(x, component = component), "plot_data")
    }
    expect_identical(plot_data(x, "table"), table)
  }
  pca <- mfrm_pca(mfrm_features(data.frame(ID = letters[1:6],
    Experience = c(1,2,4,7,9,10), Training = c(2,6,3,8,4,10)), "ID", c("Experience", "Training")))
  expect_error(as_ggplot(pca, type = "scores", component = "table"), "optional ggplot2 dependency required")
  expect_named(plot_data(pca, type = "scores")$table, c("ID", "PC1", "PC2", "Cluster"))
  pooled <- structure(list(table = transform(table, Target = ID, Status = "available"),
    settings = list(facet = "Rater", imputations = 20, ci_level = .95)), class = "mfrm_pooled")
  expect_error(as_ggplot(pooled, component = "table"), "optional ggplot2 dependency required")
  expect_identical(plot_data(pooled)$table$Lower, table$Lower)
  custom <- new_mfrm_plot_data("custom", list(table = table))
  expect_error(as_ggplot(custom), "optional ggplot2 dependency required")
})

test_that("D-study conversion keeps both panels and availability with explicit series", {
  skip_if_not_installed("ggplot2")
  d <- expand.grid(Person = paste0("P",1:8), Rater = paste0("R",1:3), Task = paste0("T",1:2))
  d$Score <- withr::with_seed(264, rnorm(8)[match(d$Person,unique(d$Person))] +
    rnorm(3)[match(d$Rater,unique(d$Rater))] + rnorm(nrow(d)))
  g <- mfrm_g_study(d,"Score")
  plans <- mfrm_d_study(g,data.frame(Raters=c(1,2,3),Tasks=2))
  for (type in c("coefficients", "sem")) {
    payload <- plot(plans,type=type,x_var="Raters",draw=FALSE)
    base <- as_ggplot(payload)
    explicit <- as_ggplot(payload,component="series")
    expect_identical(explicit$data,base$data)
    expect_identical(explicit$labels,base$labels)
    expect_equal(nlevels(explicit$data$Panel),2)
    expect_identical(explicit$data$X,payload$data$series$X)
    expect_identical(is.na(explicit$data$Value), payload$data$series$Status != "Available" |
      !is.finite(payload$data$series$Value))
    expect_error(as_ggplot(payload,component="table"), "component = 'series'", fixed=TRUE)
    expect_error(as_ggplot(plans,type=type,component="table"), "plot_data")
    expect_error(as_ggplot(payload,level=.8), "empty")
    expect_silent(ggplot2::ggplotGrob(explicit))
  }
  # An unclassified custom tabular payload still has its documented fallback.
  custom <- new_mfrm_plot_data("custom",list(table=data.frame(Label="A",Value=1)))
  expect_s3_class(as_ggplot(custom,component="table"),"ggplot")
})

test_that("the output guide connects dedicated analysis families without claiming parity", {
  for (scope in c("features", "imputation", "gtheory")) {
    row <- mfrmr_output_guide(scope)
    expect_identical(row$Scope,scope)
    expect_true(row$Question %in% mfrmr_output_guide()$Question)
    expect_match(row$NextStep,"saveRDS",fixed=TRUE)
    expect_match(row$DecisionBoundary,"not mfrm_results()",fixed=TRUE)
    expect_match(row$DecisionBoundary,"as_ggplot",fixed=TRUE)
  }
  expect_match(mfrmr_output_guide("features")$MainFunction,"mfrm_pca",fixed=TRUE)
  expect_match(mfrmr_output_guide("imputation")$MainFunction,"pool_mfrm_imputed",fixed=TRUE)
  expect_match(mfrmr_output_guide("gtheory")$MainFunction,"mfrm_d_study()",fixed=TRUE)
  expect_false(any(mfrmr_output_guide("gpcm")$Scope %in% c("features","gtheory")))
  expect_equal(nrow(mfrmr_output_guide("beginner")),6)
})


test_that("PCA grouped scores and scree views retain non-colour distinctions", {
  features <- mfrm_features(data.frame(ID=letters[1:8],
    Experience=c(1,2,3,5,9,10,12,14), Training=c(2,5,4,6,9,8,13,12),
    Workload=c(20,40,25,55,60,45,80,95)), "ID", c("Experience","Training","Workload"))
  pca <- mfrm_pca(features,components=2)
  groups <- mfrm_cluster_kmeans(pca,3,seed=17,silhouette=FALSE)
  for (preset in c("standard","monochrome")) {
    payload <- plot(pca,type="scores",groups=groups,labels=FALSE,preset=preset,draw=FALSE)$data
    expect_equal(length(unique(payload$encoding$Shape)),3)
    expect_setequal(payload$encoding$Cluster,groups$membership$Cluster)
    expect_identical(payload$table$Cluster,groups$membership$Cluster)
    expect_identical(payload$table$PC1,pca$scores$PC1)
    expect_identical(payload$table$PC2,pca$scores$PC2)
  }
  # Observe only mfrmr's graphics imports; preserve the standard S3 registry.
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  plots <- points <- legends <- list()
  local_mocked_bindings(
    plot = function(...) {
      args <- list(...)
      if (length(args) && is.numeric(args[[1]])) plots[[length(plots)+1L]] <<- args
      graphics::plot(...)
    },
    points = function(...) {
      points[[length(points)+1L]] <<- list(...)
      graphics::points(...)
    }, .package = "mfrmr")
  original_legend <- graphics::legend
  local_mocked_bindings(
    legend=function(...) {legends[[length(legends)+1L]] <<- list(...); original_legend(...)},
    .package="mfrmr")
  drawn <- plot(pca,type="scores",groups=groups,preset="monochrome",labels=FALSE)
  expect_equal(plots[[1]]$pch,drawn$data$encoding$Shape[match(drawn$data$table$Cluster,drawn$data$encoding$Cluster)])
  expect_identical(legends[[1]]$pch,drawn$data$encoding$Shape)
  points <- list() # Separate score-plot legend points from the scree renderer.
  plot(pca,preset="monochrome")
  expect_identical(plots[[2]]$type,"l")
  expect_identical(points[[1]]$pch,ifelse(pca$variance$Retained,16,1))
  expect_identical(legends[[2]]$pch,c(16,1))
})
