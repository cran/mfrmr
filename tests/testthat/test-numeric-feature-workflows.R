numeric_feature_fixture <- function() {
  mfrm_features(data.frame(ID = c("001", "1", "NA", "é", "R|5", "f", "g", "excluded"),
    Years = c(1, 2, 4, 9, 11, 13, 16, NA),
    Hours = c(5, 14, 9, 24, 18, 32, 28, 20),
    Ratings = c(60, 100, 90, 260, 310, 290, 410, 250)),
    "ID", c("Years", "Hours", "Ratings"),
    data.frame(ID = "excluded", Feature = "Years", Reason = "Not recorded"))
}

test_that("PCA preserves weighted geometry, original IDs and variance accounting", {
  x <- numeric_feature_fixture()
  w <- c(Ratings = 2, Years = 4, Hours = 1)
  p <- mfrm_pca(x, weights = w, missing = "omit")
  values <- as.matrix(x$data[1:7, x$features])
  # Independent construction of the stated geometry.
  centered <- sweep(values, 2, colMeans(values))
  divisor <- apply(values, 2, stats::sd)
  z <- sweep(sweep(centered, 2, divisor, "/"), 2, sqrt(w[x$features] / 4), "*")
  scores <- as.matrix(p$scores[1:7, -1])
  expect_equal(unname(stats::dist(scores)), unname(stats::dist(z)), tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_equal(unname(scores %*% t(p$loadings)), unname(z), tolerance = 1e-12)
  expect_equal(p$variance$Variance, sort(eigen(crossprod(z) / 6, symmetric = TRUE)$values, decreasing = TRUE))
  expect_equal(sum(p$variance$Proportion), 1)
  expect_identical(p$scores$ID, x$row_summary$ID)
  expect_true(all(is.na(p$scores[8, -1])))
  expect_identical(p$feature_data$missing, x$missing)
  expect_equal(p$transformation$Divisor, unname(divisor))
  expect_equal(p$settings$weights, w[x$features])
  reduced <- mfrm_pca(x, 1, weights = w, missing = "omit")
  expect_equal(reduced$scores$PC1, p$scores$PC1)
  expect_identical(reduced$variance$Retained, c(TRUE, FALSE, FALSE))
  expect_true(any(abs(as.vector(dist(reduced$scores[1:7, 'PC1', drop = FALSE])) -
    as.vector(dist(scores))) > 1e-6))
  unscaled <- mfrm_pca(x, scale = FALSE, missing = "omit")
  expect_equal(as.vector(dist(unscaled$scores[1:7, -1])), as.vector(dist(values)), tolerance = 1e-10)
  changed_units <- x$data
  changed_units$Years <- 12 * changed_units$Years
  converted <- mfrm_pca(mfrm_features(changed_units, "ID", x$features), weights = w, missing = "omit")
  expect_equal(as.vector(dist(converted$scores[1:7, -1])), as.vector(dist(scores)), tolerance = 1e-12)
  expect_equal(summary(p), p$variance)
  expect_output(print(p), "not ability")
})

test_that("k-means agrees with a reference fit and preserves saved comparison and profiles", {
  skip_if_not_installed("cluster")
  x <- numeric_feature_fixture()
  set.seed(1234)
  before <- .Random.seed
  out <- mfrm_cluster_kmeans(x, 2, missing = "omit", seed = 7, nstart = 30)
  expect_identical(.Random.seed, before)
  z <- scale(x$data[1:7, x$features])
  reference <- with_preserved_rng_seed(7, stats::kmeans(z, 2, nstart = 30, iter.max = 100))
  labels <- out$membership$Cluster[1:7]
  expect_equal(labels, unname(reference$cluster))
  expect_equal(out$centers, reference$centers)
  expect_equal(out$tot.withinss, sum(vapply(seq_len(7), function(i)
    sum((z[i, ] - out$centers[labels[i], ])^2), numeric(1))))
  expect_equal(out$membership$Silhouette[1:7], unname(cluster::silhouette(labels, dist(z))[, 'sil_width']))
  expect_identical(out$membership$ID, x$row_summary$ID)
  expect_true(is.na(out$membership$Cluster[8]))
  expect_null(out$medoids)
  expect_equal(out$profiles$numeric$Mean[out$profiles$numeric$Feature == 'Years'],
    as.numeric(tapply(x$data$Years[1:7], labels, mean)))
  pca <- mfrm_pca(x, missing = "omit")
  pc_groups <- mfrm_cluster_kmeans(pca, 2, seed = 7, nstart = 30)
  file <- tempfile(fileext = '.rds'); on.exit(unlink(file), add = TRUE)
  saveRDS(list(pca = pca, groups = pc_groups), file)
  saved <- readRDS(file)
  expect_identical(saved$pca, pca)
  expect_identical(saved$groups, pc_groups)
  comparison <- mfrm_cluster_compare(list(Features = out, FullPCA = saved$groups,
    PAM = mfrm_cluster(x, 2, missing = 'omit')))
  expect_equal(comparison$comparisons$AdjustedRand[1], 1)
  expect_equal(out$tot.withinss, pc_groups$tot.withinss)
  expect_identical(names(comparison$analysis_summary)[seq_len(11)],
    c('Analysis','Method','Linkage','Imputation','K','Features','Included','Excluded',
      'MinGroupSize','MaxGroupSize','MeanSilhouette'))
  expect_identical(comparison$analysis_summary$Distance, c('Euclidean','Euclidean','Gower'))
  expect_identical(comparison$analysis_summary$Space, c('Numeric features','Principal components','Mixed features'))
  expect_output(print(out), "Euclidean distance, k-means")
  expect_match(plot_data(plot(out, draw = FALSE))$subtitle, 'Euclidean')
  expect_equal(plot_data(plot(out, type='profile', feature='Years', draw=FALSE))$table,
    subset(out$profiles$numeric, Feature == 'Years'))
})

test_that("numeric workflows refuse invalid geometry and optimization results", {
  x <- numeric_feature_fixture()
  expect_error(mfrm_pca(x), 'missing values')
  expect_error(mfrm_cluster_kmeans(x, 2), 'missing values')
  complete <- mfrm_features(x$data[1:7, ], 'ID', x$features)
  factor_input <- complete$data; factor_input$Years <- factor(factor_input$Years)
  expect_error(mfrm_pca(mfrm_features(factor_input, 'ID', x$features)), 'numeric external')
  expect_error(mfrm_cluster_kmeans(mfrm_features(factor_input, 'ID', x$features), 2), 'numeric external')
  constant <- data.frame(ID=1:4, X=rep(1,4))
  expect_error(mfrm_pca(mfrm_features(constant,'ID','X')), 'must vary')
  linear <- mfrm_features(data.frame(ID=1:6, X=1:6, Y=2*(1:6)), 'ID',c('X','Y'))
  expect_equal(mfrm_pca(linear)$settings$rank, 1)
  expect_error(mfrm_pca(linear,2), 'components')
  for (bad in list(0,1.5,Inf,NA,TRUE,c(1,2))) expect_error(mfrm_pca(complete,bad),'components')
  for (bad in list(1,7,Inf,NA)) expect_error(mfrm_cluster_kmeans(complete,bad), '`k`')
  expect_error(mfrm_cluster_kmeans(complete,2,seed=-1), 'seed')
  expect_error(mfrm_cluster_kmeans(complete,2,nstart=0), 'nstart')
  expect_error(mfrm_cluster_kmeans(complete,2,iter.max=1.5), 'iter.max')
  expect_error(mfrm_pca(complete,scale=NA), 'scale')
  expect_error(mfrm_pca(complete,weights=c(Years=1,Hours=0,Ratings=1)), 'positive finite')
  expect_error(mfrm_pca(complete,weights=c(Years=1e308,Hours=1e-308,Ratings=1)), 'underflow')
  huge <- mfrm_features(data.frame(ID=1:4,X=c(-1e308,-1,1,1e308)), 'ID','X')
  expect_error(mfrm_pca(huge), 'not representable')
  expect_error(mfrm_cluster_kmeans(mfrm_pca(complete),2,scale=TRUE), 'when fitting')
  expect_error(mfrm_cluster_kmeans(complete,2,silhouette=NA), 'silhouette')
  local_mocked_bindings(kmeans = function(...) list(ifault=2), .package="mfrmr")
  expect_error(mfrm_cluster_kmeans(complete,2,silhouette=FALSE), 'converged')
})

test_that("explicitly omitted silhouettes stay unavailable and avoid quadratic limits", {
  x <- mfrm_features(data.frame(ID=seq_len(5001), X=seq_len(5001)), 'ID','X')
  expect_error(mfrm_cluster_kmeans(x,2), '5,000')
  a <- mfrm_cluster_kmeans(x,2,silhouette=FALSE,nstart=2)
  b <- mfrm_cluster_kmeans(x,3,silhouette=FALSE,nstart=2)
  expect_equal(a$settings$included,5001)
  expect_true(all(is.na(a$membership$Silhouette)))
  expect_true(all(is.na(a$cluster_summary$MeanSilhouette)))
  expect_true(all(is.na(mfrm_cluster_compare(list(Two=a,Three=b))$analysis_summary$MeanSilhouette)))
  expect_false(any(is.nan(mfrm_cluster_compare(list(Two=a,Three=b))$analysis_summary$MeanSilhouette)))
  expect_error(plot(a,draw=FALSE), 'unavailable')
})

test_that("numeric imputation comparisons pair completions and keep their own PCA", {
  skip_if_not_installed('mice')
  original <- data.frame(ID=letters[1:8], X=c(0,1,NA,3,10,11,12,NA), Y=c(3,1,2,5,12,10,14,8))
  first <- second <- original; first$X[3] <- 2; second$X[3] <- 9
  where <- is.na(original); where[8,'X'] <- FALSE
  model <- mice::as.mids(rbind(cbind(.imp=0L,original), cbind(.imp=1L,first), cbind(.imp=2L,second)), where=where)
  x <- mfrm_features(original,'ID',c('X','Y')); cells <- subset(x$missing, ID=='c')
  direct <- mfrm_cluster_imputed(x,model,cells,2,missing='omit',method='kmeans',seed=11,silhouette=FALSE)
  reduced <- mfrm_cluster_imputed(x,model,cells,2,missing='omit',method='kmeans',components=1,seed=11,silhouette=FALSE)
  expect_identical(direct$imputation_model,model)
  expect_true(all(is.na(direct$co_membership['h',])))
  expect_true(all(is.na(direct$analysis_summary$MeanSilhouette)))
  expect_identical(reduced$settings$components,1)
  expected <- matrix(0,8,8,dimnames=list(letters[1:8],letters[1:8]))
  for (i in 1:2) {
    p <- reduced$analyses[[i]]
    expect_s3_class(p$pca,'mfrm_pca')
    completion <- mice::complete(model,i)
    expect_equal(p$feature_data$data$X,completion$X)
    expected <- expected + outer(p$membership$Cluster,p$membership$Cluster,'==')/2
  }
  expect_equal(reduced$co_membership,expected)
  cmp <- mfrm_cluster_compare(list(Direct=direct,OnePC=reduced))
  expect_equal(nrow(cmp$comparisons),2)
  expect_identical(cmp$analysis_summary$Components,c(NA_integer_,NA_integer_,1L,1L))
  expect_error(mfrm_cluster_imputed(x,model,cells,2,method='pam',seed=3),'require method')
  expect_error(mfrm_cluster_imputed(x,model,cells,2,method='kmeans',linkage='average'),'no linkage')
  expect_error(mfrm_cluster_imputed(x,model,cells,2,missing='omit',method='kmeans',components=3),'Imputation 1')
})

test_that("PCA views retain stored numbers, selected geometry and matching group identities", {
  x <- numeric_feature_fixture(); p <- mfrm_pca(x,missing='omit')
  g <- mfrm_cluster_kmeans(p,2,silhouette=FALSE)
  device <- grDevices::dev.cur()
  scores <- plot_data(plot(p,type='scores',components=c(2,1),groups=g,draw=FALSE))
  expect_identical(grDevices::dev.cur(),device)
  expect_identical(scores$excluded_ids,'excluded')
  expect_equal(scores$table$PC2,p$scores$PC2[1:7])
  expect_identical(scores$table$Cluster,g$membership$Cluster[1:7])
  expect_match(scores$xlab,'PC2')
  expect_equal(plot_data(plot(p,type='scree',draw=FALSE))$table,p$variance)
  expect_equal(plot_data(plot(p,type='loadings',components=2,draw=FALSE))$table$Loading,unname(p$loadings[,2]))
  expect_error(plot(p,type='scores',components=c(1,1)), 'distinct')
  expect_error(plot(p,type='scores',components=c(1,NA)), 'distinct')
  expect_error(plot(p,type='loadings',components=4), 'retained')
  expect_error(plot(p,type='scree',components=1), 'every')
  expect_error(plot(p,groups=g), 'only used')
  g$feature_data$data$Years[1] <- 100
  expect_error(plot(p,type='scores',groups=g), 'same original')
  local_mocked_bindings(mfrm_pca=function(...) stop('Unexpected refit'))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  old <- graphics::par(c('mar','bg','fg','mfrow'))
  expect_no_warning(plot(p))
  expect_no_warning(plot(p,type='scores',labels=FALSE))
  expect_no_warning(plot(p,type='loadings'))
  expect_equal(graphics::par(names(old)),old)
})

test_that("wide PCA retains centered rank and initialization failures preserve RNG", {
  set.seed(194)
  original_seed <- .Random.seed
  on.exit(assign('.Random.seed',original_seed,envir=.GlobalEnv),add=TRUE)
  values <- matrix(rnorm(12*30),12,30)
  data <- data.frame(ID=seq_len(12),values)
  x <- mfrm_features(data,'ID',names(data)[-1])
  p <- mfrm_pca(x)
  expect_equal(p$settings$rank,11)
  expect_equal(ncol(p$loadings),11)
  expect_equal(sum(p$variance$Proportion),1)
  expect_equal(as.vector(dist(p$scores[-1])),as.vector(dist(scale(values))),tolerance=1e-12)
  rm('.Random.seed',envir=.GlobalEnv)
  a <- mfrm_cluster_kmeans(x,2,silhouette=FALSE)
  expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
  b <- mfrm_cluster_kmeans(x,2,silhouette=FALSE)
  expect_identical(a$membership,b$membership)
  expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
  assign('.Random.seed',original_seed,envir=.GlobalEnv)
  local_mocked_bindings(kmeans=function(...) { runif(1); warning('Iteration limit') }, .package="mfrmr")
  expect_error(mfrm_cluster_kmeans(x,2,silhouette=FALSE),'did not complete cleanly')
  expect_identical(.Random.seed,original_seed)
})
