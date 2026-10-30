weighted_long <- function(x) {
  d <- expand.grid(Unit = seq_len(nrow(x)), Judge = LETTERS[seq_len(ncol(x))])
  d$Score <- as.vector(x)
  d
}

# Deliberately enumerate rating pairs rather than category-count products.
weighted_oracle <- function(x, categories, type) {
  points <- if (is.numeric(categories)) categories else seq_along(categories)
  w <- function(a, b) {
    if (type == "unweighted") return(as.numeric(a == b))
    d <- abs(points[match(a, categories)] - points[match(b, categories)]) / diff(range(points))
    1 - if (type == "quadratic") d^2 else d
  }
  rates <- apply(x, 1, function(z) z[!is.na(z)], simplify = FALSE)
  rates <- Filter(length, rates)
  pa <- mean(vapply(Filter(function(z) length(z) >= 2, rates), function(z) {
    pairs <- utils::combn(seq_along(z), 2)
    mean(vapply(seq_len(ncol(pairs)), function(i) w(z[pairs[1, i]], z[pairs[2, i]]), numeric(1)))
  }, numeric(1)))
  pi <- vapply(categories, function(k) mean(vapply(rates, function(z) mean(z == k), numeric(1))), numeric(1))
  weight_sum <- sum(vapply(categories, function(a) sum(vapply(categories, function(b) w(a, b), numeric(1))), numeric(1)))
  q <- length(categories)
  pe <- weight_sum * sum(pi * (1 - pi)) / (q * (q - 1))
  pair <- x[stats::complete.cases(x[, 1:2, drop = FALSE]), 1:2, drop = FALSE]
  kp <- mean(mapply(w, pair[, 1], pair[, 2]))
  ke <- mean(vapply(pair[, 1], function(a) mean(vapply(pair[, 2], function(b) w(a, b), numeric(1))), numeric(1)))
  list(ac = (pa - pe) / (1 - pe), pa = pa, pe = pe, pi = pi,
    kappa = (kp - ke) / (1 - ke), kp = kp, ke = ke)
}

test_that("all three weightings agree with independent observation-pair enumeration", {
  x <- rbind(c(1, 2, 4, 2), c(2, 2, NA, 1), c(4, 4, 4, NA),
    c(NA, 1, NA, NA), c(1, 1, NA, NA), c(4, 2, 1, 4), rep(NA_real_, 4))
  d <- weighted_long(x)
  for (w in c("quadratic", "linear", "unweighted")) {
    ref <- weighted_oracle(x, c(1, 2, 4, 7), w)
    ac <- gwet_ac2_table(d, "Judge", "Unit", "Score", c(1, 2, 4, 7), w)
    kp <- pairwise_kappa_table(d, "Judge", "Unit", "Score", c(1, 2, 4, 7), w)
    expect_equal(ac$table$AC, ref$ac, tolerance = 1e-14)
    expect_equal(ac$table$ObservedAgreement, ref$pa)
    expect_equal(ac$table$ExpectedAgreement, ref$pe)
    expect_equal(ac$counts$UnitMeanProportion, unname(ref$pi))
    expect_equal(kp$table$Kappa[1], ref$kappa)
    expect_equal(kp$table$ObservedAgreement[1], ref$kp)
    expect_equal(kp$table$ExpectedAgreement[1], ref$ke)
    expect_equal(nrow(kp$table), 6L)
    expect_equal(ac$table$MarginalUnits, 6L)
    expect_equal(ac$table$PairableUnits, 5L)
    expect_equal(sum(ac$cases$IncludedObserved), 16L)
    expect_equal(sum(ac$cases$IncludedMargins), 17L)
  }
})

test_that("Gwet's chapter 4 complete binary example reproduces kappa", {
  x <- cbind(c(1,1,1,1,1,0,0,1,0,1), c(1,1,0,1,0,1,0,1,0,1))
  for (w in c("quadratic", "linear", "unweighted")) {
    r <- pairwise_kappa_table(weighted_long(x), "Judge", "Unit", "Score", 0:1, w)
    expect_equal(r$table$ObservedAgreement, .7)
    expect_equal(r$table$ExpectedAgreement, .54)
    expect_equal(r$table$Kappa, 8 / 23)
  }
})

test_that("QWK agrees with the established numeric-distance comparison", {
  x <- cbind(c(1,2,4,4,2,1,4), c(2,2,4,1,1,1,2))
  r <- pairwise_kappa_table(weighted_long(x), "Judge", "Unit", "Score", c(1,2,4,9))
  expect_equal(r$table$Kappa, mfrm_score_qwk(x[,1], x[,2], c(1,2,4,9))$value)
  expect_equal(r$table$ScaledObservedDisagreement / r$table$ScaledExpectedDisagreement, 1-r$table$Kappa)
  expect_equal(sum(r$counts$N), nrow(x))
  expect_equal(unname(diag(r$weights)), rep(1, 4))
})

test_that("complete-case kappa and singleton-inclusive AC retain different scopes", {
  # Gwet chapter 4, Table 4.5; that book's missing-margin extension is different.
  x <- cbind(c("A","B","C","C","B","B","A","A","B","B",NA),
             c(NA,"C","C","C","B",NA,"A","B","B","B","C"))
  d <- weighted_long(x)
  k <- pairwise_kappa_table(d, "Judge", "Unit", "Score", LETTERS[1:3])
  a <- gwet_ac2_table(d, "Judge", "Unit", "Score", LETTERS[1:3])
  ref <- weighted_oracle(x, LETTERS[1:3], "quadratic")
  expect_equal(k$table$N, 8L)
  expect_equal(k$table$Kappa, ref$kappa)
  expect_equal(k$table$ObservedAgreement, .9375)
  expect_gt(abs(k$table$Kappa - .7772), .001)
  expect_equal(a$table$AC, ref$ac)
  expect_equal(a$table$SingletonUnits, 3L)
  expect_equal(a$table$MarginalUnits, 11L)
  changed <- d; changed$Score[c(1,6)] <- "C"
  expect_identical(pairwise_kappa_table(changed, "Judge", "Unit", "Score", LETTERS[1:3])$table, k$table)
  expect_gt(abs(gwet_ac2_table(changed, "Judge", "Unit", "Score", LETTERS[1:3])$table$AC - a$table$AC), .001)
  expect_equal(a$counts$UnitMeanProportion, unname(ref$pi))
})

test_that("explicit category spacing, labels and unused categories are respected", {
  x <- cbind(c(1,2,4,4,2), c(2,2,4,1,1))
  d <- weighted_long(x)
  numeric_k <- pairwise_kappa_table(d, "Judge", "Unit", "Score", c(1,2,4))
  d$Score <- factor(d$Score, levels = c(4,1,2))
  chars <- pairwise_kappa_table(d, "Judge", "Unit", "Score", c("1","2","4"))
  expect_gt(abs(numeric_k$table$Kappa - chars$table$Kappa), .001)
  expect_identical(chars$settings$spacing, "declared_category_positions")
  expect_equal(chars$categories$Category, c("1","2","4"))
  expect_equal(unname(chars$weights), 1-(outer(1:3,1:3,"-")/2)^2)
  expect_equal(pairwise_kappa_table(d, "Judge", "Unit", "Score", c("4","2","1"))$table$Kappa, chars$table$Kappa)
  d <- weighted_long(x)
  ac <- gwet_ac2_table(d, "Judge", "Unit", "Score", c(1,2,4))
  more <- gwet_ac2_table(d, "Judge", "Unit", "Score", c(1,2,4,9))
  expect_gt(abs(ac$table$AC - more$table$AC), .001)
  expect_equal(more$counts$UnitMeanProportion[4], 0)
  expect_equal(pairwise_kappa_table(d, "Judge", "Unit", "Score", c(1,2,4,9))$table$Kappa, numeric_k$table$Kappa)
})

test_that("negative, constant, empty and one-unit results retain their meaning", {
  opposite <- weighted_long(cbind(c(0,1), c(1,0)))
  expect_equal(pairwise_kappa_table(opposite, "Judge", "Unit", "Score", 0:1)$table$Kappa, -1)
  expect_equal(gwet_ac2_table(opposite, "Judge", "Unit", "Score", 0:1)$table$AC, -1)
  constant <- opposite; constant$Score <- 0
  kp <- pairwise_kappa_table(constant, "Judge", "Unit", "Score", 0:1)$table
  expect_true(is.na(kp$Kappa)); expect_equal(kp$Status, "zero_expected_disagreement")
  expect_equal(gwet_ac2_table(constant, "Judge", "Unit", "Score", 0:1)$table$AC, 1)
  empty <- constant; empty$Score <- NA_real_
  ac <- gwet_ac2_table(empty, "Judge", "Unit", "Score", 0:1)
  expect_equal(ac$table$Status, "no_pairable_ratings")
  expect_true(is.na(ac$table$AC)); expect_true(all(is.na(ac$counts$UnitMeanProportion)))
  expect_equal(nrow(pairwise_kappa_table(empty, "Judge", "Unit", "Score", 0:1)$counts), 0L)
  single <- weighted_long(matrix(c(0,1), 1))
  expect_equal(gwet_ac2_table(single, "Judge", "Unit", "Score", 0:1)$table$AC, -1)
  expect_equal(pairwise_kappa_table(single, "Judge", "Unit", "Score", 0:1)$table$Status, "fewer_than_two_complete_units")
  ac1 <- gwet_ac2_table(opposite, "Judge", "Unit", "Score", 0:1, "unweighted")
  expect_equal(ac1$table$Coefficient, "AC1")
})

test_that("matching audits preserve missing versus absent scores and composite keys", {
  d <- weighted_long(cbind(c(1,2,NA,2), c(1,NA,2,2), c(2,1,2,1)))
  d$Task <- rep(c("t1","t2"),6); d$Person <- rep(c(1,1,2,2),3)
  d <- d[-2,]
  for (fun in list(pairwise_kappa_table, gwet_ac2_table)) {
    a <- fun(d, "Judge", c("Person","Task"), "Score", 1:2)
    b <- fun(d[nrow(d):1,], "Judge", c("Person","Task"), "Score", 1:2)
    expect_equal(a$units, b$units)
    expect_error(fun(rbind(d,d[1,]), "Judge", c("Person","Task"), "Score", 1:2), "[Rr]epeat|[Dd]uplicate")
  }
  k <- pairwise_kappa_table(d, "Judge", c("Person","Task"), "Score", 1:2)
  expect_false(k$cases$Assigned1[2]); expect_false(k$cases$Missing1[2])
  expect_true(k$cases$Missing1[3]); expect_true(k$cases$Assigned1[3])
  a <- gwet_ac2_table(d, "Judge", c("Person","Task"), "Score", 1:2)
  expect_equal(a$unit_summary$AbsentAssignments, c(0L,1L,0L,0L))
  expect_equal(a$unit_summary$MissingScores, c(0L,1L,1L,0L))
  expect_equal(nrow(a$cases), nrow(d))
})

test_that("numeric matching and distances survive extreme scales and unused extremes", {
  x <- cbind(c(0,1,2,2,1),c(1,1,2,0,1))
  ref <- pairwise_kappa_table(weighted_long(x), "Judge", "Unit", "Score", 0:2)$table$Kappa
  for (factor in c(1e-200,1e200)) {
    d <- weighted_long(x*factor)
    expect_equal(pairwise_kappa_table(d, "Judge", "Unit", "Score", (0:2)*factor)$table$Kappa, ref)
    expect_equal(gwet_ac2_table(d, "Judge", "Unit", "Score", (0:2)*factor)$table$AC,
      gwet_ac2_table(weighted_long(x), "Judge", "Unit", "Score", 0:2)$table$AC)
  }
  expect_equal(pairwise_kappa_table(weighted_long(x*1e-200), "Judge", "Unit", "Score", c((0:2)*1e-200,1e200))$table$Kappa, ref)
  y <- (x-1)*1e308
  expect_equal(pairwise_kappa_table(weighted_long(y), "Judge", "Unit", "Score", c(-1e308,0,1e308))$table$Kappa, ref)
  z <- matrix(c(-2000000000L,2000000000L,2000000000L,-2000000000L),2)
  expect_equal(pairwise_kappa_table(weighted_long(z), "Judge", "Unit", "Score", c(-2000000000L,2000000000L))$table$Kappa,-1)
  close <- weighted_long(cbind(c(1,1+1e-14),c(1+1e-14,1)))
  expect_equal(pairwise_kappa_table(close, "Judge", "Unit", "Score", c(1,1+1e-14))$table$Kappa,-1)
})

test_that("invalid category universes and unsupported score types are rejected", {
  d <- weighted_long(cbind(1:3,3:1))
  for (fun in list(pairwise_kappa_table, gwet_ac2_table)) {
    for (categories in list(NULL, 1, c(1,1,2), c(1,2,NA), c(1,2,Inf), c("1","2","3"))) {
      expect_error(fun(d,"Judge","Unit","Score",categories), "categories")
    }
    expect_error(fun(d,"Judge","Unit","Score",1:2), "outside")
    invalid <- d; invalid$Score[1] <- Inf
    expect_error(fun(invalid,"Judge","Unit","Score",1:3), "finite")
    expect_error(fun(d,"Judge","Unit","Score",1:3,"custom"), "arg")
  }
})
