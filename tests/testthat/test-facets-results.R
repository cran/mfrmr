facets_result_example <- function(...) {
  read_facets_results(system.file("extdata", "facets-example", package = "mfrmr"),
    facet_map = c("Person", "Rater", "Task"), ...)
}

facets_graph_example <- function(...) {
  read_facets_results(graph_file = system.file("extdata", "facets-graph-example",
    "graph.txt", package = "mfrmr"), ...)
}

facets_graph_test_file <- function(lines) {
  path <- tempfile(fileext = ".txt")
  withr::defer(unlink(path), envir = parent.frame())
  writeLines(lines, path, useBytes = TRUE)
  path
}

test_that("native Graphfile coordinates and category coding survive import", {
  x <- facets_graph_example(source_version = "4.5.0")
  expect_s3_class(x, "mfrm_facets_results")
  expect_false(inherits(x, "mfrm_fit"))
  expect_equal(nrow(x$graph$table), 4004L)
  expect_equal(summary(x, component = "curves")$GridRows, rep(1001L, 4))
  expect_equal(summary(x, component = "curves")$ProbabilityCheckIssues, rep(0L, 4))
  d <- as.data.frame(x, component = "curves")
  expect_identical(d, x$graph$curves)
  p <- plot(x, type = "ccc", scales = 1, draw = FALSE)$data
  expect_equal(p$plot_data$Y[p$plot_data$SourceRow == 2], c(.9519, .0474, .0007, 0))
  expect_equal(p$plot_data$ValueRaw[p$plot_data$SourceRow == 2], c(".9519", ".0474", ".0007", ".0000"))
  expect_equal(p$plot_data$ValueResolution, rep(.0001, 4004))
  expect_equal(p$plot_data$X, rep(seq(-4, 4, length.out = 1001), 4))
  expect_false(p$display$show_title)
  expect_false(p$display$show_notes)
  expect_equal(x$sources$MD5, unname(tools::md5sum(x$sources$Path)))
  expect_error(plot(x, draw = FALSE), "No element results")
  expect_equal(nrow(summary(x)), 0L)
  expect_output(print(x), "Graphfile curves")
  binary <- read_facets_results(graph_file = test_path("fixtures", "facets-graph-binary_csv.txt"))
  expect_identical(binary$graph$categories, c("1", "2"))
  expect_equal(plot(binary, type = "expected", draw = FALSE)$data$plot_data$Y[1], 1.0474)
  expect_equal(binary$graph$table$ExpCat[1], "1")
})

test_that("fixed and tab Graphfiles retain source values and precision", {
  x <- read_facets_results(graph_file = test_path("fixtures", "facets-graph-pcm_fixed.txt"))
  csv <- facets_graph_example()$graph$curves
  fixed <- x$graph$curves
  key <- function(d) paste(d$Scale, d$Measure, d$Quantity, d$Category)
  rows <- match(key(fixed), key(csv))
  expect_false(anyNA(rows))
  expect_equal(fixed$Value, csv$Value[rows])
  expect_equal(fixed$ValueResolution, csv$ValueResolution[rows])
  expect_identical(x$graph$layout, "fixed_width")
  path <- facets_graph_test_file(c('"Scale"\t"Measure"\t"Expected"\t"ExpCat"\t"Prob:10"\t"Prob:20"',
    '7\t-1.000\t12.500\t10\t.7500\t.2500', '7\t1.000\t18.000\t20\t.2000\t.8000'))
  z <- read_facets_results(graph_file = path)
  expect_equal(z$graph$categories, c("10", "20"))
  expect_equal(z$graph$layout, "tab")
  expect_equal(plot(z, type = "expected", draw = FALSE)$data$plot_data$Y, c(12.5, 18))
  expect_equal(z$graph$curves$MeasureResolution, rep(.001, 6))
})

test_that("probability checks use printed precision without renormalizing", {
  path <- facets_graph_test_file(c('Scale,Measure,Expected,ExpCat,Prob:0,Prob:1,Prob:2',
    '1,-1.00,1.00,1,.33,.33,.33', '1,0.00,1.00,1,.30,.30,.30',
    '1,1.00,1.00,1,.3300,.3300,.3300', '1,2.00,1.00,1,3.3e-1,3.3e-1,3.3e-1'))
  x <- read_facets_results(graph_file = path)
  expect_equal(x$graph$diagnostics$ProbabilitySum, c(.99, .9, .99, .99))
  expect_equal(x$graph$diagnostics$RoundingAllowance, c(.015, .015, .00015, .015))
  p <- plot(x, type = "ccc", draw = FALSE)$data$plot_data
  expect_equal(p$Available, rep(c(TRUE, FALSE, FALSE, TRUE), 3))
  expect_equal(p$Value, rep(c(.33, .30, .33, .33), 3))
  expect_equal(p$Y, rep(c(.33, NA, NA, .33), 3))
  expect_true(all(plot(x, type = "expected", draw = FALSE)$data$plot_data$Available))
  expect_equal(summary(x, component = "curves")$ProbabilityCheckIssues, 2L)
})

test_that("missing and invalid graph rows are retained and break curve segments", {
  path <- facets_graph_test_file(c('Scale,Measure,Expected,ExpCat,Prob:low,Prob:high',
    '1,-2,.2,0,.8,.2', '1,-1,NA,0,,.3', '1,NA,.5,1,.5,.5',
    '1,1,<.8,1,1.1,-.1', '1,2,.9,1,.1,.9'))
  x <- read_facets_results(graph_file = path)
  p <- plot(x, type = "ccc", draw = FALSE)$data$plot_data
  low <- p[p$Category == "low", ]
  expect_equal(low$Available, c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(low$Segment, c(0L, NA_integer_, NA_integer_, NA_integer_, 3L))
  expect_equal(low$Value[4], 1.1)
  expect_match(low$Reason[4], "out_of_range")
  expect_match(low$Reason[3], "measure:")
  expect_equal(x$graph$diagnostics$ProbabilityStatus[2], "incomplete_or_invalid")
  expect_true(p$Available[p$Category == "high" & p$SourceRow == 3])
  e <- plot(x, type = "expected", draw = FALSE)$data$plot_data
  expect_equal(e$Available, c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(e$ValueRaw[4], "<.8")
  expect_equal(nrow(x$graph$table), 5L)
})

test_that("Graphfile rejects ambiguous identity, schema and grid order", {
  header <- 'Scale,Measure,Expected,ExpCat,Prob:0,Prob:1'
  good <- '1,0,.5,0,.5,.5'
  read <- function(lines) read_facets_results(graph_file = facets_graph_test_file(lines))
  expect_error(read(c(good, good)), "headings")
  expect_error(read(c(sub('Prob:1', 'Prob:0', header), good)), "unique")
  expect_error(read(c(header, good, good)), "Duplicate")
  expect_error(read(c(header, good, '01,0.0,.5,0,.5,.5')), "Duplicate")
  expect_error(read(c(header, '0,0,.5,0,.5,.5')), "positive integers")
  expect_error(read(c(header, good, '1,2,.9,1,.1,.9', '1,1,.7,1,.3,.7')), "monotone")
  expect_error(read(c(header, good, '1,2,.9,1,.1')), "same fields")
  expect_error(read_facets_results(), "Supply")
  expect_error(read_facets_results(graph_file = c('one', 'two')), "one existing")
  decreasing <- read(c(header, '1,2,.9,1,.1,.9', good))
  expect_equal(plot(decreasing, type = "expected", draw = FALSE)$data$plot_data$X, c(2, 0))
  single <- read(c(header, good))
  expect_true(all(single$graph$curves$Available))
})

test_that("curves coexist with elements without implicit scale alignment", {
  native <- facets_graph_example()
  x <- read_facets_results(data.frame(Facet = "Task", Level = "one", Measure = 100),
    graph_file = native$sources$Path, measure_units = "transformed units", common_scale = TRUE)
  expect_equal(x$sources$SourceID, c(1L, 2L))
  expect_true(all(x$graph$curves$SourceID == 2L))
  expect_equal(plot(x, draw = FALSE)$data$plot_data$X, 100)
  expect_equal(plot(x, type = "expected", scales = 2, draw = FALSE)$data$xlab, "Graphfile measure")
  expect_true(is.na(plot(x, type = "expected", scales = 2, draw = FALSE)$data$metadata$measure_units))
  expect_equal(unique(plot_data(x, "plot_data", type = "ccc", scales = c(4, 2))$Scale), c("2", "4"))
  expect_error(plot(x, type = "ccc", facets = "Task", draw = FALSE), "scales")
  expect_error(plot(x, scales = 1, draw = FALSE), "only")
  expect_error(plot(x, type = "expected", scales = 9, draw = FALSE), "distinct")
  expect_error(plot(facets_result_example(), type = "ccc", draw = FALSE), "No Graphfile")
  expect_error(summary(facets_result_example(), component = "curves"), "No Graphfile")
  path <- tempfile(fileext = ".rds")
  withr::defer(unlink(path))
  saveRDS(x, path)
  expect_identical(readRDS(path), x)
})

test_that("Graphfile category labels retain declared Japanese encodings", {
  labels <- c("\u4f4e", "\u9ad8")
  lines <- c(paste0('Scale,Measure,Expected,ExpCat,Prob:', labels[1], ',Prob:', labels[2]),
    '1,0,.5,0,.5,.5')
  path <- facets_graph_test_file(lines)
  x <- read_facets_results(graph_file = path)
  expect_identical(x$graph$categories, labels)
  encoded <- iconv(lines, from = "UTF-8", to = "CP932")
  skip_if(anyNA(encoded), "CP932 unavailable")
  writeLines(encoded, path, useBytes = TRUE)
  y <- read_facets_results(graph_file = path, encoding = "CP932")
  expect_identical(y$graph$categories, labels)
  expect_identical(y$graph$curves, x$graph$curves)
})

test_that("FACETS native exports preserve coordinates, states and file identity", {
  x <- facets_result_example(source_version = "4.5.0")
  expect_s3_class(x, "mfrm_facets_results")
  expect_false(inherits(x, "mfrm_fit"))
  expect_equal(nrow(x$table), 11L)
  expect_equal(summary(x)$Rows, c(6L, 2L, 3L))
  expect_equal(x$table$ElementID, c(as.character(1:6), "1", "2", "1", "2", "3"))
  expect_equal(x$table$Estimate[1:6], c(-1.799409, -3.095785, 5.487158, -.062337, -1, 1.302712))
  expect_equal(x$table$EstimateRaw[1], "-1.799409")
  expect_equal(x$table$SourceStatus[1:6], c("extreme_low", "extreme_low", "extreme_high", "estimated", "anchored", "estimated"))
  expect_equal(x$table$Sign, c(rep(1, 6), rep(-1, 5)))
  expect_equal(x$table$N[4], 6.25)
  expect_equal(x$sources$MD5, unname(tools::md5sum(x$sources$Path)))
  expect_true(all(x$table$EstimateAvailable))
  expect_false(any(x$table$InfitAvailable[1:3]))
  expect_equal(x$table$Infit[1:3], rep(1, 3))
  expect_equal(as.data.frame(x)$SourceVersion, rep("4.5.0", 11))
  expect_error(fit_mfrm(x), "data|Data|argument")
})

test_that("FACETS plots preserve all rows and require explicit shared coordinates", {
  x <- facets_result_example()
  before <- x
  p <- plot(x, type = "measure", facets = "Person", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_equal(p$data$plot_data$X, x$table$Estimate[1:6])
  expect_equal(p$data$plot_data$Symbol, c(25L, 25L, 24L, 19L, 15L, 19L))
  expect_false(p$data$display$show_title)
  expect_false(p$data$display$show_notes)
  expect_error(plot(x, type = "wright", draw = FALSE), "common_scale")
  y <- facets_result_example(common_scale = TRUE, measure_units = "logits")
  w <- plot(y, type = "wright", draw = FALSE)
  expect_equal(w$data$plot_data$Y, y$table$Estimate)
  expect_equal(w$data$plot_data$X, c(rep(1, 6), rep(2, 2), rep(3, 3)))
  for (type in c("fit", "infit_outfit")) {
    p <- plot(x, type = type, draw = FALSE)
    expect_equal(nrow(p$data$plot_data), 11L)
    expect_false(any(p$data$plot_data$Available[1:3]))
    expect_true(all(is.na(p$data$plot_data$X[1:3])))
    expect_match(p$data$plot_data$Reason[1], "source_fit_not_interpretable")
    expect_equal(p$data$reference_lines$value, rep(1, if (type == "fit") 1 else 2))
  }
  expect_identical(x, before)
  expect_equal(plot_data(x, "plot_data", facets = "Person")$X, x$table$Estimate[1:6])
  expect_equal(plot_data(x, "notes"), x$notes)
  expect_error(plot(x, facets = "Missing", draw = FALSE), "facet")
  expect_error(plot(x, draw = NA), "TRUE or FALSE")
  expect_error(plot(x, show_ci = TRUE), "must be empty")
})

test_that("native group-anchor flags do not imply individually fixed measures", {
  for (variant in c("origin", "restricted")) {
    path <- test_path("fixtures", paste0("facets-group-", variant, "-score.txt"))
    x <- read_facets_results(path, facet_map = c("2" = "Rater"), source_version = "4.5.0")
    before <- serialize(x, NULL)
    expected_flags <- if (variant == "origin") c("1", "1", "-1") else rep("1", 3)
    expect_identical(x$table$StatusRaw, expected_flags)
    expect_identical(x$table$GroupRaw, c("1", "1", if (variant == "origin") "0" else "2"))
    expect_equal(mean(x$table$Estimate[1:2]), .6, tolerance = 1e-12)
    expect_true(all(is.finite(x$table$SE) & x$table$SE > 0))
    expect_true(all(x$table$EstimateAvailable))
    sm <- summary(x)
    expect_equal(sm$AnchorFlagged, if (variant == "origin") 2L else 3L)
    expect_false("Anchored" %in% names(sm))
    expect_output(print(x), "individual measure need not be fixed", fixed = TRUE)
    for (type in c("measure", "fit", "infit_outfit")) {
      p <- plot(x, type = type, draw = FALSE)$data
      expect_equal(p$plot_data$Symbol[1:2], c(15L, 15L))
      expect_true("Direct/group anchor" %in% p$legend$label)
      expect_false("Anchored" %in% p$legend$label)
      expect_identical(p$data, as.data.frame(x$table))
      note <- p$notes$Text[p$notes$Type == "anchor_status"]
      expect_length(note, 1L)
      expect_match(note, "individual measure need not be fixed", fixed = TRUE)
    }
    # A saved import predating the explanation must receive it on replotting.
    old <- unserialize(before)
    old$notes <- old$notes[old$notes$Type != "anchor_status", , drop = FALSE]
    old_bytes <- serialize(old, NULL)
    updated <- plot(old, draw = FALSE)$data
    expect_equal(sum(updated$notes$Type == "anchor_status"), 1L)
    expect_identical(updated$plot_data, plot(x, draw = FALSE)$data$plot_data)
    expect_identical(serialize(old, NULL), old_bytes)
    expect_identical(serialize(x, NULL), before)

    # The same source fields have the same meaning in named delimited input.
    tab <- data.frame(Facet = "Rater", ElementID = x$table$ElementID,
      Label = x$table$Label, Measure = x$table$EstimateRaw, SE = x$table$SE,
      Infit = x$table$InfitRaw, Outfit = x$table$OutfitRaw,
      Status = x$table$StatusRaw, Group = x$table$GroupRaw)
    csv <- tempfile(fileext = ".csv")
    withr::defer(unlink(csv))
    utils::write.csv(tab, csv, row.names = FALSE)
    y <- read_facets_results(csv)
    for (field in c("StatusRaw", "SourceStatus", "GroupRaw", "EstimateRaw",
                    "Estimate", "SE", "EstimateAvailable", "InfitAvailable", "OutfitAvailable"))
      expect_identical(y$table[[field]], x$table[[field]])
    expect_identical(summary(y), summary(x))
  }
})

test_that("native G/X extreme group members retain unavailable fit placeholders", {
  for (boundary in c("low", "high")) for (mode in c("g", "x")) {
    path <- test_path("fixtures", paste0("facets-group-extreme-", boundary,
      "-", mode, "-score.txt"))
    x <- read_facets_results(path, facet_map = c("1" = "Person"), source_version = "4.5.0")
    idx <- match(c("3", "24"), x$table$ElementID)
    d <- x$table[idx, ]
    expect_identical(d$GroupRaw, c("1", "1"))
    expect_identical(d$StatusRaw, c(if (boundary == "low") "-2" else "-3", "1"))
    expect_identical(d$SourceStatus, c(paste0("extreme_", boundary), "anchored"))
    expect_identical(d$InfitRaw[1], "1.00")
    expect_identical(d$OutfitRaw[1], "1.00")
    expect_identical(d$InfitAvailable, c(FALSE, TRUE))
    expect_identical(d$OutfitAvailable, c(FALSE, TRUE))
    expect_identical(d$InfitReason[1], "source_fit_not_interpretable")
    expect_true(all(d$EstimateAvailable))
    expect_true(all(is.finite(d$SE) & d$SE > 0))
    if (mode == "g") {
      expect_equal(mean(d$Estimate), -.4, tolerance = 1e-12)
      expect_gt(abs(diff(d$Estimate)), 1)
    } else {
      # Retain the native locations; do not invent a tail adjustment for X.
      expect_equal(d$Estimate, c(-.4, -.4), tolerance = 1e-12)
    }
    expect_equal(summary(x)$AnchorFlagged, 1L)
    expect_equal(summary(x)$Extreme, 1L)
    raw <- serialize(x, NULL)
    for (type in c("measure", "fit", "infit_outfit")) {
      p <- plot(x, type = type, draw = FALSE)$data
      z <- p$plot_data[match(c("3", "24"), p$plot_data$ElementID), ]
      expect_identical(z$Available, c(type == "measure", TRUE))
      expect_identical(z$Symbol, c(if (boundary == "low") 25L else 24L, 15L))
      expect_identical(p$data, as.data.frame(x$table))
    }
    expect_identical(serialize(x, NULL), raw)
    expect_match(x$notes$Text[x$notes$Type == "anchor_status"],
      "does not count every member", fixed = TRUE)
  }
})

test_that("reporting groups and SEs cannot silently determine an anchor subtype", {
  d <- data.frame(Facet = "Rater", ElementID = as.character(1:5), Measure = 0,
    Status = 1, Group = c("0", "1", "1", NA, "A#report"), SE = c(0, .2, NA, 0, .4))
  x <- read_facets_results(d)
  expect_identical(x$table$SourceStatus, rep("anchored", 5))
  expect_identical(x$table$GroupRaw, d$Group)
  expect_identical(x$table$SE, d$SE)
  expect_equal(summary(x)$AnchorFlagged, 5L)
  expect_identical(plot(x, draw = FALSE)$data$plot_data$Symbol, rep(15L, 5))
  expect_match(x$notes$Text[x$notes$Type == "anchor_status"], "GroupRaw, imported group size")
  # Curves have no element anchor flags, but retain the same summary schema.
  empty <- summary(facets_graph_example())
  expect_identical(empty$AnchorFlagged, integer())
  expect_false("Anchored" %in% names(empty))
})

test_that("FACETS missing, censored and unmeasurable values remain distinct", {
  x <- read_facets_results(data.frame(Facet = "Rater", ElementID = as.character(1:7),
    Label = letters[1:7], Measure = c("0", "1", NA, "<2.0", "0", "1", "2"),
    Infit = c("1", ".9", NA, "1.1", "1", "1", "-1"),
    Outfit = "1.0", Status = c(0, -4, -1, -1, -5, NA, 99)))
  expect_equal(x$table$SourceStatus, c("no_data", "unmeasurable", "estimated",
    "estimated", "single_observation", "unknown", "unknown"))
  expect_equal(x$table$EstimateAvailable, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(x$table$EstimateRaw[4], "<2.0")
  expect_equal(x$table$InfitReason[7], "invalid_mean_square")
  p <- plot(x, type = "fit", draw = FALSE)$data$plot_data
  expect_equal(p$Available, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(nrow(p), 7L)
  expect_equal(summary(x)$UnknownStatus, 2L)
  x <- read_facets_results(data.frame(Facet = "Rater", Level = "A", Measure = NA_real_))
  expect_false(plot(x, draw = FALSE)$data$plot_data$Available)
})

test_that("FACETS duplicate labels retain IDs while ambiguous keys fail", {
  d <- data.frame(Facet = "Rater", ElementID = c("01", "02"),
    Label = c("same", "same"), Measure = c(.5, -.5))
  x <- read_facets_results(d)
  expect_equal(x$table$Label, c("same", "same"))
  expect_equal(plot(x, draw = FALSE)$data$plot_data$DisplayLabel, c("same [01]", "same [02]"))
  expect_error(read_facets_results(d[c(1, 1), ]), "Duplicate")
  d$ElementID <- NULL
  expect_error(read_facets_results(d), "distinct element IDs")
  d$Label[2] <- NA_character_
  expect_error(read_facets_results(d), "no rows were discarded")
  expect_error(read_facets_results(d[FALSE, ]), "No FACETS result rows")
  d$ElementID <- c("01", "02")
  expect_equal(read_facets_results(d)$table$Label, c("same", "02"))
})

test_that("native CSV and TSV imports preserve numeric headers and nonconsecutive facet maps", {
  path <- tempfile(fileext = ".csv")
  withr::defer(unlink(path))
  lines <- c('"Measure","S.E.","Status","Sign","7","Judge","F-Number","F-Label"',
    '0.1000,0.20,-1,-1,01,"Same label",7,"Judge"',
    '-0.2000,0.30,1,-1,02,"Same label",7,"Judge"')
  for (sep in c(",", "\t")) {
    writeLines(gsub(",", sep, lines, fixed = TRUE), path)
    x <- read_facets_results(path, facet_map = c("7" = "Rater"))
    expect_equal(x$table$Facet, c("Rater", "Rater"))
    expect_equal(x$table$ElementID, c("01", "02"))
    expect_equal(x$table$Label, rep("Same label", 2))
    expect_equal(x$table$EstimateRaw, c("0.1000", "-0.2000"))
    expect_equal(x$table$SourceStatus, c("estimated", "anchored"))
  }
  expect_error(read_facets_results(c(path, path)), "more than once")
})

test_that("Japanese labels survive UTF-8 and Windows code-page delimited imports", {
  path <- tempfile(fileext = ".csv")
  withr::defer(unlink(path))
  label <- "\u8a55\u5b9a\u8005\u7532"
  text <- paste0('Facet,ElementID,Label,Measure,Status\nRater,01,', label, ',0.1250,-1\n')
  for (encoding in c("UTF-8", "CP932")) {
    encoded <- iconv(text, from = "UTF-8", to = encoding)
    if (is.na(encoded)) next
    writeBin(charToRaw(encoded), path)
    x <- read_facets_results(path, encoding = encoding)
    expect_equal(enc2utf8(x$table$Label), enc2utf8(label))
    expect_identical(x$table$EstimateRaw, "0.1250")
    saved <- tempfile(fileext = ".rds")
    withr::defer(unlink(saved))
    saveRDS(x, saved)
    expect_identical(readRDS(saved), x)
  }
})

test_that("fixed-width no-data rows remain and conflicting facet identity fails", {
  path <- tempfile(fileext = ".txt")
  withr::defer(unlink(path))
  source <- system.file("extdata", "facets-example", "score.2.txt", package = "mfrmr")
  lines <- readLines(source)
  row <- lines[3]
  substr(row, 1, 240) <- paste(rep(" ", 240), collapse = "")
  substr(row, 141, 150) <- sprintf("%10s", "0")
  writeLines(c(lines[1:2], row), path)
  x <- read_facets_results(path, facet = "Rater")
  expect_equal(nrow(x$table), 1L)
  expect_equal(x$table$ElementID, "1")
  expect_equal(x$table$SourceStatus, "no_data")
  expect_false(x$table$EstimateAvailable)
  bad <- tempfile()
  dir.create(bad)
  withr::defer(unlink(bad, recursive = TRUE))
  file.copy(source, file.path(bad, "score.7.txt"))
  expect_error(read_facets_results(bad), "disagrees")
})

test_that("display-label spacing is retained and duplicate columns cannot choose a value silently", {
  d <- data.frame(Facet = "Rater", ElementID = "01", Label = " A ", Measure = .2)
  expect_identical(read_facets_results(d)$table$Label, " A ")
  d$Conflicting <- .8
  names(d)[5] <- "Measure"
  expect_error(read_facets_results(d), "Duplicate FACETS source column names")
})

test_that("named delimited score extracts are not mistaken for fixed-width files", {
  path <- tempfile(fileext = ".txt")
  withr::defer(unlink(path))
  lines <- c('Facet,ElementID,Label,T.Score,T.Count,Measure,InfitMS',
    'Rater,01,A,12,10,0.2500,0.9')
  for (sep in c(",", "\t", ";")) {
    writeLines(gsub(",", sep, lines, fixed = TRUE), path)
    x <- read_facets_results(path)
    expect_identical(x$table$EstimateRaw, "0.2500")
    expect_equal(x$table$N, 10)
  }
})
