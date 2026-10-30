test_that("native comma and semicolon Scorefiles preserve the same source fields", {
  comma <- test_path("fixtures", "facets-native-comma-score.txt")
  semicolon <- test_path("fixtures", "facets-native-semicolon-score.txt")
  reference <- read_facets_results(comma)
  # Values are from the actual native records, including frequency counts.
  expect_identical(reference$table$ElementID, as.character(1:6))
  expect_identical(reference$table$Label, c("????", "????", "????", "003", "??? ?", "????"))
  expect_identical(reference$table$EstimateRaw,
    c("-1.799409", "-3.095785", "5.487158", "-0.062337", "-1.000000", "1.302712"))
  expect_equal(reference$table$N, c(4, 9, 9, 6.25, 4, 9))
  expect_identical(reference$table$StatusRaw, rep("1", 6))
  expect_identical(reference$table$GroupRaw, rep("0", 6))
  expect_true(all(reference$table$EstimateAvailable))
  expect_identical(reference$table$KeyType, rep("element_id", 6))
  for (format in c("auto", "scorefile")) {
    actual <- read_facets_results(semicolon, format = format)
    fields <- setdiff(names(reference$table), "Source")
    expect_identical(actual$table[fields], reference$table[fields])
    expect_identical(summary(actual), summary(reference))
    expect_identical(actual$sources$MD5, unname(tools::md5sum(semicolon)))
    for (type in c("measure", "fit", "infit_outfit")) {
      p <- plot(actual, type = type, draw = FALSE)$data$plot_data
      q <- plot(reference, type = type, draw = FALSE)$data$plot_data
      expect_identical(p[setdiff(names(p), "Source")], q[setdiff(names(q), "Source")])
      expect_false(anyDuplicated(p$DisplayLabel) > 0)
    }
    fit_table <- read_facets_fit_table(semicolon, format = format)
    expect_identical(fit_table$EstimateRaw, reference$table$EstimateRaw)
    expect_identical(fit_table$NRaw, reference$table$NRaw)
    expect_identical(fit_table$Level, reference$table$Label)
  }
})

test_that("semicolon Scorefiles work by header, score filename and directory", {
  folder <- tempfile("facets-delimiters-")
  dir.create(folder)
  withr::defer(unlink(folder, recursive = TRUE))
  source <- test_path("fixtures", "facets-native-semicolon-score.txt")
  path <- file.path(folder, "score.1.txt")
  file.copy(source, path)
  one <- read_facets_results(path, facet_map = c("1" = "Person"))
  directory <- read_facets_results(folder, facet_map = c("1" = "Person"))
  expect_identical(directory, one)
  generic <- file.path(folder, "export.txt")
  file.copy(source, generic)
  by_header <- read_facets_results(generic, facet_map = c("1" = "Person"))
  expect_identical(by_header$table[setdiff(names(one$table), "Source")],
    one$table[setdiff(names(one$table), "Source")])
  old_reader <- read_facets_fit_table(folder, facet_map = c("1" = "Person"))
  expect_identical(old_reader$EstimateRaw, one$table$EstimateRaw)
  expect_identical(old_reader$RawFacetNumber, rep("1", 6))
})

test_that("quoted header punctuation cannot select the wrong separator", {
  # These are encoding/parser controls, not native Japanese FACETS executions.
  label <- "\u8a55\u5b9a\u8005, A; B"
  header_label <- "\u6c0f\u540d,,,,,,,,,,\t\"extra\""
  path <- tempfile(fileext = ".txt")
  withr::defer(unlink(path))
  for (separator in c(",", "\t", ";")) for (encoding in c("UTF-8", "CP932")) {
    quoted <- function(x) paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
    header <- c("Measure", "S.E.", "Status", "1", header_label, "F-Number", "F-Label")
    row <- c(".2500", ".1000", "1", "007", label, "1", "Facet name")
    lines <- c(paste(quoted(header), collapse = separator), paste(quoted(row), collapse = separator))
    writeBin(charToRaw(iconv(paste0(paste(lines, collapse = "\n"), "\n"),
      from = "UTF-8", to = encoding)), path)
    x <- read_facets_results(path, encoding = encoding)
    expect_identical(x$table$ElementID, "007")
    expect_identical(enc2utf8(x$table$Label), enc2utf8(label))
    expect_identical(x$table$EstimateRaw, ".2500")
    expect_identical(x$table$StatusRaw, "1")
    # Generic named delimited tables use the same quote-aware detection.
    # Explicit bytes avoid write.table() converting Unicode through an ASCII
    # native locale before the requested file encoding is applied.
    generic <- c(paste(quoted(c("Facet", "Level", "Measure", header_label)), collapse = separator),
      paste(quoted(c("Rater", label, ".25", "preserved")), collapse = separator))
    writeBin(charToRaw(iconv(paste0(paste(generic, collapse = "\n"), "\n"),
      from = "UTF-8", to = encoding)), path)
    y <- read_facets_results(path, encoding = encoding)
    expect_identical(enc2utf8(y$table$Label), enc2utf8(label))
    expect_equal(y$table$Estimate, .25)
  }
})
