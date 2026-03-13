
library(testthat)

# normalize_preds ---------------------------------------------------------

test_that("normalize_preds data frame rows sum to 1", {
  df <- data.frame(
    A = c(0.3, 0.4, 0.5),
    B = c(0.5, 0.3, 0.2),
    C = c(0.4, 0.2, 0.1)
  )

  result <- normalize_preds(df, vars = list(c("A", "B", "C")))
  row_sums <- rowSums(result[, c("A", "B", "C")])
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-10)
})


test_that("normalize_preds with non-unit constraint", {
  df <- data.frame(
    A = c(0.3, 0.4),
    B = c(0.5, 0.6)
  )

  result <- normalize_preds(df, vars = list(c("A", "B")), constraints = 2)
  row_sums <- rowSums(result[, c("A", "B")])
  expect_equal(row_sums, rep(2, 2), tolerance = 1e-10)
})


test_that("normalize_preds preserves proportions within group", {
  df <- data.frame(
    A = c(0.3, 0.6),
    B = c(0.9, 0.6)
  )

  result <- normalize_preds(df, vars = list(c("A", "B")))

  # Ratios should be preserved
  expect_equal(result$A / result$B, df$A / df$B, tolerance = 1e-10)
})


test_that("normalize_preds multiple variable groups with different constraints", {
  df <- data.frame(
    A = c(0.3, 0.4),
    B = c(0.5, 0.6),
    C = c(0.2, 0.3),
    D = c(0.8, 0.7)
  )

  result <- normalize_preds(
    df,
    vars = list(c("A", "B"), c("C", "D")),
    constraints = c(1, 2)
  )

  # Group 1 sums to 1
  expect_equal(rowSums(result[, c("A", "B")]), rep(1, 2), tolerance = 1e-10)
  # Group 2 sums to 2
  expect_equal(rowSums(result[, c("C", "D")]), rep(2, 2), tolerance = 1e-10)
})


test_that("normalize_preds 3D array slices sum to constraint", {
  # Create a 3D array: [draws, rows, columns]
  arr <- array(
    runif(2 * 4 * 3, 0.1, 0.9),
    dim = c(2, 4, 3),
    dimnames = list(NULL, NULL, c("A", "B", "C"))
  )

  result <- normalize_preds(arr, vars = list(c("A", "B", "C")))

  # For each draw and row, A+B+C should sum to 1
  for (d in 1:2) {
    for (r in 1:4) {
      expect_equal(
        sum(result[d, r, c("A", "B", "C")]),
        1,
        tolerance = 1e-10
      )
    }
  }
})


test_that("normalize_preds errors for mismatched vars and constraints lengths", {
  df <- data.frame(A = 0.3, B = 0.7)
  expect_error(
    normalize_preds(df, vars = list(c("A", "B")), constraints = c(1, 2)),
    "vars and constraints must be same length"
  )
})


test_that("normalize_preds errors for non-numeric constraints", {
  df <- data.frame(A = 0.3, B = 0.7)
  expect_error(
    normalize_preds(df, vars = list(c("A", "B")), constraints = "one"),
    "constraints must be numeric"
  )
})


test_that("normalize_preds errors for 2D array (matrix)", {
  mat <- matrix(c(0.3, 0.7, 0.4, 0.6), nrow = 2)
  colnames(mat) <- c("A", "B")
  # A 2D matrix has class "array" in R 4.x, so it enters the array branch

  # and either fails on dimnames[[3]] or the dimension check — either way, error
  expect_error(
    normalize_preds(mat, vars = list(c("A", "B")))
  )
})


test_that("normalize_preds errors for non-array non-data.frame", {
  expect_error(
    normalize_preds(list(A = 0.3, B = 0.7), vars = list(c("A", "B"))),
    "must be an array or data frame"
  )
})


test_that("normalize_preds errors when vars not found in data frame", {
  df <- data.frame(A = 0.3, B = 0.7)
  expect_error(
    normalize_preds(df, vars = list(c("A", "X"))),
    "not in preds"
  )
})
