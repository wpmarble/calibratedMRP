
library(testthat)

# poststratify ------------------------------------------------------------

test_that("poststratify basic weighted mean matches manual computation", {
  ps <- tibble::tibble(
    county = c("A", "A", "B", "B"),
    estimate = c(0.4, 0.6, 0.3, 0.7),
    weight = c(100, 200, 150, 50)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight, by = county)

  # Manual weighted means
  wm_A <- weighted.mean(c(0.4, 0.6), c(100, 200))
  wm_B <- weighted.mean(c(0.3, 0.7), c(150, 50))

  result_A <- result$estimate[result$county == "A"]
  result_B <- result$estimate[result$county == "B"]

  expect_equal(result_A, wm_A, tolerance = 1e-10)
  expect_equal(result_B, wm_B, tolerance = 1e-10)
})


test_that("poststratify handles multiple outcomes", {
  ps <- tibble::tibble(
    county = c("A", "A", "B", "B"),
    voteshare = c(0.4, 0.6, 0.3, 0.7),
    turnout = c(0.5, 0.8, 0.6, 0.4),
    weight = c(100, 200, 150, 50)
  )

  result <- poststratify(ps, outcomes = c(voteshare, turnout),
                         weight = weight, by = county)

  expect_true("voteshare" %in% names(result))
  expect_true("turnout" %in% names(result))
  expect_equal(nrow(result), 2)
})


test_that("poststratify uniform weights equals unweighted means", {
  ps <- tibble::tibble(
    county = c("A", "A", "A"),
    estimate = c(0.2, 0.4, 0.6),
    weight = c(1, 1, 1)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight, by = county)
  expect_equal(result$estimate, mean(c(0.2, 0.4, 0.6)), tolerance = 1e-10)
})


test_that("poststratify n_out column has correct weight sums", {
  ps <- tibble::tibble(
    county = c("A", "A", "B", "B"),
    estimate = c(0.4, 0.6, 0.3, 0.7),
    weight = c(100, 200, 150, 50)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight, by = county)

  expect_true("n" %in% names(result))
  expect_equal(result$n[result$county == "A"], 300)
  expect_equal(result$n[result$county == "B"], 200)
})


test_that("poststratify SE computation matches formula", {
  # Reset the warning flag so SE message fires
  .calibratedMRP_env$poststratify_se_warned <- FALSE

  ps <- tibble::tibble(
    county = c("A", "A"),
    estimate = c(0.4, 0.6),
    estimate_se = c(0.05, 0.03),
    weight = c(100, 200)
  )

  result <- suppressMessages(
    poststratify(ps, outcomes = estimate, ses = TRUE,
                 weight = weight, by = county)
  )

  # Manual SE: sqrt(sum(w^2 * se^2) / sum(w)^2)
  w <- c(100, 200)
  se <- c(0.05, 0.03)
  expected_se <- sqrt(sum(w^2 * se^2) / sum(w)^2)

  expect_equal(result$estimate_se, expected_se, tolerance = 1e-10)
})


test_that("poststratify SE column ordering is correct", {
  .calibratedMRP_env$poststratify_se_warned <- FALSE

  ps <- tibble::tibble(
    county = c("A", "A"),
    voteshare = c(0.4, 0.6),
    voteshare_se = c(0.05, 0.03),
    turnout = c(0.5, 0.7),
    turnout_se = c(0.04, 0.02),
    weight = c(100, 200)
  )

  result <- suppressMessages(
    poststratify(ps, outcomes = c(voteshare, turnout), ses = TRUE,
                 weight = weight, by = county)
  )

  expected_order <- c("county", "voteshare", "voteshare_se",
                       "turnout", "turnout_se", "n")
  expect_equal(names(result), expected_order)
})


test_that("poststratify SE warning fires once per session", {
  .calibratedMRP_env$poststratify_se_warned <- FALSE

  ps <- tibble::tibble(
    county = "A",
    estimate = 0.5,
    estimate_se = 0.05,
    weight = 100
  )

  # First call: should produce message
  expect_message(
    poststratify(ps, outcomes = estimate, ses = TRUE,
                 weight = weight, by = county),
    "independence across cells"
  )

  # Second call: should NOT produce message (already warned)
  expect_no_message(
    poststratify(ps, outcomes = estimate, ses = TRUE,
                 weight = weight, by = county)
  )

  # Reset for other tests
  .calibratedMRP_env$poststratify_se_warned <- FALSE
})


test_that("poststratify handles multiple grouping variables", {
  ps <- tibble::tibble(
    state = c("PA", "PA", "PA", "PA"),
    county = c("A", "A", "B", "B"),
    estimate = c(0.4, 0.6, 0.3, 0.7),
    weight = c(100, 200, 150, 50)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight,
                         by = c(state, county))

  expect_equal(nrow(result), 2)
  expect_true("state" %in% names(result))
  expect_true("county" %in% names(result))
})


test_that("poststratify na.rm=TRUE ignores NAs in outcomes", {
  ps <- tibble::tibble(
    county = c("A", "A", "A"),
    estimate = c(0.4, NA, 0.6),
    weight = c(100, 200, 100)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight,
                         by = county, na.rm = TRUE)

  expect_false(is.na(result$estimate))
  expected <- weighted.mean(c(0.4, 0.6), c(100, 100), na.rm = TRUE)
  expect_equal(result$estimate, expected, tolerance = 1e-10)
})


test_that("poststratify na.rm=FALSE propagates NAs", {
  ps <- tibble::tibble(
    county = c("A", "A"),
    estimate = c(0.4, NA),
    weight = c(100, 200)
  )

  result <- poststratify(ps, outcomes = estimate, weight = weight,
                         by = county, na.rm = FALSE)

  expect_true(is.na(result$estimate))
})


test_that("poststratify ses=TRUE works with multi-variable by", {
  .calibratedMRP_env$poststratify_se_warned <- FALSE

  ps <- tibble::tibble(
    state  = c("PA", "PA", "NY", "NY"),
    county = c("A",  "B",  "C",  "D"),
    estimate    = c(0.4, 0.6, 0.3, 0.7),
    estimate_se = c(0.05, 0.03, 0.04, 0.02),
    weight = c(100, 200, 150, 50)
  )

  result <- suppressMessages(
    poststratify(ps, outcomes = estimate, ses = TRUE,
                 weight = weight, by = c(state, county))
  )

  # Each (state, county) pair is a unique group => 4 rows
  expect_equal(nrow(result), 4L)
  expect_true("state"       %in% names(result))
  expect_true("county"      %in% names(result))
  expect_true("estimate"    %in% names(result))
  expect_true("estimate_se" %in% names(result))

  # For each single-cell group, SE = sqrt(w^2 * se^2) / sum(w)^2 = se (trivial)
  # PA/A: sqrt(100^2 * 0.05^2 / 100^2) = 0.05
  expect_equal(result$estimate_se[result$county == "A" & result$state == "PA"],
               0.05, tolerance = 1e-10)
  # NY/D: sqrt(50^2 * 0.02^2 / 50^2) = 0.02
  expect_equal(result$estimate_se[result$county == "D" & result$state == "NY"],
               0.02, tolerance = 1e-10)

  # Reset flag
  .calibratedMRP_env$poststratify_se_warned <- FALSE
})
