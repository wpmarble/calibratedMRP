
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


# poststratify.calibrated_summary ------------------------------------------

# Helper: build a minimal calibrated_summary mock
make_cs_obj <- function() {
  results <- tibble::tibble(
    .cellid        = 1:4,
    state          = c("PA", "PA", "OH", "OH"),
    est_n          = c(100L, 200L, 150L, 50L),
    voteshare_mean = c(0.55, 0.60, 0.48, 0.52),
    voteshare_se   = c(0.03, 0.02, 0.04, 0.03),
    turnout_mean   = c(0.70, 0.65, 0.60, 0.55),
    turnout_se     = c(0.05, 0.04, 0.06, 0.07)
  )
  structure(list(results = results),
            class = c("calibratedMRP", "calibrated_summary", "list"))
}

test_that("poststratify.calibrated_summary auto-detects outcomes from _mean/_se pairs", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, weight = est_n, by = state)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)
  expect_true("voteshare_mean" %in% names(result))
  expect_true("voteshare_se"   %in% names(result))
  expect_true("turnout_mean"   %in% names(result))
  expect_true("turnout_se"     %in% names(result))
  expect_true("n"              %in% names(result))
})


test_that("poststratify.calibrated_summary column order is correct", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, weight = est_n, by = state)

  expected_cols <- c("state", "voteshare_mean", "voteshare_se",
                     "turnout_mean", "turnout_se", "n")
  expect_equal(names(result), expected_cols)
})


test_that("poststratify.calibrated_summary weighted mean matches manual computation", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, weight = est_n, by = state)

  # PA: cells with est_n 100 and 200
  pa_wm <- weighted.mean(c(0.55, 0.60), c(100, 200))
  oh_wm <- weighted.mean(c(0.48, 0.52), c(150, 50))

  expect_equal(result$voteshare_mean[result$state == "PA"], pa_wm, tolerance = 1e-12)
  expect_equal(result$voteshare_mean[result$state == "OH"], oh_wm, tolerance = 1e-12)
})


test_that("poststratify.calibrated_summary SE propagation matches formula", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, weight = est_n, by = state)

  # PA: sqrt(sum(w^2 * se^2) / sum(w)^2)
  w_pa  <- c(100, 200)
  se_pa <- c(0.03, 0.02)
  expected_pa_se <- sqrt(sum(w_pa^2 * se_pa^2) / sum(w_pa)^2)

  expect_equal(result$voteshare_se[result$state == "PA"], expected_pa_se, tolerance = 1e-12)
})


test_that("poststratify.calibrated_summary with explicit outcomes subsets correctly", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, outcomes = "voteshare", weight = est_n, by = state)

  expect_true("voteshare_mean" %in% names(result))
  expect_true("voteshare_se"   %in% names(result))
  expect_false("turnout_mean"  %in% names(result))
  expect_false("turnout_se"    %in% names(result))
})


test_that("poststratify.calibrated_summary errors on invalid explicit outcomes", {
  cs_obj <- make_cs_obj()
  expect_error(
    poststratify(cs_obj, outcomes = "nonexistent", weight = est_n, by = state),
    "missing from"
  )
})


test_that("poststratify.calibrated_summary is numerically equivalent to manual default call", {
  # Manual equivalent: call default on results with explicit _mean cols and ses=TRUE
  cs_obj <- make_cs_obj()
  auto_result <- poststratify(cs_obj, outcomes = "voteshare", weight = est_n, by = state)

  # Manual SE formula applied directly
  res <- cs_obj$results
  pa  <- res[res$state == "PA", ]
  oh  <- res[res$state == "OH", ]

  pa_mean <- weighted.mean(pa$voteshare_mean, pa$est_n)
  oh_mean <- weighted.mean(oh$voteshare_mean, oh$est_n)
  pa_se   <- sqrt(sum(pa$est_n^2 * pa$voteshare_se^2) / sum(pa$est_n)^2)
  oh_se   <- sqrt(sum(oh$est_n^2 * oh$voteshare_se^2) / sum(oh$est_n)^2)

  expect_equal(auto_result$voteshare_mean[auto_result$state == "PA"], pa_mean, tolerance = 1e-12)
  expect_equal(auto_result$voteshare_mean[auto_result$state == "OH"], oh_mean, tolerance = 1e-12)
  expect_equal(auto_result$voteshare_se[auto_result$state == "PA"],   pa_se,   tolerance = 1e-12)
  expect_equal(auto_result$voteshare_se[auto_result$state == "OH"],   oh_se,   tolerance = 1e-12)
})


test_that("poststratify.calibrated_summary weight sums are correct", {
  cs_obj <- make_cs_obj()
  result <- poststratify(cs_obj, weight = est_n, by = state)

  expect_equal(result$n[result$state == "PA"], 300)
  expect_equal(result$n[result$state == "OH"], 200)
})


# poststratify.calibrated_draws --------------------------------------------

# Helper: build a minimal calibrated_draws mock
make_cd_obj <- function() {
  results <- tibble::tibble(
    .cellid   = rep(1:4, 4),
    .draw     = rep(1:4, each = 4),
    state     = rep(c("PA", "PA", "OH", "OH"), 4),
    est_n     = rep(c(100L, 200L, 150L, 50L), 4),
    voteshare = c(
      0.55, 0.60, 0.48, 0.52,
      0.57, 0.62, 0.50, 0.54,
      0.53, 0.58, 0.46, 0.50,
      0.56, 0.61, 0.49, 0.53
    ),
    turnout   = c(
      0.70, 0.65, 0.60, 0.55,
      0.72, 0.67, 0.62, 0.57,
      0.68, 0.63, 0.58, 0.53,
      0.71, 0.66, 0.61, 0.56
    )
  )
  structure(list(results = results),
            class = c("calibratedMRP", "calibrated_draws", "list"))
}


test_that("poststratify.calibrated_draws summarized returns mean and SD columns", {
  cd_obj <- make_cd_obj()
  result <- poststratify(cd_obj, outcomes = c("voteshare", "turnout"),
                          weight = est_n, by = state)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)
  expect_true("voteshare_mean" %in% names(result))
  expect_true("voteshare_sd"   %in% names(result))
  expect_true("turnout_mean"   %in% names(result))
  expect_true("turnout_sd"     %in% names(result))
  expect_true("n"              %in% names(result))
})


test_that("poststratify.calibrated_draws summarized column order is correct", {
  cd_obj <- make_cd_obj()
  result <- poststratify(cd_obj, outcomes = c("voteshare", "turnout"),
                          weight = est_n, by = state)

  expected_cols <- c("state", "voteshare_mean", "voteshare_sd",
                     "turnout_mean", "turnout_sd", "n")
  expect_equal(names(result), expected_cols)
})


test_that("poststratify.calibrated_draws posterior_summary=FALSE returns per-draw rows", {
  cd_obj <- make_cd_obj()
  result <- poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state,
                          posterior_summary = FALSE)

  # 4 draws x 2 states = 8 rows
  expect_equal(nrow(result), 8L)
  expect_true(".draw"     %in% names(result))
  expect_true("state"     %in% names(result))
  expect_true("voteshare" %in% names(result))
  expect_true("n"         %in% names(result))
  # .draw should be first
  expect_equal(names(result)[1], ".draw")
})


test_that("poststratify.calibrated_draws per-draw means match manual weighted.mean", {
  cd_obj <- make_cd_obj()
  result <- poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state,
                          posterior_summary = FALSE)

  # Manual: draw 1, PA
  res <- cd_obj$results
  d1_pa <- res[res$.draw == 1 & res$state == "PA", ]
  expected_d1_pa <- weighted.mean(d1_pa$voteshare, d1_pa$est_n)

  actual_d1_pa <- result$voteshare[result$.draw == 1 & result$state == "PA"]
  expect_equal(actual_d1_pa, expected_d1_pa, tolerance = 1e-12)
})


test_that("poststratify.calibrated_draws summarized mean matches mean of per-draw means", {
  cd_obj <- make_cd_obj()
  draws_result <- poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state,
                                posterior_summary = FALSE)
  summary_result <- poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state,
                                  posterior_summary = TRUE)

  # Mean of per-draw means for PA
  pa_draw_means <- draws_result$voteshare[draws_result$state == "PA"]
  expected_pa_mean <- mean(pa_draw_means)
  expected_pa_sd   <- stats::sd(pa_draw_means)

  actual_pa_mean <- summary_result$voteshare_mean[summary_result$state == "PA"]
  actual_pa_sd   <- summary_result$voteshare_sd[summary_result$state == "PA"]

  expect_equal(actual_pa_mean, expected_pa_mean, tolerance = 1e-12)
  expect_equal(actual_pa_sd,   expected_pa_sd,   tolerance = 1e-12)
})


test_that("poststratify.calibrated_draws errors if ses=TRUE", {
  cd_obj <- make_cd_obj()
  expect_error(
    poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state, ses = TRUE),
    "ses.*argument is not supported"
  )
})


test_that("poststratify.calibrated_draws errors if outcomes is NULL or missing", {
  cd_obj <- make_cd_obj()
  expect_error(
    poststratify(cd_obj, weight = est_n, by = state),
    "'outcomes' must be specified"
  )
})


test_that("poststratify.calibrated_draws errors if outcomes not in results", {
  cd_obj <- make_cd_obj()
  expect_error(
    poststratify(cd_obj, outcomes = "nonexistent", weight = est_n, by = state),
    "missing from"
  )
})


# TC-19 and TC-20 — added in Run 2 amendment to close audit gap ---------------

test_that("TC-19: NAMESPACE contains S3method entries for new methods", {
  # Path resolution: testthat::test_path() resolves relative to the test file;
  # NAMESPACE is two levels up (tests/testthat/ -> package root).
  ns_path <- testthat::test_path("../../NAMESPACE")
  ns_lines <- readLines(ns_path)
  expect_true(any(grepl("S3method\\(poststratify,calibrated_summary\\)", ns_lines)))
  expect_true(any(grepl("S3method\\(poststratify,calibrated_draws\\)", ns_lines)))
})

test_that("TC-20: README calibrated_draws example runs without error", {
  # Construct the same mock object as the README example
  cd_results_readme <- tibble::tibble(
    .cellid   = rep(1:4, 3),
    .draw     = rep(1:3, each = 4),
    state     = rep(c("PA", "PA", "OH", "OH"), 3),
    est_n     = rep(c(100, 200, 150, 50), 3),
    voteshare = c(0.55, 0.60, 0.48, 0.52,
                  0.57, 0.62, 0.50, 0.54,
                  0.53, 0.58, 0.46, 0.50)
  )
  cd_obj_readme <- structure(
    list(results = cd_results_readme),
    class = c("calibratedMRP", "calibrated_draws", "list")
  )

  # Direct call (posterior_summary = TRUE)
  expect_no_error({
    out <- poststratify(cd_obj_readme,
                        outcomes = "voteshare",
                        weight   = est_n,
                        by       = state)
  })

  # Per-draw call + manual quantile
  expect_no_error({
    out_draws <- poststratify(cd_obj_readme,
                               outcomes = "voteshare",
                               weight   = est_n,
                               by       = state,
                               posterior_summary = FALSE)
    out_quantiles <- out_draws |>
      dplyr::summarise(
        voteshare_mean = mean(voteshare),
        voteshare_q5   = quantile(voteshare, .05),
        voteshare_q95  = quantile(voteshare, .95),
        .by = state
      )
  })
})
