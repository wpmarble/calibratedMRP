
library(testthat)

# calibrate_preds ---------------------------------------------------------

test_that("calibrate_preds produces correct calibrated column", {
  ps <- make_test_ps(n_geo = 2, n_cells_per_geo = 20, outcomes = "voteshare")
  shifts <- make_test_shifts(
    geographies = unique(ps$geography),
    outcomes = "voteshare",
    shift_vals = c(0.3, -0.2)
  )

  result <- calibrate_preds(
    ps_table = ps,
    shifts = shifts,
    preds = "voteshare",
    geography = geography
  )

  expect_true("voteshare_calib" %in% names(result))

  # Verify calibration matches manual calculation
  merged <- dplyr::left_join(ps, shifts, by = "geography")
  expected <- plogis(qlogis(merged$voteshare) + merged$voteshare_shift)
  expect_equal(result$voteshare_calib, expected, tolerance = 1e-10)
})


test_that("calibrate_preds handles multiple predictions", {
  ps <- make_test_ps(n_geo = 2, n_cells_per_geo = 20,
                     outcomes = c("voteshare", "turnout"))
  shifts <- make_test_shifts(
    geographies = unique(ps$geography),
    outcomes = c("voteshare", "turnout"),
    shift_vals = c(0.3, -0.2, 0.1, 0.4)
  )

  result <- calibrate_preds(
    ps_table = ps,
    shifts = shifts,
    preds = c("voteshare", "turnout"),
    geography = geography
  )

  expect_true("voteshare_calib" %in% names(result))
  expect_true("turnout_calib" %in% names(result))
})


test_that("calibrate_preds removes shift columns from output", {
  ps <- make_test_ps(n_geo = 2, n_cells_per_geo = 20, outcomes = "voteshare")
  shifts <- make_test_shifts(
    geographies = unique(ps$geography),
    outcomes = "voteshare",
    shift_vals = 0.5
  )

  result <- calibrate_preds(
    ps_table = ps,
    shifts = shifts,
    preds = "voteshare",
    geography = geography
  )

  expect_false("voteshare_shift" %in% names(result))
})


test_that("calibrate_preds keep_orig=FALSE drops original prediction columns", {
  ps <- make_test_ps(n_geo = 2, n_cells_per_geo = 20, outcomes = "voteshare")
  shifts <- make_test_shifts(
    geographies = unique(ps$geography),
    outcomes = "voteshare",
    shift_vals = 0.5
  )

  result <- calibrate_preds(
    ps_table = ps,
    shifts = shifts,
    preds = "voteshare",
    geography = geography,
    keep_orig = FALSE
  )

  expect_false("voteshare" %in% names(result))
  expect_true("voteshare_calib" %in% names(result))
})


test_that("calibrate_preds handles geography mismatch with NA", {
  ps <- make_test_ps(n_geo = 3, n_cells_per_geo = 10, outcomes = "voteshare")
  # Only provide shifts for 2 of 3 geographies
  shifts <- make_test_shifts(
    geographies = c("geo_1", "geo_2"),
    outcomes = "voteshare",
    shift_vals = c(0.3, -0.2)
  )

  # geo_3 rows should get NA for the calibrated column since no shift exists
  result <- suppressWarnings(calibrate_preds(
    ps_table = ps,
    shifts = shifts,
    preds = "voteshare",
    geography = geography
  ))

  geo3_rows <- result$geography == "geo_3"
  expect_true(all(is.na(result$voteshare_calib[geo3_rows])))
  expect_true(all(!is.na(result$voteshare_calib[!geo3_rows])))
})
