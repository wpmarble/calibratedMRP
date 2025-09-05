
library(testthat)

# logit_shift_internal ----------------------------------------------------

test_that("logit_shift_internal shifts probabilities correctly to match target", {
  # Example: input probs near 0.3 on average
  set.seed(123)
  x <- plogis(rnorm(100, qlogis(0.3), 0.5))
  w <- runif(length(x), .2, .8)

  # Target: 0.5 (should require positive shift)
  target <- 0.5

  shift <- logit_shift_internal(x, w, target)
  shifted <- plogis(qlogis(x) + shift)
  est <- weighted.mean(shifted, w)

  expect_equal(est, target, tolerance = 1e-4)
})

test_that("logit_shift_internal warns for extreme target values", {
  x <- plogis(rnorm(100, qlogis(0.5), 0.5))

  expect_warning(logit_shift_internal(x, target = 0))
  expect_warning(logit_shift_internal(x, target = 1))
})

test_that("logit_shift_internal errors for out-of-bounds target", {
  x <- plogis(rnorm(100, qlogis(0.5), 0.5))

  expect_error(logit_shift_internal(x, target = -0.1))
  expect_error(logit_shift_internal(x, target = 1.1))
})




# logit_shift_single ------------------------------------------------------

test_that("logit_shift_single computes correct shifts per geography", {
  # Example poststrat table: 2 counties
  set.seed(123)
  ps <- tibble::tibble(
    county = rep(c("A", "B"), each = 100),
    pred = plogis(rnorm(200, qlogis(0.3), 0.4)),
    weight = runif(length(pred), .2, .8)
  )

  # True targets: slightly higher than raw means
  calib <- tibble::tibble(
    county = c("A", "B"),
    target = c(0.5, 0.4)
  )

  shifts <- logit_shift_single(
    ps_table = ps,
    outcome = "pred",
    weight = "weight",
    geography = "county",
    calib_target = calib,
    calib_var = "target"
  )

  expect_named(shifts, c("county", "pred_shift"))
  expect_equal(nrow(shifts), 2)

  # Check that the shifts actually adjust to match target
  for (i in 1:nrow(shifts)) {
    g <- shifts$county[i]
    shift_val <- shifts$pred_shift[i]
    ps_g <- dplyr::filter(ps, county == g)
    adjusted <- logit_shift_intercept(ps_g$pred, shift_val)
    avg <- weighted.mean(adjusted, ps_g$weight)
    target <- dplyr::filter(calib, county == g)$target
    expect_equal(avg, target, tolerance = 1e-4)
  }
})

test_that("logit_shift_single handles missing targets gracefully", {
  ps <- tibble::tibble(
    county = rep("A", 50),
    pred = plogis(rnorm(50, qlogis(0.3), 0.4)),
    weight = 1
  )
  # Missing target for "A"
  calib <- tibble::tibble(
    county = "A",
    target = NA
  )
  expect_warning(
    shifts <- logit_shift_single(
      ps_table = ps,
      outcome = "pred",
      weight = "weight",
      geography = "county",
      calib_target = calib,
      calib_var = "target"
    )
  )
  expect_equal(shifts$pred_shift, 0)
})
