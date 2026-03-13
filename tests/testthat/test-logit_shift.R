
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

test_that("logit_shift_single warns once when predictions contain exact 0 or 1", {
  ps <- tibble::tibble(
    county = rep("A", 5),
    pred = c(0, 0.3, 0.5, 0.7, 1),
    weight = 1
  )
  calib <- tibble::tibble(county = "A", target = 0.5)

  expect_warning(
    shifts <- logit_shift_single(
      ps_table = ps,
      outcome = "pred",
      weight = "weight",
      geography = "county",
      calib_target = calib,
      calib_var = "target"
    ),
    "exact 0 or 1"
  )

  # Should still produce a valid shift
  expect_true(is.finite(shifts$pred_shift))
})


# logit_shift_intercept ---------------------------------------------------

test_that("logit_shift_intercept normal operation matches manual calculation", {
  # plogis(qlogis(0.5) + 1) ≈ 0.7310586
  result <- logit_shift_intercept(0.5, 1)
  expected <- plogis(qlogis(0.5) + 1)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("logit_shift_intercept is vectorized", {
  x <- c(0.2, 0.5, 0.8)
  a <- 0.5
  result <- logit_shift_intercept(x, a)
  expected <- plogis(qlogis(x) + a)
  expect_equal(result, expected, tolerance = 1e-6)
  expect_length(result, 3)
})

test_that("logit_shift_intercept with shift of zero returns input unchanged", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  result <- logit_shift_intercept(x, 0)
  expect_equal(result, x, tolerance = 1e-10)
})

test_that("logit_shift_intercept propagates NAs", {
  x <- c(0.2, NA, 0.8)
  result <- logit_shift_intercept(x, 1)
  expect_false(is.na(result[1]))
  expect_true(is.na(result[2]))
  expect_false(is.na(result[3]))
})

test_that("logit_shift_intercept clamps x=0 silently and returns finite value", {
  result <- logit_shift_intercept(0, 1)
  expect_true(is.finite(result))
  expect_gt(result, 0)
})

test_that("logit_shift_intercept clamps x=1 silently and returns finite value", {
  result <- logit_shift_intercept(1, 1)
  expect_true(is.finite(result))
  expect_lt(result, 1)
})

test_that("logit_shift_intercept errors for x out of [0, 1]", {
  expect_error(logit_shift_intercept(-0.1, 1))
  expect_error(logit_shift_intercept(1.1, 1))
  expect_error(logit_shift_intercept(c(0.5, -0.01), 1))
  expect_error(logit_shift_intercept(c(0.5, 1.5), 1))
})



# logit_shift (exported wrapper) ------------------------------------------

test_that("logit_shift handles multiple outcomes", {
  set.seed(42)
  ps <- tibble::tibble(
    county = rep(c("A", "B"), each = 100),
    voteshare = plogis(rnorm(200, qlogis(0.4), 0.5)),
    turnout = plogis(rnorm(200, qlogis(0.3), 0.4)),
    weight = runif(200, 0.2, 0.8)
  )

  targets <- tibble::tibble(
    county = c("A", "B"),
    voteshare = c(0.55, 0.45),
    turnout = c(0.65, 0.55)
  )

  shifts <- logit_shift(
    ps_table = ps,
    outcomes = c("voteshare", "turnout"),
    targets = targets,
    weight = weight,
    geography = county
  )

  # Should have geography column plus one shift column per outcome
  expect_true("county" %in% names(shifts))
  expect_true("voteshare_shift" %in% names(shifts))
  expect_true("turnout_shift" %in% names(shifts))
  expect_equal(nrow(shifts), 2)
})

test_that("logit_shift infers outcomes from column overlap with message", {
  set.seed(42)
  ps <- tibble::tibble(
    county = rep(c("A", "B"), each = 50),
    voteshare = plogis(rnorm(100, qlogis(0.4), 0.5)),
    weight = runif(100, 0.2, 0.8)
  )

  targets <- tibble::tibble(
    county = c("A", "B"),
    voteshare = c(0.55, 0.45)
  )

  expect_message(
    shifts <- logit_shift(
      ps_table = ps,
      outcomes = NULL,
      targets = targets,
      weight = weight,
      geography = county
    ),
    "inferring outcomes"
  )

  expect_true("voteshare_shift" %in% names(shifts))
})

test_that("logit_shift errors when outcomes=NULL and no column overlap", {
  ps <- tibble::tibble(
    county = rep(c("A", "B"), each = 50),
    pred_vote = plogis(rnorm(100, qlogis(0.4), 0.5)),
    weight = runif(100, 0.2, 0.8)
  )

  # targets has different column names (no overlap with ps beyond geography)
  targets <- tibble::tibble(
    county = c("A", "B"),
    actual_vote = c(0.55, 0.45)
  )

  expect_error(
    logit_shift(
      ps_table = ps,
      outcomes = NULL,
      targets = targets,
      weight = weight,
      geography = county
    ),
    "No `outcomes` provided and no overlap"
  )
})
