
library(testthat)

# logit_shift_aux ---------------------------------------------------------

test_that("logit_shift_aux hand-computed example matches conditional expectation formula", {
  # 3 outcomes: A, B (observed), C (unobserved)
  # Sigma_uo %*% solve(Sigma_oo) %*% delta_o should give the imputed shift
  outcomes <- c("A", "B", "C")
  cov_mat <- make_test_cov(outcomes, corr = 0.5, sds = c(1, 1, 1))

  # Observed shifts for A and B
  shifts <- tibble::tibble(
    geography = c("geo_1", "geo_2"),
    A_shift = c(0.5, -0.3),
    B_shift = c(0.2, 0.4)
  )

  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = c("A", "B"),
    cov = cov_mat
  )

  expect_true("C_shift" %in% names(result))

  # Manual calculation: Sigma_uo %*% solve(Sigma_oo) %*% delta_o
  cov_uo <- cov_mat["C", c("A", "B"), drop = FALSE]
  cov_oo <- cov_mat[c("A", "B"), c("A", "B")]
  premat <- cov_uo %*% MASS::ginv(cov_oo)

  for (i in 1:2) {
    delta_o <- c(shifts$A_shift[i], shifts$B_shift[i])
    expected <- as.numeric(premat %*% delta_o)
    expect_equal(result$C_shift[i], expected, tolerance = 1e-10)
  }
})


test_that("logit_shift_aux with perfect correlation imputes equal shift", {
  outcomes <- c("obs", "unobs")
  cov_mat <- make_test_cov(outcomes, corr = 1.0, sds = c(1, 1))

  shifts <- tibble::tibble(
    geography = "geo_1",
    obs_shift = 0.7
  )

  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = "obs",
    cov = cov_mat
  )

  expect_equal(result$unobs_shift, 0.7, tolerance = 1e-10)
})


test_that("logit_shift_aux with zero correlation imputes zero shift", {
  outcomes <- c("obs", "unobs")
  cov_mat <- make_test_cov(outcomes, corr = 0.0, sds = c(1, 1))

  shifts <- tibble::tibble(
    geography = "geo_1",
    obs_shift = 0.7
  )

  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = "obs",
    cov = cov_mat
  )

  expect_equal(result$unobs_shift, 0, tolerance = 1e-10)
})


test_that("logit_shift_aux with multiple observed, one unobserved", {
  outcomes <- c("A", "B", "C")
  # Custom correlation matrix
  cor_mat <- matrix(c(1, 0.3, 0.7,
                       0.3, 1, 0.5,
                       0.7, 0.5, 1), nrow = 3, byrow = TRUE)
  cov_mat <- make_test_cov(outcomes, corr = cor_mat, sds = c(1.5, 0.8, 1.2))

  shifts <- tibble::tibble(
    geography = c("geo_1"),
    A_shift = 0.4,
    B_shift = -0.2
  )

  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = c("A", "B"),
    cov = cov_mat
  )

  # Verify against manual calculation
  cov_uo <- cov_mat["C", c("A", "B"), drop = FALSE]
  cov_oo <- cov_mat[c("A", "B"), c("A", "B")]
  expected <- as.numeric(cov_uo %*% MASS::ginv(cov_oo) %*% c(0.4, -0.2))
  expect_equal(result$C_shift, expected, tolerance = 1e-10)
})


test_that("logit_shift_aux with one observed, multiple unobserved", {
  outcomes <- c("obs", "unobs1", "unobs2")
  cor_mat <- matrix(c(1, 0.6, 0.3,
                       0.6, 1, 0.5,
                       0.3, 0.5, 1), nrow = 3, byrow = TRUE)
  cov_mat <- make_test_cov(outcomes, corr = cor_mat, sds = c(1, 1, 1))

  shifts <- tibble::tibble(
    geography = "geo_1",
    obs_shift = 1.0
  )

  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = "obs",
    cov = cov_mat
  )

  expect_true("unobs1_shift" %in% names(result))
  expect_true("unobs2_shift" %in% names(result))

  # unobs1 should get a bigger imputed shift (higher corr with obs)
  expect_gt(abs(result$unobs1_shift), abs(result$unobs2_shift))
})


test_that("logit_shift_aux handles singular covariance via ginv", {
  # Create a singular covariance (outcome C = A + B perfectly)
  outcomes <- c("A", "B", "C")
  # Rank-deficient: C is perfect linear combo
  cov_mat <- matrix(c(1, 0.5, 1.5,
                       0.5, 1, 1.5,
                       1.5, 1.5, 3), nrow = 3)
  dimnames(cov_mat) <- list(outcomes, outcomes)

  shifts <- tibble::tibble(
    geography = "geo_1",
    A_shift = 0.5,
    B_shift = 0.3
  )

  # Should not error -- ginv handles singular matrices
  result <- logit_shift_aux(
    shift = shifts,
    shift_vars = c("A", "B"),
    cov = cov_mat
  )

  expect_true("C_shift" %in% names(result))
  expect_true(is.finite(result$C_shift))
})
