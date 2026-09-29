
library(testthat)

# Evaluate expr, returning its value and the number of missing-target warnings
count_missing_warnings <- function(expr) {
  n <- 0
  value <- withCallingHandlers(
    expr,
    calibratedMRP_missing_target = function(cnd) {
      n <<- n + 1
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, n = n)
}

# Geographies in ps_table without calibration targets ---------------------

test_that("logit_shift: geography absent from targets gets shift 0 with a warning", {
  ps <- make_test_ps(n_geo = 3)
  targets <- make_test_targets(ps) |> dplyr::filter(geography != "geo_3")

  out <- count_missing_warnings(
    logit_shift(ps, outcomes = c(voteshare, turnout), targets = targets,
                weight = weight, geography = geography)
  )
  expect_equal(out$n, 2)  # one per outcome
  shifts <- out$value
  g3 <- dplyr::filter(shifts, geography == "geo_3")
  expect_equal(nrow(shifts), 3)
  expect_equal(g3$voteshare_shift, 0)
  expect_equal(g3$turnout_shift, 0)
  expect_true(all(dplyr::filter(shifts, geography != "geo_3")$voteshare_shift != 0))
})

test_that("logit_shift: geography with NA target gets shift 0 with a warning", {
  ps <- make_test_ps(n_geo = 3)
  targets <- make_test_targets(ps)
  targets$voteshare[targets$geography == "geo_2"] <- NA

  expect_warning(
    shifts <- logit_shift(ps, outcomes = c(voteshare, turnout), targets = targets,
                          weight = weight, geography = geography),
    class = "calibratedMRP_missing_target"
  )
  g2 <- dplyr::filter(shifts, geography == "geo_2")
  expect_equal(g2$voteshare_shift, 0)
  expect_true(g2$turnout_shift != 0)
})

test_that("logit_shift_single: one aggregated warning naming all missing geographies", {
  ps <- make_test_ps(n_geo = 4, outcomes = "voteshare")
  targets <- make_test_targets(ps, outcomes = "voteshare") |>
    dplyr::filter(geography %in% c("geo_1", "geo_2"))

  w <- list()
  withCallingHandlers(
    logit_shift_single(ps, outcome = "voteshare", weight = "weight",
                       geography = "geography", calib_target = targets,
                       calib_var = "voteshare"),
    calibratedMRP_missing_target = function(cnd) {
      w[[length(w) + 1]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1)
  msg <- conditionMessage(w[[1]])
  expect_match(msg, "geo_3")
  expect_match(msg, "geo_4")
})

test_that("logit_shift: duplicated geographies in targets abort with informative message", {
  ps <- make_test_ps(n_geo = 3)
  targets <- make_test_targets(ps)
  targets <- dplyr::bind_rows(targets, targets[2, ])

  expect_error(
    logit_shift(ps, outcomes = c(voteshare, turnout), targets = targets,
                weight = weight, geography = geography),
    "geo_2"
  )
})


# calibrate_mrp end to end (model-dependent helpers mocked) ----------------

# Outcomes A and B have targets (anchors); C is auxiliary (imputed).
# geo_3 has no row in targets.
make_mock_setup <- function(n_draws = 3) {
  outcomes <- c("A", "B", "C")
  ps <- make_test_ps(n_geo = 3, n_cells_per_geo = 20, outcomes = character(0))
  n <- nrow(ps)
  set.seed(1)
  draws <- array(plogis(rnorm(n_draws * n * 3, qlogis(0.4), 0.6)),
                 dim = c(n_draws, n, 3),
                 dimnames = list(seq_len(n_draws), seq_len(n), outcomes))
  targets <- tibble::tibble(geography = c("geo_1", "geo_2"),
                            A = c(0.55, 0.35), B = c(0.60, 0.30))
  cov <- make_test_cov(outcomes, corr = 0.5, sds = c(1, 1, 1))
  covs <- aperm(simplify2array(rep(list(cov), n_draws)), c(3, 1, 2))
  list(ps = ps, draws = draws, targets = targets, covs = covs, outcomes = outcomes)
}

mock_model_helpers <- function(setup, env = parent.frame()) {
  local_mocked_bindings(
    generate_cell_estimates = function(model, ps_table, outcomes, draw_ids, summarize, ...) {
      d <- setup$draws[draw_ids, , outcomes, drop = FALSE]
      if (!summarize) return(d)
      for (k in outcomes) ps_table[[k]] <- colMeans(d[, , k])
      ps_table
    },
    get_re_covariance = function(model, group, tidy, draw_ids, ...) {
      setup$covs[draw_ids, , , drop = FALSE]
    },
    .env = env
  )
}

fake_model <- structure(list(), class = "brmsfit")

test_that("calibrate_mrp plugin: untargeted geography is left uncalibrated for anchors and aux", {
  setup <- make_mock_setup()
  mock_model_helpers(setup)

  out <- count_missing_warnings(
    calibrate_mrp(fake_model, ps_table = setup$ps, weight = weight,
                  targets = setup$targets, geography = geography,
                  outcomes = setup$outcomes, method = "plugin",
                  draw_ids = 1:3, keep_uncalib = TRUE)
  )
  expect_equal(out$n, 2)  # one per anchor outcome
  res <- out$value

  g3_shift <- dplyr::filter(res$logit_shifts, geography == "geo_3")
  expect_equal(unlist(g3_shift[c("A_shift", "B_shift", "C_shift")]),
               c(A_shift = 0, B_shift = 0, C_shift = 0))

  g3 <- dplyr::filter(res$results, geography == "geo_3")
  for (k in setup$outcomes) expect_equal(g3[[paste0(k, "_calib")]], g3[[k]])

  # targeted geographies are still calibrated
  g1 <- dplyr::filter(res$results, geography == "geo_1")
  expect_equal(weighted.mean(g1$A_calib, g1$weight), 0.55, tolerance = 1e-4)
})

test_that("calibrate_mrp bayes: untargeted geography uncalibrated in every draw, warned once", {
  setup <- make_mock_setup(n_draws = 3)
  mock_model_helpers(setup)

  out <- count_missing_warnings(
    calibrate_mrp(fake_model, ps_table = setup$ps, weight = weight,
                  targets = setup$targets, geography = geography,
                  outcomes = setup$outcomes, method = "bayes",
                  draw_ids = 1:3, keep_uncalib = TRUE, keep_all_ps_vars = TRUE)
  )
  # one aggregated warning per anchor outcome, not per geography x draw
  expect_equal(out$n, 2)
  res <- out$value

  g3_shift <- dplyr::filter(res$logit_shifts, geography == "geo_3")
  expect_equal(nrow(g3_shift), 3)
  expect_true(all(g3_shift$A_shift == 0 & g3_shift$B_shift == 0 & g3_shift$C_shift == 0))

  g3 <- dplyr::filter(res$results, geography == "geo_3")
  for (k in setup$outcomes) expect_equal(g3[[paste0(k, "_calib")]], g3[[k]])

  g1 <- dplyr::filter(res$results, geography == "geo_1")
  for (d in 1:3) {
    gd <- dplyr::filter(g1, .draw == d)
    expect_equal(weighted.mean(gd$B_calib, gd$weight), 0.60, tolerance = 1e-4)
  }
})
