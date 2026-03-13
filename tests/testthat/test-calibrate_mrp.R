
library(testthat)

# calibrate_mrp (integration tests) --------------------------------------
# These tests require a brms model fixture. Skip if not available.

fixture_path <- file.path(
  rprojroot::find_package_root_file(),
  "test-top-level.rdata"
)

skip_if_no_fixture <- function() {
  skip_if_not(
    file.exists(fixture_path),
    message = sprintf("Fixture file not found: %s", fixture_path)
  )
  # Also skip if the saved brms model is incompatible with current brms version
  env <- new.env()
  load(fixture_path, envir = env)
  compat <- tryCatch(
    {
      brms::posterior_epred(
        env$mod,
        newdata = env$ps_table[1:2, ],
        resp = "presvote2020twoparty",
        allow_new_levels = TRUE,
        draw_ids = 1:2
      )
      TRUE
    },
    error = function(e) FALSE
  )
  skip_if_not(compat, message = "brms model fixture incompatible with current brms version")
}


test_that("calibrate_mrp plugin: calibrated estimates match targets", {
  skip_if_no_fixture()
  load(fixture_path)  # loads: mod, univ_mod, ps_table, targets

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 50)

  plugin <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "plugin"
  )

  # Poststratify calibrated estimates and compare to targets
  calib_res <- poststratify(
    plugin$results,
    weight = est_n,
    outcomes = "presvote2020twoparty_calib",
    ses = FALSE,
    by = countyfips
  )
  calib_res <- dplyr::left_join(calib_res, targets, by = "countyfips")
  diffs <- abs(calib_res$presvote2020twoparty_calib - calib_res$presvote2020twoparty)

  expect_true(all(diffs < 1e-4))
})


test_that("calibrate_mrp plugin: output structure is correct", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 50)

  plugin <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "plugin"
  )

  expect_true(is.list(plugin))
  expect_true("results" %in% names(plugin))
  expect_true("logit_shifts" %in% names(plugin))
  expect_true("method" %in% names(plugin))
  expect_equal(plugin$method, "plugin")
})


test_that("calibrate_mrp plugin: class is calibrated_summary", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 50)

  plugin <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "plugin"
  )

  expect_s3_class(plugin, "calibrated_summary")
  expect_s3_class(plugin, "calibratedMRP")
})


test_that("calibrate_mrp plugin: keep_all_ps_vars=FALSE drops non-essential columns", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 50)

  plugin <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "plugin",
    keep_all_ps_vars = FALSE
  )

  # Should have .cellid, geography, weight, outcome columns -- not all ps_table cols
  expect_true(".cellid" %in% names(plugin$results))
  expect_true("countyfips" %in% names(plugin$results))
  expect_true("est_n" %in% names(plugin$results))
  # Should have fewer columns than original ps_table
  expect_lt(ncol(plugin$results), ncol(ps_table))
})


test_that("calibrate_mrp bayes: per-draw calibrated estimates match targets", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 10)  # few draws for speed

  bayes <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "bayes",
    posterior_summary = FALSE
  )

  # Check each draw
  maxdiffs <- numeric(length(draw_ids))
  for (i in seq_along(draw_ids)) {
    calib_res <- poststratify(
      bayes$results |> dplyr::filter(.draw == draw_ids[i]),
      weight = est_n,
      outcomes = "presvote2020twoparty_calib",
      ses = FALSE,
      by = countyfips
    )
    calib_res <- dplyr::left_join(calib_res, targets, by = "countyfips")
    maxdiffs[i] <- max(abs(
      calib_res$presvote2020twoparty_calib - calib_res$presvote2020twoparty
    ))
  }

  expect_true(all(maxdiffs < 1e-4))
})


test_that("calibrate_mrp bayes: class is calibrated_draws when posterior_summary=FALSE", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 5)

  bayes <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "bayes",
    posterior_summary = FALSE
  )

  expect_s3_class(bayes, "calibrated_draws")
  expect_s3_class(bayes, "calibratedMRP")
})


test_that("calibrate_mrp bayes: posterior_summary=TRUE returns calibrated_summary", {
  skip_if_no_fixture()
  load(fixture_path)

  set.seed(1234)
  draw_ids <- sample(seq_len(brms::ndraws(mod)), 5)

  bayes_sum <- calibrate_mrp(
    model = mod,
    ps_table = ps_table,
    outcomes = "presvote2020twoparty",
    weight = "est_n",
    targets = targets,
    geography = "countyfips",
    draw_ids = draw_ids,
    method = "bayes",
    posterior_summary = TRUE
  )

  expect_s3_class(bayes_sum, "calibrated_summary")
})
