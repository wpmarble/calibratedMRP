



# model <- mod
# weight <- "est_n"
# geography <- "countyfips"
# method <- "plugin"
# uncertainty <- "approximate"
#
#
# set.seed(1234); ps_table <- ps_cty %>% filter(state == "PA") %>% slice_sample(n = 2e3)
# targets <- targets %>% rename(presvote2020twoparty = pres2020_2pty)
#
#
# univ_form <- formula$forms$presvote2020twoparty
# univ_mod <- brm(formula = univ_form,
#                 data = surv,
#                 family = bernoulli,
#                 prior = prior(normal(0, 5), class = b),
#                 chains = 4,
#                 cores = 4,
#                 iter = 600,
#                 backend = "cmdstanr",
#                 adapt_delta = .99,
#                 max_treedepth = 12)
# save(mod, univ_mod, ps_table, targets, file = "test-top-level.rdata")


load("test-top-level.rdata")

library(dplyr)
# library(calibratedMRP)
devtools::load_all()

# parallel
options(mc.cores = parallel::detectCores())
library(future)
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 3e3 * 1024^2) # 3GB memory for `future`


set.seed(1234)
draw_ids <- sample(1:brms::ndraws(mod), 100)

# Multivariate Outcomes ---------------------------------------------------


## Test plugin -------------------------------------------------------------

plugin <- calibrate_mrp(model = mod,
                        ps_table = ps_table,
                        # outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
                        outcomes = "presvote2020twoparty",
                        weight = "est_n",
                        targets = targets,
                        geography = "countyfips",
                        draw_ids = draw_ids,
                        method = "plugin")



# Make sure calibration worked
calib_res <- poststratify(plugin$results, weight = est_n,
                          outcomes = "presvote2020twoparty_calib",
                          ses = FALSE, by = countyfips)
calib_res <- left_join(calib_res, targets, by = "countyfips")
summary(calib_res$presvote2020twoparty_calib - calib_res$presvote2020twoparty)



## Test bayes ---------------------------------------------------------------

# full bayes, return draws
bayes <- calibrate_mrp(model = mod,
                      ps_table = ps_table,
                      # outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
                      outcomes = "presvote2020twoparty",
                      weight = "est_n",
                      targets = targets,
                      geography = "countyfips",
                      draw_ids = draw_ids,
                      method = "bayes",
                      posterior_summary = FALSE)
class(bayes)
attributes(bayes)

bayes_sum <- calibrate_mrp(model = mod,
                           ps_table = ps_table,
                           # outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
                           weight = "est_n",
                           targets = targets,
                           geography = "countyfips",
                           draw_ids = draw_ids,
                           method = "bayes",
                           posterior_summary = TRUE)

# Unit test - each draw should be calibrated
test_that(desc = "Test that draw-by-draw calibration worked properly", code = {
  maxdiffs <- numeric(length(unique(bayes$results$.draw)))

  for (i in seq_along(draw_ids)) {
    calib_res <- poststratify(
      bayes$results %>% filter(.draw == draw_ids[i]),
      weight = est_n,
      outcomes = "presvote2020twoparty_calib",
      ses = FALSE,
      by = countyfips
    )
    calib_res <- left_join(calib_res, targets, by = "countyfips")
    maxdiffs[i] <- max(abs(
      calib_res$presvote2020twoparty_calib - calib_res$presvote2020twoparty
    ))
  }
  testthat::expect_true(all(maxdiffs < 1e-5))
})





# Univariate outcomes -----------------------------------------------------


## Bayes estimator -------------------------------------------------------

draw_id2 <- sample(1:brms::ndraws(univ_mod), 100)
univ_bayes <- calibrate_mrp(model = univ_mod,
                           ps_table = ps_table,
                           outcomes = "presvote2020twoparty",
                           weight = "est_n",
                           targets = targets,
                           draw_ids = draw_id2,
                           geography = "countyfips",
                           method = "bayes",
                           posterior_summary = TRUE)

# Make sure calibration worked
res <- poststratify(univ_bayes$results, weight = est_n,
                    outcomes = "presvote2020twoparty_calib_mean",
                    ses = FALSE, by = countyfips)
res <- left_join(res, targets, by = "countyfips")
summary(res$presvote2020twoparty_calib_mean - res$presvote2020twoparty)
test_that(desc = "Test univariate bayes calibration worked properly", code = {
  maxdiff <- max(abs(res$presvote2020twoparty_calib_mean - res$presvote2020twoparty))
  testthat::expect_true(maxdiff < 1e-5)
})




## Plugin estimator ---------------------------------------------------------

univ_plugin <- calibrate_mrp(model = univ_mod,
                             ps_table = ps_table,
                             outcomes = "presvote2020twoparty",
                             weight = "est_n",
                             targets = targets,
                             geography = "countyfips",
                             method = "plugin")

# Make sure calibration worked
univ_plugin <- poststratify(univ_plugin$results, weight = est_n,
                     outcomes = "presvote2020twoparty_calib",
                     ses = FALSE, by = countyfips)
univ_plugin <- left_join(univ_plugin, targets, by = "countyfips")
summary(univ_plugin$presvote2020twoparty_calib - univ_plugin$presvote2020twoparty)
test_that(desc = "Test univariate plugin calibration worked properly", code = {
  maxdiff <- max(abs(univ_plugin$presvote2020twoparty_calib - univ_plugin$presvote2020twoparty))
  testthat::expect_true(maxdiff < 1e-5)
})

