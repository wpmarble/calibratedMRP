



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
library(calibratedMRP)

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
                           draw_ids = 1:50,
                           method = "bayes",
                           posterior_summary = TRUE)

# Unit test - each draw should be calibrated
test_that(desc = "Test that draw-by-draw calibration worked properly", code = {
  maxdiffs <- numeric(length(unique(bayes$results$.draw)))
  for (i in seq_along(unique(bayes$results$.draw))) {
    calib_res <- poststratify(
      bayes$results %>% filter(.draw == i),
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


## Plug-in estimator -------------------------------------------------------


blah <- generate_cell_estimates(univ_mod, ps_table, summarize = TRUE)
shfits <- logit_shift(blah,
                     # outcomes = "presvote2020twoparty",
                     targets = targets,
                     weight = "est_n",
                     geography = "countyfips")
calib <- calibrate()


## Bayes estimator ---------------------------------------------------------

blah <- generate_cell_estimates(univ_mod, ps_table, summarize = FALSE)
debug(generate_cell_estimates)


## Test calibrate_mrp() with single outcome
what <- calibrate_mrp(model = univ_mod,
                      ps_table = ps_table,
                      outcomes = "presvote2020twoparty",
                      weight = "est_n",
                      targets = targets,
                      geography = "countyfips",
                      draw_ids = 1:100,
                      method = "plugin")

