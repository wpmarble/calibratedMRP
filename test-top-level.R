



# model <- mod
# weight <- "est_n"
# geography <- "countyfips"
# method <- "plugin"
# uncertainty <- "approximate"


set.seed(1234); ps_table <- ps_cty %>% filter(state == "PA") %>% slice_sample(n = 2e3)
targets <- targets %>% rename(presvote2020twoparty = pres2020_2pty)


# # Test plugin -------------------------------------------------------------

plugin <- calibrate_mrp(model = mod,
                        ps_table = ps_table,
                        outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
                        weight = "est_n",
                        targets = targets,
                        geography = "countyfips",
                        draw_ids = 1:100,
                        method = "plugin")



# Make sure calibration worked
calib_res <- poststratify(plugin$results, weight = est_n,
                          outcomes = "presvote2020twoparty_calib",
                          ses = FALSE, by = countyfips)
calib_res <- left_join(calib_res, targets, by = "countyfips")
summary(calib_res$presvote2020twoparty_calib - calib_res$presvote2020twoparty)



# Test bayes ---------------------------------------------------------------

# full bayes, return draws
bayes <- calibrate_mrp(model = mod,
                      ps_table = ps_table,
                      # outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
                      weight = "est_n",
                      targets = targets,
                      geography = "countyfips",
                      draw_ids = NULL,
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
