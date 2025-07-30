

calibrate_mrp <- function(model,
                          ps_table,
                          weight,
                          targets,
                          geography, # must be string, in ps_table and REs for this in `model`
                          outcomes = NULL, # defaults to outcomes in `mod`
                          method = "plugin", # or "bayes"
                          uncertainty = "approximate", # or "bayes" or "none"
                          draw_ids = NULL

){
  if (class(mod) != "brmsfit") rlang::abort("`mod` must be a `brmsfit` object")
  if (!method %in% c("plugin", "bayes")) rlang::abort("`method` must be either 'plugin' or 'bayes'")
  if (is.null(draw_ids)) draw_ids <- seq_len(ndraws(mod))
  if (is.null(outcomes)) {
    outcomes <- mod$formula[[2]]
    rlang::inform(c("No `outcomes` provided, defaulting to outcome variables from the model formula: ",
                    "i" = paste(outcomes, collapse = ", ")))
  }
  if (!any(outcomes %in% names(targets))) {
      rlang::abort(c("At least one variable in `outcomes` must be present in `targets` data frame to perform calibration.",
                     "*" = sprintf("Outcome variables: %s", paste(outcomes, collapse = ", ")),
                     "*" = sprintf("Columns in `targets`: %s", paste(names(targets), collapse = ", ")),
                     "i" = "brms automatically removes underscores from variable names; you may need to rename variables in `targets` to match."))
  } else {
    calib_vars <- intersect(outcomes, names(targets))
    rlang::inform(c("Using the following variables for calibration: ",
                    "i" = paste(calib_vars, collapse = ", ")))
  }

  # extract covariance from model
  covs <- get_re_covs(model = mod, group = geography, tidy = FALSE, draw_ids = draw_ids)

  # calculate predictions
  if (method == "plugin") {
    ps_table <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        draw_ids = draw_ids,
                                        summarize = TRUE)
    covs <- apply(covs, c(2,3), mean, simplify = TRUE)

    # calculate logit shifts for observed variables
    shifts <- logit_shift(ps_table, outcomes = calib_vars,
                          calib_target = targets,
                          weight = !!weight,
                          geography = !!geography,
                          calib_vars = calib_vars)

    # impute logit shifts for unobserved variables
    shifts <- logit_shift_aux(shifts, shift_vars = calib_vars, cov = covs)

    # generate calibrated probs
    ps_table <- calibrate_preds(ps_table = ps_table,
                                shifts = shifts,
                                preds = outcomes,
                                geography = !!geography)

    ps_table
  }
}
