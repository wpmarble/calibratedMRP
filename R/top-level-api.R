
#
# model <- mod
# ps_table <- ps_cty %>% slice(1:1000)
# weight <- "est_n"
# # targets <- targets %>% rename(presvote2020twoparty = pres2020_2pty)
# geography <- "countyfips"
# method <- "plugin"
# uncertainty <- "approximate"
#
#
#
# blah <- calibrate_mrp(model = mod,
#                       ps_table = ps_table,
#                       # outcomes = c("bidenlegitimate", "presvote2020twoparty", "bidenappr"),
#                       weight = "est_n",
#                       targets = targets,
#                       geography = "countyfips",
#                       draw_ids = 1:100,
#                       method = "plugin")



#' Calibrate Multilevel Regression and Poststratification Estimates
#'
#' This is the top-level function for generating calibrated MRP estimates. It takes
#' a fitted multivariate `brms` model, a poststratification table, and calibration
#' targets, and returns a data frame with cell-level calibrated predictions.
#'
#' Currently, the function supports a plug-in estimation approach where the posterior
#' mean and standard deviation are computed for each outcome in each cell, and
#' calibration is applied via logit intercept shifts to match known targets.
#' Full Bayesian calibration is not yet implemented.
#'
#' @param model A fitted `brmsfit` object. This will typically be a model with
#'   with multivariate binary outcomes. The model must include random intercepts
#'   for `geography`.
#' @param ps_table A data frame representing the poststratification table. Each row corresponds to a population cell.
#' @param weight A **string** giving the name of the column in `ps_table` containing population weights.
#' @param targets A data frame containing calibration targets. Must include `geography` and all outcome variables to be calibrated.
#' @param geography A **string** giving the name of the geography variable.
#'  Must be present in `ps_table` and `targets`, and used as a random intercept in `model`.
#' @param outcomes Optional character vector of outcome variables to calibrate.
#'  If `NULL`, the outcome names are inferred from the model formula.
#' @param method Calibration method, either `"plugin"` or `"bayes"` (currently not implemented).
#'  Plug-in estimates use posterior means of predictions and correlations across outcomes
#'  to compute logit shifts for calibration. Bayesian estimates compute the logit shifts
#'  separately for each posterior draw, which are then summarized.
#' @param uncertainty Method for uncertainty estimation. One of `"approximate"` (default), `"bayes"`, or `"none"`. Currently ignored.
#' @param draw_ids Optional vector of posterior draw indices to use. Defaults to all posterior draws.
#'
#' @return A data frame with the same rows as `ps_table` and additional columns for each calibrated outcome.
#'
#' @export
#'
#' @examples
#' # See vignettes.

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


  # capture NSE inputs
  weight_quo <- rlang::enquo(weight)
  geo_quo <- rlang::enquo(geography)
  outcomes_quo <- rlang::enquo(outcomes)

  # resolve to strings
  weight_var <- rlang::as_name(weight_quo)
  geo_var <- rlang::as_name(geo_quo)


  if (is.null(outcomes)) {
    outcomes <- mod$formula[[2]]
    rlang::inform(c("No `outcomes` provided, defaulting to outcome variables from the model formula: ",
                    "i" = paste(outcomes, collapse = ", ")))
  }
  if (!any(outcomes %in% names(targets))) {
      rlang::abort(c("At least one variable in `outcomes` must be present in `targets` data frame to perform calibration.",
                     "*" = sprintf("Outcome variables: %s", paste(outcomes, collapse = ", ")),
                     "*" = sprintf("`targets` variables: %s", paste(names(targets), collapse = ", ")),
                     "i" = "brms automatically removes underscores from variable names; maybe you need to rename variables in `targets` to match?"))
  } else {
    calib_vars <- intersect(outcomes, names(targets))
    rlang::inform(c("Using the following variables for calibration: ",
                    "i" = paste(calib_vars, collapse = ", ")))
  }

  # extract covariance from model
  covs <- get_re_covariance(model = mod, group = geography, tidy = FALSE, draw_ids = draw_ids)

  # calculate predictions
  if (method == "plugin") {
    ps_table <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        draw_ids = draw_ids,
                                        summarize = TRUE)
    covs <- apply(covs, c(2,3), mean, simplify = TRUE)

    # calculate logit shifts for observed variables
    shifts <- logit_shift(ps_table,
                          outcomes = calib_vars,
                          calib_target = targets,
                          weight = !!weight_var,
                          geography = !!geo_var,
                          calib_vars = calib_vars)

    # impute logit shifts for unobserved variables
    shifts <- logit_shift_aux(shifts, shift_vars = calib_vars, cov = covs)

    # generate calibrated probs
    ps_table <- calibrate_preds(ps_table = ps_table,
                                shifts = shifts,
                                preds = outcomes,
                                geography = !!geo_var)

    ps_table
  }

  # TODO: full bayes version
}
