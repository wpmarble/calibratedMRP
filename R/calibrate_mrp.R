
#' Calibrate Multilevel Regression and Poststratification Estimates
#'
#' This is the top-level function for generating calibrated MRP estimates. It takes
#' a fitted multivariate `brms` model, a poststratification table, and calibration
#' targets, and returns a data frame with calibrated estimates of
#'
#' @details
#' This function supports two estimation approaches, a plug-in estimator and a
#' full Bayesian estimator. The plug-in approach involves computing the posterior
#' mean and standard deviation for each outcome in each cell, then applying the
#' calibration procedure to these cell-level summaries. The plug-in approach uses
#' the posterior mean covariance matrix across outcomes to perform the calibration.
#' The result of the plug-in estimation is a poststratification table with additional
#' calibrated outcome variables, with the suffix "_calib" added to the original
#' outcome names.
#'
#' The full Bayesian approach performs the calibration procedure separately for
#' each draw from the posterior distribution in `model`. The result is a data
#' frame with the number of rows equal to the number of rows in `ps_table` times the
#' number of posterior draws. The data frame includes a `.draw` column indicating
#' which draw from the posterior each estimate comes from. Within each draw,
#' model-based estimates of the target variables match the ground truth exactly.
#' Researchers should summarize across `.draws` in downstream analyses to obtain
#' point estimates and uncertainty intervals.
#'
#' Full Bayesian estimation is more computationally and memory intensive than the
#' plug-in approach. To limit the resources required, calibration can be performed
#' in parallel by splitting `ps_table` into groups defined by rows in `targets`.
#' For example, if you are calibrating to county-level election results stored
#' within `targets`, you can process one county in `ps_table` at at time.
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
#'  If `NULL`, the outcome names are inferred from the `model` formula.
#' @param method Calibration method, either `"plugin"` or `"bayes"`.
#'  Plug-in estimates use posterior means of predictions and correlations across outcomes
#'  to compute logit shifts for calibration. Bayesian estimates compute the logit shifts
#'  separately for each posterior draw, which are then summarized. Defaults to `"plugin"`.
#' @param posterior_summary If `method = "bayes"`, should the function return all
#'  draws from the posterior or summarize and return the posterior mean and SD?
#' @param uncertainty Method for uncertainty estimation. One of `"approximate"` (default),
#' `"bayes"`, or `"none"`. Currently ignored.
#' @param draw_ids Optional vector of posterior draw indices to use. Defaults to all posterior draws.
#'
#' @return
#' A list containing two objects: 1) a data frame called `results` with the same
#' variables `ps_table` plus additional columns for each calibrated outcome, and
#' 2) a data frame called `logit_shifts` storing the estimated intercept-shift
#' parameters that were used to perform the calibration. If `method = "plugin"`,
#' the outputs will be posterior means and posterior standard deviations.
#' If `method = "bayes"` and `posterior_summary = FALSE`, the outputs will include
#' all posterior draws, with a `.draw` column indicating the draw index. In this
#' case value returned has class `calibrated_draws`; otherwise it has class
#' `calibrated_summary`.
#'
#' @importFrom dplyr summarise group_by mutate select across bind_rows full_join row_number
#' @export
#'
#' @examples
#' # See vignettes.
#' @references
#' Marble, William and Joshua Clinton. 2025. "Improving Small-Area Estimates of Public
#' Opinion by Calibrating to Known Population Quantities."
#' \url{https://osf.io/preprints/socarxiv/u3ekq_v1}

calibrate_mrp <- function(model,
                          ps_table,
                          weight,
                          targets,
                          geography, # must be string, in ps_table and REs for this in `model`
                          outcomes = NULL, # defaults to outcomes in `mod`
                          method = "plugin", # or "bayes"
                          posterior_summary = FALSE,
                          uncertainty = "approximate", # or "bayes" or "none"
                          draw_ids = NULL
                          ) {
  if (class(mod) != "brmsfit") rlang::abort("`mod` must be a `brmsfit` object")
  if (!method %in% c("plugin", "bayes")) rlang::abort("`method` must be either 'plugin' or 'bayes'")
  if (is.null(draw_ids)) draw_ids <- seq_len(posterior::ndraws(mod))


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
                                        outcomes = outcomes,
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


    out <- list(results = ps_table, logit_shifts = shifts)
    class(out) <- c("calibratedMRP", "calibrated_summary", class(out))
  }


  # full bayes version
  if (method == "bayes") {

    ps_table <- ps_table |>
      mutate(.rowid = row_number())

    # get posterior draws of cell-level estimates
    ps_draws <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        outcomes = outcomes,
                                        draw_ids = draw_ids,
                                        summarize = FALSE)

    ps_table_clean <- ps_table |>
      select(.rowid, all_of(geo_var), all_of(weight_var))

    res <- list()
    res_shift <- list()

    # initialize progress bar
    # rlang::inform("Calibrating posterior draws...\n")
    pb <- progress::progress_bar$new(
      format = "calibrating posterior draws [:bar] :percent",
      show_after = 0,
      total = length(draw_ids),
      clear = FALSE,
      width = 60
    )
    pb$tick(0)

    for (i in seq_along(draw_ids)){

      ps_draw_i <- as.data.frame(ps_draws[i, , ])
      ps_table_i <- cbind(ps_table_clean, ps_draw_i)
      covs_i <- covs[i , ,]


      shifts <- logit_shift(ps_table_i,
                            outcomes = calib_vars,
                            calib_target = targets,
                            weight = !!weight_var,
                            geography = !!geo_var,
                            calib_vars = calib_vars)

      # impute logit shifts for unobserved variables
      shifts <- logit_shift_aux(shifts, shift_vars = calib_vars, cov = covs_i)

      # generate calibrated probs
      ps_table_i <- calibrate_preds(ps_table = ps_table_i,
                                    shifts = shifts,
                                    preds = outcomes,
                                    geography = !!geo_var)
      res[[i]] <- ps_table_i |>
        mutate(.draw = i)

      res_shift[[i]] <- shifts |> mutate(.draw = i)

      pb$tick()  # update progress bar
    }
    res <- bind_rows(res)
    res <- res |>
      select(-all_of(c(geo_var, weight_var)))

    shifts <- bind_rows(res_shift)

    # summarize posterior if requested
    if (posterior_summary) {
      res <- res |>
        group_by(.rowid) |>
        summarise(across(-.draw, list(
          mean = ~ mean(.x),
          se = ~ sd(.x)
        ), .names = "{.col}_{.fn}"), .groups = "drop")

      shifts <- shifts |>
        summarise(across(-.draw, mean))
    }
    res <- full_join(ps_table, res, by = ".rowid")
    res <- res |> select(-.rowid)

    out <- list(results = res, logit_shifts = shifts)

    if (posterior_summary) {
      class(out) <- c("calibratedMRP", "calibrated_summary", class(out))
    } else {
      class(out) <- c("calibratedMRP", "calibrated_draws", class(out))
    }
  }

  out$method <- method
  out
}

