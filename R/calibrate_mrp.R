
#' Calibrate Multilevel Regression and Poststratification Estimates
#'
#' This is the top-level function for generating calibrated MRP estimates. It takes
#' a fitted multivariate `brms` model, a poststratification table, and calibration
#' targets, and returns a data frame with calibrated estimates of the outcomes,
#' along with the geography-specific intercept shifts that were applied to each
#' outcome.
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
#'  separately for each posterior draw. Defaults to `"plugin"`, which is less computationally
#'  intensive.
#' @param posterior_summary If `method = "bayes"`, should the function return all
#'  draws from the posterior or summarize and return the posterior mean and SD?
#' @param draw_ids Optional vector of posterior draw indices to use. Defaults to all posterior draws.
#' @param keep_uncalib Keep uncalibrated outcomes?
#' @param keep_all_ps_vars Keep all variables in `ps_table` in the output, or just
#' `geography`, `weight`, outcomes, and `.cellid`? To minimize memory overhead this
#' should be `FALSE` for most uses.
#'
#'
#'
#' @return
#' A list containing two objects: 1) a `results` data frame that corresponds to
#' rows in `ps_table` with plus additional columns added for each calibrated outcome,
#' and 2) a data frame called `logit_shifts` storing the estimated intercept-shift
#' parameters that were used to perform the calibration. If `method = "plugin"`,
#' the outputs will be posterior means and posterior standard deviations.
#' If `method = "bayes"` and `posterior_summary = FALSE`, the outputs will include
#' all posterior draws, with a `.draw` column indicating the draw index. In this
#' case value returned has class `calibrated_draws`; otherwise it has class
#' `calibrated_summary`.
#'
#' If `keep_all_ps_vars = TRUE`, the `results` data frame will include all columns
#' that originally appear in `ps_table`. However, this option is not recommended
#' for large poststratification tables as it is memory intensive, especially when
#' performing full Bayesian inference. If `keep_all_ps_vars = FALSE`, the `results`
#' data frame will include a new column `.cellid` that uniquely identifies the
#' cells in the poststratification table according to their row number in `ps_table`.
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
                          draw_ids = NULL,
                          keep_uncalib = FALSE,
                          keep_all_ps_vars = FALSE
                          ) {
  if (class(model) != "brmsfit") rlang::abort("`model` must be a `brmsfit` object")
  if (!method %in% c("plugin", "bayes")) rlang::abort("`method` must be either 'plugin' or 'bayes'")
  if (is.null(draw_ids)) draw_ids <- seq_len(posterior::ndraws(model))


  # capture NSE inputs
  weight_quo <- rlang::enquo(weight)
  geo_quo <- rlang::enquo(geography)
  outcomes_quo <- rlang::enquo(outcomes)

  # resolve to strings
  weight_var <- rlang::as_name(weight_quo)
  geo_var <- rlang::as_name(geo_quo)


  if (is.null(outcomes)) {
    outcomes <- get_outcomes(model$formula)
    rlang::inform(c("No `outcomes` provided, defaulting to outcome variables from the model formula: ",
                    "*" = paste(outcomes, collapse = ", ")))
  }
  if (!any(outcomes %in% names(targets))) {
      rlang::abort(c("At least one variable in `outcomes` must be present in `targets` data frame to perform calibration.",
                     "*" = sprintf("Outcome variables: %s", paste(outcomes, collapse = ", ")),
                     "*" = sprintf("`targets` variables: %s", paste(names(targets), collapse = ", ")),
                     "*" = "brms automatically removes underscores from variable names; maybe you need to rename variables in `targets` to match?"))
  } else {
    calib_vars <- intersect(outcomes, names(targets))
    rlang::inform(c("Using the following variables for calibration: ",
                    "*" = paste(calib_vars, collapse = ", ")))
  }

  # check if there are any outcomes that aren't in the calibration targets
  auxcalib <-  length(outcomes) > length(calib_vars)

  # extract covariance from model
  if (auxcalib){
    covs <- get_re_covariance(model = model, group = geography, tidy = FALSE, draw_ids = draw_ids)
  }

  # add a row id to ps_table
  if (".cellid" %in% names(ps_table)) {
    rlang::abort("`ps_table` must not contain a column named `.cellid`")
  }
  ps_table <- ps_table |>
    mutate(.cellid = row_number())

  # calculate predictions
  if (method == "plugin") {
    ps_table <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        outcomes = outcomes,
                                        draw_ids = draw_ids,
                                        summarize = TRUE)

    # drop most columns from ps_table
    if (!keep_all_ps_vars){
      ps_table <- ps_table |>
        select(.cellid, all_of(c(geo_var, weight_var, outcomes)))
    }

    # calculate logit shifts for observed variables
    shifts <- logit_shift(ps_table,
                          outcomes = calib_vars,
                          targets = targets,
                          weight = !!weight_var,
                          geography = !!geo_var)

    # impute logit shifts for unobserved variables
    if (auxcalib) {
      covs <- apply(covs, c(2,3), mean, simplify = TRUE)
      shifts <- logit_shift_aux(shifts, shift_vars = calib_vars, cov = covs)
    }

    # generate calibrated probs
    ps_table <- calibrate_preds(ps_table = ps_table,
                                shifts = shifts,
                                preds = outcomes,
                                geography = !!geo_var,
                                keep_orig = keep_uncalib)


    out <- list(results = ps_table, logit_shifts = shifts)
    class(out) <- c("calibratedMRP", "calibrated_summary", class(out))
  }


  # full bayes version
  if (method == "bayes") {

    # get posterior draws of cell-level estimates
    ps_draws <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        outcomes = outcomes,
                                        draw_ids = draw_ids,
                                        summarize = FALSE)

    if (!keep_all_ps_vars){
      ps_table <- ps_table |>
        select(.cellid, all_of(c(geo_var, weight_var)))
    }

    ps_table_clean <- ps_table |>
      select(.cellid, all_of(c(geo_var, weight_var)))

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

      # clunky slicing necessary to handle case with one outcome
      ps_draw_i <- as.data.frame(matrix(ps_draws[i, , , drop = FALSE],
                                        nrow = dim(ps_draws)[2],
                                        dimnames = dimnames(ps_draws)[2:3]))
      ps_table_i <- cbind(ps_table_clean, ps_draw_i)


      shifts <- logit_shift(ps_table_i,
                            outcomes = calib_vars,
                            targets = targets,
                            weight = !!weight_var,
                            geography = !!geo_var)

      # impute logit shifts for unobserved variables
      if (auxcalib) {
        covs_i <- covs[i , ,]
        shifts <- logit_shift_aux(shifts, shift_vars = calib_vars, cov = covs_i)
      }

      # generate calibrated probs
      ps_table_i <- calibrate_preds(ps_table = ps_table_i,
                                    shifts = shifts,
                                    preds = outcomes,
                                    geography = !!geo_var,
                                    keep_orig = keep_uncalib)
      res[[i]] <- ps_table_i |>
        mutate(.draw = draw_ids[i])

      res_shift[[i]] <- shifts |> mutate(.draw = draw_ids[i])

      pb$tick()  # update progress bar
    }

    res <- bind_rows(res)
    res <- res |>
      select(-all_of(c(geo_var, weight_var)))

    shifts <- bind_rows(res_shift)

    # summarize posterior if requested
    if (posterior_summary) {
      res <- res |>
        group_by(.cellid) |>
        summarise(across(-.draw, list(
          mean = ~ mean(.x),
          se = ~ sd(.x)
        ), .names = "{.col}_{.fn}"), .groups = "drop")

      shifts <- shifts |>
        group_by(!!rlang::sym(geo_var)) |>
        summarise(across(-.draw, mean))
    }
    res <- full_join(ps_table, res, by = ".cellid")
    if (!posterior_summary){
      res <- res |> relocate(.cellid, .draw)
    }

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

