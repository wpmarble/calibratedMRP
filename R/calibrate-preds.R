
#' Recalibrate probabilities after computing logit shifts
#'
#' Takes the poststratification table and a table of logit shift parameters and
#' returns the calibrated cell-level probabilities.
#'
#' @param ps_table A data frame storing poststratification cells with (uncalibrated)
#' model-based predictions.
#' @param shifts A data frame of logit shift parameters computed by [logit_shift()].
#' Each shift must be named `{pred}_{shift_suffix}`, e.g., `voteshare_shift` if
#' the prediction variable in `ps_table` is `voteshare`.
#' @param preds One or more prediction variables in `ps_table` to recalibrate.
#'   Accepts tidyselect syntax (e.g., `var`, `c(var1, var2)`, `starts_with()`).
#' @param geography Group variable variable linking `ps_table` and `shifts`.
#' This should be a single variable name (unquoted) that exists in both data frames.
#' @param shift_suffix String suffix for the shift variables, e.g. if the first outcome
#'   variable is called `presvote` then the shift variable will be called
#'   `presvote_shift`.
#' @param calib_suffix String suffix for calibrated vote probabilities. E.g., if
#' original prediction is valled `voteshare` then the calibrated prediction
#' will be called `voteshare_calib`.
#'
#' @return Returns `ps_table` with additional columns for calibrated predictions.
#'
#' @export
#'
#' @examples
#' ## Example poststratification table with predictions for voteshare and turnout
#' ps <- tibble::tibble(county = rep(c("A", "B"), each = 100),
#'                      demo_group = rep(1:100, 2),
#'                      pred_vote = plogis(rnorm(200, qlogis(0.4), 0.5)),
#'                      pred_turnout = plogis(rnorm(200, qlogis(0.3), 0.4)),
#'                      weight = runif(200, 0.2, 0.8)
#'                      )
#'
#' ## Calibration targets by county
#' targets <- tibble::tibble(
#'   county = c("A", "B"),
#'   vote_target = c(0.55, 0.45),
#'   turnout_target = c(0.65, 0.55)
#' )
#'
#' ## Compute logit shifts for both outcomes
#' shifts <- logit_shift(
#'   ps_table = ps,
#'   outcomes = c(pred_vote, pred_turnout),
#'   weight = weight,
#'   geography = county,
#'   calib_target = targets,
#'   calib_vars = c(vote_target, turnout_target)
#' )
#'
#' ## Calibrate predictions in poststratification table
#' ps <- calibrate_preds(ps_table = ps, shifts = shifts,
#'                       preds = c(pred_vote, pred_turnout),
#'                       geography = county)
#' head(ps)
calibrate_preds <- function(ps_table,
                            shifts,
                            preds,
                            geography,
                            shift_suffix = "shift",
                            calib_suffix = "calib") {

  # Resolve tidyselect
  pred_vars <- tidyselect::eval_select(rlang::enquo(preds), ps_table) |> names()
  geo_var  <- rlang::as_name(rlang::ensym(geography))
  out <- dplyr::left_join(ps_table, shifts, by = geo_var)

  # Recalibrate each prediction column manually
  for (pred in pred_vars) {
    shift_col <- paste0(pred, "_", shift_suffix)
    calib_col <- paste0(pred, "_", calib_suffix)

    out[[calib_col]] <- logit_shift_intercept(
      x = out[[pred]],
      a = out[[shift_col]]
    )
  }

  # Drop original shift columns
  shift_cols <- paste0(pred_vars, "_", shift_suffix)
  out <- dplyr::select(out, -dplyr::any_of(shift_cols))
  out
}

