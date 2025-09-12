
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
#' @param keep_orig keep original uncalibrated columns?
#'
#' @return Returns `ps_table` with additional columns for calibrated predictions.
#'
#' @export
#'
#' @examples
#' ## See examples in ?logit_shift()
calibrate_preds <- function(ps_table,
                            shifts,
                            preds,
                            geography,
                            shift_suffix = "shift",
                            keep_orig = TRUE) {

  # Resolve tidyselect
  pred_vars <- names(tidyselect::eval_select(rlang::enquo(preds), ps_table))
  geo_var  <- rlang::as_name(rlang::ensym(geography))
  out <- dplyr::left_join(ps_table, shifts, by = geo_var)

  # Recalibrate each prediction
  for (pred in pred_vars) {
    shift_col <- paste0(pred, "_", shift_suffix)
    calib_col <- paste0(pred, "_calib")

    out[[calib_col]] <- logit_shift_intercept(
      x = out[[pred]],
      a = out[[shift_col]]
    )
  }

  # Drop original shift columns
  shift_cols <- paste0(pred_vars, "_", shift_suffix)
  out <- dplyr::select(out, -dplyr::any_of(shift_cols))

  if (!keep_orig){
    out <- dplyr::select(out, -dplyr::any_of(pred_vars))
  }
  out
}

