
## Code for calibration functions.


## The "logit shift" is used to calibrate turnout/vote share estimates.
## It takes individual-level turnout probabilities, then applies a constant
## additive shift on the logit scale such that the average transformed turnout
## matches the target. Formally, let x be a vector of turnout probs. The logit
## shift finds variable a such that
##    a.star = argmin_a | target - mean(inv.logit(logit(x) + a)) |
## The logit_shift() function returns the transformed vote choice probabilities,
## namely inv.logit( logit(x) + a.star ).


# Package-level environment variables
.calibratedMRP_env <- new.env(parent = emptyenv())
.calibratedMRP_env$postratify_se_warned <- FALSE




#' Compute logit shifts for multiple outcomes given target marginals
#'
#'  Calibrate multiple sets of model-based predictions to match known
#'  geographic targets by computing logit intercept shifts for each outcome.
#'
#' @param ps_table A data frame of poststratification cells with predictions.
#' @param outcomes One or more outcome variables in `ps_table` to calibrate.
#'   Accepts tidyselect syntax (e.g., `voteshare_gov`, `c(voteshare_gov, voteshare_pres)`, `starts_with("voteshare")`).
#' @param weight Weighting variable in `ps_table` storing population counts.
#' @param geography Grouping variable found in `ps_table` and `calib_data` used
#' for calibration.
#' @param calib_target A data frame with true targets by `geography`.
#' @param calib_vars Variables in `calib_target` storing true targets.
#'   `outcomes` and `calib_vars` must be in the same order. Accepts tidyselect syntax.
#'
#' @return A data frame with logit shifts for each prediction variable and geography.
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr mutate select filter summarise rename left_join across any_of `%>%` bind_rows bind_cols
#' @importFrom rlang sym enquo enquos `!!`
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
#' shifts # data frame with intercept shift needed for each geography
#'
#' ## Calibrate predictions for each cell in poststratification table
#' ps <- calibrate_preds(ps_table = ps, shifts = shifts,
#'                       preds = c(pred_vote, pred_turnout),
#'                       geography = county)
#' head(ps)
logit_shift <- function(ps_table,
                        outcomes,
                        weight,
                        geography,
                        calib_target,
                        calib_vars){

  # Resolve tidyselect
  weight_var <- rlang::as_name(rlang::enquo(weight))
  geo_var <- rlang::as_name(rlang::enquo(geography))

  outcomes <- tidyselect::eval_select(rlang::enquo(outcomes), ps_table)
  calib_vars <- tidyselect::eval_select(rlang::enquo(calib_vars), calib_target)

  var_names <- names(outcomes)
  calib_names <- names(calib_vars)


  if (length(var_names) != length(calib_names)) {
    rlang::abort("Number of `outcomes` and `calib_vars` must match.")
  }

  # repeatedly call logit_shift_single() then combine results
  purrr::map2(var_names, calib_names,
              \(x, y)  logit_shift_single(ps_table = ps_table,
                                          outcome = x,
                                          weight = weight_var,
                                          geography = geo_var,
                                          calib_target = calib_target,
                                          calib_var = y),
              .progress = TRUE) |>
    purrr::reduce(dplyr::full_join, by = geo_var)

}



#' Find the logit-scale intercept shift such that weighted.mean(x,w) = target
#' @keywords internal
logit_shift_internal <- function(x, w = rep(1, length(x)), target, na.rm = TRUE, tol = 1e-5){
  if (target %in% c(0, 1)) {
    rlang::warn(sprintf(
      "Target is %f but logit shift requires target strictly between 0 and 1. Trimming to 0.1%% or 99.0%%.",
      target
    ))
    target <- ifelse(target == 0, 0.001, 0.999)
  }
  if (target < 0 || target > 1) {
    rlang::abort(sprintf("Target must be in [0,1]. Did you pass a percent instead of a proportion?"))
  }


  f <- function(a) (target - stats::weighted.mean(logit_shift_intercept(x, a), w = w, na.rm = na.rm))^2


  opt <- stats::optimize(f, interval = c(-6, 6))

  if (abs(abs(opt$minimum) - 6) < tol) {
    warning(
      "Boundary solution found: intercept shift = ",
      opt$minimum,
      ". Do not trust these results."
    )
  }

  opt$minimum
}



#' Shift probabilities on the logit scale
#' @keywords internal
#' @importFrom stats plogis qlogis
logit_shift_intercept <- function(x, a) {
  if (any(x[!is.na(x)] < 0) || any(x[!is.na(x)] > 1)) {
    rlang::abort(sprintf("x must be in (0, 1); got: %s", paste0(utils::head(x), collapse = ", ")))
  }
  plogis(qlogis(x) + a)
}




#' Compute logit shift for a single outcome
#'
#' Internal helper to calibrate model-based probabilities to match known
#' geographic targets by applying a constant logit-scale intercept shift.
#' Called by [logit_shift()] for each outcome.
#'
#' @param ps_table A data frame of poststratification cells with predictions.
#' @param outcome CharVariable in `ps_table` storing outcome to calibrate.
#' @param weight Weighting variable for population counts.
#' @param geography Grouping variable for calibration (e.g., county).
#' @param calib_target A data frame with true targets by geography.
#' @param calib_var Variable in `calib_target` storing true targets.
#'
#' @return A data frame with one row per geography and the estimated logit shift.
#' @importFrom tidyselect all_of
#' @importFrom dplyr select filter pull bind_rows
#' @importFrom rlang `:=`
#' @keywords internal
logit_shift_single = function(ps_table,
                              outcome,
                              weight,
                              geography,
                              calib_target,
                              calib_var) {

  # var        <- rlang::as_name(rlang::ensym(outcome))
  # print(sprintf("`outcome` argument evaluates to: %s", var))
  # weight     <- rlang::as_name(rlang::ensym(weight))
  # geography  <- rlang::as_name(rlang::ensym(geography))
  # calib_var  <- rlang::as_name(rlang::ensym(calib_var))


  ps <- ps_table |>
    dplyr::select(pred = all_of(outcome),
                  weight = all_of(weight),
                  geography = all_of(geography))

  calib_target <- calib_target |>
    dplyr::select(geography = all_of(geography),
                  calib_var = all_of(calib_var))


  # calculate logit shift for each geography
  shifts <- unique(ps$geography) %>%
    furrr::future_map(\(g) {

      # subset PS table and target to this geography
      tmp_ps <- dplyr::filter(ps, geography == g)
      tmp_targ <- dplyr::filter(calib_target, geography == g) |> dplyr::pull(calib_var)

      # calculate logit shift
      if (is.na(tmp_targ) || length(tmp_targ) == 0){
        rlang::warn(c("Calibration target missing; returning logit shift = 0.",
                      "i" = sprintf("Geography: %s ", g),
                      "i" = "This behavior is not optimal. See https://github.com/wpmarble/mrp/issues/3"))
        shift <- 0
      } else {
        shift <- logit_shift_internal(x = tmp_ps$pred,
                                      w = tmp_ps$weight,
                                      target = tmp_targ)
      }
      tibble::tibble(!!geography := g, !!paste0(outcome, "_shift") := shift)
    })
  bind_rows(shifts)
}



#' Estimate logit shift for auxiliary variables given calibrated outcomes and covariance
#'
#' @param shift A data frame storing original logit shifts for each level of geography.
#' @param shift_vars Variables in `shift` data frame that stores the logit shifts.
#' @param cov A covariance matrix for the logit shift across outcomes. Row
#'   and column names should contain the variables in `shift_vars`.
#' @param suffix The suffix for the logit shift parameters in `shift`. The
#'   prefixes are given by `shift_vars`. For example if the outcome is `dem` then
#'   the logit shift parameter might be stored in a variable called `dem_shift`.
#' @export
logit_shift_aux <- function(shift,
                            shift_vars,
                            cov,
                            suffix = "_shift"){


  # The update is \Sigma_uo \Sigma_oo^{-1} shift_o
  # where Sigma_uo is the cross-covariance of uncalibrated and calibrated outcomes,
  # Sigma_oo is the covariance of calibrated outcomes, and shift_o is the
  cov_uo <- cov[!(rownames(cov) %in% shift_vars), shift_vars, drop = FALSE]
  cov_oo <- cov[shift_vars, shift_vars]
  cov_oo_inv <- MASS::ginv(cov_oo)
  premat <- cov_uo %*% cov_oo_inv

  shift_mat <- shift |>
    dplyr::select(tidyselect::all_of(paste0(shift_vars, suffix))) |>
    as.matrix()

  updated_shift <- shift_mat %*% t(premat)
  colnames(updated_shift) <- paste0(colnames(updated_shift), suffix)
  shift <- dplyr::bind_cols(shift, updated_shift, .name_repair = "unique")

  shift
}



#' Normalize predictions to sum to a constraint within rows
#'
#' Rescales predictions of categorical variables to sum to a particular
#' constraint. This is useful when the predictions for each category are modeled
#' as independent --- for example, when three vote choice categories (Dem, Rep,
#' Other) are modeled in separate equations --- but they must sum to a specific
#' value.
#'
#'
#' @param preds A \eqn{D \times N \times K} array or an \eqn{N \times K} data frame
#'   storing predictions from the poststratification table to be normalized.
#'   If an array, the first index stores the draw number (e.g. from the posterior),
#'   the second indexes rows in the poststratification table, and the third indexes
#'   columns.
#' @param vars A list of sets of variables that must sum to one. See details.
#' @param constraints A numeric vector, of the same length as `vars`, that lists
#'   what the parameters must sum to. Alternatively, a list of the same length
#'   as `vars`, each element of which is a numeric vector that specifies what
#'   each row must sum to. Each element of this list must be the same length as
#'   the number of rows in `preds` (i.e. \eqn{N}). Defaults to 1.
#' @details For each set of variables in `vars`, this function rescales the
#'   variables so that their sum equals their corresponding constraint. For
#'   example, if `vars = list(c(v1, v2), c(v3, v4))` and `constraints = c(1, 2)`,
#'   this function will enforce the following constraints for each row in `preds`:
#'   `v1 + v2 = 1` and `v3 + v4 = 2`.

normalize_preds <- function(preds, vars, constraints = rep(1, length(vars))) {

  force(constraints)

  ## error checking
  if (!(length(vars) == length(constraints))) {
    rlang::abort("vars and constraints must be same length")
  }

  if (!is.numeric(constraints)) {
    rlang::abort("constraints must be numeric")
  }


  ## get class of preds, behavior differs based on data.frame summary or array
  ## of draws.
  arr <- "array" %in% class(preds)
  df <- "data.frame" %in% class(preds)

  if (!(arr || df)) {
    rlang::abort("preds must be an array or data frame")
  }

  # check data format and that all vars are in preds
  if (arr) {
    if (!all(unlist(vars) %in% dimnames(preds)[[3]])) {
      diffs <- setdiff(unlist(vars), dimnames(preds)[[3]])
      rlang::abort(sprintf("the following variables are not in preds: %s",
                           paste0(diffs, collapse = ", ")))
    }

    if (length(dim(preds)) != 3) {
      rlang::abort("if preds is an array, it must have three dimensions: [draws, rows, columns]")
    }
  }
  if (df) {
    if (!all(unlist(vars) %in% names(preds))) {
      diffs <- setdiff(unlist(vars), names(preds))
      rlang::abort(sprintf("the following variables are not in preds: %s",
                           paste0(diffs, collapse = ", ")))

    }
  }


  # rescale array of draws
  if (arr) {
    out <- apply(preds, 1, function(x) {
      for (i in 1:length(vars)) {
        x[,vars[[i]]] <- constraints[[i]] * (x[,vars[[i]]] / rowSums(x[,vars[[i]]]))
      }
      return(x)
    },
    simplify = FALSE)
    out <- simplify2array(out)
    out <- aperm(out, c(3, 1,2))

  }

  # rescale data frame
  if (df) {
    for (i in 1:length(vars)) {
      preds[, vars[[i]]] <- constraints[[i]] * (preds[, vars[[i]]] / rowSums(preds[, vars[[i]]]))
    }
    out <- preds
  }

  return(out)

}



