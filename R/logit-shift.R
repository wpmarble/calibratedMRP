
## Code for calibration functions.
# TODO: turn into an R package
# library(tidybayes)


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


#' Shift probabilities on the logit scale
#' @keywords internal
#' @importFrom stats plogis qlogis
logit_shift_intercept <- function(x, a) {
  if (any(x[!is.na(x)] < 0) || any(x[!is.na(x)] > 1)) {
    rlang::abort(sprintf("x must be in (0, 1); got: %s", paste0(head(x), collapse = ", ")))
  }
  plogis(qlogis(x) + a)
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


  f <- function(a) {
    (target - weighted.mean(logit_shift_intercept(x, a), w = w, na.rm = na.rm))^2
  }


  opt <- optimize(f, interval = c(-6, 6))

  if (abs(abs(opt$minimum) - 6) < tol) {
    warning(
      "Boundary solution found: intercept shift = ",
      opt$minimum,
      ". Do not trust these results."
    )
  }

  opt$minimum
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
      tmp_ps <- filter(ps, geography == g)
      tmp_targ <- filter(calib_target, geography == g) |> dplyr::pull(calib_var)

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
#'   vars = c(pred_vote, pred_turnout),
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
#'   vars = c(pred_vote, pred_turnout),
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



#' Extract covariance of random effects across outcomes
#'
#' `get_re_covs` returns the covariance parameters for random effects that are
#' correlated across outcomes. It takes a `brmsfit` object and extracts the parameters
#' that start with `sd_` and `cor_`, end with `_Intercept`, and include the pattern
#' given in the `group` argument. If `tidy = TRUE`, the output is a tidy data frame
#' in the format of `tidybayes::spread_draws`, except with columns given informative
#' names. If `tidy = FALSE`, an array of covariance matrices is returned (one
#' for each draw from the posterior).
#'
#' @param model A `brmsfit` object with random intercepts by `group`.
#' @param group Grouping variable to extract random intercept.
#' @param tidy Return a tidy data frame? If `FALSE`, returns an array of
#'   covariance matrices where the first dimension indexes draws from the
#'   posterior.
#' @param outcome_order An order for the outcome variables to return. Only used
#'   if tidy = FALSE.
#' @param draw_ids An integer vector specifying the posterior draws to be used.
#'   If NULL (the default), all draws are used.
#' @importFrom tidybayes spread_draws
#' @importFrom dplyr filter select rename_with mutate arrange group_by group_split bind_cols `%>%`
#' @export
get_re_covs <- function(model, group,
                        tidy = FALSE,
                        outcome_order = NULL,
                        draw_ids = NULL){

  if (!"brmsfit" %in% class(model)) {
    rlang::abort("model must be a brmsfit object.")
  }

  if (is.null(draw_ids)) {
    draw_ids <- 1:ndraws(model)
  }

  # get pattern to extract
  pattern <- sym(sprintf("sd_%s__.*_Intercept$|cor_%s__.*_Intercept_.*_Intercept$", group, group))

  out <- model |>
    spread_draws(!!pattern, regex = TRUE) |>
    filter(.draw %in% draw_ids ) |>
    select(c(starts_with("."), starts_with("sd_"), starts_with("cor_"))) |>
    rename_with(.cols = everything(),
                .fn = \(x) {
                  x %>%
                    gsub(sprintf("%s|Intercept", group), "", .) %>%
                    gsub("__", "_", ., fixed = TRUE) %>%
                    gsub("__", "_", ., fixed = TRUE) %>%
                    gsub("\\_$", "", .)
                })

  out <- out |>
    pivot_longer(-starts_with("."),
                 names_to = "param",
                 values_to = "est")

  out <- out |>
    mutate(
      param_og = param,
      type = case_when(
        grepl("^cor", param) ~ "cor",
        grepl("^sd", param) ~ "sd"
      ))

  todrop <- "sd_|cor_"
  out <- out %>%
    mutate(param = stringr::str_remove(param, todrop))


  out <- out %>%
    mutate(param1 = stringr::str_split(param, "_", simplify = TRUE)[,1],
           param2 = stringr::str_split(param, "_", simplify = TRUE)[,2]) |>
    mutate(param2 = ifelse(param2 == "", param1, param2)) |>
    select(c(starts_with("."), type, param1, param2, est)) |>
    arrange(param1, param2)


  if (!tidy) {
    if (is.null(outcome_order)) {
      outcome_order <- unique(out$param1)
    } else if (!setequal(out$param1, outcome_order)) {
      rlang::warn("Levels of outcome_order do not match levels found in mod. Ignoring outcome_order.")
      outcome_order <- unique(out$param1)
    }

    # Reshape to array.
    out <- out |>
      group_by(.draw) |>
      group_split() |>
      furrr::future_map(.progress = TRUE,
          .f = \(x) {
        tmp <- x |>
          select(param1, param2, est) |>
          pivot_wider(values_from = est,
                      names_from = param2) |>
          as.data.frame()
        rownames(tmp) <- tmp$param1
        tmp$param1 <- NULL
        tmp <- as.matrix(tmp)
        tmp <- tmp[outcome_order, outcome_order]
        tmp[is.na(tmp)] <- t(tmp)[is.na(tmp)]

        sds <- diag(tmp)
        diag(tmp) <- 1
        covout <- diag(sds) %*% tmp %*% diag(sds)
        dimnames(covout) <- dimnames(tmp)
        covout
      })
    out <- simplify2array(out)
    out <- aperm(out, c(3, 1,2))
    dimnames(out)[[1]] <- draw_ids
    names(dimnames(out)) <- c("draw", "par1", "par2")
  }
  out
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
    dplyr::select(all_of(paste0(shift_vars, suffix))) |>
    as.matrix()

  updated_shift <- shift_mat %*% t(premat)
  colnames(updated_shift) <- paste0(colnames(updated_shift), suffix)
  shift <- dplyr::bind_cols(shift, updated_shift, .name_repair = "unique")

  shift
}


#' Generate estimates for each poststratification cell
#'
#' Uses a fitted multivariate \code{brms} model to generate predictions for each row
#' in a poststratification table. Internally, the function uses \link[brms:posterior_epred]{brms::posterior_epred}
#' to generate draws from the posterior for each cell. There are two return options. If
#' \code{summarize = TRUE} (the default), the function summarizes the posterior draws
#' and returns the poststratification
#' table with added columns for posterior mean and SD of each outcome.
#' If \code{summarize = FALSE}, the function returns a 3D array that
#' contains each draw from the posterior for each cell and outcome. In this case,
#' the function is essentially a wrapper around \link[brms:posterior_epred]{brms::posterior_epred}.
#'
#' @param model A `brmsfit` object with multivariate binary outcomes.
#' @param ps_table A data frame of poststratification cells (must match model variables).
#' @param draw_ids Optional vector of posterior draw indices to use. If NULL, uses all draws.
#' @param outcome_suffix Optional suffix to append to posterior mean columns. Defaults to no suffix.
#' @param se_suffix Suffix to append to posterior SD columns. Defaults to "_se".
#' @param summarize Should the function return a summary data frame with posterior
#'    means and SDs? If `FALSE`, returns an array of draws from the posterior.
#' @param control A list of control parameters. See details for list of supported
#' controls.
#'
#' @details
#' This function can be memory and computationally intensive, especially when
#' working with large poststratification table and many posterior draws. To reduce
#' memory usage, the function processes the posterior samples in batches, which is
#' controlled by the `batch_size` item in the `control` list.
#'
#'
#' @return If \code{summarize = TRUE}, returns `ps_table` with additional columns
#' for outcome posterior means and posterior standard deviations.
#'  For example, if the outcomes in the `brms` model are `voteshare` and `turnout`
#'  and default suffix arguments are used, the addition columns will be called
#'  `voteshare`,`voteshare_se`, `turnout`, and `turnout_se`.
#'
#'   If \code{summarize = FALSE}, a 3D array with dimensions
#'   \code{N draws} \eqn{\times} \code{N cells} \eqn{\times} \code{N outcomes}.
#'
#' @export
#'
#' @importFrom progress progress_bar
#'
#' @examples
#' ## See vignettes.
#'
generate_cell_estimates <- function(model,
                                    ps_table,
                                    draw_ids = NULL,
                                    outcome_suffix = "",
                                    se_suffix = "_se",
                                    summarize = TRUE,
                                    control = list(batch_size = 50)) {

  if (!inherits(model, "brmsfit")) rlang::abort("`model` must be a brmsfit object.")

  # Get outcome names
  outcome_names <- formula(model)[2][[1]]

  # check that ps_table doesn't contain outcome columns (when summarizing)
  if (summarize && any(outcome_names %in% names(ps_table))) {
    name_conflicts <- outcome_names[outcome_names %in% names(ps_table)]
    rlang::abort(c("`ps_table` already contains outcome columns:",
                   "*" = paste0("Found: ", paste(name_conflicts, collapse = ", ")),
                   "i" = "Rename or remove these columns from `ps_table`"))
  } else {
    rlang::inform(c("Generating estimates for the following outcomes: ",
                    "*" = paste0(outcome_names, collapse = ", ")))
  }



  # summarize posterior
  if (summarize) {



    # process posterior samples in batches to reduce memory requirement
    batch_size <- control$batch_size
    all_draws <- if (is.null(draw_ids)) seq_len(brms::ndraws(model)) else draw_ids
    n_batches <- ceiling(length(all_draws) / batch_size)

    running_mean <- NULL
    running_M2   <- NULL
    n_total      <- 0

    # initialize progress bar
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent",
      show_after = 0,
      total = n_batches,
      clear = FALSE,
      width = 60
    )
    pb$tick(0)

    for (i in seq_len(n_batches)) {
      batch_ids <- all_draws[(((i - 1) * batch_size) + 1):min(i * batch_size, length(all_draws))]

      # Get draws: (draw x cell x outcome)
      pred_array <- brms::posterior_epred(
        model,
        newdata = ps_table,
        allow_new_levels = TRUE,
        draw_ids = batch_ids
      )

      # Summarize batch
      batch_mean <- apply(pred_array, c(2, 3), mean)
      batch_var  <- apply(pred_array, c(2, 3), var)

      n_batch    <- length(batch_ids)
      M2_batch   <- batch_var * (n_batch - 1)

      if (is.null(running_mean)) {
        # First batch
        running_mean <- batch_mean
        running_M2   <- M2_batch
        n_total      <- n_batch
      } else {
        delta <- batch_mean - running_mean
        new_n <- n_total + n_batch

        running_mean <- (n_total * running_mean + n_batch * batch_mean) / new_n
        running_M2   <- running_M2 + M2_batch + (delta^2) * n_total * n_batch / new_n
        n_total <- new_n
      }


      pb$tick()  # update progress bar
    }

    # Final posterior SD
    cell_sd <- sqrt(running_M2 / (n_total - 1))
    cell_mean <- running_mean



    # Convert to tibble with cell index
    pred_df <- tibble::as_tibble(cell_mean) %>%
      dplyr::rename_with(\(x) paste0(x, outcome_suffix)) %>%
      dplyr::mutate(.cell_id = dplyr::row_number())

    se_df <- tibble::as_tibble(cell_sd) %>%
      dplyr::rename_with(\(x) paste0(x, se_suffix)) %>%
      dplyr::mutate(.cell_id = dplyr::row_number())

    # Join and bind back to ps_table
    preds <- dplyr::left_join(pred_df, se_df, by = ".cell_id") %>%
      dplyr::select(-.cell_id)


    # reorder columns
    ordered_names <- as.vector(rbind(
      paste0(outcome_names, outcome_suffix),
      paste0(outcome_names, se_suffix)
    ))
    preds <- preds[, ordered_names]
    dplyr::bind_cols(ps_table, preds)
  } else {

    # Get draws: (draw x cell x outcome)
    pred_array <- brms::posterior_epred(
      model,
      newdata = ps_table,
      allow_new_levels = TRUE,
      draw_ids = draw_ids
    )

    dimnames(pred_array)[[1]] <- draw_ids
    dimnames(pred_array)[[2]] <- seq_len(dim(pred_array)[2])
    names(dimnames(pred_array)) <- c("draw", ".cell_id", "outcome")
    pred_array
  }
}


#' Poststratify estimates
#'
#' Compute weighted means for one or more outcome variables using
#' poststratification weights, grouped by specified variables (typically geographies
#' or demographic groups). Returns a tidy data frame with one row per group
#' and one column per poststratified variable, plus the total population weight
#' for each group. Optionally estimates uncertainty if cell-level standard errors
#' are provided.
#'
#' @param ps_table A data frame where each row represents a poststratification cell.
#'   The table must include the variables specified in `outcomes`, `weight`, `by`,
#'   and (optionally) `ses`.
#' @param outcomes One or more outcome variables to poststratify.
#'   Accepts tidyselect syntax (e.g., `var`, `c(var1, var2)`, `starts_with("prefix")`).
#'   These variables should include cell-level summaries of the outcomes of interest.
#' @param ses Should standard errors be computed for the postratified estimates under
#'   the assumption of independence across cells? See details.
#' @param se_suffix Character suffix for the columns containing standard errors
#'   for the cell-level estimates of `outcomes`. Defaults to "_se", so if the
#'   outcome is `voteshare`, the cell-level standard errors will be stored in
#'   a column called `voteshare_se`.
#' @param weight A variable giving the population weights in `ps_table`.
#' @param by One or more grouping variables. Accepts tidyselect syntax.
#' @param n_out Name of the output column that stores the total weight per group.
#' @param na.rm Should missing values in `vars` be ignored? Defaults to `TRUE`.
#'
#' @returns
#' A data frame with one row per group in `by` and columns for each variable in
#' `vars` containing the postratified estimates.
#'
#' @details
#' When the `outcomes` in `ps_table` are posterior means from a multilevel model,
#' the resulting postratified estimates are also posterior means. The uncertainty
#' estimation assumes that there is 0 covariance between cell-level posterior means
#' and computes the standard errors as \eqn{\sqrt{\sum_c w_c^2 \sigma_c^2 / (\sum w_c)^2}},
#' where \eqn{w_c} is the weight for cell \eqn{c} and \eqn{\sigma_c^2} is the
#' standard error of the posterior mean in cell \eqn{c}.
#' This feature is useful when first summarize the posterior of the multilevel
#' model at the cell level, then use this summary for downstream analysis. This
#' two-step procedure is more efficient in terms of storage and computation than
#' full Bayesian inference, which involves computing the final quantity of interest
#' separately for each draw from the posterior. However, the assumption of
#' independence across cells is typically unrealistic, and so full Bayesian inference
#' should typically yield more accurate uncertainty estimates.
#'
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#'
#' @examples
#' # Example poststrat table
#' ps <- data.frame(
#'   county = rep(c("A", "B"), each = 3),
#'   demo_group = rep(c(1, 2), 3),
#'   estimate = c(0.4, 0.5, 0.6, 0.3, 0.35, 0.4),
#'   estimate_se = c(0.05, 0.04, 0.06, 0.03, 0.02, 0.05),
#'   N = c(100, 200, 300, 150, 250, 100)
#' )
#'
#' # Poststratify cell_est to counties
#' poststratify(
#'   ps_table = ps,
#'   outcomes = estimate,
#'   ses = TRUE,
#'   weight = N,
#'   by = county
#' )
#'

poststratify <- function(ps_table, outcomes, ses = FALSE, se_suffix = "_se",
                         weight, by, n_out = "n", na.rm = TRUE) {


  group_vars <- rlang::enquo(by)
  weight_var <- rlang::enquo(weight)
  mean_vars <- tidyselect::eval_select(rlang::enquo(outcomes), ps_table) |> names()
  se_vars <- paste0(mean_vars, se_suffix)


  out <- ps_table %>%
    dplyr::summarise(dplyr::across({{ outcomes }},
                                   \(x) weighted.mean(x, {{ weight }}, na.rm = na.rm)),
              {{ n_out }} := sum({{ weight }}, na.rm = na.rm),
              .by = {{ by }})

  if (ses) {
    if (!.calibratedMRP_env$postratify_se_warned) {
      rlang::inform(c("Computing standard errors for postratified estimates under the assumption of independence across cells.",
                     "x" = "True uncertainty may be larger or smaller depending on covariance of posterior means of `outcomes` across cells",
                     "i" = "For full Bayesian inference, do XYZ",
                     "i" = "This message will only be displayed once per session"))
      .calibratedMRP_env$warned_se <- TRUE
    }

    out_se <- ps_table %>%
      dplyr::summarise(dplyr::across(all_of(se_vars),
                                     \(x) sqrt(sum({{ weight }}^2 * x^2, na.rm = na.rm)
                                               / sum({{ weight }}, na.rm = na.rm)^2 )),
                       .by = {{ by }})
    out <- left_join(out, out_se, by = rlang::englue("{{ by }}"))
  }

  # Reorder columns: [grouping vars], outcome1, outcome1_se, outcome2, outcome2_se, ..., n_out
  group_names <- names(out[tidyselect::eval_select(rlang::enquo(by), out)])


  mean_cols <- mean_vars
  se_cols <- character(0)
  if (ses && !is.null(se_suffix)) {
    se_cols <- paste0(mean_vars, se_suffix)
  }
  remainder <- setdiff(names(out), c(group_names, mean_cols, se_cols))
  col_order <- c(group_names, as.vector(rbind(mean_cols, se_cols)), remainder)
  out <- dplyr::select(out, dplyr::all_of(col_order))
  out

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



