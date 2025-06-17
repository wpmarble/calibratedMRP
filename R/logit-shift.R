
## Code for calibration functions.
# TODO: turn into an R package
library(tidybayes)


## The "logit shift" is used to calibrate turnout/vote share estimates.
## It takes individual-level turnout probabilities, then applies a constant
## additive shift on the logit scale such that the average transformed turnout
## matches the target. Formally, let x be a vector of turnout probs. The logit
## shift finds variable a such that
##    a.star = argmin_a | target - mean(inv.logit(logit(x) + a)) |
## The logit_shift() function returns the transformed vote choice probabilities,
## namely inv.logit( logit(x) + a.star ).



#' Shift probabilities on the logit scale
#' @keywords internal
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
#' @param var Variable in `ps_table` storing predictions to calibrate.
#' @param weight Weighting variable for population counts.
#' @param geography Grouping variable for calibration (e.g., county).
#' @param calib_target A data frame with true targets by geography.
#' @param calib_var Variable in `calib_target` storing true targets.
#'
#' @return A data frame with one row per geography and the estimated logit shift.
#' @keywords internal
logit_shift_single = function(ps_table,
                              var,
                              weight,
                              geography,
                              calib_target,
                              calib_var) {


  # standardizing names internally
  ps <- ps_table %>%
    rename(pred = {{ var }}, weight = {{ weight }}, geography = {{ geography }}) %>%
    select(pred, weight, geography)

  calib_target <- calib_target %>%
    rename(geography = {{ geography }}, calib_var = {{ calib_var }})


  # calculate logit shift for each geography
  shifts <- unique(ps$geography) %>%
    map( ~ {

      # subset PS table and target to this geography
      tmp_ps <- filter(ps, geography == .x)
      tmp_targ <- filter(calib_target, geography == .x) %>% pull(calib_var)

      # calculate logit shift
      if (is.na(tmp_targ) || length(tmp_targ) == 0){
        rlang::warn(c("Calibration target missing; returning logit shift = 0.",
                      "i" = sprintf("Geography: %s ", .x),
                      "i" = "This behavior is not optimal. See https://github.com/wpmarble/mrp/issues/3"))
        data.frame(geography = .x, shift = 0)
      } else {
        shift <- logit_shift_internal(x = tmp_ps$pred,
                                      w = tmp_ps$weight,
                                      target = tmp_targ)
        tibble(geography = .x, shift = shift)
      }
    })
  shifts <- bind_rows(shifts) %>%
    rename({{ geography }} := geography,
           "{{ var }}_shift" := shift)

  shifts
}


#' Compute logit shifts for multiple outcomes
#'
#'  Calibrate multiple sets of model-based predictions to match known
#'  geographic targets by computing logit intercept shifts for each outcome.
#'
#' @param ps_table A data frame of poststratification cells with predictions.
#' @param vars One or more outcome variables in `ps_table` to calibrate.
#'   Accepts tidyselect syntax (e.g., `voteshare_gov`, `c(voteshare_gov, voteshare_pres)`, `starts_with("voteshare")`).
#' @param weight Weighting variable in `ps_table` storing population counts.
#' @param geography Grouping variable found in `ps_table` and `calib_data` used
#' for calibration.
#' @param calib_target A data frame with true targets by `geography`.
#' @param calib_vars Variables in `calib_target` storing true targets.
#'   `vars` and `calib_vars` should be in the same order. Accepts tidyselect syntax.
#'
#' @return A data frame with logit shifts for each prediction variable and geography.
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr mutate select filter summarise rename left_join across any_of `%>%` bind_rows bind_cols
#'
#' @examples
#' # Example poststratification table
#' ps <- tibble::tibble(
#'   county = rep(c("A", "B"), each = 100),
#'   pred_vote = plogis(rnorm(200, qlogis(0.4), 0.5)),
#'   pred_turnout = plogis(rnorm(200, qlogis(0.3), 0.4)),
#'   weight = runif(200, 0.2, 0.8)
#' )
#'
#' # Calibration targets by county
#' calib <- tibble::tibble(
#'   county = c("A", "B"),
#'   vote_target = c(0.55, 0.45),
#'   turnout_target = c(0.65, 0.55)
#' )
#'
#' # Compute logit shifts for both outcomes
#' logit_shift(
#'   ps_table = ps,
#'   vars = c(pred_vote, pred_turnout),
#'   weight = weight,
#'   geography = county,
#'   calib_target = calib,
#'   calib_vars = c(vote_target, turnout_target)
#' )
#'
# TODO: add a final bit of example showing calibration worked right
logit_shift <- function(ps_table,
                        vars,
                        weight,
                        geography,
                        calib_target,
                        calib_vars){

  # Resolve tidyselect
  vars <- tidyselect::eval_select(rlang::enquo(vars), ps_table)
  var_names <- names(vars)

  calib_vars <- tidyselect::eval_select(rlang::enquo(calib_vars), calib_target)
  calib_names <- names(calib_vars)


  if (length(var_names) != length(calib_names)) {
    rlang::abort("Number of `vars` and `calib_vars` must match.")
  }

  # repeatedly call logit_shift_single() then combine results
  purrr::map2(var_names, calib_names, ~
                logit_shift_single(
                  ps_table = ps_table,
                  var = !!sym(.x),
                  weight = {{ weight }},
                  geography = {{ geography }},
                  calib_target = calib_target,
                  calib_var = !!sym(.y)
                )
  ) %>%
    purrr::reduce(full_join, by = rlang::englue("{{ geography }}"))

}


#' Recalibrate probabilities by applying a user-supplied logit shift
#'
#' Takes the poststratification table and a table of logit shifts and returns
#' the recalibrated probabilities.
#'
#' @param ps_table A data frame storing poststratification cells with model-based
#'   predictions attached.
#' @param shifts A data frame of logit shift parameters.
#' @param preds Variables in `ps_table` that store model-based predictions
#'   for two outomes and prefix for shift parameters in `shifts`
#' @param geography Geography variable linking `ps_table` and `shifts`.
#' @param shift_suffix Suffix for the shift variables, e.g. if the first outcome
#'   variable is called `presvote` then the shift variable will be called
#'   `presvote_shift`.
#' @param calib_suffix Suffix for recalibrated vote probabilities.


calibrate_preds <- function(ps_table, shifts, preds, geography,
                            shift_suffix = "shift",
                            calib_suffix = "calib"){

  # englue turns it into a character vector needed for _join
  geography <- rlang::englue("{{ geography }}")
  out <- left_join(ps_table, shifts, by = geography)

  # recalibrate probabilities, store in pred_calib variable
  out <- out %>%
    mutate(across(all_of(preds),
                  ~ logit_shift_intercept(x = .x,
                                          a = get(sprintf("%s_%s", cur_column(), shift_suffix))),
                  .names = "{.col}_{calib_suffix}"))

  # drop "shift" columns - just store updated probabilities
  todrop <- shifts %>%
    select(-{{ geography }}) %>%
    names() %>%
    intersect(names(out))
  out %>%
    select(-any_of(todrop))
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
#' @param mod A `brmsfit` object with random intercepts by `group`.
#' @param group Grouping variable to extract random intercept.
#' @param tidy Return a tidy data frame? If `FALSE`, returns an array of
#'   correlation matrices where the first dimension indexes draws from the
#'   posterior.
#' @param outcome_order An order for the outcome variables to return. Only used
#'   if tidy = FALSE.
#' @param draw_ids An integer vector specifying the posterior draws to be used.
#'   If NULL (the default), all draws are used.
#' @importFrom tidybayes spread_draws
get_re_covs <- function(mod, group,
                        tidy = FALSE,
                        outcome_order = NULL,
                        draw_ids = NULL){

  if (!"brmsfit" %in% class(mod)) {
    rlang::abort("mod must be a brmsfit object.")
  }

  if (is.null(draw_ids)) {
    draw_ids <- 1:ndraws(mod)
  }

  # get pattern to extract
  pattern <- sym(sprintf("sd_%s__.*_Intercept$|cor_%s__.*_Intercept_.*_Intercept$", group, group))

  out <- mod %>%
    spread_draws(!!pattern, regex = TRUE) %>%
    filter(.draw %in% draw_ids ) %>%
    select(c(starts_with("."), starts_with("sd_"), starts_with("cor_"))) %>%
    rename_with(.cols = everything(),
                .fn = ~ {
                  .x %>%
                    gsub(sprintf("%s|Intercept", group), "", .) %>%
                    gsub("__", "_", ., fixed = TRUE) %>%
                    gsub("__", "_", ., fixed = TRUE) %>%
                    gsub("\\_$", "", .)
                })

  out <- out %>%
    pivot_longer(-starts_with("."),
                 names_to = "param",
                 values_to = "est")

  out <- out %>%
    mutate(
      param_og = param,
      type = case_when(
        grepl("^cor", param) ~ "cor",
        grepl("^sd", param) ~ "sd"
      ))

  # 7/22/24 -- already got rid of `group` above.
  # todrop <- sprintf("sd_%s_|cor_%s_", group, group)
  todrop <- "sd_|cor_"
  out <- out %>%
    mutate(param = str_remove(param, todrop))


  out <- out %>%
    mutate(param1 = str_split(param, "_", simplify = TRUE)[,1],
           param2 = str_split(param, "_", simplify = TRUE)[,2]) %>%
    mutate(param2 = ifelse(param2 == "", param1, param2)) %>%
    select(c(starts_with("."), type, param1, param2, est)) %>%
    arrange(param1, param2)


  if (!tidy) {
    if (is.null(outcome_order)) {
      outcome_order <- unique(out$param1)
    } else if (!setequal(out$param1, outcome_order)) {
      rlang::warn("Levels of outcome_order do not match levels found in mod. Ignoring outcome_order.")
      outcome_order <- unique(out$param1)
    }

    # Reshape to array.
    out <- out %>%
      group_by(.draw) %>%
      group_split() %>%
      map(.progress = "Transforming covariance params to array",
          .f = ~ {
        tmp <- .x %>%
          select(param1, param2, est) %>%
          pivot_wider(values_from = est,
                      names_from = param2) %>%
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
logit_shift_aux <- function(shift,
                            shift_vars,
                            cov,
                            suffix = "_shift"){


  # The update is \Sigma_uo \Sigma_oo^{-1} shift_o
  # where Sigma_uo is the cross-covariance of uncalibrated and calibrated outcomes,
  # Sigma_oo is the covariance of calibrated outcomes, and shift_o is the
  cov_uo <- cov[!(rownames(cov) %in% shift_vars), shift_vars,drop = FALSE]
  cov_oo <- cov[shift_vars, shift_vars]
  cov_oo_inv <- MASS::ginv(cov_oo)
  premat <- cov_uo %*% cov_oo_inv

  shift_mat <- shift %>%
    select(all_of(paste0(shift_vars, suffix))) %>%
    as.matrix()

  updated_shift <- shift_mat %*% t(premat)
  colnames(updated_shift) <- paste0(colnames(updated_shift), suffix)
  shift <- bind_cols(shift, updated_shift, .name_repair = "unique")

  shift
}





#' Poststratify estimates
#'
#' Compute weighted means for one or more outcome variables using
#' poststratification weights, grouped by specified variables (typically geographies
#' or demographic groups). Returns a tidy data frame with one row per group
#' and one column per poststratified variable, plus the total population weight
#' for each group.
#'
#' @param ps_table A data frame where each row represents a poststratification cell.
#' @param vars One or more outcome variables to poststratify.
#'   Accepts tidyselect syntax (e.g., `var`, `c(var1, var2)`, `starts_with("prefix")`).
#' @param weight A variable giving the population weights in `ps_table`.
#' @param by One or more grouping variables. Accepts tidyselect syntax.
#' @param n_out Name of the output column that stores the total weight per group.
#' @param na.rm Should missing values in `vars` be ignored? Defaults to `TRUE`.
#'
#' @returns
#' A data frame with one row per group in `by` and columns for each variable in
#' `vars` containing the postratified estimates.
#'
#' @examples
#' # Example poststrat table
#' ps <- data.frame(
#'   county = rep(c("A", "B"), each = 3),
#'   demo_group = rep(c(1, 2), 3),
#'   cell_est = c(0.4, 0.5, 0.6, 0.3, 0.35, 0.4),
#'   N = c(100, 200, 300, 150, 250, 100)
#' )
#'
#' # Poststratify cell_est to counties
#' poststratify(
#'   ps_table = ps,
#'   vars = "cell_est",
#'   weight = N,
#'   by = county
#' )
poststratify = function(ps_table, vars, weight, by, n_out = "n", na.rm = TRUE) {

  ps_table %>%
    summarise(across({{ vars }},
                     ~ weighted.mean(.x, {{ weight }},
                                     na.rm = na.rm)),
              {{ n_out }} := sum({{ weight }}, na.rm = na.rm),
              .by = {{ by }})

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
#' @param preds A $D \times N \times K$ array or an $N \times K$ data frame
#'   storing predictions from the poststratification table to be normalized.
#'   If an array, the first index stores the draw number (e.g. from the posterior),
#'   the second indexes rows in the poststratification table, and the third indexes
#'   columns.
#' @param vars A list of sets of variables that must sum to one. See details.
#' @param constraints A numeric vector, of the same length as `vars`, that lists
#'   what the parameters must sum to. Alternatively, a list of the same length
#'   as `vars`, each element of which is a numeric vector that specifies what
#'   each row must sum to. Each element of this list must be the same length as
#'   the number of rows in `preds` (i.e. $N$). Defaults to 1.
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



