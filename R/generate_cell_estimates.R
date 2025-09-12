
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
#' @param outcomes Optional character vector of outcome variables to calibrate.
#'  If `NULL`, the outcome names are inferred from the `model` formula.
#' @param summarize Should the function return a summary data frame with posterior
#'    means and SDs? If `FALSE`, returns an array of draws from the posterior.
#' @param draw_ids Optional vector of posterior draw indices to use. If NULL
#'    (the default), use all draws.
#' @param outcome_suffix Optional suffix to append to posterior mean columns. Defaults to no suffix.
#' @param se_suffix Suffix to append to posterior SD columns. Defaults to "_se".
#' @param control A list of control parameters. See details for list of supported
#' controls.
#'
#' @details
#' This function generates estimates of the average outcome in each cell of the
#' poststratification table.
#'
#' This function is memory and computationally intensive, especially when
#' working with large poststratification tables. To reduce computational burden when
#' developing code, you can reduce the number of draws via `draw_ids` or work with
#' a subset of the poststratification table. However, all draws should be used for
#' final inference.
#'
#' When `summarize = TRUE`, it is possible to reduce memory usage by processing
#' the posterior samples in batches, which is controlled by the `batch_size`
#' item in the `control` list. If you have a large amount of memory available,
#' increasing the batch size may speed up computation.
#'
#' When `summarize = FALSE`, the function returns a array of posterior draws,
#' which may have large dimensions in some cases.
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
#' ## See vignette for usage and workflow.
#'
generate_cell_estimates <- function(model,
                                    ps_table,
                                    outcomes = NULL,
                                    summarize = TRUE,
                                    draw_ids = NULL,
                                    outcome_suffix = "",
                                    se_suffix = "_se",
                                    control = list(batch_size = 50)) {

  if (!inherits(model, "brmsfit")) rlang::abort("`model` must be a brmsfit object.")


  outcomes_quo <- rlang::enquo(outcomes)
  if (is.null(outcomes)) {
    outcomes <- get_outcomes(model$formula)
    rlang::inform("No `outcomes` provided, inferring outcome variables from model formula")
  }
  n_outcomes <- length(outcomes)

  if (is.null(draw_ids)) {
    draw_ids <- seq_len(brms::ndraws(model))
  }

  # check that ps_table doesn't contain outcome columns (when summarizing)
  if (summarize && any(outcomes %in% names(ps_table))) {
    name_conflicts <- outcomes[outcomes %in% names(ps_table)]
    rlang::abort(c("`ps_table` already contains outcome columns:",
                   "*" = paste0("Found: ", paste(name_conflicts, collapse = ", ")),
                   "i" = "Rename or remove these columns from `ps_table`"))
  } else {
    rlang::inform(c("Generating estimates for the following outcomes: ",
                    "*" = paste0(outcomes, collapse = ", ")))
  }



  # summarize posterior
  if (summarize) {



    # process posterior samples in batches to reduce memory requirement
    batch_size <- control$batch_size
    all_draws  <- draw_ids
    n_batches  <- ceiling(length(all_draws) / batch_size)

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
        resp = outcomes,
        newdata = ps_table,
        allow_new_levels = TRUE,
        draw_ids = batch_ids
      )

      # Summarize batch
      if (n_outcomes == 1) {
        dim(pred_array) <- c(dim(pred_array), 1)
        dimnames(pred_array)[3] <- list(outcomes)
      }
      batch_mean <- apply(pred_array, c(2, 3), mean)
      batch_var  <- apply(pred_array, c(2, 3), stats::var)

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
      paste0(outcomes, outcome_suffix),
      paste0(outcomes, se_suffix)
    ))
    preds <- preds[, ordered_names]
    dplyr::bind_cols(ps_table, preds)
  } else {

    # Get draws: (draw x cell x outcome)
    pred_array <- brms::posterior_epred(
      model,
      newdata = ps_table,
      resp = outcomes,
      allow_new_levels = TRUE,
      draw_ids = draw_ids
    )
    if (n_outcomes == 1) {
      dim(pred_array) <- c(dim(pred_array), 1)
    }
    dimnames(pred_array)[[1]] <- draw_ids
    dimnames(pred_array)[[2]] <- seq_len(dim(pred_array)[2])
    dimnames(pred_array)[[3]] <- outcomes
    names(dimnames(pred_array)) <- c("draw", ".cell_id", "outcome")
    pred_array
  }
}





#' Infer outcomes from a model formula
#'
#' @param formula a formula
#' @return character vector of outcome variable(s)
get_outcomes <- function(formula){
  if (inherits(formula, "bform")) {
    if (inherits(formula, "mvbrmsformula")) {
      outcomes <- formula$responses
    } else {
      outcomes <- formula$resp
    }
  } else {
    outcomes <- formula.tools::lhs(formula)
  }
}
