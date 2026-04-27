

#' Poststratify estimates
#'
#' Compute weighted means for one or more outcome variables using
#' poststratification weights, grouped by specified variables. Accepts a
#' plain data frame (default method), a \code{\link{calibrate_mrp}} result
#' summarized as posterior means and SEs (\code{calibrated_summary}), or a
#' full posterior draws result (\code{calibrated_draws}).
#'
#' @param ps_table A data frame (default method), or a \code{calibrated_summary}
#'   or \code{calibrated_draws} list returned by \code{\link{calibrate_mrp}}.
#' @param outcomes (Default and calibrated_summary methods) One or more outcome
#'   variables to poststratify. For the default method, accepts tidyselect
#'   syntax (e.g., \code{var}, \code{c(var1, var2)}, \code{starts_with("prefix")}).
#'   For \code{calibrated_summary}, a character vector of base outcome names
#'   (e.g., \code{c("voteshare", "turnout")}), or \code{NULL} (default) to
#'   auto-detect from columns ending in \code{_mean} that have paired \code{_se}
#'   columns. For \code{calibrated_draws}, a required character vector of
#'   outcome column names in \code{x$results} (no auto-detection).
#' @param ses (Default method only) Logical. Should standard errors be computed
#'   for the poststratified estimates under the assumption of independence across
#'   cells? See Details. Not valid for classed methods; use
#'   \code{posterior_summary = TRUE} for \code{calibrated_draws}.
#' @param se_suffix (Default method only) Character suffix for the SE columns.
#'   Defaults to \code{"_se"}.
#' @param posterior_summary (\code{calibrated_draws} method only) Logical.
#'   If \code{TRUE} (default), collapse across posterior draws and return
#'   \code{<outcome>_mean} and \code{<outcome>_sd} columns. If \code{FALSE},
#'   return one row per draw x by-group combination with bare \code{<outcome>}
#'   columns.
#' @param weight Population weight column (bare name, NSE).
#' @param by Grouping variable(s) (bare name(s), tidyselect syntax).
#' @param n_out Name of the column storing total weight per group. Default
#'   \code{"n"}.
#' @param na.rm Logical. Ignore NAs? Default \code{TRUE}.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' **Default method**: Accepts any data frame where each row is a
#' poststratification cell. Poststratifies using \code{weighted.mean()}.
#' When \code{ses = TRUE}, propagates uncertainty under the
#' independence-across-cells assumption: the SE for each by-group is
#' \eqn{\sqrt{\sum_c w_c^2 \sigma_c^2 / (\sum_c w_c)^2}}, where \eqn{w_c}
#' is the cell weight and \eqn{\sigma_c^2} is the cell-level variance. This
#' two-step procedure is more efficient in terms of storage and computation than
#' full Bayesian inference, but the independence assumption is typically
#' unrealistic.
#'
#' **\code{calibrated_summary} method**: Accepts a \code{calibrated_summary}
#' list from \code{\link{calibrate_mrp}} (produced with \code{method = "plugin"}
#' or \code{method = "bayes"} with \code{posterior_summary = TRUE}). Extracts
#' \code{x$results}, identifies outcome columns by the \code{_mean} / \code{_se}
#' suffix pairs, and returns a tibble with \code{<outcome>_mean} and
#' \code{<outcome>_se} columns per by-group.
#'
#' **\code{calibrated_draws} method**: Accepts a \code{calibrated_draws} list
#' from \code{\link{calibrate_mrp}} (produced with
#' \code{method = "bayes", posterior_summary = FALSE}). Poststratifies each
#' posterior draw separately. When \code{posterior_summary = TRUE} (default),
#' returns \code{<outcome>_mean} and \code{<outcome>_sd} columns per by-group.
#' Note that \code{_sd} is the posterior SD across draws, which is conceptually
#' distinct from the cell-level-independence-formula SE returned by the default
#' method's \code{ses = TRUE} path. When \code{posterior_summary = FALSE},
#' returns one row per draw x by-group with bare \code{<outcome>} columns and
#' a \code{.draw} column.
#'
#' @returns
#' A tibble. Shape depends on method and arguments:
#' \itemize{
#'   \item Default: one row per by-group; columns are \code{by}, outcomes,
#'     (SE columns if \code{ses = TRUE}), and \code{n_out}.
#'   \item \code{calibrated_summary}: one row per by-group; columns are
#'     \code{by}, interleaved \code{<outcome>_mean} / \code{<outcome>_se}
#'     pairs, and \code{n_out}.
#'   \item \code{calibrated_draws} with \code{posterior_summary = TRUE}: one
#'     row per by-group; columns are \code{by}, interleaved
#'     \code{<outcome>_mean} / \code{<outcome>_sd} pairs, and \code{n_out}.
#'   \item \code{calibrated_draws} with \code{posterior_summary = FALSE}: one
#'     row per draw x by-group; columns are \code{.draw}, \code{by},
#'     outcomes, and \code{n_out}.
#' }
#'
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#'
#' @examples
#' # --- Default method ---
#' ps <- data.frame(
#'   county = rep(c("A", "B"), each = 3),
#'   demo_group = rep(c(1, 2), 3),
#'   estimate = c(0.4, 0.5, 0.6, 0.3, 0.35, 0.4),
#'   estimate_se = c(0.05, 0.04, 0.06, 0.03, 0.02, 0.05),
#'   N = c(100, 200, 300, 150, 250, 100)
#' )
#' poststratify(ps_table = ps, outcomes = estimate, ses = TRUE,
#'              weight = N, by = county)
#'
#' # --- calibrated_summary method (mock object) ---
#' cs_results <- tibble::tibble(
#'   .cellid        = 1:4,
#'   state          = c("PA", "PA", "OH", "OH"),
#'   est_n          = c(100, 200, 150, 50),
#'   voteshare_mean = c(0.55, 0.60, 0.48, 0.52),
#'   voteshare_se   = c(0.03, 0.02, 0.04, 0.03)
#' )
#' cs_obj <- structure(list(results = cs_results),
#'                     class = c("calibratedMRP", "calibrated_summary", "list"))
#' poststratify(cs_obj, weight = est_n, by = state)
#'
#' # --- calibrated_draws method (mock object) ---
#' cd_results <- tibble::tibble(
#'   .cellid   = rep(1:4, 3),
#'   .draw     = rep(1:3, each = 4),
#'   state     = rep(c("PA", "PA", "OH", "OH"), 3),
#'   est_n     = rep(c(100, 200, 150, 50), 3),
#'   voteshare = c(0.55, 0.60, 0.48, 0.52,
#'                 0.57, 0.62, 0.50, 0.54,
#'                 0.53, 0.58, 0.46, 0.50)
#' )
#' cd_obj <- structure(list(results = cd_results),
#'                     class = c("calibratedMRP", "calibrated_draws", "list"))
#' # Returns posterior mean and SD per state:
#' poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state)
#' # Returns full per-draw estimates:
#' poststratify(cd_obj, outcomes = "voteshare", weight = est_n, by = state,
#'              posterior_summary = FALSE)
poststratify <- function(ps_table, ...) UseMethod("poststratify")


#' @rdname poststratify
#' @method poststratify default
#' @export
poststratify.default <- function(ps_table, outcomes, ses = FALSE, se_suffix = "_se",
                                  weight, by, n_out = "n", na.rm = TRUE) {

  weight_var <- rlang::enquo(weight)
  mean_vars <- tidyselect::eval_select(rlang::enquo(outcomes), ps_table) |> names()
  by_chr <- names(tidyselect::eval_select(rlang::enquo(by), ps_table))
  se_vars <- paste0(mean_vars, se_suffix)


  out <- ps_table |>
    dplyr::summarise(dplyr::across({{ outcomes }},
                                   \(x) stats::weighted.mean(x, {{ weight }}, na.rm = na.rm)),
                     {{ n_out }} := sum({{ weight }}, na.rm = na.rm),
                     .by = dplyr::all_of(by_chr))

  if (ses) {
    if (!.calibratedMRP_env$poststratify_se_warned) {
      rlang::inform(c("Computing standard errors for poststratified estimates under the assumption of independence across cells.",
                      "x" = "True uncertainty may be larger or smaller depending on covariance of posterior means of `outcomes` across cells",
                      "i" = "For full Bayesian inference, use `calibrate_mrp(..., method = 'bayes')`",
                      "i" = "This message will only be displayed once per session"))
      .calibratedMRP_env$poststratify_se_warned <- TRUE
    }

    out_se <- ps_table |>
      dplyr::summarise(dplyr::across(all_of(se_vars),
                                     \(x) sqrt(sum({{ weight }}^2 * x^2, na.rm = na.rm)
                                               / sum({{ weight }}, na.rm = na.rm)^2 )),
                       .by = dplyr::all_of(by_chr))
    out <- dplyr::left_join(out, out_se, by = by_chr)
  }

  # Reorder columns: [grouping vars], outcome1, outcome1_se, outcome2, outcome2_se, ..., n_out
  group_names <- by_chr


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


#' @rdname poststratify
#' @method poststratify calibrated_summary
#' @export
poststratify.calibrated_summary <- function(ps_table, outcomes = NULL,
                                             weight, by,
                                             n_out = "n", na.rm = TRUE, ...) {

  # Validate that $results exists
  if (is.null(ps_table$results)) {
    rlang::abort("'ps_table$results' is NULL or missing. The 'calibrated_summary' object must have a '$results' tibble.")
  }
  results_tbl <- ps_table$results

  # Capture weight and by quosures
  weight_quo <- rlang::enquo(weight)
  by_chr <- names(tidyselect::eval_select(rlang::enquo(by), results_tbl))

  # Resolve outcome names
  if (is.null(outcomes)) {
    # Auto-detect: find _mean columns that have paired _se columns
    all_cols <- names(results_tbl)
    mean_cols_all <- all_cols[grepl("_mean$", all_cols)]
    paired <- vapply(mean_cols_all, function(col) {
      se_col <- sub("_mean$", "_se", col)
      se_col %in% all_cols
    }, logical(1))
    mean_cols <- mean_cols_all[paired]
    if (length(mean_cols) == 0L) {
      rlang::abort("No outcome columns found. Expected columns ending in '_mean' with paired '_se' columns in '$results'.")
    }
    outcome_names <- sub("_mean$", "", mean_cols)
  } else {
    # Validate supplied base names
    outcome_names <- outcomes
    missing_mean <- paste0(outcome_names, "_mean")[!paste0(outcome_names, "_mean") %in% names(results_tbl)]
    missing_se   <- paste0(outcome_names, "_se")[!paste0(outcome_names, "_se") %in% names(results_tbl)]
    missing_all  <- c(missing_mean, missing_se)
    if (length(missing_all) > 0L) {
      rlang::abort(paste0(
        "The following expected columns are missing from '$results': ",
        paste(missing_all, collapse = ", "), "."
      ))
    }
    mean_cols <- paste0(outcome_names, "_mean")
  }

  se_cols    <- paste0(outcome_names, "_se")
  n_out_name <- as.character(n_out)

  # Compute weighted means for _mean columns
  out <- results_tbl |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(mean_cols),
                    \(x) stats::weighted.mean(x, !!weight_quo, na.rm = na.rm)),
      !!n_out_name := sum(!!weight_quo, na.rm = na.rm),
      .by = dplyr::all_of(by_chr)
    )

  # Compute propagated SEs for _se columns
  out_se <- results_tbl |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(se_cols),
                    \(x) sqrt(sum((!!weight_quo)^2 * x^2, na.rm = na.rm) /
                                sum(!!weight_quo, na.rm = na.rm)^2)),
      .by = dplyr::all_of(by_chr)
    )

  # Join SEs back
  out <- dplyr::left_join(out, out_se, by = by_chr)

  # Reorder columns: [by_chr, outcome1_mean, outcome1_se, ..., n_out]
  col_order <- c(by_chr, as.vector(rbind(mean_cols, se_cols)), n_out_name)
  dplyr::select(out, dplyr::all_of(col_order))
}


#' @rdname poststratify
#' @method poststratify calibrated_draws
#' @export
poststratify.calibrated_draws <- function(ps_table, outcomes,
                                           weight, by,
                                           posterior_summary = TRUE,
                                           ses = FALSE,
                                           n_out = "n", na.rm = TRUE, ...) {

  # Check: ses argument not supported for draws
  if (!missing(ses) && isTRUE(ses)) {
    rlang::abort(
      "The 'ses' argument is not supported for 'calibrated_draws' objects. To obtain posterior uncertainty, use 'posterior_summary = TRUE' which returns posterior mean and SD across draws."
    )
  }

  # Check: outcomes must be provided
  if (missing(outcomes) || is.null(outcomes)) {
    rlang::abort(
      paste0("'outcomes' must be specified for 'calibrated_draws' objects. ",
             "Pass a character vector of outcome column names ",
             "(e.g. outcomes = c(\"voteshare\", \"turnout\")).")
    )
  }

  # Validate that $results exists and has .draw
  if (is.null(ps_table$results)) {
    rlang::abort("'ps_table$results' is NULL or missing. The 'calibrated_draws' object must have a '$results' tibble.")
  }
  results_tbl <- ps_table$results

  if (!".draw" %in% names(results_tbl)) {
    rlang::abort("'$results' does not contain a '.draw' column. Is this really a 'calibrated_draws' object?")
  }

  # Validate outcomes exist in results
  missing_outcomes <- outcomes[!outcomes %in% names(results_tbl)]
  if (length(missing_outcomes) > 0L) {
    rlang::abort(paste0(
      "The following outcome columns are missing from '$results': ",
      paste(missing_outcomes, collapse = ", "), "."
    ))
  }

  # Capture weight and by
  weight_quo <- rlang::enquo(weight)
  by_chr     <- names(tidyselect::eval_select(rlang::enquo(by), results_tbl))
  n_out_name <- as.character(n_out)

  # Per-draw poststratification
  draw_ps <- results_tbl |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(outcomes),
                    \(x) stats::weighted.mean(x, !!weight_quo, na.rm = na.rm)),
      !!n_out_name := sum(!!weight_quo, na.rm = na.rm),
      .by = c(".draw", dplyr::all_of(by_chr))
    )

  if (!posterior_summary) {
    # Return full per-draw tibble: .draw, by_chr, outcomes, n_out
    col_order <- c(".draw", by_chr, outcomes, n_out_name)
    return(dplyr::select(draw_ps, dplyr::all_of(col_order)))
  }

  # Summarize across draws: mean and SD per by-group
  out <- draw_ps |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(outcomes),
                    list(mean = mean, sd = stats::sd),
                    .names = "{.col}_{.fn}"),
      !!n_out_name := mean(!!rlang::sym(n_out_name)),
      .by = dplyr::all_of(by_chr)
    )

  # Reorder: [by_chr, outcome1_mean, outcome1_sd, ..., n_out]
  mean_cols <- paste0(outcomes, "_mean")
  sd_cols   <- paste0(outcomes, "_sd")
  col_order <- c(by_chr, as.vector(rbind(mean_cols, sd_cols)), n_out_name)
  dplyr::select(out, dplyr::all_of(col_order))
}
