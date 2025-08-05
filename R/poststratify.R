

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
