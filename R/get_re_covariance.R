

#' Extract covariance of random effects across outcomes
#'
#' `get_re_covariance` returns the covariance parameters for random effects that are
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
#' @param show_progress Show a progress bar while constructing array? Applies if `tidy = FALSE`
#' @importFrom tidybayes spread_draws
#' @importFrom dplyr filter select rename_with mutate arrange group_by group_split bind_cols `%>%`
#' @export

get_re_covariance <- function(model,
                              group,
                              tidy = FALSE,
                              outcome_order = NULL,
                              draw_ids = NULL,
                              show_progress = FALSE) {


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
      furrr::future_map(.progress = show_progress,
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
