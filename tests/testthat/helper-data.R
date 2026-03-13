
# Shared test data constructors
# Automatically sourced by testthat before tests run.

#' Create a synthetic poststratification table
#'
#' @param n_geo Number of geographies
#' @param n_cells_per_geo Number of cells per geography
#' @param outcomes Character vector of outcome column names
#' @param seed Random seed
#' @return A tibble with geography, outcome columns, and weight column
make_test_ps <- function(n_geo = 3, n_cells_per_geo = 50,
                         outcomes = c("voteshare", "turnout"),
                         seed = 42) {
  set.seed(seed)
  n <- n_geo * n_cells_per_geo
  geo_labels <- paste0("geo_", seq_len(n_geo))

  ps <- tibble::tibble(
    geography = rep(geo_labels, each = n_cells_per_geo),
    weight = runif(n, 0.2, 1.0)
  )

  for (outcome in outcomes) {
    ps[[outcome]] <- plogis(rnorm(n, qlogis(0.4), 0.5))
  }
  ps
}


#' Create targets from a PS table (weighted means + optional perturbation)
#'
#' @param ps A PS table from make_test_ps
#' @param geography String name of geography column
#' @param outcomes Character vector of outcome column names
#' @param perturb Numeric perturbation to add to true weighted means
#' @return A tibble with one row per geography and target values
make_test_targets <- function(ps, geography = "geography",
                              outcomes = c("voteshare", "turnout"),
                              perturb = 0.05) {
  geos <- unique(ps[[geography]])
  targets <- tibble::tibble(!!geography := geos)

  for (outcome in outcomes) {
    vals <- numeric(length(geos))
    for (i in seq_along(geos)) {
      mask <- ps[[geography]] == geos[i]
      wm <- weighted.mean(ps[[outcome]][mask], ps[["weight"]][mask])
      # Perturb and clamp to (0,1)
      vals[i] <- pmin(pmax(wm + perturb, 0.01), 0.99)
    }
    targets[[outcome]] <- vals
  }
  targets
}


#' Create a shift data frame
#'
#' @param geographies Character vector of geography labels
#' @param outcomes Character vector of outcome names
#' @param shift_vals Named list or numeric vector of shift values
#' @param geography_col Name of geography column
#' @param suffix Shift column suffix
#' @return A tibble with geography and shift columns
make_test_shifts <- function(geographies = paste0("geo_", 1:3),
                             outcomes = c("voteshare"),
                             shift_vals = 0.5,
                             geography_col = "geography",
                             suffix = "_shift") {
  shifts <- tibble::tibble(!!geography_col := geographies)
  shift_vals <- rep_len(shift_vals, length(outcomes) * length(geographies))
  idx <- 1
  for (outcome in outcomes) {
    col_name <- paste0(outcome, suffix)
    vals <- shift_vals[idx:(idx + length(geographies) - 1)]
    shifts[[col_name]] <- vals
    idx <- idx + length(geographies)
  }
  shifts
}


#' Create a covariance matrix with specified correlation structure
#'
#' @param outcomes Character vector of outcome names (used as dimnames)
#' @param corr Correlation between outcomes (scalar or matrix)
#' @param sds Standard deviations for each outcome
#' @return A covariance matrix
make_test_cov <- function(outcomes = c("voteshare", "turnout", "approval"),
                          corr = 0.6,
                          sds = rep(1, length(outcomes))) {
  k <- length(outcomes)
  if (is.matrix(corr)) {
    cor_mat <- corr
  } else {
    cor_mat <- matrix(corr, nrow = k, ncol = k)
    diag(cor_mat) <- 1
  }
  cov_mat <- diag(sds) %*% cor_mat %*% diag(sds)
  dimnames(cov_mat) <- list(outcomes, outcomes)
  cov_mat
}
