

calibrate_mrp <- function(model,
                          ps_table,
                          targets,
                          geography,
                          outcomes = NULL, # defaults to outcomes in `mod`
                          method = "plugin", # or "bayes"
                          uncertainty = "approximate", # or "bayes" or "none"
                          draw_ids = NULL,

){
  if (class(mod) != "brmsfit") rlang::abort("`mod` must be a `brmsfit` object")
  if (is.null(draw_ids)) draw_ids <- seq_len(ndraws(mod))
  if (is.null(outcomes)) {
    calibration_variables <- mod$formula[[2]]
    rlang::inform(c("No `outcomes` provided, defaulting to outcome variables from the model formula: ",
                    "i" = paste(calibration_variables, collapse = ", ")))
  }
  if (!any(calibration_variables %in% names(targets))) {
      rlang::abort(c("At least one variable in `outcomes` must be present in `targets` data frame to perform calibration.",
                     "*" = sprintf("Outcome variables: %s", paste(calibration_variables, collapse = ", ")),
                     "*" = sprintf("Columns in `targets`: %s", paste(names(targets), collapse = ", ")),
                     "i" = "brms automatically removes underscores from variable names; you may need to rename variables in `targets` to match."))
  }
  if (!method %in% c("plugin", "bayes")) rlang::abort("`method` must be either 'plugin' or 'bayes'")

  # calculate predictions
  if (method == "plugin") {
    ps_table <- generate_cell_estimates(model = model,
                                        ps_table = ps_table,
                                        draw_ids = draw_ids,
                                        summarize = TRUE)
    ps_table %>% head
  }
}
