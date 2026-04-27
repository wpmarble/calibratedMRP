# Package-level environment for session state.
# Initialized once when the package namespace is loaded.
.calibratedMRP_env <- new.env(parent = emptyenv())
.calibratedMRP_env$poststratify_se_warned <- FALSE
