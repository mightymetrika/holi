## R CMD check results

0 errors | 0 warnings | 0 notes

## Update Summary

Enhanced run_sim_rstar_glm() to gracefully handle instances where model fitting
fails. This update addresses a CRAN Package Check Results error.

Added a "trace = FALSE" parameter to rstar_glm(). This argument is passed to
likelihoodAsy::rstar() as "trace = if (interactive()) trace else FALSE" to
suppress error messages and verbose output during batch processing.

## revdepcheck

This package has no reverse dependencies.

