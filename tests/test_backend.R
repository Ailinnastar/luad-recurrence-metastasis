# Run from the repository root with Rscript tests/test_backend.R.
`%||%` <- function(x, y) if (is.null(x)) y else x
source("Pipeline and Shiny/shiny_recurrence_app/model_backend.R")
ar <- load_deployed_artifacts()
genes <- required_signature_genes(ar)
stopifnot(length(genes) == 15L)
raw <- data.frame(SampleID = c("synthetic_a", "synthetic_b"), Age = c(50, 75), Stage_num = c(1, 4))
for (gene in genes) raw[[gene]] <- c(0, 0)
pred <- compute_predictions(ar, raw, 365)
stopifnot(length(pred$risks) == 2L, all(is.finite(pred$risks)))
stopifnot(all(pred$risks >= 0 & pred$risks <= 1), pred$risks[2] > pred$risks[1])
stopifnot(pred$matched_recur == 15L)
incomplete <- raw[, setdiff(names(raw), genes[[1]]), drop = FALSE]
err <- tryCatch(compute_predictions(ar, incomplete, 365), error = identity)
stopifnot(inherits(err, "error"), grepl("Missing:", conditionMessage(err)))
# Record the fallback calibration contract without implying horizon-specific validity.
stopifnot(identical(pred$risks, compute_predictions(ar, raw, 730)$risks))
cat("Backend checks passed: bounded scores, missing-gene validation, calibration contract.\n")
