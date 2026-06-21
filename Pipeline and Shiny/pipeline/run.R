#!/usr/bin/env Rscript
## Fit all engines, then build the evaluation and combined tables.
## Usage: Rscript pipeline/run.R [recurrence_top_k]   (default 5,10,15)

args <- commandArgs(trailingOnly = TRUE)
rec_k <- if (length(args)) args[[1]] else "5,10,15"

cmd <- commandArgs(trailingOnly = FALSE)
f <- gsub("~\\+~", " ", sub("^--file=", "", cmd[grep("^--file=", cmd)]))
pipeline_dir <- normalizePath(dirname(f), mustWork = TRUE)
ROOT <- normalizePath(file.path(pipeline_dir, ".."), mustWork = TRUE)
setwd(ROOT)

Sys.setenv(
  COMBO_PROJECT_ROOT = ROOT,
  COMBO_EVAL_ROOT = file.path("pipeline", "model result"),
  COMBO_SUMMARY_DIR = file.path("pipeline", "model result", "evaluation"),
  COMBO_RANDOM_SEED = "1",
  COMBO_REC_TOP_K = rec_k,
  COMBO_MET_TOP_K = Sys.getenv("COMBO_MET_TOP_K", rec_k)
)

steps <- file.path(pipeline_dir, "steps")
source(file.path(steps, "helpers.R"))

t0 <- Sys.time()
message("=== LUAD recurrence pipeline (7 engines) ===")
message("Project: ", ROOT, " | Top-K: ", rec_k)

message("\n--- Step 1: fit each engine (nested K + LODO + external) ---")
source(file.path(steps, "step1_run_all_engines.R"), local = FALSE)

message("\n--- Step 2: build evaluation tables ---")
source(file.path(steps, "step2_evaluate.R"), local = FALSE)

message("\n--- Step 3: merge all results into one TSV ---")
source(file.path(steps, "step3_combine.R"), local = FALSE)

message("\nFinished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " minutes")
message("Per-engine outputs: ", file.path(ROOT, results_root()))
message("Combined tables:    ", file.path(ROOT, results_root(), "evaluation"))
