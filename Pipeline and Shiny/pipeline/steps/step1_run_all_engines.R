#!/usr/bin/env Rscript
## Step 1: fit every engine. For each engine we source the matching engine
## body, which runs nested LODO tuning + locked external evaluation.

if (!exists("ROOT", inherits = TRUE)) {
  ROOT <- normalizePath(Sys.getenv("COMBO_PROJECT_ROOT", getwd()), mustWork = TRUE)
}
steps <- file.path(ROOT, "pipeline", "steps")
if (!exists("combo_engine_names")) source(file.path(steps, "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
known_engines <- combo_engine_names()
# Only treat arguments that name known engines as an engine selection. This
# avoids misreading the top-K argument (e.g. "5,10,15") passed through from
# run.R as an engine name.
arg_engines <- args[args %in% known_engines]
engines <- if (length(arg_engines)) arg_engines else known_engines
message("Engines: ", paste(engines, collapse = ", "))

for (engine in engines) {
  message("\n>> ", engine)
  cfg <- combo_engine_configs()[[engine]]
  if (cfg$kind == "logistic") {
    source(file.path(steps, "step1_engine_logistic.R"), local = FALSE)
  } else {
    combo_apply_engine_config(engine)
    source(file.path(steps, "step1_engine_survival.R"), local = FALSE)
  }
}

message("\nStep 1 done -> ", file.path(ROOT, results_root()))
