#!/usr/bin/env Rscript

ROOT <- normalizePath(Sys.getenv("COMBO_PROJECT_ROOT", getwd()), mustWork = TRUE)

if (!exists("combo_engine_names")) source(file.path(ROOT, "pipeline", "steps", "helpers.R"))
eval_root <- file.path(ROOT, results_root())
out_dir <- file.path(eval_root, "evaluation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

engine_label <- c(
  penalized_cox = "Penalized Cox",
  elastic_net_cox = "Elastic-net Cox",
  ridge_cox = "Ridge Cox",
  lasso_cox = "Lasso Cox",
  random_survival_forest = "Random survival forest",
  xgboost_survival_cox = "XGBoost survival Cox",
  logistic_regression = "Logistic regression"
)

read_if <- function(p) {
  if (!file.exists(p)) return(NULL)
  read.delim(p, stringsAsFactors = FALSE)
}

as_flag <- as_flag_logical

pieces <- list()

for (eng in combo_engine_names()) {
  sub <- combo_engine_eval_subdir(eng)
  base <- file.path(eval_root, sub)
  lab <- unname(engine_label[eng])

  off <- read_if(file.path(base, "combo_official_nested_summary.tsv"))
  if (!is.null(off) && nrow(off)) {
    hold <- read_if(file.path(base, "combo_lodocv_nested_tuned_by_holdout.tsv"))
    bl_mean <- NA_real_
    full_mean <- NA_real_
    auroc_full_mean <- NA_real_
    if (!is.null(hold) && nrow(hold)) {
      lodo <- lodo_means_from_holdout(
        hold,
        as.integer(off$top_k_recurrence[1]),
        as.integer(off$top_k_metastasis[1])
      )
      bl_mean <- lodo$baseline
      full_mean <- lodo$full
      auroc_full_mean <- lodo$auroc_full
    }
    lodo_c <- if (is.finite(full_mean)) full_mean else as.numeric(off$mean_cindex_lodo_nested[1])
    lodo_a <- if (is.finite(auroc_full_mean)) auroc_full_mean else as.numeric(off$mean_auroc_lodo_nested[1])
    pieces[[length(pieces) + 1L]] <- data.frame(
      record_type = "official_locked",
      engine = eng,
      engine_label = lab,
      holdout_dataset = NA_character_,
      top_k_recurrence = as.integer(off$top_k_recurrence[1]),
      top_k_metastasis = as.integer(off$top_k_metastasis[1]),
      is_baseline = NA,
      n_recurrence_genes = as.integer(off$n_recurrence_genes[1]),
      inner_cv_cindex = as.numeric(off$selection_score[1]),
      outer_cindex_lodo = lodo_c,
      outer_auroc_lodo = lodo_a,
      cindex_external = as.numeric(off$cindex_external[1]),
      auroc_external = as.numeric(off$auroc_external[1]),
      cindex_external_baseline = if ("cindex_external_baseline" %in% names(off)) as.numeric(off$cindex_external_baseline[1]) else NA_real_,
      delta_cindex_external_full_minus_baseline = if ("delta_cindex_external_full_minus_baseline" %in% names(off)) as.numeric(off$delta_cindex_external_full_minus_baseline[1]) else NA_real_,
      mean_cindex_lodo_baseline = bl_mean,
      mean_cindex_lodo_full = full_mean,
      delta_lodo_full_minus_baseline = full_mean - bl_mean,
      model_id = as.character(off$model_id[1]),
      stringsAsFactors = FALSE
    )
  }

  hold <- read_if(file.path(base, "combo_lodocv_nested_tuned_by_holdout.tsv"))
  if (!is.null(hold) && nrow(hold)) {
    if ("is_baseline" %in% names(hold)) hold$is_baseline <- as_flag(hold$is_baseline)
    pieces[[length(pieces) + 1L]] <- data.frame(
      record_type = "lodo_holdout_fold",
      engine = eng,
      engine_label = lab,
      holdout_dataset = hold$holdout_dataset,
      top_k_recurrence = as.integer(hold$top_k_recurrence),
      top_k_metastasis = as.integer(hold$top_k_metastasis),
      is_baseline = hold$is_baseline,
      n_recurrence_genes = as.integer(hold$n_recurrence_genes),
      inner_cv_cindex = as.numeric(hold$inner_cv_cindex_at_winner),
      outer_cindex_lodo = as.numeric(hold$outer_cindex_nested),
      outer_auroc_lodo = as.numeric(hold$outer_auroc_nested),
      cindex_external = NA_real_,
      auroc_external = NA_real_,
      cindex_external_baseline = NA_real_,
      delta_cindex_external_full_minus_baseline = NA_real_,
      mean_cindex_lodo_baseline = NA_real_,
      mean_cindex_lodo_full = NA_real_,
      delta_lodo_full_minus_baseline = NA_real_,
      model_id = hold$model_id,
      stringsAsFactors = FALSE
    )
  }

  inner <- read_if(file.path(base, "combo_inner_cv_k_by_holdout.tsv"))
  if (!is.null(inner) && nrow(inner)) {
    pieces[[length(pieces) + 1L]] <- data.frame(
      record_type = "inner_cv_kr_km",
      engine = eng,
      engine_label = lab,
      holdout_dataset = inner$holdout_dataset,
      top_k_recurrence = as.integer(inner$top_k_recurrence),
      top_k_metastasis = as.integer(inner$top_k_metastasis),
      is_baseline = NA,
      n_recurrence_genes = NA_integer_,
      inner_cv_cindex = as.numeric(inner$inner_cv_cindex),
      outer_cindex_lodo = NA_real_,
      outer_auroc_lodo = as.numeric(inner$inner_cv_auroc),
      cindex_external = NA_real_,
      auroc_external = NA_real_,
      cindex_external_baseline = NA_real_,
      delta_cindex_external_full_minus_baseline = NA_real_,
      mean_cindex_lodo_baseline = NA_real_,
      mean_cindex_lodo_full = NA_real_,
      delta_lodo_full_minus_baseline = NA_real_,
      model_id = NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

all_out <- do.call(rbind, pieces)
out_path <- file.path(out_dir, "all_engines_all_results_combined.tsv")
write.table(all_out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(all_out), " rows -> ", out_path)
message("  official_locked: ", sum(all_out$record_type == "official_locked"))
message("  lodo_holdout_fold: ", sum(all_out$record_type == "lodo_holdout_fold"))
message("  inner_cv_kr_km_variance: ", sum(all_out$record_type == "inner_cv_kr_km_variance"))
