#!/usr/bin/env Rscript

ROOT <- normalizePath(Sys.getenv("COMBO_PROJECT_ROOT", getwd()), mustWork = TRUE)
eval_root <- file.path(ROOT, results_root())
out_dir <- file.path(eval_root, "evaluation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

if (!exists("combo_engine_names")) source(file.path(ROOT, "pipeline", "steps", "helpers.R"))
engines <- setNames(
  lapply(combo_engine_names(), combo_engine_eval_subdir),
  combo_engine_names()
)

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

rbind_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (!length(dfs)) return(NULL)
  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  out <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    if (length(miss)) {
      for (m in miss) d[[m]] <- NA
    }
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, out)
}

as_flag <- as_flag_logical

normalize_grid_metrics <- function(d, engine, source_tag) {
  if (is.null(d) || !nrow(d)) return(NULL)
  d$engine <- engine
  d$engine_label <- unname(engine_label[engine])
  d$data_source <- source_tag
  if ("is_baseline" %in% names(d)) d$is_baseline <- as_flag(d$is_baseline)
  if ("mean_cindex_lodocv" %in% names(d)) {
    d$mean_cindex_lodo_nested <- d$mean_cindex_lodocv
    d$mean_auroc_lodo_nested <- d$mean_auc_lodocv
  }
  if ("mean_cindex_lodo" %in% names(d) && !"mean_cindex_lodo_nested" %in% names(d)) {
    d$mean_cindex_lodo_nested <- d$mean_cindex_lodo
    d$mean_auroc_lodo_nested <- d$mean_auroc_lodo
  }
  if ("auc_external" %in% names(d) && !"auroc_external" %in% names(d)) {
    d$auroc_external <- d$auc_external
  }
  if (!"is_baseline" %in% names(d)) {
    d$is_baseline <- grepl("^Baseline", d$model_id %||% "", ignore.case = TRUE) |
      (isFALSE(d$uses_metastasis %||% TRUE) & isTRUE(d$uses_recurrence %||% FALSE))
  }
  if (!"model_type" %in% names(d)) {
    d$model_type <- ifelse(isTRUE(d$is_baseline), "Baseline: Clinical+Recurrence", "Clinical+Recurrence+Metastasis")
  }
  d
}

read_engine_grid <- function(engine, subdir) {
  base <- file.path(eval_root, subdir)
  g1 <- normalize_grid_metrics(
    read_if(file.path(base, "combo_lodocv_external_metrics.tsv")),
    engine, "combo_lodocv_external_metrics"
  )
  g2 <- normalize_grid_metrics(
    read_if(file.path(base, "combo_lodo_grid_compare.tsv")),
    engine, "combo_lodo_grid_compare"
  )
  if (!is.null(g1) && nrow(g1)) return(g1)
  g2
}

read_engine_official <- function(engine, subdir) {
  base <- file.path(eval_root, subdir)
  d <- read_if(file.path(base, "combo_official_nested_summary.tsv"))
  if (is.null(d) || !nrow(d)) {
    d <- read_if(file.path(base, "combo_external_nested_final_locked.tsv"))
  }
  if (is.null(d) || !nrow(d)) return(NULL)
  d$engine <- engine
  d$engine_label <- unname(engine_label[engine])
  if ("model_engine" %in% names(d) && !"engine_display" %in% names(d)) {
    d$engine_display <- d$model_engine
  }
  d
}

read_holdout_nested <- function(subdir) {
  d <- read_if(file.path(eval_root, subdir, "combo_lodocv_nested_tuned_by_holdout.tsv"))
  if (is.null(d) || !nrow(d)) return(NULL)
  if ("is_baseline" %in% names(d)) d$is_baseline <- as_flag(d$is_baseline)
  d
}

lodo_from_holdout_fallback <- function(subdir, kr = 15L, km = 10L) {
  d <- read_holdout_nested(subdir)
  if (is.null(d)) {
    return(list(
      baseline = NA_real_, full = NA_real_,
      auroc_baseline = NA_real_, auroc_full = NA_real_,
      n_baseline_folds = 0L, n_full_folds = 0L
    ))
  }
  lodo_means_from_holdout(d, kr, km)
}

grid_from_holdout <- function(engine, subdir) {
  d <- read_holdout_nested(subdir)
  if (is.null(d)) return(NULL)
  key <- paste(
    d$is_baseline,
    d$top_k_recurrence,
    ifelse(is.na(d$top_k_metastasis), "NA", d$top_k_metastasis),
    sep = "|"
  )
  agg <- lapply(split(seq_len(nrow(d)), key, drop = TRUE), function(ii) {
    r <- d[ii[1], , drop = FALSE]
    data.frame(
      engine = engine,
      engine_label = unname(engine_label[engine]),
      data_source = "combo_lodocv_nested_tuned_by_holdout",
      model_type = if (r$is_baseline[1]) "Baseline: Clinical+Recurrence" else "Clinical+Recurrence+Metastasis",
      is_baseline = r$is_baseline[1],
      top_k_recurrence = as.integer(r$top_k_recurrence[1]),
      top_k_metastasis = if (r$is_baseline[1]) NA_integer_ else as.integer(r$top_k_metastasis[1]),
      uses_clinical = TRUE,
      uses_recurrence = TRUE,
      uses_metastasis = !r$is_baseline[1],
      mean_cindex_lodo_nested = mean(d$outer_cindex_nested[ii], na.rm = TRUE),
      mean_auroc_lodo_nested = mean(d$outer_auroc_nested[ii], na.rm = TRUE),
      cindex_external = NA_real_,
      auroc_external = NA_real_,
      model_id = r$model_id[1],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, agg)
}

read_all_official <- function() {
  pieces <- lapply(names(engines), function(nm) read_engine_official(nm, engines[[nm]]))
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(NULL)
  rbind_fill(pieces)
}

grid_pieces <- list()
compare_rows <- list()

official_by_engine <- read_all_official()

for (nm in names(engines)) {
  sub <- engines[[nm]]
  g <- read_engine_grid(nm, sub)
  if (is.null(g) || !nrow(g)) {
    g <- grid_from_holdout(nm, sub)
  }
  if (!is.null(g) && nrow(g)) grid_pieces[[length(grid_pieces) + 1L]] <- g

  off <- if (!is.null(official_by_engine) && nm %in% official_by_engine$engine) {
    official_by_engine[official_by_engine$engine == nm, , drop = FALSE]
  } else {
    read_engine_official(nm, sub)
  }

  kr_f <- if (!is.null(off) && nrow(off)) as.integer(off$top_k_recurrence[1]) else 15L
  km_f <- if (!is.null(off) && nrow(off)) as.integer(off$top_k_metastasis[1]) else 10L

  lodo_bl <- NA_real_
  lodo_full <- NA_real_
  ext_bl <- NA_real_
  ext_full <- NA_real_

  if (!is.null(g) && nrow(g)) {
    bl <- g[isTRUE(g$is_baseline), , drop = FALSE]
    if (!nrow(bl)) {
      bl <- g[grepl("Baseline", g$model_type %||% "", ignore.case = TRUE), , drop = FALSE]
    }
    if (nrow(bl)) {
      lodo_bl <- bl$mean_cindex_lodo_nested[1]
      if ("cindex_external" %in% names(bl)) ext_bl <- bl$cindex_external[1]
    }
    full_rows <- g[!isTRUE(g$is_baseline), , drop = FALSE]
    if (nrow(full_rows)) {
      hit <- full_rows[
        as.integer(full_rows$top_k_recurrence) == kr_f &
          as.integer(full_rows$top_k_metastasis) == km_f,
        ,
        drop = FALSE
      ]
      if (nrow(hit)) {
        o <- order(-hit$mean_cindex_lodo_nested, na.last = TRUE)
        hit <- hit[o[1], , drop = FALSE]
      } else {
        o <- order(-full_rows$mean_cindex_lodo_nested, na.last = TRUE)
        hit <- full_rows[o[1], , drop = FALSE]
      }
      lodo_full <- hit$mean_cindex_lodo_nested[1]
      if ("cindex_external" %in% names(hit)) ext_full <- hit$cindex_external[1]
    }
  }

  fb <- lodo_from_holdout_fallback(sub, kr = kr_f, km = km_f)
  if (fb$n_full_folds > 0L) {
    lodo_full <- fb$full
  }
  if (fb$n_baseline_folds > 0L) {
    lodo_bl <- fb$baseline
  }

  if (!is.null(off) && nrow(off)) {
    if (is.na(lodo_full)) lodo_full <- off$mean_cindex_lodo_nested[1]
    if (is.na(ext_full)) ext_full <- off$cindex_external[1]
    if (is.na(ext_bl) && "cindex_external_baseline" %in% names(off)) {
      ext_bl <- as.numeric(off$cindex_external_baseline[1])
    }
  }

  compare_rows[[length(compare_rows) + 1L]] <- data.frame(
    engine = nm,
    engine_label = unname(engine_label[nm]),
    top_k_recurrence_full = kr_f,
    top_k_metastasis_full = km_f,
    mean_cindex_lodo_nested_baseline = lodo_bl,
    mean_cindex_lodo_nested_full = lodo_full,
    delta_lodo_nested_full_minus_baseline = lodo_full - lodo_bl,
    cindex_external_baseline = ext_bl,
    cindex_external_full = ext_full,
    delta_cindex_external_full_minus_baseline = ext_full - ext_bl,
    stringsAsFactors = FALSE
  )
}

grid_all <- if (length(grid_pieces)) do.call(rbind, grid_pieces) else NULL
if (!is.null(grid_all) && nrow(grid_all) && !is.null(official_by_engine)) {
  for (eng in unique(grid_all$engine)) {
    off1 <- official_by_engine[official_by_engine$engine == eng, , drop = FALSE]
    if (!nrow(off1)) next
    ix <- grid_all$engine == eng & !grid_all$is_baseline
    if (!any(ix)) next
    kr1 <- as.integer(off1$top_k_recurrence[1])
    km1 <- as.integer(off1$top_k_metastasis[1])
    hit <- ix &
      as.integer(grid_all$top_k_recurrence) == kr1 &
      as.integer(grid_all$top_k_metastasis) == km1
    if (any(hit)) {
      grid_all$cindex_external[hit] <- off1$cindex_external[1]
      grid_all$auroc_external[hit] <- off1$auroc_external[1]
    }
  }
}

official_all <- official_by_engine
compare_all <- if (length(compare_rows)) do.call(rbind, compare_rows) else NULL

grid_cols <- c(
  "engine", "engine_label", "data_source", "model_type", "is_baseline",
  "top_k_recurrence", "top_k_metastasis",
  "uses_clinical", "uses_recurrence", "uses_metastasis",
  "mean_cindex_lodo_nested", "mean_auroc_lodo_nested",
  "cindex_external", "auroc_external", "model_id"
)

if (!is.null(grid_all) && nrow(grid_all)) {
  grid_all <- grid_all[order(
    grid_all$engine,
    grid_all$is_baseline,
    grid_all$top_k_recurrence,
    grid_all$top_k_metastasis
  ), , drop = FALSE]
  write.table(
    grid_all[, intersect(grid_cols, names(grid_all)), drop = FALSE],
    file.path(out_dir, "grid_all_combinations.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

if (!is.null(official_all) && nrow(official_all)) {
  write.table(
    official_all,
    file.path(out_dir, "official_nested_locked_per_engine.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

if (!is.null(compare_all) && nrow(compare_all)) {
  write.table(
    compare_all,
    file.path(out_dir, "baseline_vs_full_nested_lodo_external.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    compare_all[, c(
      "engine", "top_k_recurrence_full",
      "cindex_external_baseline", "cindex_external_full", "delta_cindex_external_full_minus_baseline"
    )],
    file.path(out_dir, "baseline_vs_full_external_cindex.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    compare_all[, c(
      "engine", "top_k_recurrence_full", "top_k_metastasis_full",
      "mean_cindex_lodo_nested_baseline", "mean_cindex_lodo_nested_full",
      "delta_lodo_nested_full_minus_baseline"
    )],
    file.path(out_dir, "baseline_aligned_lodo_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

if (!is.null(grid_all) && nrow(grid_all)) {
  full_only <- grid_all[!isTRUE(grid_all$is_baseline), , drop = FALSE]
  if (nrow(full_only)) {
    wide_rows <- lapply(split(full_only, full_only$engine, drop = TRUE), function(d) {
      rec_k <- sort(unique(as.integer(d$top_k_recurrence)))
      met_k <- sort(unique(as.integer(d$top_k_metastasis)))
      rec_k <- rec_k[is.finite(rec_k)]
      met_k <- met_k[is.finite(met_k)]
      if (!length(rec_k) || !length(met_k)) return(NULL)
      mat <- matrix(NA_real_, nrow = length(rec_k), ncol = length(met_k),
                    dimnames = list(paste0("kr", rec_k), paste0("km", met_k)))
      for (i in seq_len(nrow(d))) {
        r <- as.character(d$top_k_recurrence[i])
        m <- as.character(d$top_k_metastasis[i])
        if (paste0("kr", r) %in% rownames(mat) && paste0("km", m) %in% colnames(mat)) {
          mat[paste0("kr", r), paste0("km", m)] <- d$mean_cindex_lodo_nested[i]
        }
      }
      cbind(
        engine = unique(d$engine)[1],
        as.data.frame(mat, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE
      )
    })
    wide_rows <- wide_rows[!vapply(wide_rows, is.null, logical(1))]
    if (length(wide_rows)) {
      write.table(
        rbind_fill(wide_rows),
        file.path(out_dir, "grid_wide_lodo_nested_cindex.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
}

message("\nWrote evaluation tables to: ", out_dir)
if (!is.null(grid_all)) message("  grid_all_combinations.tsv (", nrow(grid_all), " rows)")
if (!is.null(compare_all)) message("  baseline_vs_full_nested_lodo_external.tsv (", nrow(compare_all), " engines)")
if (is.null(grid_all) || !nrow(grid_all)) {
  message(
    "\nNo grid files found. Run one of:\n",
    "  bash pipeline/run.sh                         # fit all engines + evaluation tables\n",
    "  bash \"model result/evaluation/run_all_evaluation.sh\"  # recompute tables from existing results\n",
    "Then re-run: Rscript \"model result/evaluation/combine_all_evaluation.R\""
  )
}

if (!is.null(compare_all)) {
  message("\nBaseline vs full (nested LODO C-index | external C-index):")
  for (i in seq_len(nrow(compare_all))) {
    message(sprintf(
      "  %s: LODO baseline=%.3f full=%.3f (delta %+.3f) | external baseline=%.3f full=%.3f (delta %+.3f)",
      compare_all$engine[i],
      compare_all$mean_cindex_lodo_nested_baseline[i],
      compare_all$mean_cindex_lodo_nested_full[i],
      compare_all$delta_lodo_nested_full_minus_baseline[i],
      compare_all$cindex_external_baseline[i],
      compare_all$cindex_external_full[i],
      compare_all$delta_cindex_external_full_minus_baseline[i]
    ))
  }
}
