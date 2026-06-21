## Shared helper functions for the LUAD recurrence pipeline.
## Sourced by every step. Holds project paths, engine configs, metrics,
## gene ranking / scoring, and the nested leave-one-dataset-out machinery.

# ---- Project paths ----

project_root <- function() {
  normalizePath(Sys.getenv("COMBO_PROJECT_ROOT", getwd()), mustWork = TRUE)
}
combo_project_root <- project_root

pipeline_root <- function(root = project_root()) {
  file.path(root, "pipeline")
}
combo_pipeline_root <- pipeline_root

geo_recurrence_dir <- function(root = project_root()) {
  file.path(pipeline_root(root), "geo", "recurrence")
}

geo_external_dir <- function(root = project_root()) {
  file.path(pipeline_root(root), "geo", "external")
}

training_dir <- function(root = project_root()) {
  env <- Sys.getenv("COMBO_TRAIN_DIR", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  file.path(pipeline_root(root), "data", "training")
}
combo_train_dir <- training_dir

external_dir <- function(root = project_root()) {
  env <- Sys.getenv("COMBO_EXT_DIR", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  file.path(pipeline_root(root), "data", "external")
}
combo_ext_dir <- external_dir

frozen_metastasis_dir <- function(root = project_root()) {
  env <- Sys.getenv("COMBO_MATA_FROZEN_DIR", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  normalizePath(file.path(pipeline_root(root), "frozen"), mustWork = FALSE)
}
combo_mata_frozen_dir <- frozen_metastasis_dir

results_root <- function() {
  Sys.getenv("COMBO_EVAL_ROOT", file.path("pipeline", "model result"))
}
combo_eval_root <- results_root

# ---- Engine registry ----

combo_engine_configs <- function() {
  list(
    penalized_cox = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "cox_glmnet", COMBO_COX_PRESET = "penalized",
      COMBO_TUNE_ALPHA = TRUE, COMBO_ALPHA_GRID_OVERRIDE = c(0.75, 1),
      COMBO_COX_ALPHA = 1, COMBO_COX_LAMBDA = "min", COMBO_RANDOM_SEED = 1L
    ),
    lasso_cox = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "cox_glmnet", COMBO_COX_PRESET = "lasso",
      COMBO_TUNE_ALPHA = TRUE, COMBO_ALPHA_GRID_OVERRIDE = c(0.75, 1),
      COMBO_COX_ALPHA = 1, COMBO_COX_LAMBDA = "1se", COMBO_RANDOM_SEED = 1L
    ),
    ridge_cox = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "cox_glmnet", COMBO_COX_PRESET = "ridge",
      COMBO_TUNE_ALPHA = TRUE, COMBO_ALPHA_GRID_OVERRIDE = c(0, 0.25),
      COMBO_COX_ALPHA = 0, COMBO_COX_LAMBDA = "min", COMBO_RANDOM_SEED = 1L
    ),
    elastic_net_cox = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "cox_glmnet", COMBO_COX_PRESET = "elastic_net",
      COMBO_TUNE_ALPHA = TRUE, COMBO_ALPHA_GRID_OVERRIDE = c(0.25, 0.5, 0.75),
      COMBO_COX_ALPHA = 0.5, COMBO_COX_LAMBDA = "min", COMBO_RANDOM_SEED = 1L
    ),
    random_survival_forest = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "rsf", COMBO_RANDOM_SEED = 1L
    ),
    xgboost_survival_cox = list(
      kind = "survival", COMBO_SURVIVAL_MODEL = "xgboost", COMBO_RANDOM_SEED = 1L,
      COMBO_XGB_NROUNDS = 150L, COMBO_XGB_NFOLD = 3L
    ),
    logistic_regression = list(kind = "logistic")
  )
}

combo_engine_names <- function() names(combo_engine_configs())

combo_engine_eval_subdir <- function(engine) {
  cfg <- combo_engine_configs()[[engine]]
  if (identical(cfg$kind, "logistic")) return("logistic_regression")
  preset <- cfg$COMBO_COX_PRESET
  if (isTRUE(cfg$COMBO_TUNE_ALPHA) && nzchar(preset)) return(paste0(preset, "_cox_alpha_tuned"))
  if (cfg$COMBO_SURVIVAL_MODEL == "rsf") return("random_survival_forest")
  if (cfg$COMBO_SURVIVAL_MODEL == "xgboost") return("xgboost_survival_cox")
  paste0(preset, "_cox")
}

combo_apply_engine_config <- function(engine) {
  cfg <- combo_engine_configs()[[engine]]
  all_keys <- unique(unlist(
    lapply(combo_engine_configs(), function(x) setdiff(names(x), "kind")),
    use.names = FALSE
  ))
  for (k in all_keys) {
    if (exists(k, envir = globalenv(), inherits = FALSE)) rm(list = k, envir = globalenv())
    Sys.unsetenv(k)
  }
  for (nm in setdiff(names(cfg), "kind")) {
    val <- cfg[[nm]]
    assign(nm, val, envir = globalenv())
    if (is.atomic(val) && length(val) == 1L) {
      do.call(Sys.setenv, setNames(list(as.character(val)), nm))
    }
  }
  invisible(cfg)
}

combo_alpha_grid_for_preset <- function(preset) {
  switch(preset,
    ridge = c(0, 0.25),
    elastic_net = c(0.25, 0.5, 0.75),
    penalized = c(0.75, 1),
    lasso = c(0.75, 1),
    sort(unique(as.numeric(trimws(strsplit(Sys.getenv("COMBO_ALPHA_GRID", "0,0.25,0.5,0.75,1"), ",")[[1]]))))
  )
}

# ---- Metrics ----

auc_binary <- function(y, x) {
  y <- as.integer(y)
  ok <- is.finite(x) & !is.na(y)
  y <- y[ok]
  x <- as.numeric(x[ok])
  if (length(unique(y)) < 2L) return(NA_real_)
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE))))
  }
  r <- rank(x, ties.method = "average")
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

cindex_harrell <- function(time, event, risk) {
  ok <- is.finite(time) & time > 0 & event %in% c(0L, 1L) & is.finite(risk)
  if (sum(ok) < 5L) return(NA_real_)
  d <- data.frame(time = time[ok], event = event[ok], risk = risk[ok])
  con <- survival::concordance(survival::Surv(time, event) ~ risk, data = d, reverse = TRUE)
  unname(con$concordance)
}

na_metrics <- function() c(cindex = NA_real_, auroc = NA_real_)

named_surv <- function(time_vec, event_vec, sample_ids) {
  time_vec <- as.numeric(time_vec)
  event_vec <- as.integer(event_vec)
  names(time_vec) <- names(event_vec) <- sample_ids
  list(time = time_vec, event = event_vec)
}
combo_ensure_named_surv <- named_surv

subset_surv <- function(time_vec, event_vec, ids) {
  list(time = time_vec[ids], event = event_vec[ids])
}
combo_subset_surv <- subset_surv

as_flag_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes")
}

normalize_holdout_baseline <- function(hold) {
  if (is.null(hold) || !nrow(hold)) return(hold)
  if ("is_baseline" %in% names(hold)) {
    hold$is_baseline <- as_flag_logical(hold$is_baseline)
  } else {
    hold$is_baseline <- grepl("^Baseline", hold$model_id, ignore.case = TRUE)
  }
  hold
}

# Mean LODO C-index / AUROC from the per-holdout table. Baseline uses the 4
# baseline folds; the full model uses the 4 folds at the locked (kr, km).
lodo_means_from_holdout <- function(hold, kr, km) {
  empty <- list(
    baseline = NA_real_, full = NA_real_,
    auroc_baseline = NA_real_, auroc_full = NA_real_,
    n_baseline_folds = 0L, n_full_folds = 0L
  )
  if (is.null(hold) || !nrow(hold)) return(empty)
  hold <- normalize_holdout_baseline(hold)
  baseline <- hold[hold$is_baseline, , drop = FALSE]
  full <- hold[
    !hold$is_baseline &
      as.integer(hold$top_k_recurrence) == as.integer(kr) &
      as.integer(hold$top_k_metastasis) == as.integer(km),
    , drop = FALSE
  ]
  if (!nrow(full)) full <- hold[!hold$is_baseline, , drop = FALSE]
  list(
    baseline = if (nrow(baseline)) mean(baseline$outer_cindex_nested, na.rm = TRUE) else NA_real_,
    full = if (nrow(full)) mean(full$outer_cindex_nested, na.rm = TRUE) else NA_real_,
    auroc_baseline = if (nrow(baseline)) mean(baseline$outer_auroc_nested, na.rm = TRUE) else NA_real_,
    auroc_full = if (nrow(full)) mean(full$outer_auroc_nested, na.rm = TRUE) else NA_real_,
    n_baseline_folds = nrow(baseline),
    n_full_folds = nrow(full)
  )
}

# ---- Expression preprocessing & gene scoring ----

as_log2_scale <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  rng <- range(mat, na.rm = TRUE)
  if (is.finite(rng[2]) && rng[2] >= 80) log2(pmax(mat, 0) + 1) else mat
}

# ComBat on the training cohorts, then standardise train and test using
# train-only mean/sd (no leakage from the held-out cohort).
combat_train_only_preprocess <- function(expr_train, expr_test, batch_train) {
  expr_train <- as.matrix(expr_train)
  expr_test <- as.matrix(expr_test)
  if (is.null(rownames(expr_train))) rownames(expr_train) <- paste0("G", seq_len(nrow(expr_train)))
  if (is.null(rownames(expr_test))) rownames(expr_test) <- rownames(expr_train)[seq_len(nrow(expr_test))]

  log2_train <- as_log2_scale(expr_train)
  train_combat <- sva::ComBat(
    dat = log2_train, batch = as.factor(batch_train),
    mod = NULL, par.prior = TRUE, prior.plots = FALSE
  )

  train_mu <- rowMeans(train_combat, na.rm = TRUE)
  names(train_mu) <- rownames(train_combat)
  train_sd <- apply(train_combat, 1, stats::sd, na.rm = TRUE)
  names(train_sd) <- rownames(train_combat)
  train_sd[is.na(train_sd) | train_sd == 0] <- 1

  standardise <- function(mat, mu, sd) {
    shared <- intersect(rownames(mat), names(mu))
    mat <- mat[shared, , drop = FALSE]
    sweep(sweep(mat, 1, mu[shared], "-"), 1, sd[shared], "/")
  }

  train_scaled <- standardise(train_combat, train_mu, train_sd)
  test_scaled <- standardise(as_log2_scale(expr_test), train_mu, train_sd)
  shared <- intersect(rownames(train_scaled), rownames(test_scaled))
  list(
    train = train_scaled[shared, , drop = FALSE],
    test = test_scaled[shared, , drop = FALSE],
    train_mu = train_mu,
    train_sd = train_sd
  )
}

combat_train_test <- function(mat_train, mat_test, batch_train, batch_test_label = NULL) {
  prep <- combat_train_only_preprocess(mat_train, mat_test, batch_train)
  list(left = prep$train, right = prep$test)
}

# Mean of weighted, standardised gene values = a single risk score per sample.
signature_score <- function(mat, genes, w) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0L) return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  ww <- w[match(g, genes)]
  ok <- is.finite(ww)
  g <- g[ok]
  ww <- ww[ok]
  if (!length(g)) return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  out <- colSums(mat[g, , drop = FALSE] * ww) / length(g)
  if (is.null(names(out)) && !is.null(colnames(mat))) names(out) <- colnames(mat)
  out
}

# Rank recurrence genes by limma differential expression (recurred vs not).
ranked_limma_recurrence <- function(mat_tr, y_tr, fdr = 0.05) {
  y_tr <- as.integer(y_tr)
  if (length(unique(y_tr)) < 2L) {
    return(data.frame(gene = character(), logFC = numeric(), adj.P.Val = numeric()))
  }
  design <- model.matrix(~y_tr)
  colnames(design)[2] <- "Recurrence"
  fit <- limma::eBayes(limma::lmFit(mat_tr, design))
  tt <- limma::topTable(fit, coef = "Recurrence", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[is.finite(tt$logFC), , drop = FALSE]
  sig <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < fdr, , drop = FALSE]
  if (nrow(sig) > 0L) return(sig[order(-abs(sig$logFC)), , drop = FALSE])
  tt[order(-abs(tt$logFC)), , drop = FALSE]
}

# Optional alternative ranking: univariate Cox per gene (used only when
# COMBO_RECURRENCE_RANK is set to cox; default is limma above).
ranked_univariate_cox_recurrence <- function(mat_tr, time_tr, event_tr, fdr = 0.05, min_events = 4L) {
  time_tr <- as.numeric(time_tr)
  event_tr <- as.integer(event_tr)
  ok <- is.finite(time_tr) & time_tr > 0 & event_tr %in% c(0L, 1L)
  if (sum(ok) < 25L || sum(event_tr[ok] == 1L) < min_events) {
    return(data.frame(gene = character(), logFC = numeric(), adj.P.Val = numeric()))
  }
  time_tr <- time_tr[ok]
  event_tr <- event_tr[ok]
  mat_tr <- mat_tr[, ok, drop = FALSE]
  genes <- rownames(mat_tr)
  coefs <- setNames(rep(NA_real_, length(genes)), genes)
  pvals <- setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    g <- as.numeric(mat_tr[i, ])
    if (sum(is.finite(g)) < 10L || stats::sd(g, na.rm = TRUE) < 1e-10) next
    dat <- data.frame(time = time_tr, event = event_tr, gene = g)
    fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ gene, data = dat), error = function(e) NULL)
    if (is.null(fit)) next
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    if (is.null(sm)) next
    coefs[i] <- suppressWarnings(as.numeric(stats::coef(fit)[1]))
    pvals[i] <- suppressWarnings(as.numeric(sm$coefficients[1, "Pr(>|z|)"]))
  }
  res <- data.frame(gene = genes, logFC = coefs, pvalue = pvals, stringsAsFactors = FALSE)
  res <- res[is.finite(res$logFC) & is.finite(res$pvalue), , drop = FALSE]
  if (!nrow(res)) return(data.frame(gene = character(), logFC = numeric(), adj.P.Val = numeric()))
  res$adj.P.Val <- p.adjust(res$pvalue, method = "BH")
  sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < fdr, , drop = FALSE]
  if (nrow(sig) > 0L) return(sig[order(sig$adj.P.Val, -abs(sig$logFC), sig$pvalue), , drop = FALSE])
  res[order(res$adj.P.Val, -abs(res$logFC), res$pvalue), , drop = FALSE]
}

rank_recurrence_genes <- function(mat_tr, y_tr, time_tr = NULL, event_tr = NULL, fdr = 0.05) {
  mode <- tolower(Sys.getenv("COMBO_RECURRENCE_RANK", "limma"))
  if (mode %in% c("cox", "univariate_cox", "uni_cox")) {
    ranked_univariate_cox_recurrence(mat_tr, time_tr, event_tr, fdr = fdr)
  } else {
    ranked_limma_recurrence(mat_tr, y_tr, fdr = fdr)
  }
}

select_recurrence_topk <- function(ranked, kr) head(ranked, min(kr, nrow(ranked)))

frozen_metastasis_gene_path <- function(dir, km) {
  file.path(dir, paste0("frozen_A_sig_top", sprintf("%02d", as.integer(km)), "_genes.tsv"))
}

load_frozen_metastasis_gene_lists <- function(dir, km_vals = c(5L, 10L, 15L)) {
  out <- list()
  for (km in km_vals) {
    tt <- read.table(frozen_metastasis_gene_path(dir, km), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    gcol <- intersect(c("gene", "Gene", "symbol"), names(tt))[1]
    if (is.na(gcol)) gcol <- names(tt)[1]
    genes <- toupper(trimws(as.character(tt[[gcol]])))
    out[[as.character(as.integer(km))]] <- unique(genes[!is.na(genes) & nzchar(genes)])
  }
  out
}

# Re-rank the frozen metastasis genes inside the training fold to get fresh
# weights (signed logFC). Keeps the metastasis score train-only.
refit_metastasis_weights <- function(mat_tr, genes, y_tr, time_tr = NULL, event_tr = NULL,
                                      rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  if (is.numeric(genes) && !is.null(names(genes))) genes <- names(genes)
  genes <- unique(toupper(trimws(as.character(genes))))
  genes <- intersect(genes, toupper(rownames(mat_tr)))
  if (!length(genes)) return(setNames(numeric(), character()))
  tt <- rank_genes_fn(mat_tr[genes, , drop = FALSE], y_tr, time_tr = time_tr, event_tr = event_tr, fdr = fdr)
  if (!nrow(tt) || !"gene" %in% names(tt)) return(setNames(numeric(), character()))
  tt$gene <- toupper(trimws(as.character(tt$gene)))
  hit <- tt[tt$gene %in% genes, , drop = FALSE]
  if (!nrow(hit)) return(setNames(numeric(), character()))
  hit <- hit[!duplicated(hit$gene), , drop = FALSE]
  ord <- match(genes, hit$gene)
  hit <- hit[ord[is.finite(ord)], , drop = FALSE]
  setNames(as.numeric(hit$logFC), hit$gene)
}

metastasis_train_test_scores <- function(adj_tr, adj_te, meta_genes, km, y_tr, time_tr = NULL,
                                          event_tr = NULL, rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  genes <- meta_genes[[as.character(as.integer(km))]]
  if (is.null(genes) || !length(genes)) {
    na_tr <- setNames(rep(NA_real_, ncol(adj_tr)), colnames(adj_tr))
    na_te <- setNames(rep(NA_real_, ncol(adj_te)), colnames(adj_te))
    return(list(train = na_tr, test = na_te, weights = setNames(numeric(), character())))
  }
  w <- refit_metastasis_weights(adj_tr, genes, y_tr, time_tr = time_tr, event_tr = event_tr,
                                rank_genes_fn = rank_genes_fn, fdr = fdr)
  list(
    train = as.numeric(signature_score(adj_tr, names(w), w)),
    test = as.numeric(signature_score(adj_te, names(w), w)),
    weights = w
  )
}

# ---- Nested LODO tuning ----

combo_k_grid_values <- function() {
  parse_k <- function(env_nm, fallback) {
    k <- suppressWarnings(as.integer(trimws(strsplit(Sys.getenv(env_nm, fallback), ",", fixed = TRUE)[[1]])))
    k[is.finite(k) & k > 0L]
  }
  rec <- parse_k("COMBO_REC_TOP_K", "5,10,15")
  if (!length(rec)) rec <- c(5L, 10L, 15L)
  met <- parse_k("COMBO_MET_TOP_K", paste(rec, collapse = ","))
  if (!length(met)) met <- rec
  list(recurrence = rec, metastasis = met)
}

combo_inner_cv_folds <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("COMBO_INNER_CV_FOLDS", "5")))
  if (!is.finite(n) || n < 2L) 5L else n
}

# Build the clinical + recurrence-score feature table (baseline model).
build_clin_rec_baseline_features <- function(adj_tr, adj_te, ranked, kr, Xc_tr, Xc_te) {
  X_tr <- as.data.frame(Xc_tr, check.names = FALSE, stringsAsFactors = FALSE)
  X_te <- as.data.frame(Xc_te, check.names = FALSE, stringsAsFactors = FALSE)
  top <- select_recurrence_topk(ranked, kr)
  if (nrow(top) == 0L) return(list(ok = FALSE, X_tr = X_tr, X_te = X_te, n_recurrence_genes = 0L))
  rw <- setNames(top$logFC, top$gene)
  X_tr$Recurrence_score <- as.numeric(signature_score(adj_tr, top$gene, rw))
  X_te$Recurrence_score <- as.numeric(signature_score(adj_te, top$gene, rw))
  list(ok = TRUE, X_tr = X_tr, X_te = X_te, n_recurrence_genes = nrow(top))
}

# Build clinical + recurrence-score + metastasis-score features (full model).
build_clin_rec_met_features <- function(adj_tr, adj_te, ranked, meta_genes, kr, km, Xc_tr, Xc_te,
                                         y_tr, time_tr = NULL, event_tr = NULL, meta_w = NULL,
                                         rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  if (is.null(meta_genes) && !is.null(meta_w)) meta_genes <- meta_w
  X_tr <- as.data.frame(Xc_tr, check.names = FALSE, stringsAsFactors = FALSE)
  X_te <- as.data.frame(Xc_te, check.names = FALSE, stringsAsFactors = FALSE)
  ms <- metastasis_train_test_scores(adj_tr, adj_te, meta_genes, km, y_tr, time_tr = time_tr,
                                     event_tr = event_tr, rank_genes_fn = rank_genes_fn, fdr = fdr)
  if (!length(ms$weights)) return(list(ok = FALSE, X_tr = X_tr, X_te = X_te, n_recurrence_genes = 0L))
  X_tr$Metastasis_score <- ms$train
  X_te$Metastasis_score <- ms$test
  top <- select_recurrence_topk(ranked, kr)
  if (nrow(top) == 0L) return(list(ok = FALSE, X_tr = X_tr, X_te = X_te, n_recurrence_genes = 0L))
  rw <- setNames(top$logFC, top$gene)
  X_tr$Recurrence_score <- as.numeric(signature_score(adj_tr, top$gene, rw))
  X_te$Recurrence_score <- as.numeric(signature_score(adj_te, top$gene, rw))
  list(ok = TRUE, X_tr = X_tr, X_te = X_te, n_recurrence_genes = nrow(top))
}

# Pick the (kr, km) row with the best metric; ties broken toward fewer genes.
pick_best_k_row <- function(df, metric_col = "inner_cv_cindex") {
  if (is.null(df) || !nrow(df)) return(NULL)
  v <- df[[metric_col]]
  if (!length(v) || all(!is.finite(v))) return(df[1, , drop = FALSE])
  mx <- max(v, na.rm = TRUE)
  cand <- df[is.finite(v) & abs(v - mx) < 1e-9, , drop = FALSE]
  if (!nrow(cand)) cand <- df
  ord <- order(
    -cand[[metric_col]],
    cand$top_k_recurrence + cand$top_k_metastasis,
    cand$top_k_recurrence,
    cand$top_k_metastasis
  )
  cand[ord[1], , drop = FALSE]
}

# Inner 5-fold CV over the kr x km grid, for one held-out cohort.
inner_cv_k_grid_one_holdout <- function(expr_gene, ptr, pte, held_out, meta_genes, meta_w = NULL,
                                        kr_vals, km_vals, clinical_cols, n_inner, seed,
                                        time_tr_all, event_tr_all, time_te_outer, event_te_outer,
                                        fit_metrics_fn, rank_genes_fn, fdr = 0.05) {
  tr_ids <- rownames(ptr)
  surv_all <- named_surv(time_tr_all, event_tr_all, tr_ids)
  time_tr_all <- surv_all$time
  event_tr_all <- surv_all$event
  set.seed(seed)
  inner_assign <- sample(rep(seq_len(n_inner), length.out = length(tr_ids)))
  grid <- expand.grid(top_k_recurrence = kr_vals, top_k_metastasis = km_vals, stringsAsFactors = FALSE)

  fold_detail_rows <- list()
  summary_rows <- list()
  for (gi in seq_len(nrow(grid))) {
    kr <- grid$top_k_recurrence[gi]
    km <- grid$top_k_metastasis[gi]
    fold_cx <- rep(NA_real_, n_inner)
    fold_au <- rep(NA_real_, n_inner)
    for (f in seq_len(n_inner)) {
      val_ids <- tr_ids[inner_assign == f]
      in_tr_ids <- tr_ids[inner_assign != f]
      if (length(val_ids) < 5L || length(in_tr_ids) < 15L) next
      ptr_in <- ptr[in_tr_ids, , drop = FALSE]
      pte_in <- ptr[val_ids, , drop = FALSE]
      cb_in <- combat_train_test(expr_gene[, in_tr_ids, drop = FALSE], expr_gene[, val_ids, drop = FALSE], ptr_in$Dataset, held_out)
      surv_in <- subset_surv(time_tr_all, event_tr_all, in_tr_ids)
      ranked_in <- rank_genes_fn(cb_in$left, ptr_in$Recurrence, time_tr = surv_in$time, event_tr = surv_in$event, fdr = fdr)
      imp_in <- impute_clinical(
        as.data.frame(ptr_in[, clinical_cols, drop = FALSE]),
        as.data.frame(pte_in[, clinical_cols, drop = FALSE])
      )
      feat_in <- build_clin_rec_met_features(
        cb_in$left, cb_in$right, ranked_in, meta_genes, kr, km,
        imp_in$train, imp_in$test, y_tr = ptr_in$Recurrence,
        time_tr = surv_in$time, event_tr = surv_in$event, meta_w = meta_w,
        rank_genes_fn = rank_genes_fn, fdr = fdr
      )
      if (!isTRUE(feat_in$ok)) next
      surv_va <- subset_surv(time_tr_all, event_tr_all, val_ids)
      mf_in <- fit_metrics_fn(surv_in$time, surv_in$event, feat_in$X_tr, surv_va$time, surv_va$event, feat_in$X_te)
      fold_cx[f] <- unname(mf_in["cindex"])
      fold_au[f] <- unname(mf_in["auroc"])
      fold_detail_rows[[length(fold_detail_rows) + 1L]] <- data.frame(
        holdout_dataset = held_out, inner_fold = f,
        top_k_recurrence = kr, top_k_metastasis = km,
        inner_fold_cindex = fold_cx[f], inner_fold_auroc = fold_au[f],
        stringsAsFactors = FALSE
      )
    }
    summary_rows[[gi]] <- data.frame(
      holdout_dataset = held_out, top_k_recurrence = kr, top_k_metastasis = km,
      inner_cv_cindex = mean(fold_cx, na.rm = TRUE), inner_cv_auroc = mean(fold_au, na.rm = TRUE),
      inner_cv_n_folds_used = sum(is.finite(fold_cx)), stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  fold_df <- if (length(fold_detail_rows)) do.call(rbind, fold_detail_rows) else NULL
  list(summary = summary_df, fold_detail = fold_df, best = pick_best_k_row(summary_df, "inner_cv_cindex"))
}

# Cache everything an outer fold needs, so the locked (kr, km) can be scored later.
stash_nested_outer_fold <- function(held_out, adj_tr, adj_te, tt_rank, Xc_tr, Xc_te,
                                     time_tr_all, event_tr_all, time_te_outer, event_te_outer, y_rec = NULL) {
  list(
    held_out = held_out, adj_tr = adj_tr, adj_te = adj_te, tt_rank = tt_rank, y_rec = y_rec,
    Xc_tr = Xc_tr, Xc_te = Xc_te, time_tr_all = time_tr_all, event_tr_all = event_tr_all,
    time_te_outer = time_te_outer, event_te_outer = event_te_outer
  )
}

evaluate_outer_clin_rec_baseline <- function(adj_tr, adj_te, tt_rank, kr, Xc_tr, Xc_te,
                                             time_tr_all, event_tr_all, time_te_outer, event_te_outer, fit_metrics_fn) {
  feat <- build_clin_rec_baseline_features(adj_tr, adj_te, tt_rank, kr, Xc_tr, Xc_te)
  if (!isTRUE(feat$ok)) return(list(metrics = na_metrics(), n_recurrence_genes = 0L))
  mf <- fit_metrics_fn(time_tr_all, event_tr_all, feat$X_tr, time_te_outer, event_te_outer, feat$X_te)
  list(metrics = mf, n_recurrence_genes = feat$n_recurrence_genes)
}

evaluate_outer_nested_winner <- function(adj_tr, adj_te, tt_rank, meta_genes, kr, km, Xc_tr, Xc_te,
                                         y_tr, meta_w = NULL, time_tr_all = NULL, event_tr_all = NULL,
                                         time_te_outer = NULL, event_te_outer = NULL, fit_metrics_fn,
                                         rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  feat <- build_clin_rec_met_features(adj_tr, adj_te, tt_rank, meta_genes, kr, km, Xc_tr, Xc_te, y_tr = y_tr,
                                      time_tr = time_tr_all, event_tr = event_tr_all, meta_w = meta_w,
                                      rank_genes_fn = rank_genes_fn, fdr = fdr)
  if (!isTRUE(feat$ok)) return(list(metrics = na_metrics(), n_recurrence_genes = 0L))
  mf <- fit_metrics_fn(time_tr_all, event_tr_all, feat$X_tr, time_te_outer, event_te_outer, feat$X_te)
  list(metrics = mf, n_recurrence_genes = feat$n_recurrence_genes)
}

run_outer_nested_baseline_clin_rec <- function(fold_cache, kr, fit_metrics_fn,
                                               rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  if (!length(fold_cache)) return(NULL)
  rows <- list()
  for (nm in names(fold_cache)) {
    fc <- fold_cache[[nm]]
    ranked <- rank_genes_fn(fc$adj_tr, fc$y_rec, time_tr = fc$time_tr_all, event_tr = fc$event_tr_all, fdr = fdr)
    outer_bl <- evaluate_outer_clin_rec_baseline(
      adj_tr = fc$adj_tr, adj_te = fc$adj_te, tt_rank = ranked, kr = kr,
      Xc_tr = fc$Xc_tr, Xc_te = fc$Xc_te, time_tr_all = fc$time_tr_all, event_tr_all = fc$event_tr_all,
      time_te_outer = fc$time_te_outer, event_te_outer = fc$event_te_outer, fit_metrics_fn = fit_metrics_fn
    )
    rows[[length(rows) + 1L]] <- data.frame(
      holdout_dataset = fc$held_out, top_k_recurrence = kr, top_k_metastasis = NA_integer_,
      inner_cv_cindex_at_winner = NA_real_, inner_cv_auroc_at_winner = NA_real_,
      outer_cindex_nested = unname(outer_bl$metrics["cindex"]),
      outer_auroc_nested = unname(outer_bl$metrics["auroc"]),
      n_recurrence_genes = outer_bl$n_recurrence_genes,
      model_id = sprintf("Baseline: Clinical+Recurrence_top%d", kr),
      is_baseline = TRUE, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

lookup_inner_cv_at_k <- function(inner_all, holdout, kr, km) {
  if (is.null(inner_all) || !nrow(inner_all)) return(c(NA_real_, NA_real_))
  hit <- inner_all[
    inner_all$holdout_dataset == holdout &
      inner_all$top_k_recurrence == kr & inner_all$top_k_metastasis == km,
    , drop = FALSE
  ]
  if (!nrow(hit)) return(c(NA_real_, NA_real_))
  c(hit$inner_cv_cindex[1], hit$inner_cv_auroc[1])
}

run_outer_nested_global_k <- function(fold_cache, kr, km, meta_genes, meta_w = NULL, fit_metrics_fn,
                                      inner_all = NULL, rank_genes_fn = rank_recurrence_genes, fdr = 0.05) {
  if (!length(fold_cache)) return(NULL)
  rows <- list()
  for (nm in names(fold_cache)) {
    fc <- fold_cache[[nm]]
    ranked <- rank_genes_fn(fc$adj_tr, fc$y_rec, time_tr = fc$time_tr_all, event_tr = fc$event_tr_all, fdr = fdr)
    outer <- evaluate_outer_nested_winner(
      adj_tr = fc$adj_tr, adj_te = fc$adj_te, tt_rank = ranked, meta_genes = meta_genes, meta_w = meta_w,
      kr = kr, km = km, Xc_tr = fc$Xc_tr, Xc_te = fc$Xc_te, y_tr = fc$y_rec,
      time_tr_all = fc$time_tr_all, event_tr_all = fc$event_tr_all,
      time_te_outer = fc$time_te_outer, event_te_outer = fc$event_te_outer,
      fit_metrics_fn = fit_metrics_fn, rank_genes_fn = rank_genes_fn, fdr = fdr
    )
    inner_at_k <- lookup_inner_cv_at_k(inner_all, fc$held_out, kr, km)
    rows[[length(rows) + 1L]] <- data.frame(
      holdout_dataset = fc$held_out, top_k_recurrence = kr, top_k_metastasis = km,
      inner_cv_cindex_at_winner = inner_at_k[1], inner_cv_auroc_at_winner = inner_at_k[2],
      outer_cindex_nested = unname(outer$metrics["cindex"]),
      outer_auroc_nested = unname(outer$metrics["auroc"]),
      n_recurrence_genes = outer$n_recurrence_genes,
      model_id = sprintf("Clinical+Recurrence_top%d+Metastasis_top%d (global grid K)", kr, km),
      is_baseline = FALSE, stringsAsFactors = FALSE
    )
    message(sprintf("  Outer LODO [%s]: kr=%d km=%d outer C=%.3f", fc$held_out, kr, km, unname(outer$metrics["cindex"])))
  }
  do.call(rbind, rows)
}

# Average inner CV scores across the outer folds and rank kr x km cells.
summarize_inner_cv_across_holdouts <- function(inner_by_holdout) {
  if (is.null(inner_by_holdout)) return(NULL)
  inner_by_holdout <- as.data.frame(inner_by_holdout, stringsAsFactors = FALSE)
  if (!nrow(inner_by_holdout)) return(NULL)
  req <- c("top_k_recurrence", "top_k_metastasis", "inner_cv_cindex", "inner_cv_auroc")
  inner_by_holdout <- inner_by_holdout[
    is.finite(inner_by_holdout$top_k_recurrence) & is.finite(inner_by_holdout$top_k_metastasis),
    req, drop = FALSE
  ]
  if (!nrow(inner_by_holdout)) return(NULL)
  key <- paste(inner_by_holdout$top_k_recurrence, inner_by_holdout$top_k_metastasis, sep = "|")
  agg <- do.call(rbind, lapply(split(seq_len(nrow(inner_by_holdout)), key), function(ii) {
    sub <- inner_by_holdout[ii, , drop = FALSE]
    data.frame(
      top_k_recurrence = sub$top_k_recurrence[1L], top_k_metastasis = sub$top_k_metastasis[1L],
      mean_inner_cv_cindex = mean(sub$inner_cv_cindex, na.rm = TRUE),
      mean_inner_cv_auroc = mean(sub$inner_cv_auroc, na.rm = TRUE),
      sd_inner_cv_cindex = if (sum(is.finite(sub$inner_cv_cindex)) < 2L) NA_real_ else stats::sd(sub$inner_cv_cindex, na.rm = TRUE),
      n_outer_folds = sum(is.finite(sub$inner_cv_cindex)), stringsAsFactors = FALSE
    )
  }))
  if (is.null(agg) || !nrow(agg)) return(NULL)
  agg <- agg[order(-agg$mean_inner_cv_cindex, agg$top_k_recurrence + agg$top_k_metastasis,
                   agg$top_k_recurrence, agg$top_k_metastasis), , drop = FALSE]
  agg$rank_by_mean_inner_cv <- seq_len(nrow(agg))
  agg
}

# Lock the final (kr, km) from inner CV only (never outer/external — no leakage).
recommend_final_k <- function(inner_averaged, metric_col = "mean_inner_cv_cindex") {
  if (is.null(inner_averaged) || !nrow(inner_averaged)) return(NULL)
  best <- pick_best_k_row(inner_averaged, metric_col)
  if (is.null(best) || !nrow(best)) return(NULL)
  data.frame(
    recommended_top_k_recurrence = best$top_k_recurrence,
    recommended_top_k_metastasis = best$top_k_metastasis,
    mean_inner_cv_cindex = best[[metric_col]],
    mean_inner_cv_auroc = if ("mean_inner_cv_auroc" %in% names(best)) best$mean_inner_cv_auroc else NA_real_,
    selection_rule = "Inner CV grid kr,km only; highest mean inner-CV C-index across outer folds. Genes refit on all development with locked K.",
    stringsAsFactors = FALSE
  )
}

# Refit recurrence + metastasis weights on all development data at the locked K.
fit_final_locked_on_development <- function(expr_gene, pheno, meta_genes, meta_w = NULL, kr, km,
                                            clinical_cols, fdr = 0.05, rank_genes_fn = NULL) {
  if (is.null(meta_genes) && !is.null(meta_w)) meta_genes <- meta_w
  kr <- as.integer(kr)
  km <- as.integer(km)
  if (is.null(rank_genes_fn)) {
    rank_genes_fn <- function(mat, y, time_tr = NULL, event_tr = NULL, fdr = 0.05) ranked_limma_recurrence(mat, y, fdr = fdr)
  }
  prep <- combat_train_only_preprocess(expr_gene, expr_gene[, 1, drop = FALSE], as.factor(as.character(pheno$Dataset)))
  adj_dev <- prep$train
  ranked <- rank_genes_fn(adj_dev, pheno$Recurrence, fdr = fdr)
  top_rec <- select_recurrence_topk(ranked, kr)

  time_dev <- if (exists("surv_time_days_train", mode = "function")) surv_time_days_train(pheno) else NULL
  event_dev <- if (exists("surv_event", mode = "function")) surv_event(pheno) else NULL

  imp <- impute_clinical(
    as.data.frame(pheno[, clinical_cols, drop = FALSE]),
    as.data.frame(pheno[, clinical_cols, drop = FALSE])
  )
  feat <- build_clin_rec_met_features(adj_dev, adj_dev, ranked, meta_genes, kr, km, imp$train, imp$test,
                                      y_tr = pheno$Recurrence, time_tr = time_dev, event_tr = event_dev,
                                      meta_w = meta_w, rank_genes_fn = rank_genes_fn, fdr = fdr)
  w_meta <- refit_metastasis_weights(adj_dev, meta_genes[[as.character(km)]], pheno$Recurrence,
                                     time_tr = time_dev, event_tr = event_dev, rank_genes_fn = rank_genes_fn, fdr = fdr)
  rec_w <- setNames(top_rec$logFC, top_rec$gene)
  list(
    kr = kr, km = km, adj_development = adj_dev, train_mu = prep$train_mu, train_sd = prep$train_sd,
    recurrence_genes = top_rec$gene, recur_weights = rec_w, a_genes = names(w_meta), meta_weights = w_meta,
    X_train = feat$X_tr, y_train = as.integer(pheno$Recurrence), n_recurrence_genes = nrow(top_rec),
    selection_note = "kr/km from inner CV only; recurrence and metastasis weights refit on all development after K locked; metastasis genes from Block A list only"
  )
}

write_final_locked_development_table <- function(fit_obj, out_dir) {
  if (is.null(fit_obj)) return(invisible(NULL))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rec_df <- data.frame(gene = fit_obj$recurrence_genes,
                       logFC_weight = as.numeric(fit_obj$recur_weights[fit_obj$recurrence_genes]),
                       stringsAsFactors = FALSE)
  meta_df <- data.frame(gene = fit_obj$a_genes,
                        logFC_weight = as.numeric(fit_obj$meta_weights[fit_obj$a_genes]),
                        stringsAsFactors = FALSE)
  write.table(
    data.frame(
      top_k_recurrence = fit_obj$kr, top_k_metastasis = fit_obj$km,
      n_recurrence_genes = fit_obj$n_recurrence_genes, n_metastasis_genes = length(fit_obj$a_genes),
      selection_note = fit_obj$selection_note, stringsAsFactors = FALSE
    ),
    file.path(out_dir, "combo_final_locked_k_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(rec_df, file.path(out_dir, "combo_final_locked_recurrence_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(meta_df, file.path(out_dir, "combo_final_locked_metastasis_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  rec_key <- toupper(trimws(rec_df$gene))
  meta_key <- toupper(trimws(meta_df$gene))
  overlap_genes <- intersect(rec_key, meta_key)
  overlap_df <- data.frame(
    gene = overlap_genes,
    recurrence_logFC_weight = rec_df$logFC_weight[match(overlap_genes, rec_key)],
    metastasis_logFC_weight = meta_df$logFC_weight[match(overlap_genes, meta_key)],
    stringsAsFactors = FALSE
  )
  write.table(overlap_df, file.path(out_dir, "combo_final_locked_overlap_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(list(recurrence = rec_df, metastasis = meta_df, overlap = overlap_df))
}

read_nested_final_k <- function(recommendation_file, averaged_file = NULL) {
  if (file.exists(recommendation_file)) {
    rec <- read.delim(recommendation_file, stringsAsFactors = FALSE)
    if (nrow(rec) >= 1L) {
      kr <- if ("recommended_top_k_recurrence" %in% names(rec)) rec$recommended_top_k_recurrence[1]
            else if ("kr" %in% names(rec)) rec$kr[1]
            else if ("top_k_recurrence" %in% names(rec)) rec$top_k_recurrence[1] else NA_integer_
      km <- if ("recommended_top_k_metastasis" %in% names(rec)) rec$recommended_top_k_metastasis[1]
            else if ("km" %in% names(rec)) rec$km[1]
            else if ("top_k_metastasis" %in% names(rec)) rec$top_k_metastasis[1] else NA_integer_
      score_col <- intersect(c("mean_inner_cv_cindex", "mean_cv_cindex", "cindex_mean", "selection_score"), names(rec))[1]
      if (is.finite(kr) && is.finite(km)) {
        return(list(
          kr = as.integer(kr), km = as.integer(km),
          score_col = if (!is.na(score_col)) score_col else "mean_inner_cv_cindex",
          selected_score = if (!is.na(score_col)) rec[[score_col]][1] else NA_real_,
          table = rec
        ))
      }
    }
  }
  avg_path <- if (is.null(averaged_file)) file.path(dirname(recommendation_file), "combo_inner_cv_k_averaged.tsv") else averaged_file
  rec <- read.delim(avg_path, stringsAsFactors = FALSE)
  best <- pick_best_k_row(rec, "mean_inner_cv_cindex")
  list(kr = as.integer(best$top_k_recurrence[1]), km = as.integer(best$top_k_metastasis[1]),
       score_col = "mean_inner_cv_cindex", selected_score = best$mean_inner_cv_cindex[1], table = rec)
}

# Refit at locked K on all development data, then evaluate once on the external cohort.
run_locked_external_clin_rec_met <- function(engine_label, out_dir, pheno, tc, expr_gene, te_gene,
                                             meta_genes, meta_w = NULL, clinical_cols, fdr,
                                             fit_metrics_fn, rank_genes_fn, time_all_fun, time_ext_fun, event_fun,
                                             recommendation_file = file.path(out_dir, "combo_final_k_recommendation.tsv")) {
  fk <- read_nested_final_k(recommendation_file)
  kr <- fk$kr
  km <- fk$km
  message("Locked external [", engine_label, "]: kr=", kr, ", km=", km)

  final_fit <- fit_final_locked_on_development(expr_gene = expr_gene, pheno = pheno, meta_genes = meta_genes,
                                               meta_w = meta_w, kr = kr, km = km, clinical_cols = clinical_cols,
                                               fdr = fdr, rank_genes_fn = rank_genes_fn)
  write_final_locked_development_table(final_fit, out_dir)

  cb_ext <- combat_train_test(expr_gene, te_gene, as.factor(as.character(pheno$Dataset)), "External_GSE68465")
  adj_ext <- cb_ext$right
  imp_ext <- impute_clinical(
    as.data.frame(pheno[, clinical_cols, drop = FALSE]),
    as.data.frame(tc[, clinical_cols, drop = FALSE])
  )
  X_te <- imp_ext$test
  X_te$Metastasis_score <- as.numeric(signature_score(adj_ext, names(final_fit$meta_weights), final_fit$meta_weights))
  X_te$Recurrence_score <- as.numeric(signature_score(adj_ext, names(final_fit$recur_weights), final_fit$recur_weights))

  mf <- fit_metrics_fn(time_all_fun(pheno), event_fun(pheno), final_fit$X_train,
                       time_ext_fun(tc), event_fun(tc), X_te)

  # Baseline (clinical + recurrence score only): drop the metastasis score and
  # re-evaluate on the same external cohort, using the same locked recurrence
  # genes/weights, so the added value of the metastasis score can be compared
  # on external data as well as on LODO.
  drop_metastasis_score <- function(X) X[, setdiff(names(X), "Metastasis_score"), drop = FALSE]
  mf_bl <- tryCatch(
    fit_metrics_fn(time_all_fun(pheno), event_fun(pheno), drop_metastasis_score(final_fit$X_train),
                   time_ext_fun(tc), event_fun(tc), drop_metastasis_score(X_te)),
    error = function(e) c(cindex = NA_real_, auroc = NA_real_)
  )

  out <- data.frame(
    model_engine = engine_label,
    model_id = sprintf("Clinical+Recurrence_top%d+Metastasis_top%d (official locked K)", kr, km),
    top_k_recurrence = kr, top_k_metastasis = km,
    uses_clinical = TRUE, uses_recurrence = TRUE, uses_metastasis = TRUE,
    n_recurrence_genes = final_fit$n_recurrence_genes,
    cindex_external = unname(mf["cindex"]), auroc_external = unname(mf["auroc"]),
    cindex_external_baseline = unname(mf_bl["cindex"]),
    auroc_external_baseline = unname(mf_bl["auroc"]),
    delta_cindex_external_full_minus_baseline = unname(mf["cindex"]) - unname(mf_bl["cindex"]),
    n_training_samples = nrow(pheno), n_external_samples = nrow(tc), external_cohort = "GSE68465",
    selection_source = basename(recommendation_file), selection_metric = fk$score_col,
    selection_score = fk$selected_score, stringsAsFactors = FALSE
  )

  nested_file <- file.path(out_dir, "combo_lodocv_nested_tuned_by_holdout.tsv")
  if (file.exists(nested_file)) {
    lodo <- lodo_means_from_holdout(read.delim(nested_file, stringsAsFactors = FALSE), kr, km)
    out$mean_cindex_lodo_nested <- lodo$full
    out$mean_auroc_lodo_nested <- lodo$auroc_full
    out$n_lodo_folds <- lodo$n_full_folds
  } else {
    out$mean_cindex_lodo_nested <- NA_real_
    out$mean_auroc_lodo_nested <- NA_real_
    out$n_lodo_folds <- NA_integer_
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(out, file.path(out_dir, "combo_external_nested_final_locked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(out, file.path(out_dir, "combo_official_nested_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out)
}

rbind_nested_tuned_tables <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, function(x) is.null(x) || !nrow(x), logical(1))]
  if (!length(parts)) return(data.frame())
  all_cols <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, parts)
}

write_nested_k_outputs <- function(out_dir, inner_by_holdout, fold_detail, nested_tuned, inner_averaged, final_k) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(d, name) {
    if (!is.null(d) && nrow(d)) write.table(d, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write_tsv(inner_by_holdout, "combo_inner_cv_k_by_holdout.tsv")
  write_tsv(fold_detail, "combo_inner_cv_k_folds.tsv")
  write_tsv(nested_tuned, "combo_lodocv_nested_tuned_by_holdout.tsv")
  write_tsv(inner_averaged, "combo_inner_cv_k_averaged.tsv")
  write_tsv(final_k, "combo_final_k_recommendation.tsv")
}
