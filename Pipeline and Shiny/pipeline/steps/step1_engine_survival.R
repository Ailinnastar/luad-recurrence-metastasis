## Step 1 (survival engines): nested leave-one-dataset-out tuning of (kr, km),
## refit at locked K on all development data, then one external evaluation.
## Sourced once per survival engine by step1_run_all_engines.R.

survival_model <- Sys.getenv("COMBO_SURVIVAL_MODEL", "cox_glmnet")
cox_lambda <- tolower(Sys.getenv("COMBO_COX_LAMBDA", "min"))
cox_alpha <- as.numeric(Sys.getenv("COMBO_COX_ALPHA", "1"))
cox_preset <- trimws(Sys.getenv("COMBO_COX_PRESET", ""))
tune_alpha <- as.logical(Sys.getenv("COMBO_TUNE_ALPHA", "TRUE"))
xgb_nrounds <- as.integer(Sys.getenv("COMBO_XGB_NROUNDS", "150"))
xgb_nfold <- as.integer(Sys.getenv("COMBO_XGB_NFOLD", "3"))
alpha_inner_folds <- as.integer(Sys.getenv("COMBO_ALPHA_INNER_FOLDS", "5"))
alpha_glmnet_nfolds <- as.integer(Sys.getenv("COMBO_ALPHA_GLMNET_NFOLDS", "3"))
cox_nfolds <- as.integer(Sys.getenv("COMBO_COX_NFOLDS", "5"))
random_seed <- as.integer(Sys.getenv("COMBO_RANDOM_SEED", "1"))

alpha_grid <- if (isTRUE(tune_alpha)) {
  if (exists("COMBO_ALPHA_GRID_OVERRIDE", inherits = TRUE)) {
    sort(unique(COMBO_ALPHA_GRID_OVERRIDE))
  } else if (nzchar(cox_preset)) {
    combo_alpha_grid_for_preset(cox_preset)
  } else {
    sort(unique(as.numeric(trimws(strsplit(Sys.getenv("COMBO_ALPHA_GRID", "0,0.25,0.5,0.75,1"), ",")[[1]]))))
  }
} else {
  numeric(0)
}

set.seed(random_seed)
suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(survival)
  library(pROC)
  if (survival_model == "cox_glmnet") library(glmnet)
  if (survival_model == "rsf") library(randomForestSRC)
  if (survival_model == "xgboost") library(xgboost)
})

ROOT <- project_root()
train_dir <- training_dir(ROOT)
ext_dir <- external_dir(ROOT)
frozen_dir <- frozen_metastasis_dir(ROOT)

subdir <- switch(survival_model,
  cox_glmnet = if (isTRUE(tune_alpha) && nzchar(cox_preset)) paste0(cox_preset, "_cox_alpha_tuned") else paste0(cox_preset, "_cox"),
  rsf = "random_survival_forest",
  xgboost = "xgboost_survival_cox"
)
out_dir <- file.path(ROOT, results_root(), subdir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_K <- c(5L, 10L, 15L)
FDR_LIMMA <- 0.05
EXT_FU_MONTHS_COL <- "months_to_last_contact_or_death:ch1"
k_grid <- combo_k_grid_values()
inner_cv_folds <- combo_inner_cv_folds()

engine_label <- switch(survival_model,
  cox_glmnet = paste0(
    switch(cox_preset, ridge = "Ridge", elastic_net = "Elastic-net", penalized = "Penalized", lasso = "Lasso", cox_preset),
    " Cox PH (glmnet); inner CV alpha in {", paste(sprintf("%.2f", alpha_grid), collapse = ","),
    "}, lambda.", cox_lambda, "; Harrell C-index + AUROC; seed=", random_seed
  ),
  rsf = paste0("Random survival forest; Harrell C-index + AUROC; seed=", random_seed),
  xgboost = paste0("XGBoost survival Cox; Harrell C-index + AUROC; seed=", random_seed)
)
message("Survival engine: ", survival_model, " -> ", out_dir)

# ---- Clinical helpers (median impute, follow-up time, event) ----

impute_clinical <- function(df_tr, df_te) {
  out_tr <- data.frame(df_tr, check.names = FALSE, stringsAsFactors = FALSE)
  out_te <- data.frame(df_te, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(out_tr)) {
    med <- suppressWarnings(median(out_tr[[nm]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    out_tr[[nm]][is.na(out_tr[[nm]])] <- med
    out_te[[nm]][is.na(out_te[[nm]])] <- med
  }
  list(train = out_tr, test = out_te)
}

surv_time_days_train <- function(df) as.numeric(df$FollowUp_days)
surv_time_days_external <- function(df) as.numeric(df[[EXT_FU_MONTHS_COL]]) * (365.25 / 12)
surv_event <- function(df) as.integer(df$Recurrence)

sanitize_X <- function(X) {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  names(X) <- make.names(names(X), unique = TRUE)
  X
}

cox_lambda_rule <- function() if (cox_lambda == "1se") "lambda.1se" else "lambda.min"

# ---- Model fitting (shared train/test prep + per-engine fit) ----

prep_cox_matrices <- function(time_tr, event_tr, X_tr, time_te, event_te, X_te) {
  X_tr <- sanitize_X(X_tr)
  X_te <- sanitize_X(X_te)
  shared <- intersect(names(X_tr), names(X_te))
  if (length(shared) == 0L) return(NULL)
  xtr <- as.matrix(X_tr[, shared, drop = FALSE])
  xte <- as.matrix(X_te[, shared, drop = FALSE])
  impute_mat <- function(M) {
    if (!ncol(M)) return(M)
    for (j in seq_len(ncol(M))) {
      v <- M[, j]
      med <- suppressWarnings(median(v[is.finite(v)], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      v[!is.finite(v)] <- med
      M[, j] <- v
    }
    M
  }
  xtr <- impute_mat(xtr)
  xte <- impute_mat(xte)
  keep <- {
    sdv <- apply(xtr, 2L, function(z) stats::sd(z, na.rm = TRUE))
    is.finite(sdv) & sdv > 1e-10
  }
  if (!any(keep)) return(NULL)
  xtr <- xtr[, keep, drop = FALSE]
  xte <- xte[, keep, drop = FALSE]
  ok_tr <- is.finite(time_tr) & time_tr > 0 & event_tr %in% c(0L, 1L)
  ok_te <- is.finite(time_te) & time_te > 0 & event_te %in% c(0L, 1L)
  if (sum(ok_tr) < 25L || sum(ok_te) < 5L) return(NULL)
  if (sum(event_tr[ok_tr] == 1L) < 4L || sum(event_te[ok_te] == 1L) < 2L) return(NULL)
  list(xtr = xtr, xte = xte, ok_tr = ok_tr, ok_te = ok_te,
       time_tr = time_tr, event_tr = event_tr, time_te = time_te, event_te = event_te)
}

# Inner CV over the alpha grid; returns the alpha with the best mean fold C-index.
select_best_alpha_cox <- function(prep, alpha_grid, inner_folds = alpha_inner_folds,
                                   glmnet_nfolds = alpha_glmnet_nfolds, seed = random_seed,
                                   lam_rule = cox_lambda_rule()) {
  idx <- which(prep$ok_tr)
  n <- length(idx)
  if (n < 30L || length(alpha_grid) < 1L) return(alpha_grid[1])
  t_ok <- prep$time_tr[idx]
  e_ok <- prep$event_tr[idx]
  x_ok <- prep$xtr[idx, , drop = FALSE]
  if (sum(e_ok == 1L) < 4L) return(alpha_grid[1])
  set.seed(seed + n + length(alpha_grid))
  fold_id <- sample(rep(seq_len(inner_folds), length.out = n))
  score_alpha <- function(alp) {
    fold_scores <- vapply(seq_len(inner_folds), function(f) {
      tr_i <- which(fold_id != f)
      te_i <- which(fold_id == f)
      if (length(te_i) < 5L || length(unique(e_ok[tr_i])) < 2L) return(0.5)
      tryCatch({
        fit_f <- glmnet::cv.glmnet(x_ok[tr_i, , drop = FALSE], survival::Surv(t_ok[tr_i], e_ok[tr_i]),
                                   family = "cox", alpha = alp, nfolds = glmnet_nfolds, grouped = FALSE)
        lp <- as.numeric(predict(fit_f, newx = x_ok[te_i, , drop = FALSE], s = lam_rule, type = "link"))
        cx <- cindex_harrell(t_ok[te_i], e_ok[te_i], lp)
        if (!is.finite(cx)) 0.5 else cx
      }, error = function(e) 0.5)
    }, numeric(1))
    mean(fold_scores, na.rm = TRUE)
  }
  alpha_grid[which.max(vapply(alpha_grid, score_alpha, numeric(1)))]
}

fit_predict_cox_glmnet <- function(prep, alpha, lam_rule = cox_lambda_rule(), glmnet_nfolds = cox_nfolds) {
  cv <- tryCatch(
    glmnet::cv.glmnet(prep$xtr[prep$ok_tr, , drop = FALSE],
                      survival::Surv(prep$time_tr[prep$ok_tr], prep$event_tr[prep$ok_tr]),
                      family = "cox", alpha = alpha, nfolds = glmnet_nfolds, grouped = FALSE),
    error = function(e) NULL
  )
  if (is.null(cv)) return(na_metrics())
  lp <- tryCatch(as.numeric(predict(cv, newx = prep$xte[prep$ok_te, , drop = FALSE], s = lam_rule, type = "link")),
                 error = function(e) NULL)
  if (is.null(lp)) return(na_metrics())
  c(cindex = cindex_harrell(prep$time_te[prep$ok_te], prep$event_te[prep$ok_te], lp),
    auroc = auc_binary(prep$event_te[prep$ok_te], lp))
}

fit_predict_xgboost <- function(prep, seed = random_seed) {
  ok_tr <- prep$ok_tr
  ok_te <- prep$ok_te
  if (sum(ok_tr) < 10L || sum(ok_te) < 3L) return(na_metrics())
  xtr <- prep$xtr[ok_tr, , drop = FALSE]
  xte <- prep$xte[ok_te, , drop = FALSE]
  time_tr <- prep$time_tr[ok_tr]
  event_tr <- prep$event_tr[ok_tr]
  set.seed(seed)
  dtrain <- xgboost::xgb.DMatrix(data = xtr, label = time_tr, weight = event_tr)
  dte <- xgboost::xgb.DMatrix(data = xte)
  params <- list(objective = "survival:cox", eval_metric = "cox-nloglik", max_depth = 3L,
                 eta = 0.05, subsample = 0.8, colsample_bytree = 0.8)
  cv_xgb <- tryCatch(
    xgboost::xgb.cv(params = params, data = dtrain, nrounds = xgb_nrounds, nfold = xgb_nfold,
                    early_stopping_rounds = 15L, verbose = 0L, seed = seed),
    error = function(e) NULL
  )
  best_n <- if (!is.null(cv_xgb)) cv_xgb$best_iteration else 80L
  if (is.null(best_n) || length(best_n) == 0L || !is.finite(best_n) || best_n < 1L) best_n <- 80L
  set.seed(seed)
  fit <- tryCatch(xgboost::xgb.train(params = params, data = dtrain, nrounds = as.integer(best_n), verbose = 0L),
                  error = function(e) NULL)
  if (is.null(fit)) return(na_metrics())
  risk <- tryCatch(as.numeric(predict(fit, dte)), error = function(e) NULL)
  if (is.null(risk) || !length(risk)) return(na_metrics())
  c(cindex = cindex_harrell(prep$time_te[ok_te], prep$event_te[ok_te], risk),
    auroc = auc_binary(prep$event_te[ok_te], risk))
}

fit_predict_rsf <- function(prep) {
  dat_tr <- data.frame(time = prep$time_tr[prep$ok_tr], status = prep$event_tr[prep$ok_tr],
                       prep$xtr[prep$ok_tr, , drop = FALSE], check.names = FALSE)
  dat_te <- data.frame(time = prep$time_te[prep$ok_te], status = prep$event_te[prep$ok_te],
                       prep$xte[prep$ok_te, , drop = FALSE], check.names = FALSE)
  if (ncol(dat_tr) != ncol(dat_te)) return(na_metrics())
  names(dat_te) <- names(dat_tr)
  fit <- tryCatch(
    randomForestSRC::rfsrc(Surv(time, status) ~ ., data = dat_tr, ntree = 80L, seed = random_seed,
                           nodesize = max(5L, min(20L, floor(nrow(dat_tr) / 15))),
                           splitrule = "logrank", forest = TRUE, importance = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(na_metrics())
  pr <- tryCatch(predict(fit, newdata = dat_te), error = function(e) NULL)
  if (is.null(pr)) return(na_metrics())
  risk <- if (!is.null(pr$predicted) && length(pr$predicted) == nrow(dat_te)) as.numeric(pr$predicted)
          else if (!is.null(pr$chf)) rowSums(pr$chf, na.rm = TRUE) else return(na_metrics())
  c(cindex = cindex_harrell(prep$time_te[prep$ok_te], prep$event_te[prep$ok_te], risk),
    auroc = auc_binary(prep$event_te[prep$ok_te], risk))
}

fit_predict_survival_metrics <- function(time_tr, event_tr, X_tr, time_te, event_te, X_te) {
  prep <- prep_cox_matrices(time_tr, event_tr, X_tr, time_te, event_te, X_te)
  if (is.null(prep)) return(na_metrics())
  if (survival_model == "cox_glmnet") {
    alpha_use <- if (isTRUE(tune_alpha)) select_best_alpha_cox(prep, alpha_grid) else cox_alpha
    return(fit_predict_cox_glmnet(prep, alpha = alpha_use))
  }
  if (survival_model == "xgboost") return(fit_predict_xgboost(prep, seed = random_seed))
  fit_predict_rsf(prep)
}

# ---- Load training (4 cohorts) + external (GSE68465), collapse probes to genes ----

e_probe <- readRDS(file.path(train_dir, "e_merged_raw.rds"))
annot <- readRDS(file.path(train_dir, "f_merged_clean.rds"))
pheno <- readRDS(file.path(train_dir, "p_merged_clean.rds"))
annot$ID <- as.character(annot$ID)
gsym <- toupper(trimws(as.character(annot[["Gene Symbol"]])))
gsym[gsym == "" | is.na(gsym)] <- NA_character_
ok <- !is.na(gsym) & rownames(e_probe) %in% annot$ID
em <- as.matrix(e_probe[ok, , drop = FALSE])
rownames(em) <- NULL
gvec <- gsym[match(rownames(e_probe)[ok], annot$ID)]
expr_gene <- as.matrix(do.call(rbind, lapply(split(seq_len(nrow(em)), gvec), function(ix) {
  if (length(ix) == 1L) matrix(em[ix, ], nrow = 1L) else matrix(colMeans(em[ix, , drop = FALSE]), nrow = 1L)
})))
rownames(expr_gene) <- names(split(seq_len(nrow(em)), gvec))
colnames(expr_gene) <- colnames(em)

te_raw <- readRDS(file.path(ext_dir, "test_GSE68465_expression.rds"))
rownames(te_raw) <- toupper(rownames(te_raw))
tc <- readRDS(file.path(ext_dir, "test_GSE68465_clinical.rds"))
common_samples <- intersect(colnames(te_raw), rownames(tc))
if (length(common_samples) != ncol(te_raw)) {
  te_raw <- te_raw[, common_samples, drop = FALSE]
  tc <- tc[common_samples, , drop = FALSE]
}

genes_xy <- sort(intersect(rownames(expr_gene), rownames(te_raw)))
expr_gene <- expr_gene[genes_xy, , drop = FALSE]
te_gene <- te_raw[genes_xy, , drop = FALSE]

clinical_cols <- c("Stage_num", "Age", "Gender_bin")
meta_genes <- load_frozen_metastasis_gene_lists(frozen_dir)

rank_genes_for_inner <- function(mat, y, time_tr, event_tr, fdr) {
  rank_recurrence_genes(mat, y, time_tr = time_tr, event_tr = event_tr, fdr = fdr)
}

# ---- Outer LODO loop: inner CV grid + cache each fold ----

cohort_ids <- sort(unique(as.character(pheno$Dataset)))
inner_k_by_holdout <- list()
inner_k_fold_detail <- list()
nested_fold_cache <- list()

for (held_out in cohort_ids) {
  tr <- rownames(pheno)[pheno$Dataset != held_out]
  te <- rownames(pheno)[pheno$Dataset == held_out]
  if (length(te) < 5L) next
  ptr <- pheno[tr, , drop = FALSE]
  pte <- pheno[te, , drop = FALSE]

  cb <- combat_train_test(expr_gene[, tr, drop = FALSE], expr_gene[, te, drop = FALSE], ptr$Dataset, held_out)
  adj_tr <- cb$left
  adj_te <- cb$right

  ranked <- rank_recurrence_genes(adj_tr, ptr$Recurrence,
                                  time_tr = surv_time_days_train(ptr), event_tr = surv_event(ptr), fdr = FDR_LIMMA)
  imp <- impute_clinical(as.data.frame(ptr[, clinical_cols, drop = FALSE]),
                         as.data.frame(pte[, clinical_cols, drop = FALSE]))

  time_tr <- setNames(surv_time_days_train(ptr), rownames(ptr))
  event_tr <- setNames(surv_event(ptr), rownames(ptr))
  time_te <- setNames(surv_time_days_train(pte), rownames(pte))
  event_te <- setNames(surv_event(pte), rownames(pte))

  hold_i <- match(held_out, cohort_ids)
  inner_res <- inner_cv_k_grid_one_holdout(
    expr_gene = expr_gene[, tr, drop = FALSE], ptr = ptr, pte = pte, held_out = held_out,
    meta_genes = meta_genes, kr_vals = k_grid$recurrence, km_vals = k_grid$metastasis,
    clinical_cols = clinical_cols, n_inner = inner_cv_folds, seed = random_seed + hold_i,
    time_tr_all = time_tr, event_tr_all = event_tr, time_te_outer = time_te, event_te_outer = event_te,
    fit_metrics_fn = fit_predict_survival_metrics, rank_genes_fn = rank_genes_for_inner, fdr = FDR_LIMMA
  )
  inner_k_by_holdout[[length(inner_k_by_holdout) + 1L]] <- inner_res$summary
  if (!is.null(inner_res$fold_detail)) inner_k_fold_detail[[length(inner_k_fold_detail) + 1L]] <- inner_res$fold_detail

  nested_fold_cache[[held_out]] <- stash_nested_outer_fold(
    held_out = held_out, adj_tr = adj_tr, adj_te = adj_te, tt_rank = ranked,
    Xc_tr = imp$train, Xc_te = imp$test, time_tr_all = time_tr, event_tr_all = event_tr,
    time_te_outer = time_te, event_te_outer = event_te, y_rec = ptr$Recurrence
  )
  message(sprintf("  Inner CV grid [%s]: %d kr x km cells", held_out, nrow(inner_res$summary)))
}

# ---- Lock K from inner CV, score it across outer folds, refit on development ----

inner_all <- do.call(rbind, inner_k_by_holdout)
fold_all <- if (length(inner_k_fold_detail)) do.call(rbind, inner_k_fold_detail) else NULL
inner_avg <- summarize_inner_cv_across_holdouts(inner_all)
final_k <- recommend_final_k(inner_avg)

nested_all <- data.frame()
if (!is.null(final_k)) {
  kr_g <- as.integer(final_k$recommended_top_k_recurrence)
  km_g <- as.integer(final_k$recommended_top_k_metastasis)
  full_nested <- run_outer_nested_global_k(nested_fold_cache, kr = kr_g, km = km_g, meta_genes = meta_genes,
                                           fit_metrics_fn = fit_predict_survival_metrics, inner_all = inner_all,
                                           rank_genes_fn = rank_genes_for_inner, fdr = FDR_LIMMA)
  baseline_nested <- run_outer_nested_baseline_clin_rec(nested_fold_cache, kr = kr_g,
                                                        fit_metrics_fn = fit_predict_survival_metrics,
                                                        rank_genes_fn = rank_genes_for_inner, fdr = FDR_LIMMA)
  nested_all <- rbind_nested_tuned_tables(full_nested, baseline_nested)
}

write_nested_k_outputs(out_dir, inner_all, fold_all, nested_all, inner_avg, final_k)

if (!is.null(final_k)) {
  message("Recommended final K: kr=", final_k$recommended_top_k_recurrence,
          " km=", final_k$recommended_top_k_metastasis,
          " (mean inner CV C-index=", sprintf("%.3f", final_k$mean_inner_cv_cindex), ")")
  final_fit_dev <- fit_final_locked_on_development(
    expr_gene = expr_gene, pheno = pheno, meta_genes = meta_genes,
    kr = as.integer(final_k$recommended_top_k_recurrence), km = as.integer(final_k$recommended_top_k_metastasis),
    clinical_cols = clinical_cols, fdr = FDR_LIMMA, rank_genes_fn = rank_genes_for_inner
  )
  write_final_locked_development_table(final_fit_dev, out_dir)
}

# ---- One external evaluation at the locked K ----

run_locked_external_clin_rec_met(
  engine_label = engine_label, out_dir = out_dir, pheno = pheno, tc = tc,
  expr_gene = expr_gene, te_gene = te_gene, meta_genes = meta_genes, clinical_cols = clinical_cols,
  fdr = FDR_LIMMA, fit_metrics_fn = fit_predict_survival_metrics, rank_genes_fn = rank_genes_for_inner,
  time_all_fun = surv_time_days_train, time_ext_fun = surv_time_days_external, event_fun = surv_event
)

message("Done: ", out_dir)
