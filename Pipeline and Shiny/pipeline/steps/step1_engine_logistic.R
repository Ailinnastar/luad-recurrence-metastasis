## Step 1 (logistic engine): same nested LODO design as the survival engines,
## but the model is a binomial GLM on clinical + recurrence/metastasis scores.
## Sourced by step1_run_all_engines.R.

random_seed <- as.integer(Sys.getenv("COMBO_RANDOM_SEED", "1"))
set.seed(random_seed)

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(survival)
  if (requireNamespace("pROC", quietly = TRUE)) library(pROC)
})

ROOT <- project_root()
train_dir <- training_dir(ROOT)
ext_dir <- external_dir(ROOT)
frozen_dir <- frozen_metastasis_dir(ROOT)
out_dir <- file.path(ROOT, results_root(), "logistic_regression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_K <- c(5L, 10L, 15L)
FDR_LIMMA <- 0.05
EXT_FU_MONTHS_COL <- "months_to_last_contact_or_death:ch1"
k_grid <- combo_k_grid_values()
inner_cv_folds <- combo_inner_cv_folds()
message("Logistic regression -> ", out_dir)

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

glm_train_pred_metrics <- function(y_tr, X_tr, time_te, event_te, X_te) {
  if (ncol(X_tr) == 0L) return(na_metrics())
  df_tr <- data.frame(y = y_tr, X_tr, check.names = FALSE)
  df_te <- data.frame(X_te, check.names = FALSE)
  fit <- suppressWarnings(tryCatch(glm(y ~ ., data = df_tr, family = binomial()), error = function(e) NULL))
  if (is.null(fit)) return(na_metrics())
  pr <- tryCatch(as.numeric(predict(fit, newdata = df_te, type = "response")), error = function(e) NULL)
  if (is.null(pr) || length(pr) != length(event_te)) return(na_metrics())
  c(cindex = cindex_harrell(as.numeric(time_te), as.integer(event_te), pr),
    auroc = auc_binary(as.integer(event_te), pr))
}

glm_metrics_for_inner <- function(time_tr, event_tr, X_tr, time_te, event_te, X_te) {
  glm_train_pred_metrics(event_tr, X_tr, time_te, event_te, X_te)
}
rank_genes_for_inner <- function(mat, y, time_tr, event_tr, fdr) ranked_limma_recurrence(mat, y, fdr = fdr)

# ---- Load + collapse probes to genes (same as survival engine) ----

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

# ---- Outer LODO loop: inner CV grid + cache each fold ----

datasets <- sort(unique(as.character(pheno$Dataset)))
inner_k_by_holdout <- list()
inner_k_fold_detail <- list()
nested_fold_cache <- list()

for (hold in datasets) {
  tr <- rownames(pheno)[pheno$Dataset != hold]
  te <- rownames(pheno)[pheno$Dataset == hold]
  if (length(te) < 5L) next
  ptr <- pheno[tr, , drop = FALSE]
  pte <- pheno[te, , drop = FALSE]

  cb <- combat_train_test(expr_gene[, tr, drop = FALSE], expr_gene[, te, drop = FALSE], ptr$Dataset, hold)
  adj_tr <- cb$left
  adj_te <- cb$right
  ranked <- ranked_limma_recurrence(adj_tr, ptr$Recurrence, fdr = FDR_LIMMA)
  imp <- impute_clinical(as.data.frame(ptr[, clinical_cols, drop = FALSE]),
                         as.data.frame(pte[, clinical_cols, drop = FALSE]))

  time_tr <- setNames(as.numeric(ptr$FollowUp_days), rownames(ptr))
  event_tr <- setNames(as.integer(ptr$Recurrence), rownames(ptr))
  time_te <- setNames(as.numeric(pte$FollowUp_days), rownames(pte))
  event_te <- setNames(as.integer(pte$Recurrence), rownames(pte))

  hold_i <- match(hold, datasets)
  inner_res <- inner_cv_k_grid_one_holdout(
    expr_gene = expr_gene[, tr, drop = FALSE], ptr = ptr, pte = pte, held_out = hold,
    meta_genes = meta_genes, kr_vals = k_grid$recurrence, km_vals = k_grid$metastasis,
    clinical_cols = clinical_cols, n_inner = inner_cv_folds, seed = random_seed + hold_i,
    time_tr_all = time_tr, event_tr_all = event_tr, time_te_outer = time_te, event_te_outer = event_te,
    fit_metrics_fn = glm_metrics_for_inner, rank_genes_fn = rank_genes_for_inner, fdr = FDR_LIMMA
  )
  inner_k_by_holdout[[length(inner_k_by_holdout) + 1L]] <- inner_res$summary
  if (!is.null(inner_res$fold_detail)) inner_k_fold_detail[[length(inner_k_fold_detail) + 1L]] <- inner_res$fold_detail

  nested_fold_cache[[hold]] <- stash_nested_outer_fold(
    held_out = hold, adj_tr = adj_tr, adj_te = adj_te, tt_rank = ranked,
    Xc_tr = imp$train, Xc_te = imp$test, time_tr_all = time_tr, event_tr_all = event_tr,
    time_te_outer = time_te, event_te_outer = event_te, y_rec = ptr$Recurrence
  )
  message(sprintf("  Inner CV grid [%s]: %d kr x km cells", hold, nrow(inner_res$summary)))
}

# ---- Lock K, score across outer folds, refit on development ----

inner_all <- do.call(rbind, inner_k_by_holdout)
fold_all <- if (length(inner_k_fold_detail)) do.call(rbind, inner_k_fold_detail) else NULL
inner_avg <- summarize_inner_cv_across_holdouts(inner_all)
final_k <- recommend_final_k(inner_avg)

nested_all <- data.frame()
if (!is.null(final_k)) {
  kr_g <- as.integer(final_k$recommended_top_k_recurrence)
  km_g <- as.integer(final_k$recommended_top_k_metastasis)
  full_nested <- run_outer_nested_global_k(nested_fold_cache, kr = kr_g, km = km_g, meta_genes = meta_genes,
                                           fit_metrics_fn = glm_metrics_for_inner, inner_all = inner_all,
                                           rank_genes_fn = rank_genes_for_inner, fdr = FDR_LIMMA)
  baseline_nested <- run_outer_nested_baseline_clin_rec(nested_fold_cache, kr = kr_g,
                                                        fit_metrics_fn = glm_metrics_for_inner,
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
  engine_label = "Logistic regression", out_dir = out_dir, pheno = pheno, tc = tc,
  expr_gene = expr_gene, te_gene = te_gene, meta_genes = meta_genes, clinical_cols = clinical_cols,
  fdr = FDR_LIMMA, fit_metrics_fn = glm_metrics_for_inner, rank_genes_fn = rank_genes_for_inner,
  time_all_fun = function(df) as.numeric(df$FollowUp_days),
  time_ext_fun = function(df) as.numeric(df[[EXT_FU_MONTHS_COL]]) * (365.25 / 12),
  event_fun = function(df) as.integer(df$Recurrence)
)

message("Done: ", out_dir)
