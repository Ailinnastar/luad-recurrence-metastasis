## Shiny backend: deployed baseline formula (stage + age + recurrence score).
## Metastasis score is not used in the linear predictor (coefficient = 0 in the full model).

DEPLOYED_RECURRENCE_WEIGHTS <- c(
  TUBA1C = 0.6896,
  CHRDL1 = -0.6425,
  CDC45 = 0.6281,
  KRBOX4 = -0.6166,
  ALG3 = 0.6096,
  FTCD = 0.6033,
  SDPR = -0.6028,
  MYBL2 = 0.6009,
  GAPDH = 0.6006,
  UBE2C = 0.5992,
  TUBB3 = 0.5942,
  SCN7A = -0.5921,
  CENPN = 0.5891,
  PSME3 = 0.5866,
  FOXM1 = 0.5852
)

DEPLOYED_LP_COEF <- c(
  Stage_num = 0.309,
  Age = 0.0033,
  Recurrence_score = 3.008
)

## Maps Cox-style LP to horizon risk when training calibration is unavailable.
DEPLOYED_CALIB_INTERCEPT <- -2.37
DEPLOYED_CALIB_SLOPE <- 0.32

signature_score <- function(mat, genes, w) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) {
    return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  }
  ww <- w[match(g, genes)]
  ok <- is.finite(ww)
  g <- g[ok]
  ww <- ww[ok]
  if (!length(g)) {
    return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  }
  sx <- colSums(mat[g, , drop = FALSE] * ww)
  out <- sx / length(genes)
  if (is.null(names(out)) && !is.null(colnames(mat))) names(out) <- colnames(mat)
  out
}

recurrence_risk_calibrated <- function(lp, calib, horizon_days) {
  horizon_days <- as.numeric(horizon_days)
  if (!is.finite(horizon_days) || horizon_days <= 0) horizon_days <- 365

  if (!is.null(calib$intercept) && !is.null(calib$slope)) {
    p <- stats::plogis(as.numeric(calib$intercept) + as.numeric(calib$slope) * lp)
    return(pmin(pmax(p, 0), 1))
  }

  lp_tr <- calib$lp_train
  time_tr <- calib$time_train
  event_tr <- calib$event_train
  if (is.null(lp_tr) || !length(lp_tr)) {
    return(pmin(pmax(stats::plogis(lp), 0), 1))
  }

  y <- as.integer(event_tr == 1L & time_tr <= horizon_days)
  ok <- is.finite(lp_tr) & is.finite(y)
  if (sum(ok) < 20L) {
    return(pmin(pmax(stats::ecdf(lp_tr)(lp), 0), 1))
  }
  df_cal <- data.frame(y = y[ok], lp_tr = lp_tr[ok])
  cal <- tryCatch(
    stats::glm(y ~ lp_tr, data = df_cal, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(cal)) {
    return(pmin(pmax(stats::ecdf(lp_tr)(lp), 0), 1))
  }
  nd <- data.frame(lp_tr = as.numeric(lp))
  p <- as.numeric(stats::predict(cal, newdata = nd, type = "response"))
  pmin(pmax(p, 0), 1)
}

compute_linear_predictor <- function(model, newdata) {
  preds <- model$predictors
  X <- as.data.frame(newdata[, preds, drop = FALSE])
  lp <- rep(0, nrow(X))
  for (nm in preds) {
    v <- suppressWarnings(as.numeric(X[[nm]]))
    if (length(v) != nrow(X)) v <- rep(v[1], nrow(X))
    v[!is.finite(v)] <- 0
    lp <- lp + as.numeric(model$coef[[nm]]) * v
  }
  lp
}

predict_risk <- function(model, newdata, horizon_days = 365, calib = NULL) {
  if (is.list(model) && identical(model$type, "fixed_lp")) {
    lp <- compute_linear_predictor(model, newdata)
    return(recurrence_risk_calibrated(lp, calib, horizon_days))
  }
  if (is.list(model) && !is.null(model$coef) && !is.null(model$predictors)) {
    lp <- compute_linear_predictor(model, newdata)
    return(recurrence_risk_calibrated(lp, calib, horizon_days))
  }
  as.numeric(stats::predict(model, newdata = newdata, type = "response"))
}

load_deployed_artifacts <- function() {
  rec_w <- DEPLOYED_RECURRENCE_WEIGHTS
  list(
    model = list(
      type = "fixed_lp",
      coef = DEPLOYED_LP_COEF,
      predictors = names(DEPLOYED_LP_COEF)
    ),
    recur_genes = names(rec_w),
    recurrence_genes = names(rec_w),
    recur_weights = rec_w,
    training_gene_symbols = names(rec_w),
    lp_train = NULL,
    time_train = NULL,
    event_train = NULL,
    calib_intercept = DEPLOYED_CALIB_INTERCEPT,
    calib_slope = DEPLOYED_CALIB_SLOPE,
    top_k_recurrence = length(rec_w),
    model_kind = "deployed_baseline_formula",
    model_label = paste(
      "Deployed baseline: stage, age, and 15-gene recurrence score",
      "(metastasis score not included)"
    )
  )
}

training_gene_symbols <- function(ar) {
  if (!is.null(ar$training_gene_symbols)) {
    return(ar$training_gene_symbols)
  }
  if (!is.null(ar$recur_weights)) {
    return(names(ar$recur_weights))
  }
  character(0)
}

normalize_gender <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("male", "m", "1", "1.0")] <- "male"
  out[x %in% c("female", "f", "0", "0.0")] <- "female"
  out
}

patient_expr_matrix <- function(raw) {
  gene_cols <- setdiff(names(raw), c(
    "SampleID", "geo_accession", "Age", "Gender", "Gender_bin", "Stage", "Stage_num",
    "Smoking", "GeneAlt", "A_score", "Recurrence_score", "Metastasis_score"
  ))
  gene_cols <- gene_cols[vapply(raw[gene_cols], function(z) is.numeric(z) || is.integer(z), logical(1))]
  if (!length(gene_cols)) {
    return(NULL)
  }
  mat <- t(as.matrix(raw[, gene_cols, drop = FALSE]))
  rownames(mat) <- toupper(rownames(mat))
  storage.mode(mat) <- "double"
  mat
}

required_signature_genes <- function(ar) {
  unique(toupper(names(ar$recur_weights %||% DEPLOYED_RECURRENCE_WEIGHTS)))
}

align_expr_for_scoring <- function(mat_new, ar, min_genes = NULL) {
  sig <- required_signature_genes(ar)
  g <- intersect(rownames(mat_new), sig)
  need <- min_genes %||% length(sig)
  if (length(g) < need) {
    miss <- setdiff(sig, g)
    stop(
      "Need at least ", need, " of ", length(sig), " signature genes. ",
      "Missing: ", paste(head(miss, 8), collapse = ", "),
      if (length(miss) > 8L) " ..." else ""
    )
  }
  mat_new
}

score_patients <- function(ar, raw, min_genes = NULL) {
  mat <- patient_expr_matrix(raw)
  if (is.null(mat)) stop("No numeric gene columns found in upload.")

  sig <- required_signature_genes(ar)
  min_genes <- min_genes %||% length(sig)
  mat <- align_expr_for_scoring(mat, ar, min_genes = min_genes)

  rec_w <- ar$recur_weights
  rs <- as.numeric(signature_score(mat, names(rec_w), rec_w))

  n <- nrow(raw)
  if (is.null(raw$SampleID)) {
    raw$SampleID <- paste0("Sample_", seq_len(n))
  }
  out <- raw
  if (!"Age" %in% names(out)) out$Age <- 65
  if (!"Stage_num" %in% names(out)) out$Stage_num <- 3L
  if (!"Gender_bin" %in% names(out)) {
    g <- if ("Gender" %in% names(out)) normalize_gender(out$Gender) else NA
    out$Gender_bin <- ifelse(g == "male", 1L, 0L)
  }
  for (nm in c("Stage_num", "Age", "Gender_bin")) {
    v <- suppressWarnings(as.numeric(out[[nm]]))
    if (length(v) != nrow(out)) v <- rep(v[1], nrow(out))
    v[!is.finite(v)] <- if (nm == "Stage_num") 3 else if (nm == "Gender_bin") 0 else 65
    out[[nm]] <- v
  }
  out$Recurrence_score <- rs
  out
}

compute_predictions <- function(ar, raw, horizon_days, conf_level = 0.95, min_genes = NULL) {
  horizon_days <- as.numeric(horizon_days)
  if (!is.finite(horizon_days) || horizon_days <= 0) horizon_days <- 365

  sig_n <- length(required_signature_genes(ar))
  min_genes <- min_genes %||% sig_n

  scored <- score_patients(ar, raw, min_genes = min_genes)
  preds <- ar$model$predictors %||% names(DEPLOYED_LP_COEF)
  X <- scored[, preds, drop = FALSE]
  for (nm in names(X)) X[[nm]] <- suppressWarnings(as.numeric(X[[nm]]))

  calib <- list(
    lp_train = ar$lp_train,
    time_train = ar$time_train,
    event_train = ar$event_train,
    intercept = ar$calib_intercept,
    slope = ar$calib_slope
  )
  pr <- predict_risk(ar$model, X, horizon_days, calib = calib)
  pr <- pmin(pmax(pr, 0), 1)

  mat_for_match <- patient_expr_matrix(raw)
  matched <- 0L
  if (!is.null(mat_for_match)) {
    matched <- length(intersect(required_signature_genes(ar), rownames(mat_for_match)))
  }

  list(
    df = scored,
    risks = pr,
    matched_recur = matched
  )
}

as_prediction_risks <- function(risks) {
  if (is.null(risks)) {
    return(numeric(0))
  }
  if (is.list(risks)) {
    return(as.numeric(unlist(risks, use.names = FALSE)))
  }
  as.numeric(risks)
}

load_uploaded_table <- function(path, ext, sep = ",") {
  ext <- tolower(ext %||% "")
  if (ext == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }
  stop("Unsupported file type. Use CSV, TSV/TXT, or RDS.")
}

to_patient_gene_matrix <- function(expr_obj) {
  ex <- as.matrix(expr_obj)
  if (!is.numeric(ex)) stop("Expression data must be numeric.")
  rn <- rownames(ex)
  cn <- colnames(ex)
  if (is.null(rn) || is.null(cn)) stop("Expression matrix needs row and column names.")
  looks_samples_in_rows <- mean(grepl("^GSM|^GSE|patient|sample", cn, ignore.case = TRUE)) >
    mean(grepl("^GSM|^GSE|patient|sample", rn, ignore.case = TRUE))
  if (looks_samples_in_rows) {
    return(ex)
  }
  t(ex)
}

expression_to_gene_matrix <- function(expr_raw) {
  if (is.matrix(expr_raw)) {
    return(to_patient_gene_matrix(expr_raw))
  }
  if (is.data.frame(expr_raw)) {
    df <- as.data.frame(expr_raw, stringsAsFactors = FALSE)
    id_col <- intersect(c("SampleID", "geo_accession", "sample_id", "ID"), names(df))[1]
    if (!is.na(id_col)) {
      gene_cols <- setdiff(names(df), id_col)
      if (!length(gene_cols)) {
        stop("Expression table needs gene columns besides the sample ID column.")
      }
      mat <- as.matrix(df[, gene_cols, drop = FALSE])
      mode(mat) <- "numeric"
      rownames(mat) <- as.character(df[[id_col]])
      return(to_patient_gene_matrix(t(mat)))
    }
    if (!is.null(rownames(df)) && nrow(df) > 1L) {
      return(to_patient_gene_matrix(df))
    }
  }
  stop(
    "Expression file should be an RDS matrix (genes x samples), or a table with SampleID and gene columns, or a gene-by-sample matrix with row and column names."
  )
}

build_combined_upload <- function(clin_raw, expr_raw) {
  ex <- expression_to_gene_matrix(expr_raw)
  clin <- as.data.frame(clin_raw, stringsAsFactors = FALSE)
  id_col <- intersect(c("SampleID", "geo_accession", "sample_id", "ID"), names(clin))[1]
  if (is.na(id_col)) stop("Clinical file needs SampleID or geo_accession.")
  clin$SampleID <- as.character(clin[[id_col]])
  samp <- colnames(ex)
  if (!all(clin$SampleID %in% samp)) {
    miss <- setdiff(clin$SampleID, samp)
    stop("Clinical IDs not found in expression columns: ", paste(head(miss, 5), collapse = ", "))
  }
  ex_sub <- ex[, clin$SampleID, drop = FALSE]
  gene_df <- as.data.frame(t(ex_sub), stringsAsFactors = FALSE, check.names = FALSE)
  cbind(clin, gene_df)
}

cohort_prioritization <- function(sample_ids, risks) {
  risks <- as.numeric(risks)
  o <- order(-risks, sample_ids)
  rk <- seq_along(o)
  z <- (risks - mean(risks)) / (stats::sd(risks) + 1e-8)
  pct <- round(100 * stats::ecdf(risks)(risks), 1)
  lab <- ifelse(
    risks >= 0.75, "Very high",
    ifelse(risks >= 0.50, "High",
      ifelse(risks >= 0.25, "Moderate", "Lower priority")
    )
  )
  data.frame(
    SampleID = sample_ids[o],
    cohort_rank = rk,
    recurrence_risk = risks[o],
    percentile_urgency = pct[o],
    z_vs_cohort_mean = z[o],
    priority_label = lab[o],
    stringsAsFactors = FALSE
  )
}

human_gene_direction <- function(weight) {
  if (!is.finite(weight)) return("")
  if (weight >= 0) "higher expression raises score" else "higher expression lowers score"
}

`%||%` <- function(x, y) if (is.null(x)) y else x
