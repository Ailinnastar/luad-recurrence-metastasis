#!/usr/bin/env Rscript
## Step 0b — download GEO (if needed), merge recurrence cohorts + external GSE68465.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

cmd <- commandArgs(trailingOnly = FALSE)
f <- gsub("~\\+~", " ", sub("^--file=", "", cmd[grep("^--file=", cmd)]))
steps_dir <- if (length(f) && nzchar(f[[1]])) {
  normalizePath(dirname(f[[1]]), mustWork = TRUE)
} else {
  normalizePath("pipeline/steps", mustWork = FALSE)
}
root <- normalizePath(file.path(steps_dir, "..", ".."), mustWork = TRUE)
Sys.setenv(COMBO_PROJECT_ROOT = root)
source(file.path(steps_dir, "helpers.R"))

geo_rec <- geo_recurrence_dir(root)
geo_ext <- geo_external_dir(root)
train_dir <- training_dir(root)
ext_dir <- external_dir(root)
dir.create(geo_rec, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_ext, recursive = TRUE, showWarnings = FALSE)
dir.create(train_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)

load_eset <- function(gse_id, dest) {
  g <- getGEO(gse_id, destdir = dest, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  if (is.list(g)) g[[1]] else g
}

pick_col <- function(p, patterns) {
  nms <- names(p)
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  NULL
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

followup_days <- function(p, dataset) {
  # Column names and unit conventions match the GEO series matrices that
  # produced the committed integrated training set (validated by cor = 1.0
  # against p_merged_clean.rds$FollowUp_days).
  if (dataset == "GSE31210") {
    c <- pick_col(p, c("days before relapse/censor", "days before relapse"))
    if (is.null(c)) stop("GSE31210: cannot find follow-up days column")
    as_num(p[[c]])
  } else if (dataset == "GSE30219") {
    c <- pick_col(p, c("disease free survival in months", "disease free survival.*month"))
    if (is.null(c)) stop("GSE30219: cannot find follow-up (months) column")
    as_num(p[[c]]) * 30
  } else if (dataset == "GSE37745") {
    c <- pick_col(p, c("days to recurrence / to last visit", "days to recurrence"))
    if (is.null(c)) stop("GSE37745: cannot find follow-up days column")
    as_num(p[[c]])
  } else if (dataset == "GSE50081") {
    c <- pick_col(p, c("disease-free survival time", "disease.?free survival time"))
    if (is.null(c)) stop("GSE50081: cannot find follow-up column")
    as_num(p[[c]])
  } else {
    stop("Unknown dataset: ", dataset)
  }
}

maybe_log2 <- function(e) {
  rng <- range(e, na.rm = TRUE)
  if (rng[2] > 50) log2(e + 1) else e
}

gene_symbol_col <- function(fd) {
  for (nm in c("Gene Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol")) {
    if (nm %in% names(fd)) return(nm)
  }
  pick_col(fd, c("symbol", "gene.*symbol"))
}

collapse_probes_to_genes <- function(e_probe, f_annot) {
  sym_col <- gene_symbol_col(f_annot)
  if (is.null(sym_col)) stop("Annotation missing gene symbol column")
  f_annot$ID <- as.character(f_annot$ID)
  if (!"ID" %in% names(f_annot) && nrow(f_annot) == nrow(e_probe)) {
    f_annot$ID <- rownames(e_probe)
  }
  gsym <- toupper(trimws(as.character(f_annot[[sym_col]])))
  gsym[gsym == "" | is.na(gsym)] <- NA_character_
  ok <- !is.na(gsym) & f_annot$ID %in% rownames(e_probe)
  em <- as.matrix(e_probe[f_annot$ID[ok], , drop = FALSE])
  rownames(em) <- NULL
  gvec <- gsym[match(f_annot$ID[ok], f_annot$ID)]
  keys <- split(seq_len(nrow(em)), gvec)
  do.call(
    rbind,
    lapply(names(keys), function(g) {
      ix <- keys[[g]]
      if (length(ix) == 1L) {
        matrix(em[ix, ], nrow = 1L, dimnames = list(g, colnames(em)))
      } else {
        matrix(colMeans(em[ix, , drop = FALSE]), nrow = 1L, dimnames = list(g, colnames(em)))
      }
    })
  )
}

# --- Training cohorts ---
gse1 <- load_eset("GSE31210", geo_rec)
gse2 <- load_eset("GSE30219", geo_rec)
gse6 <- load_eset("GSE37745", geo_rec)
gse7 <- load_eset("GSE50081", geo_rec)

e1_raw <- exprs(gse1)
e2_raw <- exprs(gse2)
e6_raw <- exprs(gse6)
e7_raw <- exprs(gse7)
p1_raw <- pData(gse1)
p2_raw <- pData(gse2)
p6_raw <- pData(gse6)
p7_raw <- pData(gse7)

p1 <- p1_raw[
  p1_raw$`tissue:ch1` == "primary lung tumor" &
    p1_raw$`exclude for prognosis analysis due to incomplete resection or adjuvant therapy:ch1` == "none" &
    !is.na(p1_raw$`relapse:ch1`),
]
p1$Recurrence <- ifelse(p1$`relapse:ch1` == "relapsed", 1L, 0L)
p1$Dataset <- "GSE31210"

p2 <- p2_raw[
  p2_raw$`histology:ch1` == "ADC" &
    !is.na(p2_raw$`relapse (event=1; no event=0):ch1`) &
    p2_raw$`relapse (event=1; no event=0):ch1` != "na",
]
p2$Recurrence <- as.integer(p2$`relapse (event=1; no event=0):ch1`)
p2$Dataset <- "GSE30219"

p6 <- p6_raw[
  p6_raw$`histology:ch1` == "adeno" &
    p6_raw$`recurrence:ch1` %in% c("yes", "no"),
]
p6$Recurrence <- ifelse(p6$`recurrence:ch1` == "yes", 1L, 0L)
p6$Dataset <- "GSE37745"

p7 <- p7_raw[
  p7_raw$`histology:ch1` == "adenocarcinoma" &
    p7_raw$`recurrence:ch1` %in% c("Y", "N"),
]
p7$Recurrence <- ifelse(p7$`recurrence:ch1` == "Y", 1L, 0L)
p7$Dataset <- "GSE50081"

e1 <- e1_raw[, intersect(colnames(e1_raw), rownames(p1)), drop = FALSE]
e2 <- e2_raw[, intersect(colnames(e2_raw), rownames(p2)), drop = FALSE]
e6 <- e6_raw[, intersect(colnames(e6_raw), rownames(p6)), drop = FALSE]
e7 <- e7_raw[, intersect(colnames(e7_raw), rownames(p7)), drop = FALSE]
p1 <- p1[colnames(e1), , drop = FALSE]
p2 <- p2[colnames(e2), , drop = FALSE]
p6 <- p6[colnames(e6), , drop = FALSE]
p7 <- p7[colnames(e7), , drop = FALSE]

e1 <- maybe_log2(e1)

stage_map <- c(
  "IA" = 1, "1A" = 1, "1a" = 1, "IB" = 2, "1B" = 2, "1b" = 2,
  "IIA" = 3, "2A" = 3, "2a" = 3, "IIB" = 4, "2B" = 4, "2b" = 4,
  "IIIA" = 5, "3A" = 5, "3a" = 5, "IIIB" = 6, "3B" = 6, "3b" = 6,
  "II" = 3, "IV" = 7, "4" = 7, "M0" = 1
)

clin1 <- data.frame(
  Dataset = "GSE31210",
  Recurrence = p1$Recurrence,
  Age = as_num(p1$`age (years):ch1`),
  Gender_bin = ifelse(tolower(trimws(p1$`gender:ch1`)) %in% c("female", "f"), 0L, 1L),
  Stage_num = stage_map[as.character(p1$`pathological stage:ch1`)],
  FollowUp_days = followup_days(p1, "GSE31210"),
  row.names = rownames(p1),
  stringsAsFactors = FALSE
)
clin2 <- data.frame(
  Dataset = "GSE30219",
  Recurrence = p2$Recurrence,
  Age = as_num(p2$`age at surgery:ch1`),
  Gender_bin = ifelse(tolower(trimws(p2$`gender:ch1`)) %in% c("female", "f"), 0L, 1L),
  Stage_num = stage_map[as.character(p2$`pm stage:ch1`)],
  FollowUp_days = followup_days(p2, "GSE30219"),
  row.names = rownames(p2),
  stringsAsFactors = FALSE
)
clin6 <- data.frame(
  Dataset = "GSE37745",
  Recurrence = p6$Recurrence,
  Age = as_num(p6$`age:ch1`),
  Gender_bin = ifelse(tolower(trimws(p6$`gender:ch1`)) %in% c("female", "f"), 0L, 1L),
  Stage_num = stage_map[as.character(p6$`tumor stage:ch1`)],
  FollowUp_days = followup_days(p6, "GSE37745"),
  row.names = rownames(p6),
  stringsAsFactors = FALSE
)
clin7 <- data.frame(
  Dataset = "GSE50081",
  Recurrence = p7$Recurrence,
  Age = as_num(p7$`age:ch1`),
  Gender_bin = ifelse(tolower(trimws(p7$`Sex:ch1`)) %in% c("female", "f"), 0L, 1L),
  Stage_num = stage_map[as.character(p7$`Stage:ch1`)],
  FollowUp_days = followup_days(p7, "GSE50081"),
  row.names = rownames(p7),
  stringsAsFactors = FALSE
)

p_merged <- rbind(clin1, clin2, clin6, clin7)
common_probes <- Reduce(intersect, list(rownames(e1), rownames(e2), rownames(e6), rownames(e7)))
e_merged <- cbind(
  e1[common_probes, , drop = FALSE],
  e2[common_probes, , drop = FALSE],
  e6[common_probes, , drop = FALSE],
  e7[common_probes, , drop = FALSE]
)
stopifnot(all(colnames(e_merged) == rownames(p_merged)))

f1 <- fData(gse1)
f2 <- fData(gse2)
f6 <- fData(gse6)
f7 <- fData(gse7)
f1$ID <- rownames(f1)
f2$ID <- rownames(f2)
f6$ID <- rownames(f6)
f7$ID <- rownames(f7)
sym1 <- toupper(trimws(as.character(f1[[gene_symbol_col(f1)]])))
sym2 <- toupper(trimws(as.character(f2[[gene_symbol_col(f2)]])))
sym6 <- toupper(trimws(as.character(f6[[gene_symbol_col(f6)]])))
sym7 <- toupper(trimws(as.character(f7[[gene_symbol_col(f7)]])))
names(sym1) <- f1$ID
names(sym2) <- f2$ID
names(sym6) <- f6$ID
names(sym7) <- f7$ID
sym_union <- c(sym1, sym2, sym6, sym7)
sym_union <- sym_union[!is.na(sym_union) & sym_union != ""]
f_merged <- data.frame(
  ID = common_probes,
  `Gene Symbol` = sym_union[common_probes],
  row.names = common_probes,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

saveRDS(e_merged, file.path(train_dir, "e_merged_raw.rds"))
saveRDS(p_merged, file.path(train_dir, "p_merged_clean.rds"))
saveRDS(f_merged, file.path(train_dir, "f_merged_clean.rds"))

expr_gene <- collapse_probes_to_genes(e_merged, f_merged)

# --- External GSE68465 ---
gse11 <- load_eset("GSE68465", geo_ext)
p11 <- pData(gse11)
e11 <- exprs(gse11)
f11 <- fData(gse11)
p11 <- p11[p11$`disease_state:ch1` == "Lung Adenocarcinoma", , drop = FALSE]
p11 <- p11[p11$`first_progression_or_relapse:ch1` %in% c("Yes", "No"), , drop = FALSE]
p11$Recurrence <- ifelse(p11$`first_progression_or_relapse:ch1` == "Yes", 1L, 0L)
common_s <- intersect(colnames(e11), rownames(p11))
e11 <- e11[, common_s, drop = FALSE]
p11 <- p11[common_s, , drop = FALSE]
e11 <- maybe_log2(e11)

sym_col11 <- gene_symbol_col(f11)
probe_genes_11 <- toupper(trimws(as.character(f11[[sym_col11]])))
names(probe_genes_11) <- rownames(f11)
valid <- names(probe_genes_11)[!is.na(probe_genes_11) & probe_genes_11 != ""]
e11_valid <- e11[valid, , drop = FALSE]
rownames(e11_valid) <- probe_genes_11[valid]
unique_genes_11 <- unique(rownames(e11_valid))
keep_idx <- vapply(unique_genes_11, function(g) {
  idx <- which(rownames(e11_valid) == g)
  if (length(idx) == 1L) idx else idx[which.max(apply(e11_valid[idx, , drop = FALSE], 1, var))]
}, integer(1))
e11_genes <- e11_valid[keep_idx, , drop = FALSE]
rownames(e11_genes) <- unique_genes_11

common_genes <- intersect(rownames(expr_gene), rownames(e11_genes))
e_test_common <- e11_genes[common_genes, , drop = FALSE]

p11$T_stage <- gsub(".*pT([0-9]).*", "\\1", p11$`disease_stage:ch1`)
p11$N_stage <- gsub(".*pN([0-9X]).*", "\\1", p11$`disease_stage:ch1`)
p11$Stage_num <- ifelse(
  p11$T_stage == "1" & p11$N_stage == "0", 1L,
  ifelse(p11$T_stage == "2" & p11$N_stage == "0", 2L,
    ifelse(p11$T_stage %in% c("3", "4") & p11$N_stage == "0", 3L,
      ifelse(p11$N_stage == "1", 4L, ifelse(p11$N_stage == "2", 5L, NA_integer_))))
)
p11$Age <- as_num(p11$`age:ch1`)
p11$Gender_bin <- ifelse(
  tolower(trimws(p11$`Sex:ch1`)) %in% c("female", "f"), 0L,
  ifelse(tolower(trimws(p11$`Sex:ch1`)) %in% c("male", "m"), 1L, NA_integer_)
)

saveRDS(e_test_common, file.path(ext_dir, "test_GSE68465_expression.rds"))
saveRDS(p11, file.path(ext_dir, "test_GSE68465_clinical.rds"))
saveRDS(common_genes, file.path(ext_dir, "common_genes_train_test.rds"))
