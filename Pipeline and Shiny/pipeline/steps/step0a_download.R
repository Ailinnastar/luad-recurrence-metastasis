#!/usr/bin/env Rscript
## Step 0a — download metastasis GEO supplementary files (cached under pipeline/geo/metastasis).

suppressPackageStartupMessages(library(GEOquery))

cmd <- commandArgs(trailingOnly = FALSE)
f <- gsub("~\\+~", " ", sub("^--file=", "", cmd[grep("^--file=", cmd)]))
steps_dir <- if (length(f) && nzchar(f[[1]])) dirname(normalizePath(f[[1]], mustWork = TRUE)) else "pipeline/steps"
dest <- Sys.getenv("COMBO_METASTASIS_INPUTS", file.path(dirname(steps_dir), "geo", "metastasis"))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

need <- function(name) !file.exists(file.path(dest, name))

copy_match <- function(files, pattern, out_name) {
  if (!need(out_name)) return(invisible(NULL))
  hit <- grep(pattern, files, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) return(invisible(NULL))
  file.copy(hit[[1]], file.path(dest, out_name), overwrite = TRUE)
}

all_under <- function(base) {
  if (!dir.exists(base)) return(character())
  list.files(base, recursive = TRUE, full.names = TRUE)
}

dl_supp <- function(gse) {
  sub <- file.path(dest, gse)
  if (!dir.exists(sub)) getGEOSuppFiles(gse, baseDir = dest)
  all_under(file.path(dest, gse))
}

if (need("GSE161116_Raw_data.txt.gz")) {
  copy_match(dl_supp("GSE161116"), "Raw.*\\.txt\\.gz$", "GSE161116_Raw_data.txt.gz")
}
if (need("GSE248830_Raw_data.csv.gz")) {
  copy_match(dl_supp("GSE248830"), "Raw.*\\.csv\\.gz$", "GSE248830_Raw_data.csv.gz")
}
if (need("GSE271259_processed_data.xlsx")) {
  copy_match(dl_supp("GSE271259"), "\\.xlsx$", "GSE271259_processed_data.xlsx")
}
if (need("GSE200563_processed_data.txt")) {
  f <- dl_supp("GSE200563")
  copy_match(f, "processed.*\\.txt$", "GSE200563_processed_data.txt")
  if (need("GSE200563_processed_data.txt")) copy_match(f, "processed.*data.*\\.txt$", "GSE200563_processed_data.txt")
}
if (need("samples_roi_table.csv")) {
  f <- c(dl_supp("GSE200563"), all_under(dest))
  copy_match(f, "roi.*\\.csv$|samples.*roi.*\\.csv$", "samples_roi_table.csv")
}
if (need("GSE223499_RAW.tar")) {
  copy_match(dl_supp("GSE223499"), "RAW\\.tar$", "GSE223499_RAW.tar")
}
if (need("GSE223499_family.soft")) {
  soft <- list.files(dest, pattern = "GSE223499.*family\\.soft(\\.gz)?$", full.names = TRUE)
  if (!length(soft)) {
    getGEO("GSE223499", destdir = dest, GSEMatrix = FALSE)
    soft <- list.files(dest, pattern = "GSE223499.*family\\.soft(\\.gz)?$", full.names = TRUE)
  }
  if (length(soft)) {
    out <- file.path(dest, "GSE223499_family.soft")
    if (grepl("\\.gz$", soft[[1]], ignore.case = TRUE)) {
      con <- gzfile(soft[[1]], "rt")
      writeLines(readLines(con, warn = FALSE), out)
      close(con)
    } else {
      file.copy(soft[[1]], out, overwrite = TRUE)
    }
  }
}
