#!/usr/bin/env python3
"""Step 0a (metastasis, Path A): build the 5-study brain-vs-lung gene support
table, then freeze the top-K (5/10/15) metastasis gene lists used by the
recurrence models. Inputs live in pipeline/geo/metastasis."""
from __future__ import annotations

import os
import subprocess
import tarfile
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
META_IN = Path(os.environ.get("COMBO_METASTASIS_INPUTS", PIPELINE_ROOT / "geo/metastasis"))
A_OUT = Path(os.environ.get("COMBO_METASTASIS_RESULTS", PIPELINE_ROOT / "metastasis_results"))
FROZEN_OUT = Path(os.environ.get("COMBO_MATA_FROZEN_DIR", PIPELINE_ROOT / "frozen"))
A_OUT.mkdir(parents=True, exist_ok=True)


def limma_two_group(expr_gene_x_sample, groups, blocks=None):
    """Brain (B) vs lung (L) differential expression via limma, called out to R.
    Uses duplicateCorrelation blocking when paired samples are supplied."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        expr_path = td / "expr.tsv"
        meta_path = td / "meta.tsv"
        out_path = td / "out.tsv"
        expr_gene_x_sample.to_csv(expr_path, sep="\t")
        pd.DataFrame(
            {
                "sample": list(expr_gene_x_sample.columns),
                "group": groups,
                "block": (blocks if blocks is not None else ["NA"] * len(groups)),
            }
        ).to_csv(meta_path, sep="\t", index=False)
        r_code = f"""
        suppressPackageStartupMessages(library(limma))
        expr <- read.delim("{expr_path.as_posix()}", check.names=FALSE, row.names=1)
        meta <- read.delim("{meta_path.as_posix()}", stringsAsFactors=FALSE)
        expr <- as.matrix(expr[, meta$sample, drop=FALSE])
        grp <- factor(meta$group, levels=c("L","B"))
        design <- model.matrix(~0 + grp)
        colnames(design) <- c("L","B")
        blk <- meta$block
        use_block <- any(blk != "NA") && length(unique(blk[blk != "NA"])) > 1
        if (use_block) {{
          dc <- duplicateCorrelation(expr, design, block=blk)
          fit <- lmFit(expr, design, block=blk, correlation=dc$consensus)
        }} else {{
          fit <- lmFit(expr, design)
        }}
        cont <- makeContrasts(BminusL = B - L, levels=design)
        fit2 <- eBayes(contrasts.fit(fit, cont))
        tt <- topTable(fit2, coef="BminusL", number=Inf, sort.by="none")
        out <- data.frame(gene=rownames(tt), logFC=tt$logFC, t=tt$t, p=tt$P.Value, q=tt$adj.P.Val)
        write.table(out, file="{out_path.as_posix()}", sep="\\t", row.names=FALSE, quote=FALSE)
        """
        subprocess.run(["Rscript", "-e", r_code], check=True)
        return pd.read_csv(out_path, sep="\t")


def build_a_spatial_200563():
    expr = pd.read_csv(META_IN / "GSE200563_processed_data.txt", sep="\t").rename(columns={"Gene#": "gene"})
    expr["gene"] = expr["gene"].astype(str).str.upper().str.strip()
    expr = expr.groupby("gene", as_index=False).mean().set_index("gene")
    smap = pd.read_csv(META_IN / "samples_roi_table.csv")
    L = smap[smap["roi_code"].astype(str).str.startswith("L") & ~smap["roi_code"].astype(str).str.startswith("LB")]
    LB = smap[smap["roi_code"].astype(str).str.startswith("LB")]
    common = sorted(set(L["patient_id"]).intersection(set(LB["patient_id"])))
    L_rows, LB_rows = [], []
    for pid in common:
        lcols = [c for c in L[L["patient_id"] == pid]["roi_code"] if c in expr.columns]
        bcols = [c for c in LB[LB["patient_id"] == pid]["roi_code"] if c in expr.columns]
        if not lcols or not bcols:
            continue
        L_rows.append(expr[lcols].mean(axis=1))
        LB_rows.append(expr[bcols].mean(axis=1))
    Lm = pd.concat(L_rows, axis=1).T
    Bm = pd.concat(LB_rows, axis=1).T
    Bm_log = np.log1p(Bm)
    Lm_log = np.log1p(Lm)
    patient_ids = [f"PID_{i+1}" for i in range(Bm_log.shape[0])]
    paired = pd.concat(
        [Bm_log.T.set_axis([f"B_{p}" for p in patient_ids], axis=1),
         Lm_log.T.set_axis([f"L_{p}" for p in patient_ids], axis=1)],
        axis=1,
    )
    de = limma_two_group(paired, ["B"] * len(patient_ids) + ["L"] * len(patient_ids), patient_ids + patient_ids)
    out = pd.DataFrame(
        {
            "gene": de["gene"].astype(str),
            "A_spatial_logFC": pd.to_numeric(de["logFC"], errors="coerce"),
            "A_spatial_t": pd.to_numeric(de["t"], errors="coerce"),
            "A_spatial_p": pd.to_numeric(de["p"], errors="coerce"),
            "A_spatial_q": pd.to_numeric(de["q"], errors="coerce"),
            "A_spatial_n_pairs": Lm.shape[0],
        }
    )
    out.to_csv(A_OUT / "A_spatial_paired_LB_vs_L.tsv", sep="\t", index=False)
    return out


def build_a_bulk_271259():
    xlsx = META_IN / "GSE271259_processed_data.xlsx"
    meta = pd.read_excel(xlsx, sheet_name="Table1")[["Sample number", "Tumor_Type"]]
    meta["Sample number"] = meta["Sample number"].astype(str)
    meta["Tumor_Type"] = meta["Tumor_Type"].astype(str).str.strip()
    mat = pd.read_excel(xlsx, sheet_name="Table3").rename(columns={"Gene": "gene"})
    mat["gene"] = mat["gene"].astype(str).str.split("!", n=1).str[0].str.upper().str.strip()
    num_cols = [c for c in mat.columns if c != "gene"]
    for c in num_cols:
        mat[c] = pd.to_numeric(mat[c], errors="coerce")
    mat = mat.groupby("gene", as_index=False)[num_cols].mean().set_index("gene")
    mt = meta.set_index("Sample number")
    brain = [c for c in num_cols if mt.reindex([c])["Tumor_Type"].iloc[0] == "Brain"]
    lung = [c for c in num_cols if mt.reindex([c])["Tumor_Type"].iloc[0] == "Lung"]
    logx = np.log2(mat[brain + lung] + 1.0)
    de = limma_two_group(logx, ["B"] * len(brain) + ["L"] * len(lung))
    out = pd.DataFrame(
        {
            "gene": de["gene"].astype(str),
            "logFC_bulk_direct": pd.to_numeric(de["logFC"], errors="coerce"),
            "t_bulk_direct": pd.to_numeric(de["t"], errors="coerce"),
            "p_bulk_direct": pd.to_numeric(de["p"], errors="coerce"),
            "q_bulk_direct": pd.to_numeric(de["q"], errors="coerce"),
            "n_brain_bulk": len(brain),
            "n_lung_bulk": len(lung),
        }
    )
    out.to_csv(A_OUT / "A_bulk_GSE271259_direct_from_table3.tsv", sep="\t", index=False)
    return out


def build_a_nanostring_161116():
    d = pd.read_csv(META_IN / "GSE161116_Raw_data.txt.gz", sep="\t")
    d = d[d["Class Name"].astype(str).str.contains("Endogenous", case=False, na=False)].copy()
    expr = d.rename(columns={"Probe Name": "gene"})
    cols = [c for c in expr.columns if isinstance(c, str) and (c.endswith("L") or c.endswith("B"))]
    expr = expr[["gene"] + cols]
    for c in cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
    expr["gene"] = expr["gene"].astype(str).str.upper().str.strip()
    expr = expr.groupby("gene", as_index=False)[cols].mean().set_index("gene")
    logx = np.log2(expr + 1)
    B = [c for c in cols if c.endswith("B")]
    L = [c for c in cols if c.endswith("L")]
    de = limma_two_group(logx[B + L], ["B"] * len(B) + ["L"] * len(L))
    out = pd.DataFrame(
        {
            "gene": de["gene"].astype(str),
            "logFC_161116": pd.to_numeric(de["logFC"], errors="coerce"),
            "t_161116": pd.to_numeric(de["t"], errors="coerce"),
            "p_161116": pd.to_numeric(de["p"], errors="coerce"),
            "q_161116": pd.to_numeric(de["q"], errors="coerce"),
        }
    )
    out.to_csv(A_OUT / "A_validation_GSE161116_B_vs_L.tsv", sep="\t", index=False)
    return out


def build_a_nanostring_248830():
    d = pd.read_csv(META_IN / "GSE248830_Raw_data.csv.gz").rename(columns={"Probe Name": "gene"})
    p_cols = [c for c in d.columns if str(c).startswith("Primary LUAD ")]
    b_cols = [c for c in d.columns if str(c).startswith("BM-LUAD ")]
    d = d[["gene"] + p_cols + b_cols].copy()
    d["gene"] = d["gene"].astype(str).str.upper().str.strip()
    for c in p_cols + b_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.groupby("gene", as_index=False)[p_cols + b_cols].mean().set_index("gene")
    logx = np.log2(d + 1)
    pairs = [(pc, f"BM-LUAD {pc.split()[-1]}") for pc in p_cols if f"BM-LUAD {pc.split()[-1]}" in b_cols]
    pa = [x for x, _ in pairs]
    pb = [x for _, x in pairs]
    pair_ids = [f"PAIR_{i+1}" for i in range(len(pairs))]
    paired = pd.concat(
        [logx[pb].set_axis([f"B_{p}" for p in pair_ids], axis=1),
         logx[pa].set_axis([f"L_{p}" for p in pair_ids], axis=1)],
        axis=1,
    )
    de = limma_two_group(paired, ["B"] * len(pair_ids) + ["L"] * len(pair_ids), pair_ids + pair_ids)
    out = pd.DataFrame(
        {
            "gene": de["gene"].astype(str),
            "logFC_248830": pd.to_numeric(de["logFC"], errors="coerce"),
            "t_248830": pd.to_numeric(de["t"], errors="coerce"),
            "p_248830": pd.to_numeric(de["p"], errors="coerce"),
            "q_248830": pd.to_numeric(de["q"], errors="coerce"),
        }
    )
    out.to_csv(A_OUT / "A_validation_GSE248830_paired_BM_vs_Primary.tsv", sep="\t", index=False)
    return out


def build_a_sc_223499():
    tar_path = META_IN / "GSE223499_RAW.tar"
    soft_path = META_IN / "GSE223499_family.soft"
    source_by_gsm = {}
    cur = None
    for ln in soft_path.read_text(errors="ignore").splitlines():
        if ln.startswith("^SAMPLE ="):
            cur = ln.split("=", 1)[1].strip()
        elif cur and ln.startswith("!Sample_source_name_ch1 ="):
            source_by_gsm[cur] = ln.split("=", 1)[1].strip()

    expr_by_sample = {}
    records = []
    with tarfile.open(tar_path, "r") as tf:
        for m in [x for x in tf.getmembers() if x.name.endswith("_sn_raw_feature_bc_matrix.h5")]:
            gsm = m.name.split("_", 1)[0]
            src = source_by_gsm.get(gsm, "")
            grp = "BRAIN" if "brain metastasis" in src.lower() else ("PRIMARY" if "primary tumor" in src.lower() else "OTHER")
            with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
                tmp.write(tf.extractfile(m).read())
                tmppath = tmp.name
            try:
                with h5py.File(tmppath, "r") as h:
                    names = np.array(h["matrix/features/name"]).astype(str)
                    idx = np.array(h["matrix/indices"])
                    val = np.array(h["matrix/data"], dtype=np.float64)
                    sums = np.bincount(idx, weights=val, minlength=len(names))
                expr_by_sample[gsm] = pd.Series(sums, index=pd.Index(names).str.upper()).groupby(level=0).sum()
                records.append({"gsm": gsm, "source": src, "group": grp})
            finally:
                os.unlink(tmppath)

    X = pd.DataFrame(expr_by_sample).fillna(0.0)
    Xlog = np.log1p(X.divide(X.sum(axis=0), axis=1) * 1e6)
    brain = [r["gsm"] for r in records if r["group"] == "BRAIN"]
    prim = [r["gsm"] for r in records if r["group"] == "PRIMARY"]
    de = limma_two_group(Xlog[brain + prim], ["B"] * len(brain) + ["L"] * len(prim))
    out = pd.DataFrame(
        {
            "gene": de["gene"].astype(str),
            "logFC_223499": pd.to_numeric(de["logFC"], errors="coerce"),
            "t_223499": pd.to_numeric(de["t"], errors="coerce"),
            "p_223499": pd.to_numeric(de["p"], errors="coerce"),
            "q_223499": pd.to_numeric(de["q"], errors="coerce"),
            "n_brain_223499": len(brain),
            "n_primary_223499": len(prim),
        }
    )
    out.to_csv(A_OUT / "A_GSE223499_pseudobulk_BM_vs_Primary.tsv", sep="\t", index=False)
    return out


def build_gene_support(a_sp, a_bk, a161, a248, a223):
    """Merge the 5 studies, mark FDR<0.05 hits, rank by mean significant logFC,
    and write the gene support table + the top-K membership flags."""
    M = (
        a_sp[["gene", "A_spatial_logFC", "A_spatial_t", "A_spatial_p", "A_spatial_q"]]
        .rename(columns={"A_spatial_logFC": "sp_logFC", "A_spatial_t": "sp_sign_raw", "A_spatial_p": "sp_p", "A_spatial_q": "sp_q"})
        .merge(a_bk[["gene", "logFC_bulk_direct", "t_bulk_direct", "p_bulk_direct", "q_bulk_direct"]], on="gene", how="outer")
        .merge(a161[["gene", "logFC_161116", "t_161116", "p_161116", "q_161116"]], on="gene", how="outer")
        .merge(a248[["gene", "logFC_248830", "t_248830", "p_248830", "q_248830"]], on="gene", how="outer")
        .merge(a223[["gene", "logFC_223499", "t_223499", "p_223499", "q_223499"]], on="gene", how="outer")
    )
    M["gene"] = M["gene"].astype(str).str.upper()
    M["sp_sig"] = M["sp_q"] < 0.05
    M["bk_sig"] = M["q_bulk_direct"] < 0.05
    M["n161_sig"] = M["q_161116"] < 0.05
    M["n248_sig"] = M["q_248830"] < 0.05
    M["a223499_sig"] = M["q_223499"] < 0.05
    sig_cols = ["sp_sig", "bk_sig", "n161_sig", "n248_sig", "a223499_sig"]
    logfc_cols = ["sp_logFC", "logFC_bulk_direct", "logFC_161116", "logFC_248830", "logFC_223499"]
    lc = M[logfc_cols].to_numpy(dtype=float)
    sg = M[sig_cols].fillna(False).to_numpy()
    with np.errstate(all="ignore"):
        M["A5_sig_mean_logFC"] = np.nanmean(np.where(sg, lc, np.nan), axis=1)
    rankable = M[np.isfinite(M["A5_sig_mean_logFC"])].copy()
    rankable["A5_abs_sig_mean_logFC"] = rankable["A5_sig_mean_logFC"].abs()
    rankable = rankable.sort_values("A5_abs_sig_mean_logFC", ascending=False).reset_index(drop=True)
    rankable["A5_sig_logfc_rank"] = np.arange(1, len(rankable) + 1)
    M = M.merge(rankable[["gene", "A5_sig_logfc_rank", "A5_abs_sig_mean_logFC"]], on="gene", how="left")
    for k in (5, 10, 15):
        gk = set(rankable.head(min(k, len(rankable)))["gene"])
        M[f"A5_sig_top{k:02d}_logfc"] = M["gene"].isin(gk)
    M.to_csv(A_OUT / "A_multidataset_5study_gene_support.tsv", sep="\t", index=False)
    for k in (5, 10, 15):
        top = rankable.head(min(k, len(rankable)))[["gene", "A5_sig_mean_logFC", "A5_abs_sig_mean_logFC", "A5_sig_logfc_rank"]]
        top.to_csv(A_OUT / f"A_genes_5study_sig_top{k:02d}_logfc.tsv", sep="\t", index=False)
    return M


def freeze_top_k(support):
    """Write frozen top-K (5/10/15) metastasis gene lists from the support table."""
    FROZEN_OUT.mkdir(parents=True, exist_ok=True)
    km_suffix = {5: "05", 10: "10", 15: "15"}
    for k in (5, 10, 15):
        col = f"A5_sig_top{k:02d}_logfc"
        genes = support.loc[support[col].fillna(False).astype(bool), "gene"]
        genes = genes.astype(str).str.upper().str.strip()
        genes = genes[genes.ne("")].drop_duplicates().tolist()
        pd.DataFrame({"gene": genes}).to_csv(FROZEN_OUT / f"frozen_A_sig_top{km_suffix[k]}_genes.tsv", sep="\t", index=False)
        if k == 10:
            pd.DataFrame({"gene": genes}).to_csv(FROZEN_OUT / "frozen_A_signature_genes.tsv", sep="\t", index=False)


def main():
    support = build_gene_support(
        build_a_spatial_200563(),
        build_a_bulk_271259(),
        build_a_nanostring_161116(),
        build_a_nanostring_248830(),
        build_a_sc_223499(),
    )
    freeze_top_k(support)


if __name__ == "__main__":
    main()
