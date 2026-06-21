# LUAD recurrence pipeline

Self-contained under `pipeline/`: download, integrate, model, evaluate.

## One command

```bash
cd "recur inter"
bash pipeline/run.sh
```

| Step    | What                                                                     | Output location                                                       |
| ------- | ------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| **0a**  | Download metastasis GEO → build + freeze lists (only if `frozen/` empty) | `pipeline/frozen/`, `pipeline/geo/metastasis/`                        |
| **0b**  | Download recurrence + external GEO → merge RDS                           | `pipeline/data/training/`, `pipeline/data/external/`, `pipeline/geo/` |
| **1–3** | Seven engines + evaluation tables                                        | `pipeline/model result/`                                              |

Requires R (`GEOquery`, `limma`, `glmnet`, `survival`, …) and Python (`pandas`, `numpy`, `h5py`) for step 0a.

## Layout (everything under `pipeline/`)

```
pipeline/
├── frozen/                 # committed metastasis gene lists (required)
├── geo/                    # GEO download cache (recreated by 0a/0b)
├── steps/                  # all code
├── run.R
└── run.sh
```
