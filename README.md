# Lung Cancer Recurrence ML

A university research project investigating whether brain-metastasis transcriptomic signatures improve LUAD recurrence prediction. Includes a reproducible pipeline, final report, and an interactive platform exploring future clinician decision support.

**Start here:** [Analysis website](https://ailinnastar.github.io/luad-recurrence-metastasis/) | [Interactive platform](https://elenna-shine-shiny.shinyapps.io/shiny_recurrence_app/) | [Model card](docs/MODEL_CARD.md)

## What the project demonstrates

- Integration of public GEO gene-expression and clinical cohorts.
- Comparison of seven modelling engines with cross-cohort and external validation.
- A Shiny interface using age, stage, and a 15-gene recurrence score.
- A distinction between biologically plausible signatures and consistent predictive improvement: the deployed baseline excludes the metastasis score.

## Prerequisites and checks

The app requires R packages `shiny`, `bslib`, `ggplot2`, `DT`, and `scales`. The pipeline additionally needs Bioconductor packages `GEOquery`, `Biobase`, `limma`, `sva`, and CRAN packages `survival`, `glmnet`, `pROC`, `randomForestSRC`, `xgboost`. Python requirements are in `Pipeline and Shiny/pipeline/steps/requirements-modeling.txt`.

Run lightweight software checks from the repository root:

```bash
Rscript tests/test_backend.R
```

Full GEO downloads and model training are separate, resource-intensive steps. The lightweight checks do not reproduce the published evaluation.

## Research scope

This repository contains university team submission materials. Ailinna Zhang's portfolio links to this shared project; repository ownership does not imply sole authorship of every component. Existing source and report attributions are preserved.

This is a research prototype. The default fixed calibration does not change its probability mapping when the follow-up horizon changes; outputs are not validated patient-care recommendations. See the [model card](docs/MODEL_CARD.md) for details.

Integrating brain-metastasis transcriptomic signatures for recurrence prediction in **lung adenocarcinoma (LUAD)**.

## Contents

| Folder | Description |
|--------|-------------|
| [THE FINAL REPORT/](THE%20FINAL%20REPORT/) | Final report (`report.html`) and reproducible Quarto source |
| [Pipeline and Shiny/](Pipeline%20and%20Shiny/) | Analysis pipeline and Shiny recurrence-prediction app |

## Quick start

**Run the pipeline:**

```bash
cd "Pipeline and Shiny"
bash pipeline/run.sh
```

**Run the Shiny app locally:**

```r
setwd("Pipeline and Shiny/shiny_recurrence_app")
shiny::runApp(".", host = "127.0.0.1", port = 3840)
```

**Live demo:** https://elenna-shine-shiny.shinyapps.io/shiny_recurrence_app/

**Project website:** https://ailinnastar.github.io/luad-recurrence-metastasis/

See each subfolder's README for full reproduction instructions.
