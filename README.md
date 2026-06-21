# LUAD Recurrence & Metastasis Predictor

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
