# Recurrence Prediction Shiny App

## Run locally

```r
setwd("/path/to/recur inter/shiny_recurrence_app")
shiny::runApp(".", host = "127.0.0.1", port = 3840)
```

**Live app:** https://elenna-shine-shiny.shinyapps.io/shiny_recurrence_app/

## Model

Baseline formula (from the final report):

- **Linear predictor:** `0.309 × stage + 0.0033 × age + 3.008 × recurrence score`
- **Recurrence score:** weighted average of 15 genes (TUBA1C, CHRDL1, CDC45, …)
- **Not used:** sex, metastasis score (coefficient zero in the full structural model)

Model weights and coefficients are defined in `model_backend.R`.

## Files

- `app.R` — Shiny UI and server
- `model_backend.R` — scoring and prediction logic
- `demo_data/` — example uploads for local testing
