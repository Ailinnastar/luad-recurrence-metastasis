# LUAD research model card

## Purpose

Research and education about transcriptomic recurrence prediction in LUAD, with an interactive platform exploring future clinician decision support. This repository does not establish clinical utility.

## Inputs and baseline

The application scores age, stage, and expression of 15 recurrence-signature genes. Coefficients live in `Pipeline and Shiny/shiny_recurrence_app/model_backend.R`. The brain-metastasis score is excluded from the deployed linear predictor.

## Evaluation and calibration

Read the final report for cohort construction, metrics, and external validation. The app's default artifacts supply a fixed logistic intercept and slope. In this path, `recurrence_risk_calibrated()` does not use `horizon_days` to change the returned probability. Changing the follow-up window does not produce a validated horizon-specific survival estimate. Demonstration priority labels must not be used as clinical thresholds.

## Data and limitations

Analysis scripts retrieve public GEO cohorts. Use synthetic examples for demonstrations. Do not commit identifiable patient records, deployment credentials, or private clinical exports. Cohort selection, censoring, assay differences, preprocessing, and incomplete clinical variables can affect generalisation.

Backend tests verify software behaviour, not discrimination, calibration, clinical safety, or usefulness for patient care. Full model training remains a separate workflow.
