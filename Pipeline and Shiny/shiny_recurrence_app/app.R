## Lung cancer recurrence follow-up explorer (Shiny).
## Model backend: model_backend.R (deployed baseline formula; no metastasis score).

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(DT)
  library(scales)
})

app_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
  if (length(f) && nzchar(f[[1]])) {
    dirname(normalizePath(f[[1]], mustWork = FALSE))
  } else {
    getwd()
  }
})
source(file.path(app_dir, "model_backend.R"), local = FALSE)

options(shiny.maxRequestSize = 500 * 1024^2)

file_input_uploaded <- function(fi) {
  !is.null(fi) &&
    is.data.frame(fi) &&
    nrow(fi) > 0L &&
    "datapath" %in% names(fi) &&
    nzchar(fi$datapath[[1]] %||% "")
}

read_table_for_upload <- function(path, ext, sep_pref = "\t") {
  ext <- tolower(ext %||% tools::file_ext(path))
  sep_use <- if (ext == "csv") "," else sep_pref
  load_uploaded_table(path, ext, sep = sep_use)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal stage buckets 1–7; UI shows only the clinical labels below.
STAGE_NUM_CHOICES <- c(
  "IA, 1A, 1a, or M0" = "1",
  "IB, 1B, or 1b" = "2",
  "II, IIA, 2A, or 2a" = "3",
  "IIB, 2B, or 2b" = "4",
  "IIIA, 3A, or 3a" = "5",
  "IIIB, 3B, or 3b" = "6",
  "IV or 4" = "7"
)

stage_num_row_label <- function(k) {
  labs <- c(
    "IA / 1A / 1a / M0",
    "IB / 1B / 1b",
    "II / IIA / 2A / 2a",
    "IIB / 2B / 2b",
    "IIIA / 3A / 3a",
    "IIIB / 3B / 3b",
    "IV / 4"
  )
  ki <- as.integer(k)
  if (length(ki) != 1L || is.na(ki) || ki < 1L || ki > 7L) {
    return("Unknown")
  }
  labs[[ki]]
}

# Gene helpers (used with training artifacts from build_training_artifacts).
model_gene_sets <- function(ar) {
  if (is.null(ar)) {
    return(list(recurrence = character(0)))
  }
  rec <- unique(toupper(names(ar$recur_weights %||% character(0))))
  rec <- rec[!is.na(rec) & nzchar(rec)]
  list(recurrence = rec)
}

signature_gene_names <- function(ar) {
  model_gene_sets(ar)$recurrence
}

extract_recurrence_gene_names <- function(ar) {
  model_gene_sets(ar)$recurrence
}

render_model_genes_panel <- function(gene_cache, ar = NULL) {
  if (is.null(gene_cache)) {
    return(div(class = "small-help", "Loading the 15-gene panel..."))
  }
  rec <- gene_cache$recurrence %||% character(0)
  if (!length(rec)) {
    return(div(class = "small-help", "Gene list could not be loaded. Refresh the page."))
  }
  NULL
}

gene_input_id <- function(gene) {
  paste0("geneval_", gsub("[^A-Za-z0-9]+", "_", toupper(as.character(gene))))
}

signature_gene_inputs_ui <- function(genes, weights = NULL) {
  genes <- unique(toupper(genes[!is.na(genes) & nzchar(genes)]))
  if (!length(genes)) {
    return(div(class = "small-help", "Waiting for signature genes..."))
  }
  gene_field <- function(g) {
    numericInput(
      gene_input_id(g),
      label = g,
      value = NA,
      step = 0.01
    )
  }
  mid <- ceiling(length(genes) / 2)
  fluidRow(
    column(6, lapply(genes[seq_len(mid)], gene_field)),
    column(6, lapply(genes[seq(mid + 1, length(genes))], gene_field))
  )
}

collect_signature_gene_values <- function(input, genes) {
  genes <- unique(toupper(genes[!is.na(genes) & nzchar(genes)]))
  vals <- setNames(rep(NA_real_, length(genes)), genes)
  missing <- character(0)
  for (g in genes) {
    v <- suppressWarnings(as.numeric(input[[gene_input_id(g)]]))
    if (!is.finite(v)) {
      missing <- c(missing, g)
    } else {
      vals[g] <- v
    }
  }
  list(values = vals, missing = missing)
}

risk_label <- function(x) {
  ifelse(x >= 0.75, "Very high",
    ifelse(x >= 0.50, "High",
      ifelse(x >= 0.25, "Moderate", "Low")
    )
  )
}

followup_note <- function(x) {
  ifelse(x >= 0.75, "Worth a closer look in clinic",
    ifelse(x >= 0.50, "Consider more frequent follow-up",
      ifelse(x >= 0.25, "Standard follow-up, but stay alert", "Lower estimated risk for this window")
    )
  )
}

model_plain_english <- function() {
  tagList(
    p("The calculator combines three inputs:"),
    tags$ul(
      tags$li(strong("Cancer stage"), " — how advanced the tumour is"),
      tags$li(strong("Age"), " — in years"),
      tags$li(
        strong("Recurrence score"), " — a summary of 15 gene expression values ",
        "(TUBA1C, CHRDL1, CDC45, and 12 others)"
      )
    ),
    p(
      "Sex and metastasis-related genes were tested during model development but are ",
      strong("not"), " used in the deployed formula because they did not improve validation performance."
    )
  )
}

palette_values <- function(name) {
  palettes <- list(
    "Ocean teal" = c(
      primary = "#0d9488", dark = "#0f172a", accent = "#ccfbf1", bg = "#f0fdfa", nav = "#0f766e",
      gauge1 = "#5eead4", gauge2 = "#2dd4bf", gauge3 = "#f59e0b", gauge4 = "#dc2626",
      plot_line = "#0d9488", marker = "#0f172a", chart_fill = "#14b8a6",
      row_vh = "#ccfbf1", row_hi = "#d1fae5", row_mod = "#fef3c7", row_lo = "#e0f2fe"
    ),
    "Slate indigo" = c(
      primary = "#6366f1", dark = "#1e1b4b", accent = "#e0e7ff", bg = "#f5f5ff", nav = "#4338ca",
      gauge1 = "#a5b4fc", gauge2 = "#818cf8", gauge3 = "#f472b6", gauge4 = "#be123c",
      plot_line = "#6366f1", marker = "#1e1b4b", chart_fill = "#818cf8",
      row_vh = "#e0e7ff", row_hi = "#ede9fe", row_mod = "#fef9c3", row_lo = "#ecfccb"
    ),
    "Arctic cyan" = c(
      primary = "#0891b2", dark = "#164e63", accent = "#cffafe", bg = "#f0fdff", nav = "#0e7490",
      gauge1 = "#67e8f9", gauge2 = "#22d3ee", gauge3 = "#fbbf24", gauge4 = "#ea580c",
      plot_line = "#0891b2", marker = "#164e63", chart_fill = "#06b6d4",
      row_vh = "#cffafe", row_hi = "#e0f2fe", row_mod = "#fef3c7", row_lo = "#f1f5f9"
    ),
    "Graphite rose" = c(
      primary = "#be185d", dark = "#1c1917", accent = "#fce7f3", bg = "#fafaf9", nav = "#9d174d",
      gauge1 = "#86efac", gauge2 = "#fbbf24", gauge3 = "#fb923c", gauge4 = "#be185d",
      plot_line = "#be185d", marker = "#1c1917", chart_fill = "#db2777",
      row_vh = "#fce7f3", row_hi = "#ffedd5", row_mod = "#fef9c3", row_lo = "#ecfdf5"
    ),
    "Clinical red" = c(
      primary = "#b43c2f", dark = "#171713", accent = "#f7dfda", bg = "#f6f3ef", nav = "#b43c2f",
      gauge1 = "#2f7d5c", gauge2 = "#b39b2e", gauge3 = "#d9822b", gauge4 = "#b43c2f",
      plot_line = "#2f7d5c", marker = "#171713", chart_fill = "#b43c2f",
      row_vh = "#fde8e5", row_hi = "#fff0df", row_mod = "#fff8d9", row_lo = "#edf7f2"
    ),
    "Hospital blue" = c(
      primary = "#2f6f9f", dark = "#123047", accent = "#e1f0f7", bg = "#f5f8fa", nav = "#2f6f9f",
      gauge1 = "#5b9bd5", gauge2 = "#8ecae6", gauge3 = "#f4a261", gauge4 = "#e76f51",
      plot_line = "#2f6f9f", marker = "#123047", chart_fill = "#3a86ff",
      row_vh = "#dbeafe", row_hi = "#e0f2fe", row_mod = "#fef3c7", row_lo = "#ecfdf5"
    ),
    "Sage green" = c(
      primary = "#2f7d5c", dark = "#17352a", accent = "#e2f3ea", bg = "#f5f8f4", nav = "#2f7d5c",
      gauge1 = "#4ade80", gauge2 = "#a3e635", gauge3 = "#eab308", gauge4 = "#b45309",
      plot_line = "#2f7d5c", marker = "#17352a", chart_fill = "#16a34a",
      row_vh = "#dcfce7", row_hi = "#ecfccb", row_mod = "#fef9c3", row_lo = "#e0f2fe"
    ),
    "Warm neutral" = c(
      primary = "#8b5e3c", dark = "#27221d", accent = "#f1e5d7", bg = "#f8f5f1", nav = "#8b5e3c",
      gauge1 = "#a8a29e", gauge2 = "#d6d3d1", gauge3 = "#f59e0b", gauge4 = "#b45309",
      plot_line = "#78716c", marker = "#27221d", chart_fill = "#8b5e3c",
      row_vh = "#ffedd5", row_hi = "#fef3c7", row_mod = "#e7e5e4", row_lo = "#ecfdf5"
    ),
    "Royal purple" = c(
      primary = "#6950a1", dark = "#241c3f", accent = "#eee9fb", bg = "#f7f5fb", nav = "#6950a1",
      gauge1 = "#c4b5fd", gauge2 = "#a78bfa", gauge3 = "#f472b6", gauge4 = "#7c3aed",
      plot_line = "#6950a1", marker = "#241c3f", chart_fill = "#7c3aed",
      row_vh = "#ede9fe", row_hi = "#fae8ff", row_mod = "#fef9c3", row_lo = "#ecfccb"
    )
  )
  palettes[[name]] %||% palettes[["Ocean teal"]]
}

sanitize_hex <- function(x, default = "#0d9488") {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x) || !grepl("^#[0-9A-Fa-f]{6}$", x)) return(default)
  x
}

blend_white <- function(hex, t) {
  hex <- sanitize_hex(hex)
  r <- col2rgb(hex)[, 1]
  rgb((1 - t) * r[1] / 255 + t, (1 - t) * r[2] / 255 + t, (1 - t) * r[3] / 255 + t)
}

custom_palette_from_hex <- function(primary, dark, accent, bg, nav) {
  pr <- sanitize_hex(primary)
  dk <- sanitize_hex(dark, "#0f172a")
  ac <- sanitize_hex(accent, "#ccfbf1")
  bgg <- sanitize_hex(bg, "#f8fafc")
  nv <- sanitize_hex(nav, pr)
  list(
    primary = pr, dark = dk, accent = ac, bg = bgg, nav = nv,
    gauge1 = blend_white(pr, 0.55), gauge2 = blend_white(pr, 0.35),
    gauge3 = "#f59e0b", gauge4 = "#dc2626",
    plot_line = pr, marker = dk, chart_fill = pr,
    row_vh = ac, row_hi = blend_white(pr, 0.75), row_mod = "#fef3c7", row_lo = blend_white(pr, 0.88)
  )
}

empty_state <- function(title, body) {
  div(
    class = "empty-state",
    h4(title),
    p(body)
  )
}

value_card <- function(title, value, subtitle = NULL, class = "") {
  div(
    class = paste("metric-card", class),
    div(class = "metric-title", title),
    div(class = "metric-value", value),
    if (!is.null(subtitle)) div(class = "metric-subtitle", subtitle)
  )
}

theme_toolbar_ui <- div(
  class = "theme-toolbar",
  fluidRow(
    column(
      width = 4,
      selectInput(
        "palette",
        label = "Colour palette",
        choices = c(
          "Ocean teal", "Slate indigo", "Arctic cyan", "Graphite rose",
          "Clinical red", "Hospital blue", "Sage green", "Warm neutral", "Royal purple"
        ),
        selected = "Ocean teal",
        width = "100%"
      )
    ),
    column(
      width = 4,
      radioButtons(
        "color_mode",
        label = "Colour mode",
        choices = c("Preset" = "preset", "Custom hex" = "custom"),
        selected = "preset",
        inline = TRUE
      )
    ),
    column(
      width = 4,
      conditionalPanel(
        condition = "input.color_mode == 'custom'",
        helpText("Custom colours apply when this mode is selected."),
        textInput("hex_primary", "Primary", value = "#0d9488", placeholder = "#RRGGBB"),
        textInput("hex_nav", "Navbar", value = "#0f766e", placeholder = "#RRGGBB")
      )
    )
  ),
  conditionalPanel(
    condition = "input.color_mode == 'custom'",
    fluidRow(
      column(3, textInput("hex_dark", "Text / buttons", value = "#0f172a", placeholder = "#RRGGBB")),
      column(3, textInput("hex_accent", "Light accent", value = "#ccfbf1", placeholder = "#RRGGBB")),
      column(3, textInput("hex_bg", "Page background", value = "#f0fdfa", placeholder = "#RRGGBB")),
      column(
        3,
        helpText(style = "margin-top:28px; font-size:12px;", "Use full 6-digit hex with #.")
      )
    )
  )
)

ui <- page_navbar(
  title = tags$span(
    class = "app-brand-text",
    tags$span(class = "brand-long", "Lung Cancer Recurrence Follow-up Explorer"),
    tags$span(class = "brand-short", "Recurrence explorer")
  ),
  window_title = "Recurrence explorer",
  navbar_options = navbar_options(
    collapsible = TRUE,
    bg = "#0f766e",
    theme = "dark",
    underline = TRUE
  ),
  theme = bslib::bs_add_variables(
    bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = "#0d9488",
      bg = "#f0fdfa",
      fg = "#1e293b",
      base_font = font_google("Inter"),
      heading_font = font_google("Libre Baskerville")
    ),
    "enable-dark-mode" = FALSE
  ),
  header = tags$head(
    tags$meta(name = "color-scheme", content = "light"),
    tags$script(HTML("document.documentElement.setAttribute('data-bs-theme','light');")),
    tags$style(HTML("
      :root, html[data-bs-theme='light'] {
        color-scheme: light;
        --app-primary: #0d9488;
        --app-dark: #0f172a;
        --app-accent: #ccfbf1;
        --app-bg: #f0fdfa;
        --app-nav: #0f766e;
        --bs-body-bg: #f0fdfa;
        --bs-body-color: #1e293b;
      }
      html[data-bs-theme='dark'] {
        color-scheme: light;
        --bs-body-bg: #f0fdfa;
        --bs-body-color: #1e293b;
      }
      body {
        background: var(--app-bg) !important;
        color: #1e293b !important;
      }
      .tab-content, .sidebar, .bslib-page-fill > .tab-content {
        background-color: var(--app-bg) !important;
        color: #1e293b !important;
      }
      .form-control, .form-select, textarea {
        background-color: #fff !important;
        color: #1e293b !important;
      }
      /* Navbar: do not inherit page background on .container-fluid (was white-on-white tabs) */
      .navbar,
      .navbar .container-fluid,
      .navbar .navbar-collapse,
      .navbar .navbar-nav {
        background-color: var(--app-nav) !important;
        background: var(--app-nav) !important;
        border-color: transparent;
      }
      .navbar {
        border-bottom: 4px solid var(--app-primary);
        padding-top: 0.35rem;
        padding-bottom: 0.35rem;
      }
      .navbar > .container-fluid {
        flex-wrap: wrap;
        align-items: center;
        row-gap: 0.25rem;
      }
      .navbar-brand,
      .navbar-brand:hover,
      .navbar-brand:focus {
        color: #ffffff !important;
        font-weight: 700;
        max-width: min(100%, 42vw);
        white-space: normal;
        line-height: 1.2;
        font-size: clamp(0.95rem, 2.2vw, 1.15rem);
        margin-right: 0.5rem;
      }
      .app-brand-text,
      .app-brand-text .brand-long,
      .app-brand-text .brand-short {
        color: #ffffff !important;
      }
      .brand-long { display: none; }
      .brand-short { display: inline; }
      @media (min-width: 1100px) {
        .brand-long { display: inline; }
        .brand-short { display: none; }
        .navbar-brand { max-width: 38%; }
      }
      .navbar .nav-link,
      .navbar .navbar-nav .nav-link,
      .navbar-nav .nav-link,
      .nav-underline .nav-link {
        color: #ffffff !important;
        font-weight: 700;
        white-space: nowrap;
        padding: 0.45rem 0.75rem !important;
        opacity: 1 !important;
      }
      .navbar .nav-link:hover,
      .navbar .nav-link:focus {
        color: #ecfdf5 !important;
      }
      .navbar .nav-link.active,
      .navbar-nav .nav-link.active,
      .nav-underline .nav-link.active {
        color: #ffffff !important;
        border-bottom-color: #ffffff !important;
        border-bottom-width: 3px;
      }
      /* Menu (☰) always visible; tabs in collapsible panel at every screen width */
      .navbar .navbar-toggler,
      .navbar-expand .navbar-toggler,
      .navbar-expand-sm .navbar-toggler,
      .navbar-expand-md .navbar-toggler,
      .navbar-expand-lg .navbar-toggler,
      .navbar-expand-xl .navbar-toggler {
        display: block !important;
        visibility: visible !important;
        margin-left: auto;
        border-color: rgba(255,255,255,.75) !important;
      }
      .navbar-toggler-icon {
        filter: invert(1);
      }
      .navbar .navbar-collapse {
        flex-basis: 100%;
        width: 100%;
        background-color: var(--app-nav) !important;
      }
      .navbar .navbar-collapse:not(.show) {
        display: none !important;
      }
      .navbar-expand .navbar-collapse:not(.show),
      .navbar-expand-sm .navbar-collapse:not(.show),
      .navbar-expand-md .navbar-collapse:not(.show),
      .navbar-expand-lg .navbar-collapse:not(.show),
      .navbar-expand-xl .navbar-collapse:not(.show) {
        display: none !important;
      }
      .navbar .navbar-collapse.show {
        display: flex !important;
        flex-direction: column;
        align-items: flex-start;
      }
      .navbar .navbar-nav {
        flex-direction: column;
        align-items: flex-start;
        padding: 0.35rem 0 0.5rem 0;
        width: 100%;
      }
      .navbar .navbar-nav .nav-link {
        width: 100%;
        border-bottom: none !important;
      }
      .navbar .navbar-nav .nav-link.active {
        background: rgba(255,255,255,.18);
        border-radius: 8px;
      }
      .hero {
        background: linear-gradient(135deg, var(--app-nav) 0%, var(--app-primary) 100%);
        color: #fff;
        padding: 30px 38px;
        border-radius: 18px;
        margin: 18px 0 22px 0;
        box-shadow: 0 14px 32px rgba(0,0,0,.12);
      }
      .hero h1 { margin: 0; font-size: 34px; letter-spacing: .2px; color: #fff; }
      .hero p { margin: 10px 0 0 0; color: #ecfdf5; max-width: 1040px; font-size: 16px; }
      .panel-card {
        background: #fff;
        border: 1px solid #e8e1da;
        border-radius: 16px;
        padding: 22px;
        box-shadow: 0 10px 28px rgba(36,30,24,.06);
        margin-bottom: 18px;
      }
      .section-label {
        font-size: 12px;
        letter-spacing: 1.4px;
        text-transform: uppercase;
        color: #6f6860;
        font-weight: 800;
        margin-bottom: 14px;
      }
      .section-label::before {
        content: '';
        display: inline-block;
        width: 7px;
        height: 7px;
        background: var(--app-primary);
        border-radius: 50%;
        margin-right: 9px;
        vertical-align: middle;
      }
      .metric-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 14px;
      }
      .metric-card {
        background: #faf8f5;
        border: 1px solid #ebe3dc;
        border-radius: 14px;
        padding: 16px;
      }
      .metric-title {
        color: #736b63;
        font-size: 11px;
        letter-spacing: 1.2px;
        text-transform: uppercase;
        font-weight: 800;
      }
      .metric-value {
        color: #231f1c;
        font-size: 28px;
        font-weight: 850;
        line-height: 1.1;
        margin-top: 7px;
      }
      .metric-subtitle {
        color: #766f67;
        font-size: 13px;
        margin-top: 5px;
      }
      .main-risk .metric-value { font-size: 56px; color: var(--app-primary); }
      .risk-badge {
        display: inline-block;
        padding: 5px 10px;
        border-radius: 999px;
        background: var(--app-accent);
        color: var(--app-primary);
        font-size: 12px;
        font-weight: 850;
        letter-spacing: .6px;
        text-transform: uppercase;
      }
      .small-help { color: #756e67; font-size: 13px; }
      .note-box {
        background: var(--app-accent);
        color: #3f332d;
        padding: 14px 16px;
        border-radius: 12px;
        margin-top: 14px;
        font-size: 14px;
      }
      .mode-hint {
        font-size: 13px;
        color: #64748b;
        margin: 8px 0 12px 0;
        padding: 8px 12px;
        background: #f1f5f9;
        border-radius: 8px;
        border-left: 3px solid var(--app-primary);
      }
      .intro-steps ol { padding-left: 1.25rem; }
      .intro-steps li { margin-bottom: 0.45rem; }
      .data-format-pre {
        font-family: ui-monospace, monospace;
        font-size: 12px;
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 12px;
        overflow-x: auto;
        white-space: pre;
        margin: 8px 0 14px 0;
      }
      .empty-state {
        border: 1px dashed #d8cec5;
        background: #fbfaf8;
        border-radius: 14px;
        padding: 24px;
        color: #6f6860;
      }
      .btn-primary {
        background-color: var(--app-dark) !important;
        border-color: var(--app-dark) !important;
        padding: 11px 14px;
        font-weight: 800;
        letter-spacing: .7px;
        text-transform: uppercase;
      }
      .sidebar .btn-primary,
      .panel-card > .btn-primary,
      .layout-sidebar .btn-primary {
        width: 100%;
      }
      .gene-input-grid .form-group label {
        font-size: 13px;
        font-weight: 700;
        color: #334155;
      }
      .gene-input-grid .form-control {
        max-width: 100%;
      }
      .form-control, .form-select { border-radius: 10px; }
      table.dataTable tbody tr:hover { background-color: var(--app-accent) !important; }
      @media (max-width: 900px) {
        .metric-grid { grid-template-columns: 1fr; }
        .hero h1 { font-size: 26px; }
      }
    ")),
    uiOutput("palette_css")
  ),
  nav_panel(
    "Introduction",
    div(
      class = "hero",
      h1("Lung cancer recurrence risk calculator"),
      p("A simple tool to estimate LUAD recurrence risk for one patient or a whole cohort.")
    ),
    div(
      class = "panel-card",
      div(class = "section-label", "Two ways to use this tool"),
      tags$dl(
        tags$dt(strong("Single patient")),
        tags$dd(
          "Best when you are reviewing ", strong("one person"), ". Enter their age, stage, and 15 gene values, ",
          "then see an estimated risk percentage."
        ),
        tags$dt(strong("Cohort explorer")),
        tags$dd(
          "Best when you have ", strong("many patients"), ". Upload one spreadsheet and the app will score everyone ",
          "and rank them from highest to lowest estimated risk."
        )
      )
    ),
    div(
      class = "panel-card",
      fluidRow(
        column(
          6,
          div(
            class = "intro-steps",
            div(class = "section-label", "Single patient — steps"),
            tags$ol(
              tags$li("Open the ", strong("Single patient"), " tab."),
              tags$li("Choose how far ahead to look (for example 365 days)."),
              tags$li("Enter ", strong("age"), " and ", strong("cancer stage"), "."),
              tags$li("Fill in all ", strong("15 gene expression"), " boxes."),
              tags$li("Click ", strong("Estimate recurrence risk"), ".")
            )
          )
        ),
        column(
          6,
          div(
            class = "intro-steps",
            div(class = "section-label", "Cohort explorer — steps"),
            tags$ol(
              tags$li("Open the ", strong("Cohort explorer"), " tab."),
              tags$li(
                "Upload a ", strong("CSV or TSV"), " with one row per patient (ID, age, stage, and gene columns)."
              ),
              tags$li("Choose your follow-up window and click ", strong("Score the cohort"), "."),
              tags$li("Review the ranked table and charts — highest risk appears first.")
            )
          )
        )
      )
    ),
    div(
      class = "panel-card",
      div(class = "section-label", "What goes into the estimate"),
      model_plain_english(),
      hr(),
      div(class = "section-label", "How to read the result"),
      p(
        "You get an estimated chance of recurrence within your chosen time window. ",
        "A higher percentage means the model sees more risk for that period."
      ),
      p(
        "This is a decision-support tool only — always interpret it alongside the full clinical picture."
      )
    ),
    div(
      class = "panel-card",
      div(class = "section-label", "Appearance (optional)"),
      p(class = "small-help", "Change colours here. Use the ", strong("menu icon"), " (☰) in the top bar to open Introduction, Single patient, and Cohort explorer."),
      theme_toolbar_ui
    )
  ),
  nav_panel(
    "Single patient",
    div(
      class = "hero",
      h1("Estimate risk for one patient"),
      p("Enter clinical details and the 15 gene values, then review the summary on the right.")
    ),
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        div(class = "section-label", "Time window"),
        numericInput(
          "horizon",
          "How far ahead should we look? (days)",
          value = 365,
          min = 30,
          step = 30
        ),
        hr(),
        div(class = "section-label", "About the patient"),
        textInput("sp_id", "Patient ID (optional)", value = "Patient_1"),
        numericInput("sp_age", "Age (years)", value = 65, min = 0, max = 200, step = 1),
        selectInput(
          "sp_stage_num",
          "Cancer stage",
          choices = STAGE_NUM_CHOICES,
          selected = "3"
        ),
        hr(),
        div(class = "section-label", "Gene expression (15 genes)"),
        uiOutput("signature_gene_help"),
        div(class = "gene-input-grid", uiOutput("signature_gene_inputs")),
        actionButton("run_single", "Estimate recurrence risk", class = "btn-primary")
      ),
      div(
        class = "panel-card",
        div(class = "section-label", "Patient recurrence summary"),
        uiOutput("single_cards"),
        plotOutput("single_gauge", height = 115),
        plotOutput("single_survival", height = 290),
        uiOutput("single_interpretation")
      )
    )
  ),
  nav_panel(
    "Cohort explorer",
    div(
      class = "hero",
      h1("Compare and rank a cohort"),
      p(
        "Upload a spreadsheet with one row per patient. The app scores everyone and lists them from ",
        "highest to lowest estimated recurrence risk for your chosen follow-up window."
      )
    ),
    fluidRow(
      column(
        4,
        div(
          class = "panel-card",
          div(class = "section-label", "Upload cohort"),
          p(
            class = "small-help",
            "One file with SampleID, Age, Stage (or Stage_num), and columns for all 15 signature genes ",
            "(TUBA1C, CHRDL1, CDC45, etc.). Extra gene columns are fine."
          ),
          fileInput(
            "datafile",
            "Cohort file (CSV or TSV)",
            accept = c(".csv", ".tsv", ".txt")
          ),
          selectInput(
            "sep",
            "Column separator (TSV/TXT only; CSV uses comma)",
            choices = c("Tab" = "	", "Comma" = ",", "Semicolon" = ";"),
            selected = "	"
          ),
          numericInput(
            "cohort_horizon",
            "How far ahead should we look? (days)",
            value = 365,
            min = 30,
            step = 30
          ),
          actionButton("run_cohort", "Score the cohort", class = "btn-primary")
        )
      ),
      column(
        8,
        div(
          class = "panel-card",
          div(class = "section-label", "Results"),
          uiOutput("cohort_cards"),
          plotOutput("cohort_distribution", height = 280),
          plotOutput("cohort_rank_plot", height = 300),
          div(class = "section-label", style = "margin-top:18px;", "All patients (highest risk first)"),
          DTOutput("cohort_table_dt")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  artifacts <- reactiveVal(NULL)
  rv_single <- reactiveValues(res = NULL)
  cohort_result <- reactiveVal(NULL)
  recur_genes_cached <- reactiveVal(NULL)
  genes_boot_done <- reactiveVal(FALSE)
  signature_inputs_built <- reactiveVal(FALSE)

  load_model_genes <- function() {
    ar <- ensure_ar()
    recur_genes_cached(model_gene_sets(ar))
  }

  bootstrap_model_genes <- function() {
    if (genes_boot_done()) {
      return()
    }
    genes_boot_done(TRUE)
    tryCatch(
      {
        withProgress(
          session = session,
          message = "Loading model and gene lists...",
          value = 0.4,
          {
            load_model_genes()
          }
        )
      },
      error = function(e) {
        genes_boot_done(FALSE)
        showNotification(
          paste("Could not load model:", conditionMessage(e)),
          type = "error",
          duration = 10
        )
      }
    )
  }

  followup_days <- reactive({
    h <- suppressWarnings(as.numeric(input$cohort_horizon))
    if (is.finite(h) && h > 0) {
      return(h)
    }
    h2 <- suppressWarnings(as.numeric(input$horizon))
    if (is.finite(h2) && h2 > 0) {
      return(h2)
    }
    365
  })

  load_cohort_raw <- function() {
    if (!file_input_uploaded(input$datafile)) {
      stop("Upload a cohort file (CSV or TSV).")
    }
    ext <- tools::file_ext(input$datafile$name)
    read_table_for_upload(input$datafile$datapath, ext, input$sep)
  }

  output$signature_gene_help <- renderUI({
    render_model_genes_panel(recur_genes_cached())
  })

  output$signature_gene_inputs <- renderUI({
    gc <- recur_genes_cached()
    genes <- gc$recurrence %||% character(0)
    if (!length(genes)) {
      return(div(class = "small-help", "Loading gene boxes..."))
    }
    signature_inputs_built(TRUE)
    signature_gene_inputs_ui(genes, DEPLOYED_RECURRENCE_WEIGHTS)
  })

  active_palette <- reactive({
    if (is.null(input$color_mode) || input$color_mode != "custom") {
      return(palette_values(input$palette %||% "Ocean teal"))
    }
    custom_palette_from_hex(
      input$hex_primary %||% "#0d9488",
      input$hex_dark %||% "#0f172a",
      input$hex_accent %||% "#ccfbf1",
      input$hex_bg %||% "#f0fdfa",
      input$hex_nav %||% input$hex_primary %||% "#0d9488"
    )
  })

  output$palette_css <- renderUI({
    p <- active_palette()
    tags$style(HTML(sprintf("
      :root, html[data-bs-theme='light'] {
        --app-primary: %s;
        --app-dark: %s;
        --app-accent: %s;
        --app-bg: %s;
        --app-nav: %s;
        --bs-body-bg: %s;
        --bs-body-color: #1e293b;
      }
      body { background: %s !important; color: #1e293b !important; }
      .navbar, .navbar .container-fluid, .navbar .navbar-collapse, .navbar .navbar-nav {
        background-color: %s !important;
        background: %s !important;
      }
      .navbar-brand, .navbar .nav-link, .nav-underline .nav-link {
        color: #ffffff !important;
      }
    ", p[["primary"]], p[["dark"]], p[["accent"]], p[["bg"]], p[["nav"]], p[["bg"]], p[["bg"]],
       p[["nav"]], p[["nav"]])))
  })

  ensure_ar <- function() {
    if (!is.null(artifacts())) {
      return(artifacts())
    }
    ar <- load_deployed_artifacts()
    artifacts(ar)
    artifacts()
  }

  observeEvent(input$run_single, {
    tryCatch(
      {
        ar <- ensure_ar()
        gc <- recur_genes_cached()
        genes <- gc$recurrence %||% character(0)
        if (!length(genes)) {
          stop("Signature genes are still loading. Wait a moment and try again.")
        }
        collected <- collect_signature_gene_values(input, genes)
        if (length(collected$missing)) {
          stop(
            "Enter a numeric value for every signature gene. Missing: ",
            paste(head(collected$missing, 8), collapse = ", "),
            if (length(collected$missing) > 8L) " ..." else ""
          )
        }
        sid <- trimws(input$sp_id)
        if (!nzchar(sid)) sid <- "Patient_1"
        stg_num <- as.integer(input$sp_stage_num)
        stage_lbl <- stage_num_row_label(stg_num)
        raw <- data.frame(
          SampleID = sid,
          Age = as.numeric(input$sp_age),
          Stage_num = stg_num,
          Stage = stage_lbl,
          stringsAsFactors = FALSE
        )
        for (nm in names(collected$values)) {
          raw[[nm]] <- collected$values[[nm]]
        }
        min_g <- length(genes)
        rv_single$res <- compute_predictions(
          ar,
          raw,
          as.numeric(input$horizon),
          0.95,
          min_genes = min_g
        )
      },
      error = function(e) showNotification(conditionMessage(e), type = "error", duration = 8)
    )
  })

  observeEvent(input$run_cohort, {
    tryCatch(
      {
        ar <- ensure_ar()
        raw <- load_cohort_raw()
        res <- compute_predictions(ar, raw, followup_days(), 0.95)
        pri <- cohort_prioritization(res$df$SampleID, as_prediction_risks(res$risks))
        pri <- pri[order(pri$cohort_rank), ]
        cohort_result(list(pred = res, priority = pri))
      },
      error = function(e) showNotification(conditionMessage(e), type = "error", duration = 8)
    )
  })

  output$single_cards <- renderUI({
    r <- rv_single$res
    if (is.null(r)) {
      return(empty_state(
        "No estimate yet",
        "Fill in age, stage, and all 15 gene values on the left, then click Estimate recurrence risk."
      ))
    }
    risk <- r$risks[[1]]
    div(
      class = "metric-grid",
      value_card("Estimated recurrence risk", percent(risk, accuracy = 0.1), span(class = "risk-badge", risk_label(risk)), "main-risk"),
      value_card("Time window", paste0(input$horizon, " days"), "From baseline to your chosen horizon"),
      value_card("Suggested follow-up", risk_label(risk), followup_note(risk))
    )
  })

  output$single_gauge <- renderPlot({
    r <- rv_single$res
    req(r)
    pal <- active_palette()
    risk <- r$risks[[1]]
    ggplot(data.frame(x = 0:100, y = 1), aes(x, y)) +
      geom_segment(aes(x = 0, xend = 25, y = 1, yend = 1), linewidth = 9, lineend = "round", colour = pal[["gauge1"]]) +
      geom_segment(aes(x = 25, xend = 50, y = 1, yend = 1), linewidth = 9, lineend = "round", colour = pal[["gauge2"]]) +
      geom_segment(aes(x = 50, xend = 75, y = 1, yend = 1), linewidth = 9, lineend = "round", colour = pal[["gauge3"]]) +
      geom_segment(aes(x = 75, xend = 100, y = 1, yend = 1), linewidth = 9, lineend = "round", colour = pal[["gauge4"]]) +
      geom_point(aes(x = 100 * risk, y = 1), size = 6, colour = pal[["marker"]]) +
      annotate("text", x = c(0, 25, 50, 75, 100), y = .78, label = c("0%", "25%", "50%", "75%", "100%"), size = 3.5, colour = "#6f6860") +
      annotate("text", x = 100 * risk, y = 1.25, label = paste("Patient", percent(risk, accuracy = 0.1)), fontface = "bold", size = 4) +
      coord_cartesian(xlim = c(-3, 103), ylim = c(.6, 1.4), clip = "off") +
      theme_void()
  })

  output$single_survival <- renderPlot({
    r <- rv_single$res
    req(r)
    pal <- active_palette()
    risk <- r$risks[[1]]
    horizon <- as.numeric(input$horizon)
    df <- data.frame(
      days = c(0, horizon),
      survival = c(1, 1 - risk)
    )
    ggplot(df, aes(days, survival)) +
      geom_step(linewidth = 1.2, colour = pal[["plot_line"]]) +
      geom_point(size = 3, colour = pal[["marker"]]) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_x_continuous(labels = comma) +
      labs(
        title = "Estimated chance of staying recurrence-free",
        subtitle = "From baseline (day 0) to your chosen follow-up window",
        x = "Days from baseline",
        y = "Recurrence-free probability"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$single_interpretation <- renderUI({
    r <- rv_single$res
    if (is.null(r)) {
      return(NULL)
    }
    risk <- r$risks[[1]]
    stg_disp <- if ("Stage_num" %in% names(r$df)) {
      stage_num_row_label(r$df$Stage_num[[1]])
    } else {
      as.character(r$df$Stage[[1]])
    }
    div(
      h4("Summary"),
      p(
        strong("Patient: "),
        paste0(r$df$SampleID[[1]], " | Age ", r$df$Age[[1]], " | Stage ", stg_disp)
      ),
      p(
        strong("Estimated recurrence risk: "),
        percent(risk, accuracy = 0.1),
        " over the next ", input$horizon, " days."
      ),
      if ("Recurrence_score" %in% names(r$df)) {
        p(
          strong("Recurrence score: "),
          sprintf("%.3f", r$df$Recurrence_score[[1]]),
          " (from ", r$matched_recur, " of 15 signature genes)."
        )
      },
      div(
        class = "note-box",
        followup_note(risk), "."
      )
    )
  })

  output$cohort_cards <- renderUI({
    cr <- cohort_result()
    if (is.null(cr)) {
      return(empty_state("No cohort loaded yet", "Upload a file and click Score the cohort."))
    }
    r <- cr$pred
    div(
      class = "metric-grid",
      value_card("Patients", nrow(r$df), "Samples included"),
      value_card(
        "Mean estimated risk",
        percent(mean(as_prediction_risks(r$risks)), accuracy = 0.1),
        paste0("For the selected follow-up length (", followup_days(), " days)")
      ),
      value_card("Signature genes found", paste0(r$matched_recur, " / 15"), "Per patient upload")
    )
  })

  output$cohort_distribution <- renderPlot({
    cr <- cohort_result()
    req(cr)
    pal <- active_palette()
    z <- cr$priority
    ggplot(z, aes(recurrence_risk)) +
      geom_histogram(bins = 20, colour = "white", fill = pal[["chart_fill"]], alpha = .85) +
      geom_vline(xintercept = mean(z$recurrence_risk), linewidth = 1, linetype = "dashed") +
      scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(title = "Distribution of estimated recurrence risk", x = "Estimated recurrence risk", y = "Number of patients") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$cohort_rank_plot <- renderPlot({
    cr <- cohort_result()
    req(cr)
    pal <- active_palette()
    z <- cr$priority
    z$rank_order <- seq_len(nrow(z))
    ggplot(z, aes(rank_order, recurrence_risk)) +
      geom_col(fill = pal[["chart_fill"]], alpha = .85) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(z$recurrence_risk, na.rm = TRUE) * 1.08)) +
      labs(title = "Patients ranked from highest to lowest estimated risk", x = "Cohort rank", y = "Estimated recurrence risk") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$cohort_table_dt <- renderDT({
    cr <- cohort_result()
    if (is.null(cr)) {
      return(datatable(data.frame(Message = "Upload cohort data and run the model."), options = list(dom = "t")))
    }
    z <- cr$priority
    z$recurrence_risk_pct <- round(100 * z$recurrence_risk, 2)
    z <- z[, c("SampleID", "cohort_rank", "recurrence_risk_pct", "percentile_urgency", "z_vs_cohort_mean", "priority_label")]
    names(z) <- c("Patient", "Rank", "Estimated risk, %", "Relative urgency percentile", "Z-score vs cohort", "Priority group")
    pc <- active_palette()
    datatable(
      z,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 10, scrollX = TRUE)
    ) %>%
      formatStyle(
        "Priority group",
        target = "row",
        backgroundColor = styleEqual(
          c("Very high", "High", "Moderate", "Lower priority"),
          c(pc[["row_vh"]], pc[["row_hi"]], pc[["row_mod"]], pc[["row_lo"]])
        )
      )
  })

  observe({
    bootstrap_model_genes()
  })
}

shinyApp(ui, server)
