# LUAD recurrence project

```bash
cd "recur inter"
bash pipeline/run.sh
```

Step **0a** downloads metastasis GEO data, then builds frozen gene lists.  
Step **0b** downloads recurrence + GSE68465.  
Steps **1–3** write all results to **`model result/`**

Pipeline code: [pipeline/README.md](pipeline/README.md)

Shiny app: [shiny_recurrence_app/](shiny_recurrence_app/)
