# RNAi Screen — Hub Dispersion Analysis

## Overview

R script generating the RNAi screen figure: hub dispersion (centroid
distance score) for each knockdown gene vs. the LacZ control, split by
Non-expressing / Dorsal / Ventral condition, plus a summary barplot of
percent change vs. LacZ. Reads directly from
`ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx`, the final output of
`../roi-classification/`.

**`hub_dispersion_rnai_screen_summary.R`** — single vertical
multi-panel plot: Dorsal, Ventral, Non-expressing panels (each showing every
screened gene vs. LacZ), followed by a summary barplot aligned to the same
x-axis positions as the panels above it (LacZ kept as a blank placeholder so
gene positions line up between the panels and the summary).

**Genes covered** (fixed plotting order used in the script): LacZ (control),
Awh, pnr, Ubx, Blimp-1, ttk, lola, Exd, bowl, hth, pnt, Med13, CTCF, BEAF-32,
Su(Hw), Mod(mdg4).


---

## How to run

```r
Rscript hub_dispersion_rnai_screen_summary.R
```

Edit `xlsx_path` near the top of the file before running.


---

## Origin note

This script was extracted as a standalone file from a larger combined R
Markdown notebook, which also contained an earlier draft of the same
analysis (missing several genes — bowl, hth, pnt, Med13 — from the plotting
order). Only this final, complete version was kept.
