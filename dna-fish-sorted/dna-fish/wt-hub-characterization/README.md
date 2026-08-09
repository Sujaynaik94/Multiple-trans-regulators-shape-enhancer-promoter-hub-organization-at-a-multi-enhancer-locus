# WT Hub Characterization

## Overview

R scripts that generate the plots and statistics characterizing
enhancer-promoter hub organization in wild-type (WT/SXLGFP) nuclei — all
read directly from `ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx`, the
final output of `../roi-classification/`.

Three independent, standalone scripts:

1. **`inter_probe_distances_wt_ne_dorsal_ventral.R`** — half-violin/boxplot
   figure of the three pairwise inter-probe distances (Promoter-E6,
   Promoter-DG2, E6-DG2) across WT Non-expressing, Dorsal, and Ventral
   nuclei. Filters out pairwise distances < 0.01 µm. Wilcoxon rank-sum tests
   (Non-expressing vs Dorsal, Non-expressing vs Ventral) with BH correction,
   applied separately within each probe-pair group. Includes embryo median
   overlays and n/N labels.

2. **`hub_dispersion_wt_ne_vs_dorsal_ventral.R`** — hub dispersion
   (centroid distance score) comparison, WT Non-expressing vs Dorsal and
   Ventral. Nucleus/triplet-level Wilcoxon tests, BH-corrected across the 2
   comparisons. No manual embryo exclusions.

3. **`hub_dispersion_nc14_vs_wt_ne.R`** — hub dispersion comparison, nc14
   whole-embryo nuclei vs WT Non-expressing nuclei (developmental
   comparison). Nucleus/triplet-level Wilcoxon rank-sum test. No manual
   embryo exclusions.

Each script is fully self-contained — its own data loading, statistics, and
plotting/save steps.


---

## Shared conventions across all three scripts

- **Data source:** `ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx`, sheet
  `Filtered_Triplets`.
- **WT gene ID:** `SXLGFP` (the wild-type/control condition used for
  developmental- and domain-level comparisons; see `../../scatac-chip-reanalysis.md`-adjacent
  notes for the biological rationale).
- **Statistics:** Wilcoxon rank-sum tests at the nucleus/triplet level, with
  Benjamini-Hochberg correction applied within each script's set of planned
  comparisons — not across all three scripts jointly.
- **Plot style:** transparent half-violins, jittered individual
  nuclei/triplets, pooled boxplots, embryo median dots, n/N labels,
  significance brackets (including explicit "ns" where applicable).


---

## How to run

Each script can be run independently in RStudio or via `Rscript`:

```r
Rscript inter_probe_distances_wt_ne_dorsal_ventral.R
Rscript hub_dispersion_wt_ne_vs_dorsal_ventral.R
Rscript hub_dispersion_nc14_vs_wt_ne.R
```

Edit the `xlsx_path` (and, for the nc14 script, any nc14-specific input
path) near the top of each file before running.


---

## Origin note

These three scripts were extracted as standalone files from a larger
combined R Markdown notebook. The notebook also contained two earlier
drafts of the RNAi hub dispersion analysis (now in `../rnai-screen/`) that
are not included here, since only the final, complete version was kept.
