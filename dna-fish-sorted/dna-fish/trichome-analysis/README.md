# Trichome Analysis

## Overview

R script generating the trichome phenotype scoring figure for the RNAi
screen — cuticle preparation and trichome counts, Control vs. selected
transcription factor/Hox genes vs. boundary element genes, in a single
combined panel.

**`trichome_screen_control_vs_tf_boundary_genes.R`** — combined panel plot,
Control → TF/Hox genes → boundary elements, gene labels at 45°, reduced
inter-gene spacing. Wilcoxon rank-sum tests (each gene vs. Control) with BH
correction applied within each region (Dorsal, Ventral). Includes median
dotted-line overlay and summary statistics.

**Genes covered:** Control, Awh, pnr, Ubx, Blimp-1, ttk, lola, Exd (TF/Hox
group); CTCF, BEAF-32, Su(Hw), Mod(mdg4) (boundary element group).

**Regions:** Dorsal+Lateral (relabeled "Dorsal" for the plot), Ventral.


---

## How to run

```r
Rscript trichome_screen_control_vs_tf_boundary_genes.R
```

Check the data-loading section near the top of the file for the expected
input path/format before running.


---

## Origin note

This script was extracted as a standalone file from a larger combined R
Markdown notebook, which also contained an earlier draft (without BH
correction on the Wilcoxon tests) and a second copy of this same final
version (identical except for whitespace). Only this one, final version was
kept.
