# svb-enhancer-hub

Code accompanying Naik et al., "Multiple trans-regulators shape enhancer-promoter hub organization at a multi-enhancer locus" (manuscript in preparation / submitted to *Nature Communications*).

This repository contains the analysis pipelines and figure-generation scripts used in the paper. Raw sequencing and imaging data are deposited separately (see **Data availability** below).

## Repository structure

```
svb-enhancer-hub/
├── umi4c/                       # UMI-4C promoter-enhancer contact profile analysis
├── dna-fish/
│   ├── chromatic-aberration-correction/ # ZEN registration file (not a script) for channel-to-channel correction
│   ├── nuclei-segmentation/      # Cellpose-based nuclei segmentation from raw DNA-FISH images
│   ├── spot-detection/            # RS-FISH Fiji/ImageJ macro for spot detection on raw images
│   ├── triplet-matching/         # Nucleus-restricted, voxel-corrected triplet matching (full screen + SXL-GFP)
│   ├── roi-classification/       # Dorsal / Ventral / Non-expressing region assignment → final summary output
│   ├── shared/
│   │   └── plotted-stats/        # KDE contour, triangle reconstruction, Simes p-values (PLOTTED-style),
│   │                              # one script per supplementary figure
│   ├── wt-hub-characterization/  # Pairwise distance analysis and hub dispersion (CDS)
│   │                              # in wild-type nuclei (NE/Dorsal/Ventral)
│   ├── rnai-screen/               # RNAi screen hub dispersion panels + summary barplot
│   └── trichome-analysis/         # Cuticle/trichome phenotype scoring, Control vs TF/Hox vs boundary genes
├── microc-hic/                  # Micro-C / Hi-C processing: ICE normalization, obs/exp transform
├── scatac-chip-reanalysis/      # scATAC-seq and ChIP-seq reanalysis at the svb locus
├── figures/                     # Scripts that generate main and supplementary figures,
│                                 # calling into the pipelines above
├── environment.yml              # Conda environment (Python dependencies)
├── r_packages.txt               # R package list (version-pinned where relevant)
└── .gitignore
```

Each pipeline folder has its own README with parameters and run order.

## Pipelines

- **DNA-FISH** (`dna-fish/`) — 3D DNA-FISH triplet matching from imaging data. Every script here was confirmed by tracing actual input/output file paths between scripts, not by filename alone — only the scripts that actually feed into the final output are included.
  - `chromatic-aberration-correction/` — not a script. `channel_alignment_correction_2025-12.xml`, the ZEN (Zeiss) `ScanfieldTransformation` registration file used to correct channel-to-channel chromatic shift during image export — this is what "Channel Alignment (Extended)" in every embryo filename refers to. Applied in ZEN before any script in this repo touches the data. See the folder's own README.
  - `nuclei-segmentation/` — `run_segmentation_batch.py` (loops over CZI files in a folder, calls the segmentation script per file) and `nuclei_segmentation.py` (auto-selects the DAPI channel from CZI/OME-TIFF metadata, applies Gaussian smoothing, runs Cellpose 3D nuclei segmentation on a resized volume, rescales the resulting mask back to native resolution, and saves it as `{name}_mask.tif`). Runs after chromatic aberration correction, before spot detection.
  - `spot-detection/` — RS-FISH Fiji/ImageJ macro (`rsfish_spot_detection_macro.ijm`, marked "FINAL WORKING VERSION" in its own header) with `rsfish_parameters_all.csv` — per-embryo, per-channel detection parameters (SigmaDoG, ThresholdDoG, intensity range, final detection count) for 229 embryos across 15 genes (BEAF32, Blimp1, CTCF, EXD, Hth, LOLAJ, LacZ, Med13, Pnr, SUHW, Ttk, UBX, cp190, mdg4, pnt). Compiled from per-embryo `_RSFISH_LOG.txt` files (both standalone per-channel and combined multi-channel log formats), with conflicting duplicate values auto-resolved in favor of the standalone per-channel file when both existed for the same embryo/channel. `Awh` and `bowl` are not yet included — their logs weren't found in the same location; confirm where their spot-detection output lives before treating this as complete for the full screen. Detects fluorescent spots per channel, producing the CSVs `triplet-matching/` consumes.
  - `triplet-matching/` — `triplet_matching_full_screen.py` (main RNAi screen) and `triplet_matching_sxlgfp.py` (SXL-GFP male embryo branch, run in parallel). Both apply per-embryo voxel size correction from CZI metadata (fixed a ~30% distance underestimation) and output `ALL_GENES_FINAL_TRIPLETS_FIXED_VOXEL.csv`. See the folder's own README for full pipeline documentation.
  - `roi-classification/` — `update_roi_mapping.py` reads `ALL_GENES_FINAL_TRIPLETS_FIXED_VOXEL.csv` directly (confirmed by exact path match) and updates the Dorsal/Ventral/Non-expressing ROI mapping; **requires a pre-existing `ROI_TRIPLET_BATCH.csv`** as a starting input (a one-time historical mapping build, not regenerated by any script in this repo — treat it as a required data file, not something to rerun). It must have at minimum these three columns:
    | Column | Description |
    |---|---|
    | `gene_id` | Gene/condition label for the embryo (e.g. `LacZ`, `Awh`, `SXLGFP`) |
    | `triplet_path` | Path to that embryo's `..._TRIPLETS.csv` file (basename, minus the `_TRIPLETS` suffix, is used to derive the embryo name) |
    | `roi_path` | Path to that embryo's ROI `.zip` file (Dorsal/Ventral/Non-expressing polygon annotations) |

    `roi_filter_batch.py` then reads the updated mapping and produces the final `ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx`.
  - `shared/plotted-stats/` — three self-contained scripts, each reading `ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx` directly and generating one supplementary figure's difference contour plots (KDE density estimation over the 3 pairwise distances, triangle reconstruction with promoter fixed at origin, subtraction contours, and Simes-combined KS test significance): `contour_wt_nonexpressing_vs_dorsal_ventral.py` (WT Non-expressing vs Dorsal/Ventral), `contour_nc14_vs_wt_nonexpressing.py` (nc14 whole embryo vs WT Non-expressing), and `contour_rnai_screen_vs_lacz.py` (RNAi screen — each knockdown gene vs LacZ control, split by Non-expressing/Dorsal/Ventral).
  - `wt-hub-characterization/` — three self-contained R scripts, each reading `ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx` directly: `inter_probe_distances_wt_ne_dorsal_ventral.R` (half-violin plots of the 3 pairwise distances), `hub_dispersion_wt_ne_vs_dorsal_ventral.R` (hub dispersion, Non-expressing vs Dorsal/Ventral), and `hub_dispersion_nc14_vs_wt_ne.R` (hub dispersion, nc14 whole embryo vs Non-expressing). See the folder's own README.
  - `rnai-screen/` — `hub_dispersion_rnai_screen_summary.R`: hub dispersion for each of 15 screened genes (Awh, Pnr, Ubx, Blimp-1, ttk, lola, Exd, bowl, pnt, hth, Med13, CTCF, BEAF-32, Su(Hw), Mod(mdg4)) vs LacZ control, split by Non-expressing/Dorsal/Ventral, plus a summary barplot. See the folder's own README.
  - `trichome-analysis/` — `trichome_screen_control_vs_tf_boundary_genes.R`: cuticle/trichome phenotype scoring, Control vs TF/Hox genes vs boundary element genes, Dorsal and Ventral regions. See the folder's own README.
- **Micro-C / Hi-C** — ICE normalization (hicExplorer, `--filterThreshold -1.5 5`) followed by obs/exp transformation and log1p scaling; visualized with pyGenomeTracks across 8 resolutions (50–10,000 bp).
- **scATAC-seq / ChIP-seq reanalysis** — CPM-normalized pseudobulk accessibility tracks at the *svb* locus (dm6, chrX:4,935,000–5,090,000) from Calderon et al. 2022 (GSE190130), integrated with H3K27ac/H3K27me3 ChIP-seq from Gonzaga-Saavedra et al. 2025 (GSE299311).
- **UMI-4C** — Promoter-enhancer contact profiles from UMI-4C viewpoints at the *svb* promoter and known enhancers (DG2, E6, and others).

## Requirements

See `environment.yml` for Python dependencies and `r_packages.txt` for R dependencies. Both Python and R are used across pipelines (Python for DNA-FISH/scATAC processing and plotting, R for statistical figure generation).

```bash
conda env create -f environment.yml
conda activate svb-enhancer-hub
```

## Data availability

- Public sequencing datasets used for reanalysis: GSE190130 (Calderon et al. 2022), GSE299311 (Gonzaga-Saavedra et al. 2025), GSE202018, GSE83851, GSE180376, GSE41354.
- Raw DNA-FISH imaging data and newly generated sequencing data: accession numbers to be added upon deposition.

## Citation

If you use this code, please cite:
> Naik et al. (in preparation). Multiple trans-regulators shape enhancer-promoter hub organization at a multi-enhancer locus.

## License

MIT (see `LICENSE`).
