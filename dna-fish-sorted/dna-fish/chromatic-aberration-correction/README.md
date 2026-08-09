# Chromatic Aberration Correction

## Overview

Chromatic aberration correction was **not** performed by a custom script.
It was performed in **ZEN** (Zeiss's acquisition/processing software) using
its built-in **Channel Alignment** feature, applied during image
processing/export on the Zeiss LSM 900 confocal.

This is what **"Channel Alignment (Extended)"**, appearing in every embryo
filename throughout this repository (e.g.
`LacZ_15_HU-Airyscan Processing-94-Channel Alignment (Extended)-114.czi`),
refers to — the correction is already applied by the time a CZI file reaches
any of the scripts in `../nuclei-segmentation/`, `../spot-detection/`, or
downstream steps. There is no separate correction step to run in this
pipeline.


---

## `channel_alignment_correction_2025-12.xml`

This is the ZEN **`ScanfieldTransformation`** registration file used for the
correction. It defines a per-channel affine transformation (translation,
rotation, scaling) that registers each fluorescence channel to a reference
view, correcting the small channel-to-channel spatial offset introduced by
chromatic aberration in the optical path.

**Channels covered** (4-channel acquisition):
| Channel | Label |
|---|---|
| View 1 (reference) | Dy635-T1 |
| View 2 | At550-T3 |
| View 3 | AF488-T4 |
| View 4 | DAPI-T5 |

**Acquisition parameters recorded in the file:**
- Objective: Plan-Apochromat 63x/1.40 Oil DIC M27
- Camera/detector: LSM 900
- Voxel size: X/Y = 0.0706 µm, Z = 0.14 µm (`<Scaling>` block)
- Scan mode: bidirectional, frame scanning, laser blanking on

**Per-view registration** (`<Registration><ViewN>`): each non-reference view
has its own affine transformation matrix and translation vector, plus an
`RMSE` value reporting the registration error (residual misalignment in nm)
for that channel — View 2: RMSE ≈ 103 nm, View 3: RMSE ≈ 57 nm, View 4:
RMSE ≈ 62 nm.

This correction file is loaded into ZEN and applied at image export time; it
is not read or applied by any script in this repository.


---

## Practical note for reproducing the pipeline

If re-processing raw microscope files from scratch, the Channel Alignment
step must be applied in ZEN using this (or an equivalently-generated)
correction file **before** exporting to CZI/TIFF for the rest of the
pipeline. Skipping this step means the promoter/DG2/E6 spot coordinates in
different channels will carry a systematic offset relative to each other,
which is exactly the kind of error the correction file above is designed to
remove.
