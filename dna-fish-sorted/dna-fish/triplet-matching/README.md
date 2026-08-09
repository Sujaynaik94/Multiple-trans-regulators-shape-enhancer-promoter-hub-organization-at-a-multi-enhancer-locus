# Nucleus-Restricted, Voxel-Corrected Triplet Matching

## Table of Contents
1. Overview
2. Files in this pipeline
3. Folder structure assumptions
4. Parameters — full reference
5. Pipeline steps — detailed walkthrough
6. The exact 3-way matching algorithm — how it works
7. Memory management
8. Validated results
9. How to run
10. Output file formats
11. Known limitations / open items
12. Troubleshooting


---

## 1. Overview

This pipeline assigns RS-FISH-detected spots (3 imaging channels: promoter,
enhancer "E", and boundary element "DG") to Cellpose-segmented nuclei, then
finds valid 3-color "triplets" — one spot per channel, all three mutually
within a distance cutoff — restricted to spots within the same nucleus. It
runs after nuclei segmentation and spot detection, and its output feeds
directly into ROI classification (`../roi-classification/`).

Two scripts, covering two branches of the RNAi screen:

- **`triplet_matching_full_screen.py`** — the main screen. Loops over gene
  folders, applies the same validated per-nucleus matching to every embryo,
  and reads each embryo's true voxel size from its own CZI metadata rather
  than assuming a fixed value (see Section 4).

- **`triplet_matching_sxlgfp.py`** — the same pipeline, applied to the
  SXL-GFP male-embryo branch (a separate imaging/collection set), run
  alongside the main screen.

Both scripts implement the identical core per-nucleus matching logic
(`match_triplets_exact()`, defined inline in each — there is no separate
module file).


---

## 2. Files in this pipeline

| File | Purpose |
|---|---|
| `triplet_matching_full_screen.py` | Main RNAi screen, all gene folders, per-embryo voxel correction from CZI metadata |
| `triplet_matching_sxlgfp.py` | SXL-GFP male-embryo branch, same matching logic and voxel correction |

**Known caveat:** `triplet_matching_full_screen.py` currently has
`GENE_FOLDERS_FILTER` set to a specific list of 7 gene folders (the ones
confirmed to have voxel/metadata mismatches), not `None`. Set it to `None`
before running to process every gene folder. See Section 11.


---

## 3. Folder structure assumptions

- **Masks**: organized in per-gene subfolders under a common root, e.g.:
  ```
  Segmentation_masks\
    LACZ\
      LacZ_1_HU-Airyscan Processing-96-Channel Alignment (Extended)-115_mask.tif
      ...
    AWH_MASK\
      AWH_2_HU-Airyscan Processing-08-Channel Alignment (Extended)-10_mask.tif
      awh_1_HR-Airyscan Processing-07-Channel Alignment (Extended)-09_mask.tif
      ...
  ```
  Gene subfolder naming is NOT fully consistent (some have a `_MASK` suffix,
  some don't). The script does not assume a fixed subfolder naming pattern;
  it treats every subfolder under the masks root as a gene folder and strips
  a trailing `_MASK`/`_mask` suffix (if present) to derive a clean gene label.

  Filename casing is inconsistent between embryos of the same gene (e.g.
  `awh_1_HR` lowercase vs `AWH_2_HU` uppercase) — but the mask and its 3
  matching spot CSVs always use the *same* casing as each other for any given
  embryo. Matching is done by exact base-filename string match.

- **Spot-detection CSVs**: all genes' CSVs are mixed together in one flat
  folder (`AP_files\spotdetetcetion\`), NOT organized by gene. Each embryo's
  3 channel CSVs are matched to its mask via the base filename:

  ```
  Mask:  AWH_2_HU-Airyscan Processing-08-Channel Alignment (Extended)-10_mask.tif
  CH1:   AWH_2_HU-Airyscan Processing-08-Channel Alignment (Extended)-10_CH1_633nm_RSFISH.csv
  CH2:   AWH_2_HU-Airyscan Processing-08-Channel Alignment (Extended)-10_CH2_555nm_RSFISH.csv
  CH3:   AWH_2_HU-Airyscan Processing-08-Channel Alignment (Extended)-10_CH3_488nm_RSFISH.csv
  ```

  Base name = mask filename with the `_mask.tif` suffix removed.

- **Raw CZI files** (`ap_files_root`): indexed separately so each embryo's
  true voxel size can be read from its own metadata. Folders listed in
  `SKIP_FOLDERS` (e.g. output/spot-detection/QC folders) are excluded when
  indexing.


---

## 4. Parameters — full reference

| Parameter | Value | Rationale |
|---|---|---|
| `voxel_size` (per embryo) | read from each embryo's own CZI metadata (`Scaling\|Distance\|Value`) | Per-embryo correction — fixes a ~30% distance underestimation found in a subset of wild-type embryos with mismatched/downsampled voxel metadata |
| `FALLBACK_VOXEL` | X/Y = 0.0495 µm, Z = 0.14 µm | Used only if an embryo's own CZI can't be found/read; flagged clearly in `VOXEL_SIZE_USED_PER_EMBRYO_FULLSCREEN.csv` — should be rare |
| `MAX_DIST` | 1.2 µm | Distance cutoff for accepting a triplet |
| `TOP_N` | 8000 | Fixed intensity-based filter, keeps brightest N spots per channel. Fixed absolute count, not adaptive per embryo — embryo-to-embryo variability in raw spot count is expected, and an adaptive threshold risked over-fitting the cutoff to each embryo's individual noise profile |
| `MAX_TRIPLETS_PER_NUCLEUS` | 2 | Diploid cap — nuclei with more than 2 valid triplets are dropped ENTIRELY, since Drosophila somatic nuclei are diploid and >2 real triplets indicates noise crowding or ambiguous data |
| `candidate_limit` (internal, in `match_triplets_exact`) | 200 | If a nucleus has more than 200 geometrically-valid candidate triangles, exact backtracking is skipped in favor of a faster Hungarian-based approximation, and the nucleus is flagged (`exact_match = False`) for manual review |

### Channel-to-role mapping
- `CH3_488nm` → `prom` (promoter/svb)
- `CH2_555nm` → `E` (enhancer)
- `CH1_633nm` → `DG` (boundary element)


---

## 5. Pipeline steps — detailed walkthrough

### Step 1 — Load mask
The Cellpose 3D label mask (`shape = Z, Y, X`, integer labels, 0 = background)
is opened using memory-mapped reading (see Section 7) rather than loading the
full array into RAM.

### Step 2 — Determine voxel size for this embryo
The embryo's raw CZI file is located and its true voxel size read from
metadata. If it can't be found or read, `FALLBACK_VOXEL` is used and the
embryo is flagged in the voxel-size log.

### Step 3 — Load spots, apply intensity filter, assign to nuclei
For each of the 3 channels:
1. Read the RS-FISH CSV (`x, y, z, intensity` columns, pixel units)
2. Sort by `intensity` descending, keep only the top `TOP_N` (8000) spots
3. Round `(x, y, z)` to the nearest voxel index (clipped to valid mask bounds)
4. Look up `nucleus_ID = labels[z, y, x]` for every spot
5. Drop every spot with `nucleus_ID == 0` (background) — this both restricts
   matching to real nuclei and removes most false-positive detections in one
   step
6. Print spot survival counts at each stage

**Coordinate mapping validation:** the unflipped/unswapped `labels[z,y,x]`
indexing was tested against X-flip, Y-flip, X+Y-flip (180° rotation), and a
full X/Y axis swap — on two independent embryos with two independent masks.
The unmodified mapping won decisively in every test, confirming it is correct,
despite CZI metadata showing `X|AxisOrientation=-1`, `Y|AxisOrientation=-1`,
and `ImageFlip=ExchangeXY` fields that might have suggested a transform was
needed.

### Step 4 — Per-nucleus exact 3-way matching
For every nucleus with at least one surviving spot in all 3 channels:
1. Convert that nucleus's spot coordinates to microns (using this embryo's
   own voxel size from Step 2)
2. Run `match_triplets_exact()` (see Section 6)
3. Record every accepted triplet: pixel + micron coordinates for all 3 spots,
   all 3 pairwise distances, total distance, nucleus ID, and whether the match
   was found exactly or via fallback

### Step 5 — Diploid cap filter
1. Count triplets per nucleus
2. Drop any nucleus with more than `MAX_TRIPLETS_PER_NUCLEUS` (2) triplets
   entirely, since exceeding the diploid expectation is treated as a signal
   that the whole nucleus's data is unreliable rather than something to
   selectively prune
3. Log which nuclei were dropped and their triplet counts

### Step 6 — Save output
One master CSV across all processed embryos/genes, plus a skip-log CSV and a
per-embryo voxel-size log.


---

## 6. The exact 3-way matching algorithm — how it works

1. **Compute all pairwise distances** between every prom/E/DG spot within the
   nucleus (3 distance matrices: prom↔E, prom↔DG, E↔DG).

2. **Build a candidate list**: only keep `(a, b, c)` triangles where all three
   pairwise distances are already ≤ `MAX_DIST`. This is normally a very short
   list — even when raw per-channel spot counts are moderate — because the
   triple-distance constraint is restrictive.

3. **If the candidate list is small** (≤ `candidate_limit`, default 200):
   solve exactly via backtracking — search over all ways of selecting a
   maximum number of mutually disjoint candidate triangles (no spot reused
   across two different triplets), breaking ties by minimum total distance.

4. **If the candidate list is unexpectedly large** (rare — usually a crowded
   or oversegmented nucleus): falls back to an independent-Hungarian
   approximation, and flags the nucleus (`exact_match = False`) in the output
   for manual review.

**Validation performed:** tested on synthetic data with 2 known "real" allele
triplets plus random noise spots — the algorithm correctly recovered both real
triplets and ignored the noise. Also stress-tested with 15 densely-packed
random spots per channel (worst case for candidate-list size) to confirm the
fallback path engages correctly and completes quickly.


---

## 7. Memory management

Segmentation masks in this dataset are large (2.5–5 GB each as `int32`
arrays). Loading the full mask into RAM for every embryo in sequence risked
multi-GB memory allocation failures. The fix: `tifffile.memmap()`, which
returns an array-like object backed directly by the file on disk — indexing
into it (`labels[z, y, x]`) only pages in the small regions actually touched
by the lookup, rather than requiring one huge contiguous allocation upfront.
This was verified to work correctly (indexed values match a full in-memory
read) on a synthetic test file before being adopted. `del` + `gc.collect()`
between embryos is used as a light-touch additional safety measure.


---

## 8. Validated results

### Single embryo (LacZ_15_HU-114)
| Stage | Count |
|---|---|
| Raw spots (prom / E / DG) | 18,296 / 23,473 / 55,693 |
| After TOP_N=8000 filter | 8,000 / 8,000 / 8,000 |
| After nucleus assignment (background dropped) | 1,472 / 1,979 / 2,285 |
| Nuclei with signal in all 3 channels | 399 |
| Triplets before diploid cap | 222 |
| Nuclei dropped by diploid cap | 4 (had 3–4 triplets each) |
| **Final triplets** | **209** (from 193 nuclei: 177 with 1 triplet, 16 with 2) |

### Sanity checks performed
- **Distance distributions**: median pairwise distances in accepted triplets
  are ~0.5–0.6 µm, well below the 1.2 µm cutoff — confirms genuine tight
  clustering, not a loose threshold admitting weak/coincidental matches.
- **Coordinate mapping**: cross-validated against 4 flip/swap alternatives on
  2 independent embryos; unmodified mapping won every time.
- **TOP_N intensity filter impact**: reduced max triplets-per-nucleus from up
  to 9 (before filter) down to ≤4 (after filter, before the diploid cap) —
  strong evidence that most of the "crowding" was low-intensity background
  noise, not real biology or oversegmentation.
- **Visual QC in Fiji**: matched triplet points exported as point ROIs and
  overlaid on the nucleus mask. Points land inside real nuclei and form
  sensible tight clusters. (This is a basic spot-check of nucleus
  assignment/matching correctness only — separate from the downstream
  Dorsal/Ventral/Non-expressing region classification step in
  `../roi-classification/`.)


---

## 9. How to run

Edit the `USER INPUTS` block at the top of the script:

```python
output_dir = r"...\OUTPUT_FULL_SCREEN_FIXED_VOXEL"
masks_root = r"...\Segmentation_masks"
spots_root = r"...\spotdetetcetion"
ap_files_root = r"...\AP_files"
GENE_FOLDERS_FILTER = None   # None = process every gene folder
```

Then:
```cmd
python triplet_matching_full_screen.py
```

(`triplet_matching_sxlgfp.py` follows the same pattern for the SXL-GFP
branch.)

The script will:
1. List every subfolder under `masks_root` (each = one gene), unless
   restricted by `GENE_FOLDERS_FILTER`
2. Find every `*_mask.tif` in each gene folder
3. Locate the 3 matching spot CSVs by base filename in `spots_root`
4. Locate and read this embryo's own CZI file for its true voxel size
5. Skip (and log) any embryo with missing CSVs, an unexpected filename, or an
   unreadable CZI (falls back to `FALLBACK_VOXEL`, flagged in the log)
6. Run the full validated pipeline on every valid embryo
7. Save one master CSV, one skip-log CSV, and one voxel-size-used log


---

## 10. Output file formats

### Master triplets CSV
| Column | Description |
|---|---|
| `gene` | Gene label |
| `embryo` | Embryo base name |
| `nucleus_id` | Cellpose label ID |
| `exact_match` | `True` = solved by exact backtracking; `False` = fallback used, recommend manual review |
| `x/y/z_svb_px`, `x/y/z_E_px`, `x/y/z_DG_px` | Pixel coordinates per spot |
| `x/y/z_svb_um`, `x/y/z_E_um`, `x/y/z_DG_um` | Micron coordinates per spot |
| `dist_svb_E`, `dist_svb_DG`, `dist_E_DG` | Pairwise distances (µm) |
| `total_dist` | Sum of the 3 pairwise distances |

### Skip log CSV
| Column | Description |
|---|---|
| `gene` | Gene label |
| `file` | Embryo base name (or raw filename if pattern didn't match) |
| `reason` | Why it was skipped: missing CSV(s), unexpected filename casing, or a caught processing error |

### Voxel-size-used log CSV
| Column | Description |
|---|---|
| `embryo` | Embryo base name |
| `voxel_x`, `voxel_y`, `voxel_z` | Voxel size actually used (µm) |
| `source` | Whether read from this embryo's own CZI, or `FALLBACK_VOXEL` |


---

## 11. Known limitations / open items

- `triplet_matching_full_screen.py` currently has `GENE_FOLDERS_FILTER`
  set to the 7 genes with confirmed voxel mismatches (`Beaf32_MASKS,
  CTCF_MASKS, Hth_MASKS, LACZ, SUHW_MASKS, Ubx_MASKS, cp190_MASKS`), not
  `None`. Confirm whether this has since been re-run against all gene
  folders, and whether the resulting `ALL_GENES_FINAL_TRIPLETS_FIXED_VOXEL.csv`
  is the complete, final dataset for every gene in the screen, before
  treating it as such.
- `TOP_N=8000` is a fixed absolute count, applied identically regardless of
  each embryo's raw spot count. Intentional (see Section 4), but means very
  sparse embryos may have little/no filtering applied while very dense
  embryos are filtered hard.
- The diploid cap (2) assumes standard diploid somatic nuclei; would need
  reconsideration if the screen includes known polyploid/multinucleate
  regions.
- Visual QC has been done for one test embryo (LacZ_15_HU-114) only.


---

## 12. Troubleshooting

**`OSError: [Errno 22] Invalid argument` when reading a mask** — transient
network-drive read issue. The script automatically retries once after a 2
second pause. If it persists, try copying the mask locally first.

**`MemoryError: Unable to allocate X GiB`** — should be resolved by the
memory-mapped reading approach (Section 7). If it recurs, check actual
available system RAM during the run — other applications competing for
memory can still cause issues even with memory-mapped file access.

**`PermissionError` when saving to a network drive** — seen intermittently
when a file is already open elsewhere. Save to a local folder first, then
copy results to the network drive manually once confirmed correct.

**A nucleus shows `exact_match = False`** — this nucleus had more than 200
geometrically-valid candidate triangles, triggering the Hungarian-fallback
path instead of exact backtracking. Recommend manually inspecting these
nuclei (likely oversegmentation or unusually dense spot detection) before
trusting their triplets.

**An embryo's voxel size came from `FALLBACK_VOXEL`** — check the voxel-size
log for that embryo. This means its CZI file couldn't be located or read;
confirm the fallback value is actually appropriate for that embryo before
trusting its distances.
