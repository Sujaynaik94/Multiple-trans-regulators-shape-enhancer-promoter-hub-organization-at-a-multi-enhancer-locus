import pandas as pd
from shapely.geometry import Point, Polygon
from roifile import roiread
import os

# =======================================================
# USER PATHS
# =======================================================
batch_csv = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\ROI\ROI_TRIPLET_BATCH_UPDATED.csv"

output_base = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\ROI"
os.makedirs(output_base, exist_ok=True)

combined_output_xlsx = os.path.join(
    output_base,
    "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

# =======================================================
# LOAD BATCH TABLE
# =======================================================
batch_df = pd.read_csv(batch_csv)

required_cols = [
    "gene_id",
    "triplet_path",
    "roi_path"
]

for col in required_cols:
    if col not in batch_df.columns:
        raise ValueError(
            f"Missing column in CSV: {col}\n"
            f"Available columns: {batch_df.columns.tolist()}"
        )

print(f"\nLoaded batch table with {len(batch_df)} entries")

# =======================================================
# STORAGE
# =======================================================
combined_rows = []

# =======================================================
# POINT-IN-ROI FUNCTION  (IMPORTANT FIX)
# =======================================================
def in_roi(row, poly):
    A  = Point(row["x_svb_px"], row["y_svb_px"])
    E  = Point(row["x_E_px"],   row["y_E_px"])
    DG = Point(row["x_DG_px"],  row["y_DG_px"])
    # covers = inside OR on boundary (restores old behavior)
    return poly.covers(A) and poly.covers(E) and poly.covers(DG)

# =======================================================
# MAIN BATCH LOOP
# =======================================================
for idx, row in batch_df.iterrows():

    gene        = row["gene_id"]
    triplet_csv = row["triplet_path"]
    roi_zip     = row["roi_path"]

    print("\n====================================================")
    print(f"Processing {idx+1}/{len(batch_df)}")
    print(f"Gene        : {gene}")
    print(f"Triplets    : {triplet_csv}")
    print(f"ROI ZIP     : {roi_zip}")

    # ---------------------------------------------------
    # SAFETY CHECKS
    # ---------------------------------------------------
    if not isinstance(triplet_csv, str) or not os.path.isfile(triplet_csv):
        print("Triplet CSV missing -- skipping")
        continue

    # ---------------------------------------------------
    # LOAD TRIPLETS
    # ---------------------------------------------------
    df = pd.read_csv(triplet_csv)

    required_triplet_cols = [
        "x_svb_px","y_svb_px",
        "x_E_px","y_E_px",
        "x_DG_px","y_DG_px"
    ]

    for c in required_triplet_cols:
        if c not in df.columns:
            raise ValueError(f"{triplet_csv} -> missing column: {c}")

    print(f"Triplets loaded: {len(df)}")

    # ---------------------------------------------------
    # DEBUG: explicit diagnostic for SXLGFP specifically
    # ---------------------------------------------------
    if gene == "SXLGFP":
        print(f"  [DEBUG] Columns in df: {df.columns.tolist()}")
        if "voxel_x" in df.columns:
            print(f"  [DEBUG] voxel_x dtype: {df['voxel_x'].dtype}")
            print(f"  [DEBUG] voxel_x first 3 values: {df['voxel_x'].head(3).tolist()}")
            print(f"  [DEBUG] voxel_x non-null count: {df['voxel_x'].notna().sum()} / {len(df)}")
        else:
            print(f"  [DEBUG] *** voxel_x column NOT PRESENT in this df! ***")

    # ---------------------------------------------------
    # NO ROI CASE (e.g. nc14): keep all triplets as-is, unclassified,
    # instead of dropping them entirely
    # ---------------------------------------------------
    if not isinstance(roi_zip, str) or not os.path.isfile(roi_zip):
        print("ROI ZIP missing/empty (e.g. nc14, no ROI expected) -- keeping triplets unfiltered")

        df_unfiltered = df.copy()
        df_unfiltered["gene_id"] = gene
        df_unfiltered["roi"] = "no_ROI"
        df_unfiltered["triplet_source"] = os.path.basename(triplet_csv)

        combined_rows.append(df_unfiltered)
        print(f"  '{gene}' (no ROI): {len(df_unfiltered)} triplets kept as-is")
        continue

    # ---------------------------------------------------
    # LOAD ROIs
    # ---------------------------------------------------
    roi_objects = roiread(roi_zip)
    polygons = {}

    for roi in roi_objects:
        roi_name = roi.name if roi.name else "ROI"

        coords = [(float(x), float(y)) for x, y in roi.coordinates()]
        if coords[0] != coords[-1]:
            coords.append(coords[0])

        polygons[roi_name] = Polygon(coords)

    print(f"ROIs loaded: {list(polygons.keys())}")

    # ---------------------------------------------------
    # FILTER + COLLECT
    # ---------------------------------------------------
    for roi_name, poly in polygons.items():

        mask = df.apply(lambda r: in_roi(r, poly), axis=1)
        df_roi = df.loc[mask].copy()

        if df_roi.empty:
            continue

        df_roi["gene_id"] = gene
        df_roi["roi"] = roi_name
        df_roi["triplet_source"] = os.path.basename(triplet_csv)

        combined_rows.append(df_roi)

        print(f"  ROI '{roi_name}': {len(df_roi)} triplets")

# =======================================================
# WRITE EXCEL WITH 2 SHEETS
# =======================================================
if combined_rows:

    combined_df = pd.concat(combined_rows, ignore_index=True)

    # DEBUG: check voxel_x right before writing
    print(f"\n[DEBUG] Final combined_df columns: {combined_df.columns.tolist()}")
    if "voxel_x" in combined_df.columns:
        sxlgfp_mask = combined_df["gene_id"] == "SXLGFP"
        print(f"[DEBUG] voxel_x dtype in combined_df: {combined_df['voxel_x'].dtype}")
        print(f"[DEBUG] SXLGFP rows: {sxlgfp_mask.sum()}")
        print(f"[DEBUG] SXLGFP voxel_x non-null count: {combined_df.loc[sxlgfp_mask, 'voxel_x'].notna().sum()}")
        print(f"[DEBUG] SXLGFP voxel_x first 3 values: {combined_df.loc[sxlgfp_mask, 'voxel_x'].head(3).tolist()}")
    else:
        print(f"[DEBUG] *** voxel_x column NOT PRESENT in combined_df! ***")

    summary_df = (
        combined_df
        .groupby(["gene_id", "roi", "triplet_source"])
        .size()
        .reset_index(name="triplet_count")
        .sort_values(["gene_id", "roi", "triplet_source"])
    )

    with pd.ExcelWriter(combined_output_xlsx, engine="openpyxl") as writer:
        combined_df.to_excel(writer, sheet_name="Filtered_Triplets", index=False)
        summary_df.to_excel(writer, sheet_name="Summary_Counts", index=False)

    print("\nExcel file written:")
    print(combined_output_xlsx)
    print(f"Total filtered triplets: {len(combined_df)}")

else:
    print("\nNo triplets passed ROI filtering -- no output written")
