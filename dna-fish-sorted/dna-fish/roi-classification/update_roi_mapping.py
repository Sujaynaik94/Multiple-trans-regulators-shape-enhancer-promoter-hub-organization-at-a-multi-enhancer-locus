import os
import pandas as pd

# ==========================================================
# USER INPUTS
# ==========================================================

new_triplets_csv = r"F:\RSFISH\script for nuclei assignment_july2026\OUTPUT_FULL_SCREEN_FIXED_VOXEL\ALL_GENES_FINAL_TRIPLETS_FIXED_VOXEL.csv"
split_output_dir = r"F:\RSFISH\script for nuclei assignment_july2026\OUTPUT_FULL_SCREEN_FIXED_VOXEL\SPLIT_BY_EMBRYO"

old_mapping_csv = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\ROI\ROI_TRIPLET_BATCH.csv"

# Writes to a NEW file for review -- does NOT overwrite ROI_TRIPLET_BATCH.csv.
# Review the summary, then manually rename/replace once you're confident.
output_mapping_csv = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\ROI\ROI_TRIPLET_BATCH_UPDATED.csv"

os.makedirs(split_output_dir, exist_ok=True)


# ==========================================================
# STEP 1: SPLIT THE CORRECTED FULL-SCREEN TRIPLETS BY EMBRYO
# ==========================================================

df_new = pd.read_csv(new_triplets_csv)

if "embryo" not in df_new.columns or "gene" not in df_new.columns:
    raise ValueError(f"{new_triplets_csv} missing 'embryo' or 'gene' column")

print(f"Loaded {len(df_new)} corrected triplets across {df_new['gene'].nunique()} gene(s), "
      f"{df_new['embryo'].nunique()} embryo(s)\n")

embryo_to_gene = {}      # embryo -> corrected gene label
embryo_to_new_path = {}  # embryo -> new split file path

for embryo, df_embryo in df_new.groupby("embryo"):
    gene_label = df_embryo["gene"].iloc[0]
    out_path = os.path.join(split_output_dir, f"{embryo}_TRIPLETS.csv")
    df_embryo.to_csv(out_path, index=False)
    embryo_to_gene[embryo] = gene_label
    embryo_to_new_path[embryo] = out_path

print(f"Split into {len(embryo_to_gene)} per-embryo file(s) in:\n  {split_output_dir}\n")


# ==========================================================
# STEP 2: BUILD LOOKUP FROM OLD MAPPING (embryo name -> old row info)
# ==========================================================

df_old = pd.read_csv(old_mapping_csv)

required_cols = ["gene_id", "triplet_path", "roi_path"]
for col in required_cols:
    if col not in df_old.columns:
        raise ValueError(f"{old_mapping_csv} missing column: {col}")


def derive_embryo_name(triplet_path):
    if not isinstance(triplet_path, str):
        return None
    base = os.path.splitext(os.path.basename(triplet_path))[0]
    if base.endswith("_TRIPLETS"):
        base = base[: -len("_TRIPLETS")]
    return base


df_old["_derived_embryo"] = df_old["triplet_path"].apply(derive_embryo_name)


# ==========================================================
# STEP 3: REBUILD -- remove old rows for affected embryos, add corrected ones
# ==========================================================

affected_embryos = set(embryo_to_gene.keys())

# rows to KEEP unchanged: everything whose derived embryo is NOT one we just corrected
rows_kept = df_old[~df_old["_derived_embryo"].isin(affected_embryos)].drop(columns=["_derived_embryo"])

# for each affected embryo, find its OLD roi_path (must preserve this -- ROI
# polygons haven't changed, only the triplet source data has)
new_rows = []
embryos_with_no_old_mapping = []

for embryo, gene_label in embryo_to_gene.items():
    old_matches = df_old[df_old["_derived_embryo"] == embryo]

    if len(old_matches) == 0:
        embryos_with_no_old_mapping.append(embryo)
        continue

    if len(old_matches) > 1:
        print(f"  WARNING: '{embryo}' had {len(old_matches)} old rows (expected 1) -- using the first, "
              f"please verify manually")

    old_roi_path = old_matches.iloc[0]["roi_path"]

    new_rows.append({
        "gene_id": gene_label,
        "triplet_path": embryo_to_new_path[embryo],
        "roi_path": old_roi_path,
    })

df_new_rows = pd.DataFrame(new_rows)
df_updated = pd.concat([rows_kept, df_new_rows], ignore_index=True)

df_updated.to_csv(output_mapping_csv, index=False)


# ==========================================================
# SUMMARY
# ==========================================================

print("=== SUMMARY ===")
print(f"Old mapping table rows: {len(df_old)}")
print(f"Rows kept unchanged (unaffected embryos/genes): {len(rows_kept)}")
print(f"Rows removed (old, corrected embryos): {len(df_old) - len(rows_kept)}")
print(f"Rows added (new, corrected embryos): {len(df_new_rows)}")
print(f"New mapping table total rows: {len(df_updated)}")

if embryos_with_no_old_mapping:
    print(f"\n*** {len(embryos_with_no_old_mapping)} embryo(s) had NO existing row in the old mapping "
          f"(never had an ROI assigned before) -- these need MANUAL roi_path entry: ***")
    for e in embryos_with_no_old_mapping:
        print(f"  {e}")

print(f"\nSaved: {output_mapping_csv}")
print(f"\n*** This is a NEW file -- review it, then manually rename/replace ROI_TRIPLET_BATCH.csv "
      f"once you're confident it's correct. Nothing has been overwritten yet. ***")
