import os
import gc
import time
import re
import czifile
import xml.etree.ElementTree as ET
import tifffile as tiff
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

# ==========================================================
# USER INPUTS
# ==========================================================

output_dir = r"F:\RSFISH\script for nuclei assignment_july2026\OUTPUT_FULL_SCREEN_FIXED_VOXEL"

masks_root = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\Segmentation_masks"
spots_root = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\AP_files\spotdetetcetion"
ap_files_root = r"Z:\Sujay\RNAi screen\DNA-FISH_RNAi_screen\AP_file_2026\AP_files"

# folders under ap_files_root that are NOT gene folders -- skip when indexing CZIs
SKIP_FOLDERS = {"OUTPUT", "spotdetetcetion", "sxlgfpmale dna fish_15082025", "awhgfp",
                 "Set2", "set2_ROI", "spotdetetction2"}

# TESTING FILTER: restrict to specific gene mask-folder(s) only. Set to None to
# process every gene folder found under masks_root (full screen run).
GENE_FOLDERS_FILTER = ["Beaf32_MASKS", "CTCF_MASKS", "Hth_MASKS", "LACZ",
                        "SUHW_MASKS", "Ubx_MASKS", "cp190_MASKS"]
# ^ the 7 CONFIRMED genes with voxel mismatches (AWH_MASKS removed -- the
# earlier "awhgfp" mismatch hit belongs to a separate, unrelated experiment,
# not AWH_MASKS. AWH_MASKS's own raw CZI location is still unknown and has
# NOT yet been checked for voxel mismatches -- do not assume it's clean.)

CHANNEL_SUFFIXES = {
    "prom": "_CH3_488nm_RSFISH.csv",
    "E":    "_CH2_555nm_RSFISH.csv",
    "DG":   "_CH1_633nm_RSFISH.csv",
}

# Fallback voxel size, used ONLY if an embryo's own CZI can't be found/read.
# Flagged clearly in the output log -- should be rare if the CZI index works.
FALLBACK_VOXEL = {"X": 0.0495, "Y": 0.0495, "Z": 0.14}

MAX_DIST = 1.2
TOP_N = 8000
MAX_TRIPLETS_PER_NUCLEUS = 2


# ==========================================================
# DERIVED PATHS
# ==========================================================

os.makedirs(output_dir, exist_ok=True)
skip_log_csv = os.path.join(output_dir, "SKIPPED_EMBRYOS_LOG.csv")
voxel_log_csv = os.path.join(output_dir, "VOXEL_SIZE_USED_PER_EMBRYO_FULLSCREEN.csv")
output_csv = os.path.join(output_dir, "ALL_GENES_FINAL_TRIPLETS_FIXED_VOXEL.csv")


# ==========================================================
# STEP 0: BUILD A CZI INDEX ACROSS ALL OF ap_files_root, ONE TIME
# (base filename -> full path), so we never have to guess folder names
# ==========================================================

def build_czi_index(ap_files_root, skip_folders):
    index = {}
    collisions = []

    for folder in os.listdir(ap_files_root):
        if folder in skip_folders:
            continue
        folder_path = os.path.join(ap_files_root, folder)
        if not os.path.isdir(folder_path):
            continue
        for fname in os.listdir(folder_path):
            if fname.lower().endswith(".czi"):
                base_name = fname[:-4]  # strip .czi
                full_path = os.path.join(folder_path, fname)
                if base_name in index:
                    collisions.append((base_name, index[base_name], full_path))
                else:
                    index[base_name] = full_path

    return index, collisions


print("Building CZI file index across all of AP_files (one-time scan)...")
czi_index, czi_collisions = build_czi_index(ap_files_root, SKIP_FOLDERS)
print(f"Indexed {len(czi_index)} CZI file(s)")
if czi_collisions:
    print(f"WARNING: {len(czi_collisions)} base-filename collision(s) found (same name in multiple folders):")
    for base_name, path1, path2 in czi_collisions:
        print(f"  '{base_name}': kept '{path1}', ignored '{path2}'")
print()


# ==========================================================
# PER-EMBRYO VOXEL SIZE (from the CZI index built above)
# ==========================================================

def get_voxel_size_um(czi_path):
    with czifile.CziFile(czi_path) as czi:
        xml_str = czi.metadata()
    root = ET.fromstring(xml_str)
    voxel = {}
    for dist in root.iter("Distance"):
        axis = dist.get("Id")
        value_elem = dist.find("Value")
        if value_elem is not None and axis in ("X", "Y", "Z"):
            voxel[axis] = float(value_elem.text) * 1e6
    return voxel


def get_voxel_for_embryo(embryo_base_name, voxel_log_rows):
    czi_path = czi_index.get(embryo_base_name)

    if czi_path is None:
        print(f"    WARNING: no CZI found in index for '{embryo_base_name}' -- using fallback voxel size")
        voxel_log_rows.append({
            "embryo": embryo_base_name, "source": "FALLBACK (CZI not found in index)",
            "voxel_x": FALLBACK_VOXEL["X"], "voxel_y": FALLBACK_VOXEL["Y"], "voxel_z": FALLBACK_VOXEL["Z"]
        })
        return FALLBACK_VOXEL["X"], FALLBACK_VOXEL["Y"], FALLBACK_VOXEL["Z"]

    try:
        voxel = get_voxel_size_um(czi_path)
        vx, vy, vz = voxel["X"], voxel["Y"], voxel["Z"]
        voxel_log_rows.append({
            "embryo": embryo_base_name, "source": "CZI metadata",
            "voxel_x": vx, "voxel_y": vy, "voxel_z": vz
        })
        return vx, vy, vz
    except Exception as e:
        print(f"    WARNING: failed to read voxel size from {czi_path} ({e}) -- using fallback")
        voxel_log_rows.append({
            "embryo": embryo_base_name, "source": f"FALLBACK (read error: {e})",
            "voxel_x": FALLBACK_VOXEL["X"], "voxel_y": FALLBACK_VOXEL["Y"], "voxel_z": FALLBACK_VOXEL["Z"]
        })
        return FALLBACK_VOXEL["X"], FALLBACK_VOXEL["Y"], FALLBACK_VOXEL["Z"]


# ==========================================================
# EXACT 3-WAY PER-NUCLEUS TRIPLET MATCHING (same validated logic used all day)
# ==========================================================

def match_triplets_exact(A_um, B_um, C_um, max_dist, candidate_limit=200, max_backtrack_calls=50000):
    nA, nB, nC = len(A_um), len(B_um), len(C_um)
    if nA == 0 or nB == 0 or nC == 0:
        return [], True

    D_AB = cdist(A_um, B_um)
    D_AC = cdist(A_um, C_um)
    D_BC = cdist(B_um, C_um)

    candidates = []
    for a in range(nA):
        b_ok = np.where(D_AB[a] <= max_dist)[0]
        if len(b_ok) == 0:
            continue
        c_ok = np.where(D_AC[a] <= max_dist)[0]
        if len(c_ok) == 0:
            continue
        for b in b_ok:
            if D_BC[b][c_ok].min() > max_dist:
                continue
            for c in c_ok:
                if D_BC[b, c] <= max_dist:
                    total = D_AB[a, b] + D_AC[a, c] + D_BC[b, c]
                    candidates.append((a, b, c, total))

    if len(candidates) == 0:
        return [], True

    def hungarian_fallback():
        row_ind, col_ind = linear_sum_assignment(D_AB)
        map_AB = {a: b for a, b in zip(row_ind, col_ind) if D_AB[a, b] <= max_dist}
        row_ind2, col_ind2 = linear_sum_assignment(D_AC)
        map_AC = {a: c for a, c in zip(row_ind2, col_ind2) if D_AC[a, c] <= max_dist}
        results = []
        for a in map_AB:
            if a not in map_AC:
                continue
            b, c = map_AB[a], map_AC[a]
            if D_BC[b, c] <= max_dist:
                results.append((a, b, c, D_AB[a, b] + D_AC[a, c] + D_BC[b, c]))
        return results

    if len(candidates) > candidate_limit:
        return hungarian_fallback(), False

    candidates.sort(key=lambda t: t[3])
    best = {"count": -1, "total": np.inf, "chosen": []}
    call_counter = {"n": 0}

    class BacktrackAborted(Exception):
        pass

    def backtrack(idx, used_a, used_b, used_c, chosen, total):
        call_counter["n"] += 1
        if call_counter["n"] > max_backtrack_calls:
            raise BacktrackAborted()
        if len(chosen) + (len(candidates) - idx) < best["count"]:
            return
        if idx == len(candidates):
            if (len(chosen) > best["count"]) or (len(chosen) == best["count"] and total < best["total"]):
                best["count"] = len(chosen)
                best["total"] = total
                best["chosen"] = list(chosen)
            return
        a, b, c, d = candidates[idx]
        backtrack(idx + 1, used_a, used_b, used_c, chosen, total)
        if a not in used_a and b not in used_b and c not in used_c:
            chosen.append((a, b, c, d))
            backtrack(idx + 1, used_a | {a}, used_b | {b}, used_c | {c}, chosen, total + d)
            chosen.pop()

    try:
        backtrack(0, set(), set(), set(), [], 0.0)
        return best["chosen"], True
    except BacktrackAborted:
        return hungarian_fallback(), False


def pixels_to_um(px, vx, vy, vz):
    out = px.copy().astype(float)
    out[:, 0] *= vx
    out[:, 1] *= vy
    out[:, 2] *= vz
    return out


def load_filter_assign(spots_file, role, labels, Xdim, Ydim, Zdim):
    df = pd.read_csv(spots_file)
    required = {"x", "y", "z", "intensity"}
    if not required.issubset(df.columns):
        raise ValueError(f"{spots_file} missing columns {required}")

    df = df.sort_values("intensity", ascending=False).reset_index(drop=True)
    if len(df) > TOP_N:
        df = df.iloc[:TOP_N].copy()

    x_i = np.round(df["x"]).astype(int).clip(0, Xdim - 1)
    y_i = np.round(df["y"]).astype(int).clip(0, Ydim - 1)
    z_i = np.round(df["z"]).astype(int).clip(0, Zdim - 1)

    df["nucleus_ID"] = labels[z_i.to_numpy(), y_i.to_numpy(), x_i.to_numpy()]
    df = df[df["nucleus_ID"] != 0].copy()
    return df


def normalize_gene_label(folder_name):
    """Strip a trailing _MASK/_MASKS (any case) to get a clean gene label."""
    return re.sub(r'_masks?$', '', folder_name, flags=re.IGNORECASE)


def process_embryo(mask_path, spots_paths, embryo_name, gene_label, voxel_log_rows):
    vx, vy, vz = get_voxel_for_embryo(embryo_name, voxel_log_rows)

    labels = None
    last_err = None
    for attempt in range(2):
        try:
            labels = tiff.memmap(mask_path)
            break
        except OSError as e:
            last_err = e
            time.sleep(2)
    if labels is None:
        raise RuntimeError(f"Failed to read mask after retry: {last_err}")

    Zdim, Ydim, Xdim = labels.shape

    dfs = {role: load_filter_assign(path, role, labels, Xdim, Ydim, Zdim)
           for role, path in spots_paths.items()}

    del labels
    gc.collect()

    shared_nuclei = set(dfs["prom"]["nucleus_ID"]) & set(dfs["E"]["nucleus_ID"]) & set(dfs["DG"]["nucleus_ID"])

    rows = []
    n_fallback = 0

    for nid in sorted(shared_nuclei):
        df_prom = dfs["prom"][dfs["prom"]["nucleus_ID"] == nid].reset_index(drop=True)
        df_E    = dfs["E"][dfs["E"]["nucleus_ID"] == nid].reset_index(drop=True)
        df_DG   = dfs["DG"][dfs["DG"]["nucleus_ID"] == nid].reset_index(drop=True)

        A_px = df_prom[["x", "y", "z"]].to_numpy()
        B_px = df_E[["x", "y", "z"]].to_numpy()
        C_px = df_DG[["x", "y", "z"]].to_numpy()

        A_um = pixels_to_um(A_px, vx, vy, vz)
        B_um = pixels_to_um(B_px, vx, vy, vz)
        C_um = pixels_to_um(C_px, vx, vy, vz)

        results, exact = match_triplets_exact(A_um, B_um, C_um, MAX_DIST)
        if not exact:
            n_fallback += 1

        for a, b, c, total_dist in results:
            rows.append({
                "gene": gene_label,
                "embryo": embryo_name,
                "nucleus_id": nid,
                "exact_match": exact,
                "voxel_x": vx, "voxel_y": vy, "voxel_z": vz,
                "x_svb_px": A_px[a, 0], "y_svb_px": A_px[a, 1], "z_svb_px": A_px[a, 2],
                "x_E_px":   B_px[b, 0], "y_E_px":   B_px[b, 1], "z_E_px":   B_px[b, 2],
                "x_DG_px":  C_px[c, 0], "y_DG_px":  C_px[c, 1], "z_DG_px":  C_px[c, 2],
                "x_svb_um": A_um[a, 0], "y_svb_um": A_um[a, 1], "z_svb_um": A_um[a, 2],
                "x_E_um":   B_um[b, 0], "y_E_um":   B_um[b, 1], "z_E_um":   B_um[b, 2],
                "x_DG_um":  C_um[c, 0], "y_DG_um":  C_um[c, 1], "z_DG_um":  C_um[c, 2],
                "dist_svb_E":  np.linalg.norm(A_um[a] - B_um[b]),
                "dist_svb_DG": np.linalg.norm(A_um[a] - C_um[c]),
                "dist_E_DG":   np.linalg.norm(B_um[b] - C_um[c]),
                "total_dist":  total_dist,
            })

    df_all = pd.DataFrame(rows)
    n_before_cap = len(df_all)
    n_dropped_nuclei = 0

    if len(df_all) > 0:
        counts = df_all.groupby("nucleus_id").size()
        good_nuclei = counts[counts <= MAX_TRIPLETS_PER_NUCLEUS].index
        n_dropped_nuclei = (counts > MAX_TRIPLETS_PER_NUCLEUS).sum()
        df_final = df_all[df_all["nucleus_id"].isin(good_nuclei)].copy()
    else:
        df_final = df_all.copy()

    stats = {
        "gene": gene_label, "embryo": embryo_name,
        "voxel_x": vx, "voxel_y": vy, "voxel_z": vz,
        "n_nuclei_shared": len(shared_nuclei),
        "n_triplets_before_cap": n_before_cap,
        "n_nuclei_dropped_cap": n_dropped_nuclei,
        "n_fallback_nuclei": n_fallback,
        "n_final_triplets": len(df_final),
    }

    return df_final, stats


# ==========================================================
# MAIN LOOP
# ==========================================================

all_final_dfs = []
all_stats = []
skipped = []
voxel_log_rows = []

gene_folders = [d for d in os.listdir(masks_root) if os.path.isdir(os.path.join(masks_root, d))]
if GENE_FOLDERS_FILTER is not None:
    gene_folders = [d for d in gene_folders if d in GENE_FOLDERS_FILTER]
    print(f"TEST MODE: restricted to gene folder(s): {GENE_FOLDERS_FILTER}\n")

print(f"Processing {len(gene_folders)} gene folder(s): {gene_folders}\n")

for gene_folder in gene_folders:
    gene_dir = os.path.join(masks_root, gene_folder)
    gene_label = normalize_gene_label(gene_folder)

    mask_files = [f for f in os.listdir(gene_dir) if f.lower().endswith("_mask.tif")]
    print(f"=== Gene: {gene_label} ({gene_folder}) -- {len(mask_files)} mask file(s) ===")

    for mask_fname in mask_files:
        if not mask_fname.endswith("_mask.tif"):
            skipped.append({"gene": gene_label, "file": mask_fname, "reason": "unexpected filename casing"})
            continue

        base_name = mask_fname[: -len("_mask.tif")]
        mask_path = os.path.join(gene_dir, mask_fname)

        spots_paths = {role: os.path.join(spots_root, base_name + suffix)
                        for role, suffix in CHANNEL_SUFFIXES.items()}

        missing = [p for p in spots_paths.values() if not os.path.isfile(p)]
        if missing:
            skipped.append({"gene": gene_label, "file": base_name, "reason": f"missing CSV(s): {len(missing)}"})
            print(f"  SKIP {base_name} -- missing {len(missing)} CSV(s)")
            continue

        print(f"  Processing: {base_name} ...")
        try:
            df_final, stats = process_embryo(mask_path, spots_paths, base_name, gene_label, voxel_log_rows)
            all_final_dfs.append(df_final)
            all_stats.append(stats)
            print(f"    -> {stats['n_final_triplets']} final triplets  (voxel_x={stats['voxel_x']:.5f})")
        except Exception as e:
            skipped.append({"gene": gene_label, "file": base_name, "reason": f"ERROR: {e}"})
            print(f"    ERROR: {e}")
        finally:
            gc.collect()

    print()


# ==========================================================
# SAVE OUTPUT
# ==========================================================

if all_final_dfs:
    df_master = pd.concat(all_final_dfs, ignore_index=True)
else:
    df_master = pd.DataFrame()

df_master.to_csv(output_csv, index=False)
print(f"\nSaved: {output_csv}")
print(f"Total final triplets: {len(df_master)}")

if skipped:
    pd.DataFrame(skipped).to_csv(skip_log_csv, index=False)
    print(f"Saved skip log ({len(skipped)} skipped): {skip_log_csv}")

pd.DataFrame(voxel_log_rows).to_csv(voxel_log_csv, index=False)
print(f"Saved voxel log: {voxel_log_csv}")

if all_stats:
    df_stats = pd.DataFrame(all_stats)
    print(f"\n=== SUMMARY ===")
    print(f"Embryos processed: {len(df_stats)}")
    print(f"Embryos skipped: {len(skipped)}")
    print(f"\nPer-gene totals:")
    print(df_stats.groupby("gene").agg(
        n_embryos=("embryo", "count"),
        total_triplets=("n_final_triplets", "sum"),
    ))
