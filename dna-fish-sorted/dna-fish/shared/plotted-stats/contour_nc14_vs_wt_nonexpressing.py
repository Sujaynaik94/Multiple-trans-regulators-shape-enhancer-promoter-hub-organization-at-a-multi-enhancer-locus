"""
PLOTTED contour generation — nc14 whole embryo vs WT Non-expressing
Based directly on Le et al. 2026 (StathopoulosLab/find_spots notebook)

Coordinate system:
  P  (promoter) = (0, 0)       — fixed at origin
  DG2           = (d_P_DG, 0)  — on x-axis  [= d3 in their notation]
  E6            = (x_E6, y_E6) — positive y-space

Mapping to their variable names:
  d1 = dist_P_E   (E6-P distance,  their E2-PPE)
  d2 = dist_P_DG  (DG2-P distance, their E1-PPE) -- NOTE: in their code d2 is plotted on x-axis
  d3 = dist_DG_E  (DG2-E6,         their E1-E2)  -- d3 goes on x-axis = DG2 position

Wait - let us be precise. In their system:
  E1 = fixed at (0,0)
  E2 = on x-axis at (E1-E2, 0)   -> d3 = E1-E2
  PPE = in positive y            -> computed from d1=E2-PPE, d2=E1-PPE

In our system:
  P  = fixed at (0,0)
  DG2 = on x-axis at (P-DG2, 0) -> d3 = dist_P_DG
  E6  = in positive y            -> computed from d1=dist_DG_E (E6-DG2), d2=dist_P_E (E6-P)

So:
  d1 = dist_DG_E  (distance from DG2 to E6, their E2-PPE analog)
  d2 = dist_P_E   (distance from P to E6,   their E1-PPE analog)
  d3 = dist_P_DG  (distance from P to DG2,  their E1-E2 analog, goes on x-axis)

Triangle equations (from their eq 8-11):
  x_E6 = (d3^2 + d2^2 - d1^2) / (2*d3)
  y_E6 = sqrt(d2^2 - x_E6^2)

Where:
  d3 = dist_P_DG  (DG2 x-position)
  d1 = dist_DG_E  
  d2 = dist_P_E
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, ks_2samp
from scipy.interpolate import NearestNDInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import sys
import os


# ============================================================
# FONT / VECTOR OUTPUT SETTINGS
# Match the RNAi script so text remains editable in Illustrator
# ============================================================
FONT_SIZE = 8

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
    "figure.titlesize": FONT_SIZE,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# ============================================================
# INPUTS — edit these paths
# ============================================================
XLSX_PATH  = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/ROI/ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
SHEET      = "Filtered_Triplets"
WT_GENE_ID = "SXLGFP"
NC14_CSV   = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/triplets_fromnucleisegment/nc14_MASKS_FINAL_TRIPLETS.csv"
OUTPUT_DIR = "F:/Thesis/Figures/PLOTTED_NC14_VS_WT_NONEXP_NEW_DATASET"
MDT = 2.0   # maximum distance threshold µm

# ============================================================
# GRID (matches their xi/yi in plot_contour_stage)
# ============================================================
XI = np.linspace(-0.1, 0.8, 100)
YI = np.linspace(0.0,  1.0,  75)

# ============================================================
# COLORS
# ============================================================
WHITE_PURPLE_ORANGE = LinearSegmentedColormap.from_list(
    "white_purple_hot_orange",
    [
        (0.00, "white"),
        (0.15, "#D8BFD8"),
        (0.25, "purple"),
        (0.55, "#FFD580"),
        (0.80, "#FF8C00"),
        (1.00, "#FF4500"),
    ]
)

DIVERGING_RB = LinearSegmentedColormap.from_list(
    "red_white_blue",
    [
        (0.0,  "#053061"),
        (0.25, "#4393c3"),
        (0.5,  "white"),
        (0.75, "#d6604d"),
        (1.0,  "#67001f"),
    ]
)

CONDITION_COLORS = {
    "Non-expressing": "#686765",
    "Dorsal":         "#cc2829",
    "Ventral":        "#2b3f99",
}

# ============================================================
# LOAD DATA
# ============================================================
def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def compute_distances(df):
    x_DG = pick_col(df, ["DG_CoM_X", "x_DG_um", "DG_x_um", "x_DG", "DG_x"])
    y_DG = pick_col(df, ["DG_CoM_Y", "y_DG_um", "DG_y_um", "y_DG", "DG_y"])
    z_DG = pick_col(df, ["DG_CoM_Z", "z_DG_um", "DG_z_um", "z_DG", "DG_z"])
    x_E  = pick_col(df, ["E_CoM_X",  "x_E_um",  "E_x_um",  "x_E",  "E_x"])
    y_E  = pick_col(df, ["E_CoM_Y",  "y_E_um",  "E_y_um",  "y_E",  "E_y"])
    z_E  = pick_col(df, ["E_CoM_Z",  "z_E_um",  "E_z_um",  "z_E",  "E_z"])
    x_P  = pick_col(df, ["P_CoM_X",  "x_svb_um", "x_P_um", "x_promoter_um", "x_svb", "x_promoter"])
    y_P  = pick_col(df, ["P_CoM_Y",  "y_svb_um", "y_P_um", "y_promoter_um", "y_svb", "y_promoter"])
    z_P  = pick_col(df, ["P_CoM_Z",  "z_svb_um", "z_P_um", "z_promoter_um", "z_svb", "z_promoter"])

    def d3r(r, xA, yA, zA, xB, yB, zB):
        return np.sqrt((r[xA]-r[xB])**2 + (r[yA]-r[yB])**2 + (r[zA]-r[zB])**2)

    df["dist_P_DG"] = df.apply(lambda r: d3r(r, x_P, y_P, z_P, x_DG, y_DG, z_DG), axis=1)
    df["dist_P_E"]  = df.apply(lambda r: d3r(r, x_P, y_P, z_P, x_E,  y_E,  z_E),  axis=1)
    df["dist_DG_E"] = df.apply(lambda r: d3r(r, x_DG, y_DG, z_DG, x_E, y_E, z_E), axis=1)
    df = df.dropna(subset=["dist_P_DG", "dist_P_E", "dist_DG_E"])
    df = df[
        (df["dist_P_DG"] <= MDT) & (df["dist_P_DG"] > 0) &
        (df["dist_P_E"]  <= MDT) & (df["dist_P_E"]  > 0) &
        (df["dist_DG_E"] <= MDT) & (df["dist_DG_E"] > 0)
    ]
    return df[["dist_P_DG", "dist_P_E", "dist_DG_E"]].copy()

def load_distances(xlsx_path, sheet, nc14_csv):
    """
    Load two sources:
      1. WT Non-expressing from the new combined Excel dataset
      2. nc14 whole embryo from nc14_MASKS_FINAL_TRIPLETS.csv

    Both datasets already contain precomputed distances:
      dist_svb_E
      dist_svb_DG
      dist_E_DG

    These are renamed internally to:
      dist_P_E
      dist_P_DG
      dist_DG_E
    """

    results = {}

    # ========================================================
    # WT NON-EXPRESSING
    # ========================================================

    print("Loading WT Non-expressing from combined Excel dataset")
    print(f"  Excel: {xlsx_path}")
    print(f"  Sheet: {sheet}")

    df_wt = pd.read_excel(
        xlsx_path,
        sheet_name=sheet
    )

    required_wt = {
        "gene_id",
        "roi",
        "dist_svb_E",
        "dist_svb_DG",
        "dist_E_DG",
    }

    missing_wt = required_wt.difference(df_wt.columns)

    if missing_wt:
        raise KeyError(
            "Missing required WT columns: "
            + ", ".join(sorted(missing_wt))
        )

    df_wt = df_wt[
        df_wt["gene_id"]
        .astype(str)
        .str.strip()
        .str.upper()
        .eq(WT_GENE_ID)
    ].copy()

    if df_wt.empty:
        raise ValueError(
            f"No rows found for gene_id == {WT_GENE_ID!r}"
        )

    df_wt["roi_key"] = (
        df_wt["roi"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    print("\nWT ROI values before filtering:")
    print(
        sorted(
            df_wt["roi"]
            .dropna()
            .astype(str)
            .unique()
        )
    )

    df_wt = df_wt[
        df_wt["roi_key"].isin(
            [
                "nonexp",
                "non-exp",
                "non_exp",
                "nonexpressing",
                "non-expressing",
            ]
        )
    ].copy()

    if df_wt.empty:
        raise ValueError(
            "No WT Non-expressing rows were found."
        )

    df_wt["dist_P_E"] = pd.to_numeric(
        df_wt["dist_svb_E"],
        errors="coerce"
    )
    df_wt["dist_P_DG"] = pd.to_numeric(
        df_wt["dist_svb_DG"],
        errors="coerce"
    )
    df_wt["dist_DG_E"] = pd.to_numeric(
        df_wt["dist_E_DG"],
        errors="coerce"
    )

    before_wt = len(df_wt)

    df_wt = df_wt.dropna(
        subset=[
            "dist_P_E",
            "dist_P_DG",
            "dist_DG_E",
        ]
    ).copy()

    df_wt = df_wt[
        (df_wt["dist_P_DG"] > 0) &
        (df_wt["dist_P_DG"] <= MDT) &
        (df_wt["dist_P_E"] > 0) &
        (df_wt["dist_P_E"] <= MDT) &
        (df_wt["dist_DG_E"] > 0) &
        (df_wt["dist_DG_E"] <= MDT)
    ].copy()

    wt_distances = df_wt[
        [
            "dist_P_DG",
            "dist_P_E",
            "dist_DG_E",
        ]
    ].copy()

    print(
        f"  WT Non-expressing retained: "
        f"{len(wt_distances)} of {before_wt}"
    )

    results["Non-expressing"] = wt_distances

    # ========================================================
    # nc14 WHOLE EMBRYO
    # ========================================================

    print("\nLoading nc14 whole embryo")
    print(f"  CSV: {nc14_csv}")

    df_nc14 = pd.read_csv(
        nc14_csv
    )

    required_nc14 = {
        "embryo",
        "dist_svb_E",
        "dist_svb_DG",
        "dist_E_DG",
    }

    missing_nc14 = required_nc14.difference(
        df_nc14.columns
    )

    if missing_nc14:
        raise KeyError(
            "Missing required nc14 columns: "
            + ", ".join(sorted(missing_nc14))
        )

    print("\nnc14 embryos and row counts:")
    print(
        df_nc14["embryo"]
        .value_counts()
        .sort_index()
    )

    df_nc14["dist_P_E"] = pd.to_numeric(
        df_nc14["dist_svb_E"],
        errors="coerce"
    )
    df_nc14["dist_P_DG"] = pd.to_numeric(
        df_nc14["dist_svb_DG"],
        errors="coerce"
    )
    df_nc14["dist_DG_E"] = pd.to_numeric(
        df_nc14["dist_E_DG"],
        errors="coerce"
    )

    before_nc14 = len(df_nc14)

    df_nc14 = df_nc14.dropna(
        subset=[
            "dist_P_E",
            "dist_P_DG",
            "dist_DG_E",
        ]
    ).copy()

    df_nc14 = df_nc14[
        (df_nc14["dist_P_DG"] > 0) &
        (df_nc14["dist_P_DG"] <= MDT) &
        (df_nc14["dist_P_E"] > 0) &
        (df_nc14["dist_P_E"] <= MDT) &
        (df_nc14["dist_DG_E"] > 0) &
        (df_nc14["dist_DG_E"] <= MDT)
    ].copy()

    nc14_distances = df_nc14[
        [
            "dist_P_DG",
            "dist_P_E",
            "dist_DG_E",
        ]
    ].copy()

    print(
        f"  nc14 retained: "
        f"{len(nc14_distances)} of {before_nc14}"
    )

    results["nc14"] = nc14_distances

    return results


# ============================================================
# KDE (Scott's rule, normalized)
# ============================================================
def compute_kde(data, x_range=None, n_points=500):
    data = np.asarray(data, float)
    data = data[np.isfinite(data)]
    if x_range is None:
        x_range = np.linspace(0, MDT, n_points)
    kde = gaussian_kde(data, bw_method='scott')
    y = kde(x_range)
    # normalize so area = 1
    y = y / np.trapezoid(y, x_range)
    return x_range, y


# ============================================================
# TRIANGLE VERTICES (vectorized, exact from their cell 134)
# ============================================================
def compute_triangle_vertices_from_PDF(d1_arr, d2_arr, p_d1_arr, p_d2_arr, d3):
    """
    d1 = dist_DG_E  (E6-DG2, their E2-PPE)
    d2 = dist_P_E   (E6-P,   their E1-PPE)
    d3 = dist_P_DG  (DG2-P,  their E1-E2) -- fixed for this iteration

    E6 coordinates:
      x_E6 = (d3^2 + d2^2 - d1^2) / (2*d3)
      y_E6 = sqrt(d2^2 - x_E6^2)
    """
    D1, D2 = np.meshgrid(d1_arr, d2_arr, indexing='ij')
    P1, P2 = np.meshgrid(p_d1_arr, p_d2_arr, indexing='ij')
    p = P1 * P2  # joint PDF = product (independence assumption)

    # Triangle inequality
    valid = (D1 + D2 >= d3) & (D1 + d3 >= D2) & (D2 + d3 >= D1)

    x_E6 = np.where(valid, (d3**2 + D2**2 - D1**2) / (2 * d3), np.nan)
    y2   = np.where(valid, D2**2 - x_E6**2, np.nan)
    y_E6 = np.where(y2 >= 0, np.sqrt(np.maximum(y2, 0)), np.nan)

    return x_E6, y_E6, np.where(valid, p, np.nan)


# ============================================================
# WEIGHTED Z SURFACE (exact from their cell 135)
# ============================================================
def compute_weighted_Z(d1_arr, d2_arr, p_d1_arr, p_d2_arr,
                       d3_arr, p_d3_arr, xi, yi,
                       num_samples=50, subsample_stride=3):
    """
    Exact implementation of their compute_weighted_Z.
    d1 = dist_DG_E, d2 = dist_P_E, d3 = dist_P_DG (x-axis)
    """
    X_grid, Y_grid = np.meshgrid(xi, yi)

    # Evenly spaced d3 sampling (avoids bias)
    indices_even = np.linspace(0, len(d3_arr) - 1, num_samples, dtype=int)
    d3_sampled   = d3_arr[indices_even]
    p_d3_sampled = p_d3_arr[indices_even]
    p_d3_sampled = p_d3_sampled / np.sum(p_d3_sampled)

    # Riemann sum weights
    bin_edges = np.concatenate([
        [d3_sampled[0] - (d3_sampled[1] - d3_sampled[0]) / 2],
        (d3_sampled[:-1] + d3_sampled[1:]) / 2,
        [d3_sampled[-1] + (d3_sampled[-1] - d3_sampled[-2]) / 2]
    ])
    bin_widths = np.diff(bin_edges)
    weights    = p_d3_sampled * bin_widths
    weights    = weights / np.sum(weights)

    # Subsample d1 and d2 arrays for performance
    d1_sub    = d1_arr[::subsample_stride]
    p_d1_sub  = p_d1_arr[::subsample_stride]
    d2_sub    = d2_arr[::subsample_stride]
    p_d2_sub  = p_d2_arr[::subsample_stride]

    Z_weighted   = np.zeros_like(X_grid, dtype=float)
    total_weight = 0.0
    dx = xi[1] - xi[0]
    dy = yi[1] - yi[0]

    for d3, w in zip(d3_sampled, weights):
        x_coords, y_coords, p_vals = compute_triangle_vertices_from_PDF(
            d1_sub, d2_sub, p_d1_sub, p_d2_sub, d3
        )
        valid = ~np.isnan(x_coords) & ~np.isnan(y_coords) & ~np.isnan(p_vals)
        if np.count_nonzero(valid) == 0:
            continue

        interpolator = NearestNDInterpolator(
            list(zip(x_coords[valid], y_coords[valid])),
            p_vals[valid]
        )
        Z_interp = interpolator(X_grid, Y_grid)
        area = np.sum(Z_interp) * dx * dy
        if area > 0:
            Z_interp   /= area
            Z_weighted += Z_interp * w
        total_weight += w

    if total_weight > 0:
        Z_weighted /= total_weight

    print(f"  Z surface area: {np.sum(Z_weighted) * dx * dy:.4f}")
    return Z_weighted, X_grid, Y_grid


# ============================================================
# PLOT SINGLE CONTOUR
# ============================================================
def plot_contour(Z, X_grid, Y_grid,
                 d3_arr, p_d3_arr,
                 condition_name,
                 vmin=0, vmax=3.5,
                 ax=None, fig=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    # Peak of Z = MPP of E6
    peak_idx = np.unravel_index(np.argmax(Z), Z.shape)
    peak_x   = X_grid[peak_idx]
    peak_y   = Y_grid[peak_idx]

    # Modal d3 = most probable DG2 x-position
    i_mode = np.argmax(p_d3_arr)
    d3_mode = d3_arr[i_mode]

    levels = np.linspace(vmin, vmax, 15)
    cf = ax.contourf(X_grid, Y_grid, Z,
                     levels=levels,
                     cmap=WHITE_PURPLE_ORANGE,
                     vmin=vmin, vmax=vmax)

    # Triangle lines
    E1 = np.array([0.0, 0.0])
    E2 = np.array([d3_mode, 0.0])
    E6_mpp = np.array([peak_x, peak_y])

    ax.plot([E1[0], E2[0]],    [E1[1], E2[1]],    'k-', lw=1.5)
    ax.plot([E1[0], E6_mpp[0]], [E1[1], E6_mpp[1]], 'k-', lw=1.5)
    ax.plot([E2[0], E6_mpp[0]], [E2[1], E6_mpp[1]], 'k-', lw=1.5)

    # Reference dots
    ax.plot(*E1,     'o', color='#8B3A8B', ms=10, zorder=5)   # P — magenta
    ax.plot(*E2,     'o', color='#F5C842', ms=10, zorder=5)   # DG2 — yellow
    ax.plot(*E6_mpp, 'o', color='#4FB3E8', ms=10, zorder=5)   # E6 MPP — blue

    # Labels
    ax.text(E1[0]-0.05, E1[1]-0.07, 'P',    fontsize=12, fontweight='bold', ha='right')
    ax.text(E2[0]+0.03, E2[1]-0.07, 'DG2',  fontsize=12, fontweight='bold')
    ax.text(E6_mpp[0]+0.03, E6_mpp[1]+0.04, 'E6', fontsize=12, fontweight='bold')

    # Grayscale bar for DG2 (d3) positional likelihood
    ax_pos = ax.get_position()
    bar_height = 0.04
    bar_bottom = ax_pos.y0 - 0.08
    ax_bar = fig.add_axes([ax_pos.x0, bar_bottom, ax_pos.width, bar_height])

    # Map d3 values to x positions on plot
    d3_vis = d3_arr[(d3_arr >= XI[0]) & (d3_arr <= XI[-1])]
    p_vis  = p_d3_arr[(d3_arr >= XI[0]) & (d3_arr <= XI[-1])]
    p_norm = p_vis / np.max(p_vis) if np.max(p_vis) > 0 else p_vis

    ax_bar.imshow(p_norm[np.newaxis, :],
                  aspect='auto', origin='lower',
                  cmap='Greys', vmin=0, vmax=1,
                  extent=[d3_vis[0], d3_vis[-1], 0, 1])
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("DG2 distance from P (µm)", fontsize=10)
    ax_bar.set_xlim(XI[0], XI[-1])

    ax.set_xlim(XI[0], XI[-1])
    ax.set_ylim(YI[0], YI[-1])
    ax.set_xlabel("x (µm)", fontsize=12)
    ax.set_ylabel("y (µm)", fontsize=12)
    ax.set_title(condition_name,
                 fontsize=14, fontweight='bold',
                 color=CONDITION_COLORS.get(condition_name, 'black'))

    return cf, fig, ax


# ============================================================
# SUBTRACTION CONTOUR
# ============================================================
def plot_subtraction(Z1, Z2, X_grid, Y_grid,
                     peak1, peak2,
                     d3_mode1, d3_mode2,
                     d3_arr1, p_d3_arr1,
                     d3_arr2, p_d3_arr2,
                     label1, label2,
                     simes_p, driven_by,
                     vmin=-1.8, vmax=1.8,
                     ax=None, fig=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    Z_diff = Z1 - Z2
    # Mask edges where data is sparse to avoid artifacts
    x_mask = (X_grid > 0.65) | (X_grid < -0.05)
    Z_diff[x_mask] = 0.0

    cf = ax.contourf(X_grid, Y_grid, Z_diff,
                     levels=np.linspace(vmin, vmax, 20),
                     cmap=DIVERGING_RB,
                     vmin=vmin, vmax=vmax)

    # MPP markers — filled circles for condition 1, hollow squares for condition 2
    E1 = np.array([0.0, 0.0])
    E2_1 = np.array([d3_mode1, 0.0])
    E2_2 = np.array([d3_mode2, 0.0])

    # Condition 1 — filled circles
    ax.plot(*E1,      'o', color='#8B3A8B', ms=10, zorder=5)
    ax.plot(*E2_1,    'o', color='#F5C842', ms=10, zorder=5)
    ax.plot(*peak1,   'o', color='#4FB3E8', ms=10, zorder=5)

    # Condition 2 — hollow squares
    ax.plot(*E1,      's', color='#8B3A8B', ms=10, zorder=5,
            mfc='none', mew=2)
    ax.plot(*E2_2,    's', color='#F5C842', ms=10, zorder=5,
            mfc='none', mew=2)
    ax.plot(*peak2,   's', color='#4FB3E8', ms=10, zorder=5,
            mfc='none', mew=2)

    # Mode triangle lines — solid for cond1, dashed for cond2
    for peak, E2, ls in [(peak1, E2_1, '-'), (peak2, E2_2, '--')]:
        ax.plot([E1[0], E2[0]],   [E1[1], E2[1]],   color='black', ls=ls, lw=1.2)
        ax.plot([E1[0], peak[0]], [E1[1], peak[1]],  color='black', ls=ls, lw=1.2)
        ax.plot([E2[0], peak[0]], [E2[1], peak[1]],  color='black', ls=ls, lw=1.2)

    # Two grayscale bars below
    ax_pos = ax.get_position()
    for k, (d3a, p3a, bar_y) in enumerate([
        (d3_arr1, p_d3_arr1, ax_pos.y0 - 0.06),
        (d3_arr2, p_d3_arr2, ax_pos.y0 - 0.11),
    ]):
        d3_vis = d3a[(d3a >= XI[0]) & (d3a <= XI[-1])]
        p_vis  = p3a[(d3a >= XI[0]) & (d3a <= XI[-1])]
        p_norm = p_vis / np.max(p_vis) if np.max(p_vis) > 0 else p_vis
        ax_bar = fig.add_axes([ax_pos.x0, bar_y, ax_pos.width, 0.03])
        ax_bar.imshow(p_norm[np.newaxis, :], aspect='auto', origin='lower',
                      cmap='Greys', vmin=0, vmax=1,
                      extent=[d3_vis[0], d3_vis[-1], 0, 1])
        ax_bar.set_yticks([])
        ax_bar.set_xlim(0, 0.8)
        bar_label = [label1, label2][k]
        if k == 1:
            ax_bar.set_xlabel(f"DG2-P distance (µm)", fontsize=8)
            ax_bar.set_title(bar_label, fontsize=8, loc='left', pad=2)
        else:
            ax_bar.set_xticklabels([])
            ax_bar.set_title(bar_label, fontsize=8, loc='left', pad=2)

    ax.set_xlim(XI[0], XI[-1])
    ax.set_ylim(YI[0], YI[-1])
    ax.set_xlabel("x (µm)", fontsize=12)
    ax.set_ylabel("y (µm)", fontsize=12)
    ax.set_title(f"{label1}\n– {label2}", fontsize=12, fontweight='bold')

    # Global Simes p-value below
    p_str = f"Global p = {simes_p:.4f}" if simes_p >= 0.0001 else f"Global p = {simes_p:.2e}"
    if driven_by:
        annotation = f"{p_str}\ndriven by {driven_by}"
    else:
        annotation = f"{p_str}\nnot significant"
    ax.text(0.5, -0.35, annotation,
            transform=ax.transAxes, ha='center', fontsize=10,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

    return cf, fig, ax


# ============================================================
# SIMES GLOBAL P
# ============================================================
def simes_global_p(pvals):
    p = np.sort(np.clip(np.asarray(pvals, float), 0, 1))
    m = len(p)
    scaled = (m / np.arange(1, m + 1)) * p
    return float(np.min(scaled))


def run_ks_tests(df1, df2):
    names = ["E6-DG2", "E6-P", "DG2-P"]
    pvals = []
    for col in ["dist_DG_E", "dist_P_E", "dist_P_DG"]:
        p = ks_2samp(df1[col].values, df2[col].values, alternative='two-sided').pvalue
        pvals.append(p)
    m = len(pvals)
    order = np.argsort(pvals)
    p_sorted = np.array(pvals)[order]
    ranks = np.arange(1, m + 1)
    q = np.minimum.accumulate(((m / ranks) * p_sorted)[::-1])[::-1]
    padj = np.empty(m)
    padj[order] = np.minimum(1.0, q)

    sig_names = [names[i] for i in range(m) if padj[i] < 0.05]
    driven_by = ", ".join(sig_names) if sig_names else ""
    global_p  = simes_global_p(pvals)
    return pvals, padj, global_p, driven_by


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Load data ---
    data = load_distances(XLSX_PATH, SHEET, NC14_CSV)
    CONDITIONS = ["nc14", "Non-expressing"]

    # --- Compute KDEs ---
    print("\nComputing KDEs...")
    x_range = np.linspace(0, MDT, 500)
    kdes = {}
    for cond in CONDITIONS:
        df = data[cond]
        kdes[cond] = {
            "d1": compute_kde(df["dist_DG_E"].values, x_range),
            "d2": compute_kde(df["dist_P_E"].values,  x_range),
            "d3": compute_kde(df["dist_P_DG"].values, x_range),
        }

    # --- Compute Z surfaces ---
    print("\nComputing contour surfaces...")
    surfaces = {}
    for cond in CONDITIONS:
        print(f"\n{cond}:")
        d1_x, d1_y = kdes[cond]["d1"]
        d2_x, d2_y = kdes[cond]["d2"]
        d3_x, d3_y = kdes[cond]["d3"]
        Z, X_grid, Y_grid = compute_weighted_Z(
            d1_x, d2_x, d1_y, d2_y,
            d3_x, d3_y, XI, YI,
            num_samples=50, subsample_stride=3
        )
        peak_idx = np.unravel_index(np.argmax(Z), Z.shape)
        peak_xy  = np.array([X_grid[peak_idx], Y_grid[peak_idx]])
        d3_mode  = d3_x[np.argmax(d3_y)]
        surfaces[cond] = {
            "Z": Z, "X_grid": X_grid, "Y_grid": Y_grid,
            "peak": peak_xy, "d3_mode": d3_mode,
            "d3_arr": d3_x, "p_d3": d3_y
        }

    # --- Individual contour plots ---
    print("\nPlotting individual contours...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.subplots_adjust(bottom=0.2, wspace=0.35)
    vmax_global = max(s["Z"].max() for s in surfaces.values())

    for ax, cond in zip(axes, CONDITIONS):
        s = surfaces[cond]
        plot_contour(
            s["Z"], s["X_grid"], s["Y_grid"],
            s["d3_arr"], s["p_d3"],
            cond,
            vmin=0, vmax=vmax_global,
            ax=ax, fig=fig
        )

    fig.suptitle("Chromatin conformation — svb locus\nP at origin | DG2 on x-axis | E6 in positive y",
                 fontsize=14, y=1.01)

    out_contour = os.path.join(OUTPUT_DIR, "plotted_contour_nc14_vs_WT_nonexp_NEW_DATASET.pdf")
    # Save as PDF (standard) and SVG (editable text in Illustrator)
    fig.savefig(out_contour, bbox_inches='tight', dpi=150)
    fig.savefig(os.path.splitext(out_contour)[0] + '.svg', bbox_inches='tight')
    fig.savefig(out_contour.replace('.pdf', '.svg'), bbox_inches='tight')
    print(f"Saved: {out_contour}")
    plt.close(fig)

    # --- Subtraction contour: Non-expressing vs nc14 ---
    print("\nComputing subtraction contour...")
    s1, s2 = surfaces["nc14"], surfaces["Non-expressing"]
    pvals, padj, global_p, driven_by = run_ks_tests(
        data["nc14"], data["Non-expressing"]
    )
    print(f"nc14 vs Non-expressing: global p = {global_p:.4e}, driven by {driven_by}")

    Z_diff_max = np.max(np.abs(s1["Z"] - s2["Z"]))

    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 8))
    fig2.subplots_adjust(bottom=0.28)

    plot_subtraction(
        s1["Z"], s2["Z"], s1["X_grid"], s1["Y_grid"],
        s1["peak"], s2["peak"],
        s1["d3_mode"], s2["d3_mode"],
        s1["d3_arr"], s1["p_d3"],
        s2["d3_arr"], s2["p_d3"],
        "nc14", "Non-expressing",
        global_p, driven_by,
        vmin=-Z_diff_max, vmax=Z_diff_max,
        ax=ax2, fig=fig2
    )

    fig2.suptitle("Subtraction contour — nc14 vs Non-expressing\nsvb locus",
                  fontsize=14, y=1.01)

    out_sub = os.path.join(OUTPUT_DIR, "plotted_subtraction_nc14_vs_WT_nonexp_NEW_DATASET.pdf")
    # Save as PDF and SVG
    fig2.savefig(out_sub, bbox_inches='tight', dpi=150)
    fig2.savefig(os.path.splitext(out_sub)[0] + '.svg', bbox_inches='tight')
    fig2.savefig(out_sub.replace('.pdf', '.svg'), bbox_inches='tight')
    print(f"Saved: {out_sub}")
    plt.close(fig2)

    print("\nDONE.")
