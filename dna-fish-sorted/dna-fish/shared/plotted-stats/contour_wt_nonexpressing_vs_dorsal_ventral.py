"""
PLOTTED contour generation for svb locus DNA-FISH data
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
# INPUTS — NEW COMBINED WT DATASET
# ============================================================
XLSX_PATH = (
    "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/"
    "AP_file_2026/ROI/"
    "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

SHEET = "Filtered_Triplets"
WT_GENE_ID = "SXLGFP"

OUTPUT_DIR = "F:/Thesis/Figures/PLOTTED_WT_NEW_DATASET"
MDT = 2.0   # maximum distance threshold, µm

ROI_MAP = {
    "nonexp": "Non-expressing",
    "non-exp": "Non-expressing",
    "non_exp": "Non-expressing",
    "nonexpressing": "Non-expressing",
    "non-expressing": "Non-expressing",
    "dorsal": "Dorsal",
    "ventral": "Ventral",
}

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
def load_distances(xlsx_path, sheet):
    """
    Load WT SXLGFP rows from the new combined Filtered_Triplets sheet,
    split them by ROI condition, and calculate the three pairwise distances.
    """

    print(f"Loading WT data from: {xlsx_path}")
    print(f"Sheet: {sheet}")

    df = pd.read_excel(
        xlsx_path,
        sheet_name=sheet
    )

    required = {
        "gene_id",
        "roi",
        "x_svb_um",
        "y_svb_um",
        "z_svb_um",
        "x_E_um",
        "y_E_um",
        "z_E_um",
        "x_DG_um",
        "y_DG_um",
        "z_DG_um",
    }

    missing = required.difference(df.columns)

    if missing:
        raise KeyError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    # Keep WT only.
    df = df[
        df["gene_id"]
        .astype(str)
        .str.strip()
        .str.upper()
        .eq(WT_GENE_ID)
    ].copy()

    if df.empty:
        raise ValueError(
            f"No rows found for gene_id == {WT_GENE_ID!r}"
        )

    # Normalize ROI strings, including trailing spaces.
    df["roi_key"] = (
        df["roi"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    df["Condition"] = df["roi_key"].map(ROI_MAP)

    print("\nRaw WT ROI values:")
    print(
        sorted(
            df["roi"]
            .dropna()
            .astype(str)
            .unique()
        )
    )

    print("\nMapped WT condition counts before distance filtering:")
    print(
        df["Condition"].value_counts(
            dropna=False
        )
    )

    # Calculate pairwise distances directly from µm coordinates.
    df["dist_P_DG"] = np.sqrt(
        (pd.to_numeric(df["x_svb_um"], errors="coerce") -
         pd.to_numeric(df["x_DG_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["y_svb_um"], errors="coerce") -
         pd.to_numeric(df["y_DG_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["z_svb_um"], errors="coerce") -
         pd.to_numeric(df["z_DG_um"], errors="coerce")) ** 2
    )

    df["dist_P_E"] = np.sqrt(
        (pd.to_numeric(df["x_svb_um"], errors="coerce") -
         pd.to_numeric(df["x_E_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["y_svb_um"], errors="coerce") -
         pd.to_numeric(df["y_E_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["z_svb_um"], errors="coerce") -
         pd.to_numeric(df["z_E_um"], errors="coerce")) ** 2
    )

    df["dist_DG_E"] = np.sqrt(
        (pd.to_numeric(df["x_DG_um"], errors="coerce") -
         pd.to_numeric(df["x_E_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["y_DG_um"], errors="coerce") -
         pd.to_numeric(df["y_E_um"], errors="coerce")) ** 2
        +
        (pd.to_numeric(df["z_DG_um"], errors="coerce") -
         pd.to_numeric(df["z_E_um"], errors="coerce")) ** 2
    )

    before_dropna = len(df)

    df = df.dropna(
        subset=[
            "Condition",
            "dist_P_DG",
            "dist_P_E",
            "dist_DG_E",
        ]
    ).copy()

    print(
        "\nRows removed because of missing condition or distance: "
        f"{before_dropna - len(df)}"
    )

    before_mdt = len(df)

    df = df[
        (df["dist_P_DG"] > 0) &
        (df["dist_P_DG"] <= MDT) &
        (df["dist_P_E"] > 0) &
        (df["dist_P_E"] <= MDT) &
        (df["dist_DG_E"] > 0) &
        (df["dist_DG_E"] <= MDT)
    ].copy()

    print(
        "Rows removed by the 0 < distance <= MDT filter: "
        f"{before_mdt - len(df)}"
    )

    results = {}

    for condition in [
        "Non-expressing",
        "Dorsal",
        "Ventral",
    ]:
        subset = df[
            df["Condition"] == condition
        ][
            [
                "dist_P_DG",
                "dist_P_E",
                "dist_DG_E",
            ]
        ].copy()

        print(f"{condition}: n = {len(subset)}")

        if len(subset) < 2:
            raise ValueError(
                f"Insufficient WT data for {condition}: n={len(subset)}"
            )

        results[condition] = subset

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
    data = load_distances(XLSX_PATH, SHEET)

    # --- Compute KDEs ---
    print("\nComputing KDEs...")
    x_range = np.linspace(0, MDT, 500)
    kdes = {}
    for cond, df in data.items():
        kdes[cond] = {
            "d1": compute_kde(df["dist_DG_E"].values, x_range),  # E6-DG2
            "d2": compute_kde(df["dist_P_E"].values,  x_range),  # E6-P
            "d3": compute_kde(df["dist_P_DG"].values, x_range),  # DG2-P (x-axis)
        }

    # --- Compute Z surfaces ---
    print("\nComputing contour surfaces (this takes a few minutes)...")
    surfaces = {}
    for cond in ["Non-expressing", "Dorsal", "Ventral"]:
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
    fig, axes = plt.subplots(1, 3, figsize=(21, 7))
    fig.subplots_adjust(bottom=0.2, wspace=0.35)

    vmax_global = max(s["Z"].max() for s in surfaces.values())

    for ax, cond in zip(axes, ["Non-expressing", "Dorsal", "Ventral"]):
        s = surfaces[cond]
        plot_contour(
            s["Z"], s["X_grid"], s["Y_grid"],
            s["d3_arr"], s["p_d3"],
            cond,
            vmin=0, vmax=vmax_global,
            ax=ax, fig=fig
        )

    fig.suptitle("Chromatin conformation — WT svb locus\nP at origin | DG2 on x-axis | E6 in positive y",
                 fontsize=14, y=1.01)

    out_contour = os.path.join(OUTPUT_DIR, "plotted_contour_WT_NEW_DATASET.pdf")
    fig.savefig(out_contour, bbox_inches='tight', dpi=150)
    fig.savefig(out_contour.replace('.pdf', '.svg'), bbox_inches='tight')
    print(f"Saved: {out_contour}")

    # --- Subtraction contours ---
    print("\nComputing subtraction contours...")
    comparisons = [
        ("Non-expressing", "Dorsal"),
        ("Non-expressing", "Ventral"),
        ("Dorsal",         "Ventral"),
    ]

    fig2, axes2 = plt.subplots(1, 3, figsize=(21, 8))
    fig2.subplots_adjust(bottom=0.30, wspace=0.4)

    Z_diff_max = max(
        np.max(np.abs(surfaces[c1]["Z"] - surfaces[c2]["Z"]))
        for c1, c2 in comparisons
    )

    for ax, (cond1, cond2) in zip(axes2, comparisons):
        s1, s2 = surfaces[cond1], surfaces[cond2]
        pvals, padj, global_p, driven_by = run_ks_tests(data[cond1], data[cond2])
        print(f"\n{cond1} vs {cond2}: global p = {global_p:.4e}, driven by {driven_by}")

        plot_subtraction(
            s1["Z"], s2["Z"], s1["X_grid"], s1["Y_grid"],
            s1["peak"], s2["peak"],
            s1["d3_mode"], s2["d3_mode"],
            s1["d3_arr"], s1["p_d3"],
            s2["d3_arr"], s2["p_d3"],
            cond1, cond2,
            global_p, driven_by,
            vmin=-Z_diff_max, vmax=Z_diff_max,
            ax=ax, fig=fig2
        )

    fig2.suptitle("Subtraction contours — WT svb locus",
                  fontsize=14, y=1.01)

    out_sub = os.path.join(OUTPUT_DIR, "plotted_subtraction_WT_NEW_DATASET.pdf")
    fig2.savefig(out_sub, bbox_inches='tight', dpi=150)
    fig2.savefig(out_sub.replace('.pdf', '.svg'), bbox_inches='tight')
    print(f"Saved: {out_sub}")

    print("\nDONE.")
