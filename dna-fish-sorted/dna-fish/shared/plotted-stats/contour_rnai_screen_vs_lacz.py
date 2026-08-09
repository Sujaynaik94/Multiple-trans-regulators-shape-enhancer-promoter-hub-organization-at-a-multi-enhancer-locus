"""
PLOTTED contour generation for RNAi screen data
Each gene vs LacZ control, per condition (Non-expressing, Dorsal, Ventral)
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, ks_2samp
from scipy.interpolate import NearestNDInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# INPUTS
# ============================================================
XLSX_PATH  = "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/AP_file_2026/ROI/ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
SHEET      = "Filtered_Triplets"
OUTPUT_DIR = "F:/Thesis/Figures/PLOTTED_RNAi_NEW_DATASET"
MDT        = 2.0

# ============================================================
# FONT SETTINGS
# ============================================================
FONT_SIZE = 14

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
# GENE ORDER + RECODE
# ============================================================
GENE_RECODE = {
    "LACZ":    "LacZ",
    "CTCF":    "CTCF",
    "SUHW":    "Su(Hw)",
    "MDG4":    "Mod(mdg4)",
    "BEAF-32": "BEAF-32",
    "BEAF32":  "BEAF-32",
    "BLIMP1":  "Blimp-1",
    "LOLAJ":   "lola",
    "LOLA":    "lola",
    "PNR":     "pnr",
    "TTK":     "ttk",
    "UBX":     "Ubx",
    "EXD":     "Exd",
    "HTH":     "hth",
    "MED13":   "Med13",
    "AWH":     "Awh",
    "BOWL":    "bowl",
    "CP190":   "CP190",
    "PNT":     "pnt",
}

PLOT_ORDER = [
    "LacZ",
    "Awh",
    "pnr",
    "Ubx",
    "Blimp-1",
    "ttk",
    "lola",
    "Exd",
    "CTCF",
    "BEAF-32",
    "Su(Hw)",
    "Mod(mdg4)",
    "CP190",
    "bowl",
    "pnt",
    "hth",
    "Med13",
]

CONTROL_GENE = "LacZ"

ROI_MAP = {
    "nonexp":                  "Non-expressing",
    "non-exp":                 "Non-expressing",
    "non_exp":                 "Non-expressing",
    "nonexpressing":           "Non-expressing",
    "non-expressing":          "Non-expressing",
    "shavenbaby off":          "Non-expressing",

    "dorsal":                  "Dorsal",
    "shavenbaby on (dorsal)":  "Dorsal",

    "ventral":                 "Ventral",
    "shavenbaby on (ventral)": "Ventral",
}

CONDITIONS = ["Non-expressing", "Dorsal", "Ventral"]

GENE_COLORS = {
    "LacZ":      "#686765",
    "Awh":       "#008088",
    "pnr":       "#008088",
    "Ubx":       "#008088",
    "Blimp-1":   "#008088",
    "ttk":       "#008088",
    "lola":      "#008088",
    "Exd":       "#008088",
    "CTCF":      "#E07B39",
    "BEAF-32":   "#E07B39",
    "Su(Hw)":    "#E07B39",
    "Mod(mdg4)": "#E07B39",
}

XI = np.linspace(-0.2, 1.0, 120)
YI = np.linspace(0.0, 0.8, 90)

DIVERGING_RB = LinearSegmentedColormap.from_list(
    "rwb",
    [(0.0, "#053061"), (0.25, "#4393c3"), (0.5, "white"),
     (0.75, "#d6604d"), (1.0, "#67001f")]
)

KDE_RANGE = np.linspace(0.001, 1.5, 500)

def load_rnai_data(xlsx_path, sheet):
    print(f"Loading RNAi data from {sheet}...")

    df = pd.read_excel(
        xlsx_path,
        sheet_name=sheet
    )

    required_columns = {
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

    missing_columns = required_columns.difference(df.columns)

    if missing_columns:
        raise KeyError(
            "Missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    # ========================================================
    # CLEAN GENE IDS FROM THE NEW DATASET
    #
    # Examples:
    # LACZ_MASKS     -> LACZ     -> LacZ
    # AWH_MASKS      -> AWH      -> Awh
    # BEAF-32_MASKS  -> BEAF-32  -> BEAF-32
    # ========================================================

    df["gene_key"] = (
        df["gene_id"]
        .astype(str)
        .str.strip()
        .str.replace(r"_MASKS$", "", regex=True, case=False)
        .str.replace(r"_[0-9]+$", "", regex=True)
        .str.upper()
    )

    df["gene"] = (
        df["gene_key"]
        .map(GENE_RECODE)
        .fillna(df["gene_key"])
    )

    print("\nGene ID to cleaned-gene mapping:")

    print(
        df[
            ["gene_id", "gene_key", "gene"]
        ]
        .drop_duplicates()
        .sort_values("gene_id")
        .to_string(index=False)
    )

    # ========================================================
    # CLEAN ROI VALUES
    # Handles nonexp, non-exp and trailing spaces such as
    # "ventral "
    # ========================================================

    df["roi_key"] = (
        df["roi"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    df["Condition"] = df["roi_key"].map(ROI_MAP)

    print("\nRaw ROI values:")

    print(
        sorted(
            df["roi"]
            .dropna()
            .astype(str)
            .unique()
        )
    )

    print("\nMapped condition counts:")

    print(
        df["Condition"].value_counts(
            dropna=False
        )
    )

    # ========================================================
    # CALCULATE THE THREE PAIRWISE DISTANCES
    # ========================================================

    df["dist_P_DG"] = np.sqrt(
        (df["x_svb_um"] - df["x_DG_um"])**2 +
        (df["y_svb_um"] - df["y_DG_um"])**2 +
        (df["z_svb_um"] - df["z_DG_um"])**2
    )

    df["dist_P_E"] = np.sqrt(
        (df["x_svb_um"] - df["x_E_um"])**2 +
        (df["y_svb_um"] - df["y_E_um"])**2 +
        (df["z_svb_um"] - df["z_E_um"])**2
    )

    df["dist_DG_E"] = np.sqrt(
        (df["x_DG_um"] - df["x_E_um"])**2 +
        (df["y_DG_um"] - df["y_E_um"])**2 +
        (df["z_DG_um"] - df["z_E_um"])**2
    )

    # ========================================================
    # FILTER
    # ========================================================

    before_dropna = len(df)

    df = df.dropna(
        subset=[
            "dist_P_DG",
            "dist_P_E",
            "dist_DG_E",
            "Condition",
            "gene",
        ]
    )

    print(
        f"\nRows removed because of missing distance/condition/gene: "
        f"{before_dropna - len(df)}"
    )

    before_distance_filter = len(df)

    df = df[
        (df["dist_P_DG"] <= MDT) & (df["dist_P_DG"] > 0) &
        (df["dist_P_E"]  <= MDT) & (df["dist_P_E"]  > 0) &
        (df["dist_DG_E"] <= MDT) & (df["dist_DG_E"] > 0) &
        (df["gene"].isin(PLOT_ORDER))
    ].copy()

    print(
        f"Rows removed by distance/gene filtering: "
        f"{before_distance_filter - len(df)}"
    )

    print(f"\nFiltered rows: {len(df)}")

    print("\nFiltered rows by condition:")

    print(
        df.groupby(
            "Condition",
            observed=True
        ).size()
    )

    print("\nFiltered rows by gene:")

    print(
        df.groupby(
            "gene",
            observed=True
        ).size()
    )

    print("\nLacZ rows by condition:")

    print(
        df[
            df["gene"] == CONTROL_GENE
        ]
        .groupby(
            "Condition",
            observed=True
        )
        .size()
    )

    return df

def compute_kde(data, x_range):
    data = np.asarray(data, float)
    data = data[np.isfinite(data) & (data > 0)]
    if len(data) < 10:
        return x_range, np.zeros_like(x_range)

    kde = gaussian_kde(data, bw_method="scott")
    y = kde(x_range)
    area = np.trapezoid(y, x_range)
    if area > 0:
        y /= area
    return x_range, y

def compute_triangle_vertices_from_PDF(d1_arr, d2_arr, p_d1_arr, p_d2_arr, d3):
    D1, D2 = np.meshgrid(d1_arr, d2_arr, indexing="ij")
    P1, P2 = np.meshgrid(p_d1_arr, p_d2_arr, indexing="ij")
    p = P1 * P2

    valid = (D1 + D2 >= d3) & (D1 + d3 >= D2) & (D2 + d3 >= D1)

    x_E6 = np.where(
        valid & (d3 > 0),
        (d3**2 + D2**2 - D1**2) / (2 * d3),
        np.nan
    )
    y2 = np.where(valid, D2**2 - x_E6**2, np.nan)
    y_E6 = np.where(
        np.isfinite(y2) & (y2 >= 0),
        np.sqrt(np.maximum(y2, 0)),
        np.nan
    )

    return x_E6, y_E6, np.where(valid, p, np.nan)

def compute_weighted_Z(d1_arr, d2_arr, p_d1_arr, p_d2_arr,
                       d3_arr, p_d3_arr,
                       num_samples=50, subsample_stride=3):
    X_grid, Y_grid = np.meshgrid(XI, YI)

    indices_even = np.linspace(0, len(d3_arr) - 1, num_samples, dtype=int)
    d3_sampled = d3_arr[indices_even]
    p_d3_sampled = p_d3_arr[indices_even]
    p_d3_sampled = p_d3_sampled / (np.sum(p_d3_sampled) + 1e-12)

    bin_edges = np.concatenate([
        [d3_sampled[0] - (d3_sampled[1] - d3_sampled[0]) / 2],
        (d3_sampled[:-1] + d3_sampled[1:]) / 2,
        [d3_sampled[-1] + (d3_sampled[-1] - d3_sampled[-2]) / 2]
    ])

    bin_widths = np.diff(bin_edges)
    weights = p_d3_sampled * bin_widths
    weights = weights / (np.sum(weights) + 1e-12)

    d1_sub = d1_arr[::subsample_stride]
    p_d1_sub = p_d1_arr[::subsample_stride]
    d2_sub = d2_arr[::subsample_stride]
    p_d2_sub = p_d2_arr[::subsample_stride]

    Z_weighted = np.zeros_like(X_grid, dtype=float)
    total_weight = 0.0
    dx = XI[1] - XI[0]
    dy = YI[1] - YI[0]

    for d3, w in zip(d3_sampled, weights):
        if d3 <= 0:
            continue

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
            Z_interp /= area
            Z_weighted += Z_interp * w

        total_weight += w

    if total_weight > 0:
        Z_weighted /= total_weight

    return Z_weighted, X_grid, Y_grid

def compute_surface(df_sub, label):
    if len(df_sub) < 30:
        print(f"SKIP {label}: n={len(df_sub)} too few")
        return None

    _, p_d1 = compute_kde(df_sub["dist_DG_E"].values, KDE_RANGE)
    _, p_d2 = compute_kde(df_sub["dist_P_E"].values, KDE_RANGE)
    _, p_d3 = compute_kde(df_sub["dist_P_DG"].values, KDE_RANGE)

    Z, X_grid, Y_grid = compute_weighted_Z(
        KDE_RANGE, KDE_RANGE, p_d1, p_d2,
        KDE_RANGE, p_d3,
        num_samples=50,
        subsample_stride=3
    )

    peak_idx = np.unravel_index(np.argmax(Z), Z.shape)
    peak_xy = np.array([X_grid[peak_idx], Y_grid[peak_idx]])
    d3_mode = KDE_RANGE[np.argmax(p_d3)]

    return {
        "Z": Z,
        "X_grid": X_grid,
        "Y_grid": Y_grid,
        "peak": peak_xy,
        "d3_mode": d3_mode,
        "d3_arr": KDE_RANGE,
        "p_d3": p_d3,
        "n": len(df_sub)
    }

def simes_global_p(pvals):
    p = np.sort(np.clip(np.asarray(pvals, float), 0, 1))
    m = len(p)
    return float(np.min((m / np.arange(1, m + 1)) * p))

def run_ks_tests(df1, df2):
    names = ["E6-DG2", "E6-P", "DG2-P"]
    pvals = []

    for col in ["dist_DG_E", "dist_P_E", "dist_P_DG"]:
        p = ks_2samp(
            df1[col].values,
            df2[col].values,
            alternative="two-sided"
        ).pvalue
        pvals.append(p)

    pvals = np.array(pvals)
    m = len(pvals)
    order = np.argsort(pvals)
    p_sorted = pvals[order]
    q = np.minimum.accumulate(((m / np.arange(1, m + 1)) * p_sorted)[::-1])[::-1]

    padj = np.empty(m)
    padj[order] = np.minimum(1.0, q)

    sig = [names[i] for i in range(m) if padj[i] < 0.05]
    driven_by = ", ".join(sig) if sig else ""

    return simes_global_p(pvals), driven_by, pvals, padj


def plot_subtraction_panel(ax,
                           s_ctrl, s_gene,
                           ctrl_label, gene_label,
                           global_p, driven_by,
                           gene_color,
                           show_xlabel=False,
                           show_ylabel=False,
                           vmin=-1.8, vmax=1.8):
    """Draw one proportional contour panel with KDE bars aligned to the contour axes."""

    Z_diff = s_ctrl["Z"] - s_gene["Z"]
    x_mask = (s_ctrl["X_grid"] > 0.90) | (s_ctrl["X_grid"] < -0.15)
    Z_diff = Z_diff.copy()
    Z_diff[x_mask] = 0.0

    ax.contourf(
        s_ctrl["X_grid"], s_ctrl["Y_grid"], Z_diff,
        levels=np.linspace(vmin, vmax, 20),
        cmap=DIVERGING_RB, vmin=vmin, vmax=vmax
    )

    promoter = np.array([0.0, 0.0])
    for peak, d3_mode, marker, is_ctrl in [
        (s_ctrl["peak"], s_ctrl["d3_mode"], "o", True),
        (s_gene["peak"], s_gene["d3_mode"], "s", False),
    ]:
        dg2 = np.array([d3_mode, 0.0])
        linestyle = "-" if is_ctrl else "--"
        linewidth = 1.8 if is_ctrl else 1.6

        ax.plot([promoter[0], dg2[0]], [promoter[1], dg2[1]],
                color="black", ls=linestyle, lw=linewidth, zorder=4)
        ax.plot([promoter[0], peak[0]], [promoter[1], peak[1]],
                color="black", ls=linestyle, lw=linewidth, zorder=4)
        ax.plot([dg2[0], peak[0]], [dg2[1], peak[1]],
                color="black", ls=linestyle, lw=linewidth, zorder=4)

        fill_p = "#8B3A8B" if is_ctrl else "none"
        fill_dg = "#F5C842" if is_ctrl else "none"
        fill_e = "#4FB3E8" if is_ctrl else "none"
        ax.plot(*promoter, marker, color="#8B3A8B", ms=10.5,
                mfc=fill_p, mew=1.8, zorder=6)
        ax.plot(*dg2, marker, color="#F5C842", ms=10.5,
                mfc=fill_dg, mew=1.8, zorder=6)
        ax.plot(*peak, marker, color="#4FB3E8", ms=10.5,
                mfc=fill_e, mew=1.8, zorder=6)

    p_str = f"p={global_p:.4f}" if global_p >= 0.0001 else f"p={global_p:.1e}"
    p_ann = f"{p_str}\n{driven_by}" if driven_by else f"{p_str}\nns"
    ax.text(
        0.975, 0.965, p_ann,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=12, color="black" if driven_by else "grey",
        bbox=dict(facecolor="white", edgecolor="grey", alpha=0.88,
                  pad=2.2, linewidth=0.6), zorder=10
    )

    ax.set_xlim(XI[0], XI[-1])
    ax.set_ylim(YI[0], YI[-1])
    # Equal data scaling keeps all configurations geometrically proportional.
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=14, length=3, pad=2)
    ax.set_xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])

    # Keep numerical x-axis tick labels on every contour panel.
    # Only the bottom row receives the axis title.
    if show_xlabel:
        ax.set_xlabel("E6 x-position (µm)", fontsize=14, labelpad=4)
    if show_ylabel:
        ax.set_ylabel("E6 y-position (µm)", fontsize=14, labelpad=4)
    else:
        ax.set_yticklabels([])

    # KDE bars are inset relative to the FINAL contour-axis box.
    # This exactly matches their physical width to the contour panel, even
    # after equal-aspect scaling has reduced the contour width in the grid cell.
    x_bar = np.linspace(0, 0.8, 240)
    # Leave a dedicated gap below the contour x tick labels, then draw
    # slightly taller bars. This prevents the tick numbers from being hidden.
    bar_specs = [
        (-0.205, s_ctrl["d3_arr"], s_ctrl["p_d3"], ctrl_label),
        (-0.285, s_gene["d3_arr"], s_gene["p_d3"], gene_label),
    ]
    bar_axes = []
    for y0, d3_arr, p_d3, label in bar_specs:
        bar_ax = ax.inset_axes([0.0, y0, 1.0, 0.052], transform=ax.transAxes)
        bar_axes.append(bar_ax)
        p_interp = np.interp(x_bar, d3_arr, p_d3)
        p_norm = p_interp / (np.max(p_interp) + 1e-10)
        bar_ax.imshow(
            p_norm[np.newaxis, :], aspect="auto", origin="lower",
            cmap="Greys", vmin=0, vmax=1, extent=[0, 0.8, 0, 1]
        )
        bar_ax.set_xlim(0, 0.8)
        bar_ax.set_ylim(0, 1)
        bar_ax.set_yticks([])
        bar_ax.text(0.012, 0.5, label, transform=bar_ax.transAxes,
                    fontsize=11, va="center", ha="left",
                    color="white", fontweight="bold")
        for spine in bar_ax.spines.values():
            spine.set_linewidth(0.65)

    bar_axes[0].set_xticks([])
    bar_axes[1].set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    bar_axes[1].tick_params(axis="x", labelsize=14, length=2.2, pad=1)
    # Show the DG2-P numerical scale for every panel; keep the title only
    # on the bottom row to avoid repeated text and overlap.
    if show_xlabel:
        bar_axes[1].set_xlabel("DG2-P distance (µm)", fontsize=13, labelpad=3)


def precompute_surfaces(df, genes):
    all_surfaces = {}
    all_ctrl = {}
    all_ctrl_df = {}
    z_diff_maxima = []

    for condition in CONDITIONS:
        df_cond = df[df["Condition"] == condition]
        df_ctrl = df_cond[df_cond["gene"] == CONTROL_GENE]
        all_ctrl_df[condition] = df_ctrl
        s_ctrl = compute_surface(df_ctrl, f"LacZ_{condition}")
        all_ctrl[condition] = s_ctrl
        all_surfaces[condition] = {}

        if s_ctrl is None:
            continue

        for gene in genes:
            df_gene = df_cond[df_cond["gene"] == gene]
            s_gene = compute_surface(df_gene, f"{gene}_{condition}")
            all_surfaces[condition][gene] = s_gene
            if s_gene is not None:
                z_diff_maxima.append(np.max(np.abs(s_ctrl["Z"] - s_gene["Z"])))

    vmax = np.percentile(z_diff_maxima, 90) if z_diff_maxima else 1.8
    return all_ctrl, all_ctrl_df, all_surfaces, float(vmax)


def make_group_figure(df, genes, group_name, output_stem):
    """Create one 16 x 20 inch figure without panel or label overlap."""
    genes = [g for g in genes if g in df["gene"].values]
    if not genes:
        print(f"No genes available for {group_name}; figure skipped.")
        return

    all_ctrl, all_ctrl_df, all_surfaces, vmax = precompute_surfaces(df, genes)

    fig = plt.figure(figsize=(16, 20), constrained_layout=False)
    outer = fig.add_gridspec(
        nrows=len(genes), ncols=len(CONDITIONS),
        left=0.075, right=0.965, bottom=0.105, top=0.925,
        wspace=0.16, hspace=0.40
    )

    contour_axes = {}
    for row_idx, gene in enumerate(genes):
        contour_axes[gene] = []
        gene_color = GENE_COLORS.get(gene, "#008088")

        for col_idx, condition in enumerate(CONDITIONS):
            # One axis per cell. The two KDE bars are added as aligned insets
            # inside plot_subtraction_panel(), exactly as in the original code.
            ax = fig.add_subplot(outer[row_idx, col_idx])
            contour_axes[gene].append(ax)

            s_ctrl = all_ctrl.get(condition)
            s_gene = all_surfaces.get(condition, {}).get(gene)

            if s_ctrl is None or s_gene is None:
                ax.text(0.5, 0.5, "insufficient data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=17, color="grey")
                ax.set_xlim(XI[0], XI[-1]); ax.set_ylim(YI[0], YI[-1])
                ax.set_aspect("equal", adjustable="box")
            else:
                df_ctrl = all_ctrl_df[condition]
                df_gene = df[(df["Condition"] == condition) & (df["gene"] == gene)]
                global_p, driven_by, _, _ = run_ks_tests(df_ctrl, df_gene)
                plot_subtraction_panel(
                    ax,
                    s_ctrl, s_gene, CONTROL_GENE, gene,
                    global_p, driven_by, gene_color,
                    show_xlabel=(row_idx == len(genes) - 1),
                    show_ylabel=(col_idx == 0),
                    vmin=-vmax, vmax=vmax
                )

        # One clean row label, outside the plotting region.
        left_ax = contour_axes[gene][0]
        pos = left_ax.get_position()
        fig.text(0.025, pos.y0 + pos.height / 2, gene,
                 ha="center", va="center", rotation=90,
                 fontsize=20, fontweight="bold", color=gene_color)

    condition_colors = {
        "Non-expressing": "#686765",
        "Dorsal": "#cc2829",
        "Ventral": "#2b3f99",
    }
    for col_idx, condition in enumerate(CONDITIONS):
        pos = contour_axes[genes[0]][col_idx].get_position()
        fig.text(pos.x0 + pos.width / 2, 0.947, condition,
                 ha="center", va="center", fontsize=22,
                 fontweight="bold", color=condition_colors[condition])

    fig.suptitle(
        f"{group_name}: LacZ versus RNAi knockdown subtraction contours\n"
        "Red = enriched in LacZ | Blue = enriched in knockdown",
        fontsize=24, fontweight="bold", y=0.988
    )

    import matplotlib.lines as mlines
    legend_elements = [
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#4FB3E8", markeredgecolor="#4FB3E8",
                      markersize=9, label="LacZ MPP (filled circle)"),
        mlines.Line2D([0], [0], marker="s", color="w",
                      markerfacecolor="none", markeredgecolor="#4FB3E8",
                      markersize=9, markeredgewidth=1.7,
                      label="Knockdown MPP (open square)"),
        mlines.Line2D([0], [0], color="black", lw=1.8, linestyle="-",
                      label="LacZ mode triangle"),
        mlines.Line2D([0], [0], color="black", lw=1.6, linestyle="--",
                      label="Knockdown mode triangle"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
               fontsize=16, frameon=True, bbox_to_anchor=(0.5, 0.055))

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    cbar_ax = fig.add_axes([0.27, 0.025, 0.46, 0.012])
    sm = ScalarMappable(cmap=DIVERGING_RB, norm=Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks([-vmax, 0, vmax])
    cbar.set_ticklabels([f"{-vmax:.1f}", "0", f"{vmax:.1f}"])
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("LacZ − knockdown probability density", fontsize=16, labelpad=3)

    pdf_path = os.path.join(OUTPUT_DIR, output_stem + ".pdf")
    svg_path = os.path.join(OUTPUT_DIR, output_stem + ".svg")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved PDF: {pdf_path}")
    print(f"Saved SVG: {svg_path}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_rnai_data(XLSX_PATH, SHEET)

    # Split the large TF set into two pages so each contour remains large.
    tf_part1 = [
        "Awh", "pnr", "Ubx", "Blimp-1", "ttk", "lola"
    ]
    tf_part2 = [
        "Exd", "bowl", "pnt", "hth", "Med13"
    ]
    boundary_genes = [
        "CTCF", "BEAF-32", "Su(Hw)", "Mod(mdg4)", "CP190"
    ]

    make_group_figure(
        df, tf_part1,
        "Transcription-factor RNAi — part 1",
        "plotted_RNAi_TF_part1_large_contours_16x20"
    )
    make_group_figure(
        df, tf_part2,
        "Transcription-factor and regulatory-protein RNAi — part 2",
        "plotted_RNAi_TF_part2_large_contours_16x20"
    )
    make_group_figure(
        df, boundary_genes,
        "Boundary-protein RNAi",
        "plotted_RNAi_BOUNDARY_large_contours_16x20"
    )
    print("DONE.")


# Called unconditionally so reticulate::source_python() executes the analysis.
main()



