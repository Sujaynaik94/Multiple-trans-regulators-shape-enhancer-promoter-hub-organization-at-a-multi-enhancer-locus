# ============================================================
# HUB DISPERSION
# nc14 WHOLE EMBRYO vs WT NON-EXPRESSING
#
# NO MANUAL EMBRYO EXCLUSIONS
#
# Statistics:
# - Wilcoxon rank-sum test
# - nucleus/triplet-level observations

##important: manually save in pdf cairo portrait mode with 6x7
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(gghalves)
  library(tibble)
})

# ============================================================
# INPUT FILES
# ============================================================

nc14_csv <- paste0(
  "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/",
  "AP_file_2026/triplets_fromnucleisegment/",
  "nc14_MASKS_FINAL_TRIPLETS.csv"
)

xlsx_path <- paste0(
  "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/",
  "AP_file_2026/ROI/",
  "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

sheet <- "Filtered_Triplets"
WT_GENE_ID <- "SXLGFP"

NONEXP_ALIASES <- c(
  "nonexp",
  "non_exp",
  "non-exp",
  "nonexpressing",
  "non-expressing",
  "ne"
)

RADIAL_CUTOFF <- 0.5

# ============================================================
# PSEUDO-LOG AXIS
# ============================================================

SIGMA <- 0.001

y_breaks <- c(
  0.1,
  0.2,
  0.5,
  1,
  2
)

# ============================================================
# GROUP ORDER AND COLOURS
# ============================================================

group_levels <- c(
  "nc14",
  "Non-expressing"
)

group_colors <- c(
  "nc14" = "#686765",
  "Non-expressing" = "#686765"
)

group_spacing <- tibble::tibble(
  Group = factor(
    group_levels,
    levels = group_levels
  ),
  xpos = c(1, 2)
)

# ============================================================
# PLOT STYLE
# ============================================================

box_w <- 0.35
violin_w <- 1.0

dot_band_w <- box_w / 2
dot_left_shift <- dot_band_w
dot_jit_x <- dot_band_w * 0.9

FONT_FAMILY <- "Arial"
FONT_AXIS_TITLE <- 22
FONT_AXIS_TEXT <- 22
FONT_N_LABEL <- 6
FONT_STARS <- 8

# ============================================================
# HELPER FUNCTIONS
# ============================================================

d3 <- function(x1, y1, z1, x2, y2, z2) {
  sqrt(
    (x1 - x2)^2 +
      (y1 - y2)^2 +
      (z1 - z2)^2
  )
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    NA_character_
  } else {
    hit[[1]]
  }
}

find_embryo_col <- function(df, label) {
  col <- pick_col(
    df,
    c(
      "FileName_Image",
      "filename_image",
      "FileName",
      "filename",
      "EmbryoID",
      "embryo_id",
      "embryo",
      "Embryo",
      "Image",
      "image",
      "triplet_source"
    )
  )

  if (is.na(col)) {
    stop(
      paste0(
        "Could not find embryo/image ID column in ",
        label,
        ".\nColumns present:\n",
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  col
}

compute_dispersion <- function(df, label) {
  x_DG <- pick_col(
    df,
    c("DG_CoM_X", "x_DG_um", "DG_x_um", "x_DG", "DG_x")
  )
  y_DG <- pick_col(
    df,
    c("DG_CoM_Y", "y_DG_um", "DG_y_um", "y_DG", "DG_y")
  )
  z_DG <- pick_col(
    df,
    c("DG_CoM_Z", "z_DG_um", "DG_z_um", "z_DG", "DG_z")
  )

  x_E <- pick_col(
    df,
    c("E_CoM_X", "x_E_um", "E_x_um", "x_E", "E_x")
  )
  y_E <- pick_col(
    df,
    c("E_CoM_Y", "y_E_um", "E_y_um", "y_E", "E_y")
  )
  z_E <- pick_col(
    df,
    c("E_CoM_Z", "z_E_um", "E_z_um", "z_E", "E_z")
  )

  x_P <- pick_col(
    df,
    c(
      "P_CoM_X",
      "x_svb_um",
      "x_P_um",
      "x_promoter_um",
      "x_svb",
      "x_promoter"
    )
  )
  y_P <- pick_col(
    df,
    c(
      "P_CoM_Y",
      "y_svb_um",
      "y_P_um",
      "y_promoter_um",
      "y_svb",
      "y_promoter"
    )
  )
  z_P <- pick_col(
    df,
    c(
      "P_CoM_Z",
      "z_svb_um",
      "z_P_um",
      "z_promoter_um",
      "z_svb",
      "z_promoter"
    )
  )

  coordinate_cols <- c(
    x_DG, y_DG, z_DG,
    x_E, y_E, z_E,
    x_P, y_P, z_P
  )

  if (any(is.na(coordinate_cols))) {
    stop(
      paste0(
        "Missing coordinate columns in ",
        label,
        ".\nColumns present:\n",
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  message("\nCoordinate columns used for ", label, ":")

  print(
    tibble::tibble(
      coordinate = c(
        "DG x", "DG y", "DG z",
        "E x", "E y", "E z",
        "P x", "P y", "P z"
      ),
      selected_column = coordinate_cols
    ),
    n = Inf
  )

  df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      centroid_x = (
        .data[[x_DG]] +
          .data[[x_E]] +
          .data[[x_P]]
      ) / 3,

      centroid_y = (
        .data[[y_DG]] +
          .data[[y_E]] +
          .data[[y_P]]
      ) / 3,

      centroid_z = (
        .data[[z_DG]] +
          .data[[z_E]] +
          .data[[z_P]]
      ) / 3,

      d_DG = d3(
        .data[[x_DG]],
        .data[[y_DG]],
        .data[[z_DG]],
        centroid_x,
        centroid_y,
        centroid_z
      ),

      d_E = d3(
        .data[[x_E]],
        .data[[y_E]],
        .data[[z_E]],
        centroid_x,
        centroid_y,
        centroid_z
      ),

      d_P = d3(
        .data[[x_P]],
        .data[[y_P]],
        .data[[z_P]],
        centroid_x,
        centroid_y,
        centroid_z
      ),

      Value = d_DG + d_E + d_P
    ) %>%
    dplyr::ungroup()
}

# ============================================================
# LOAD nc14 WHOLE-EMBRYO DATA
# ============================================================

message("\nReading nc14 CSV:")
message(nc14_csv)

df_nc14_raw <- readr::read_csv(
  nc14_csv,
  show_col_types = FALSE
)

embryo_col_nc14 <- find_embryo_col(
  df_nc14_raw,
  "nc14 CSV"
)

df_nc14_raw <- df_nc14_raw %>%
  dplyr::mutate(
    EmbryoID = as.character(
      .data[[embryo_col_nc14]]
    )
  )

df_nc14 <- compute_dispersion(
  df_nc14_raw,
  "nc14 CSV"
) %>%
  dplyr::filter(
    is.finite(Value),
    Value > 0
  ) %>%
  dplyr::transmute(
    Group = "nc14",
    EmbryoID,
    Value
  )

message(
  "\nnc14: ",
  nrow(df_nc14),
  " triplets from ",
  dplyr::n_distinct(df_nc14$EmbryoID),
  " embryo(s)"
)

# ============================================================
# LOAD WT NON-EXPRESSING DATA
# ============================================================

message("\nReading combined Excel:")
message(xlsx_path)

df_wt_raw <- suppressMessages(
  readxl::read_excel(
    xlsx_path,
    sheet = sheet,
    guess_max = 100000
  )
)

if (!"gene_id" %in% names(df_wt_raw)) {
  stop("Column 'gene_id' was not found.", call. = FALSE)
}

if (!"roi" %in% names(df_wt_raw)) {
  stop("Column 'roi' was not found.", call. = FALSE)
}

message(
  "\ngene_id values present: ",
  paste(
    sort(unique(df_wt_raw$gene_id)),
    collapse = ", "
  )
)

message(
  "\nROI values present for gene_id == '",
  WT_GENE_ID,
  "': ",
  paste(
    sort(
      unique(
        df_wt_raw$roi[
          df_wt_raw$gene_id == WT_GENE_ID
        ]
      )
    ),
    collapse = ", "
  )
)

df_wt_filtered <- df_wt_raw %>%
  dplyr::filter(
    gene_id == WT_GENE_ID,
    tolower(
      trimws(
        as.character(roi)
      )
    ) %in% NONEXP_ALIASES
  )

if (nrow(df_wt_filtered) == 0) {
  stop(
    paste0(
      "No WT non-expressing rows were found for gene_id == '",
      WT_GENE_ID,
      "'."
    ),
    call. = FALSE
  )
}

embryo_col_wt <- find_embryo_col(
  df_wt_filtered,
  "WT Non-expressing"
)

df_wt_filtered <- df_wt_filtered %>%
  dplyr::mutate(
    EmbryoID = as.character(
      .data[[embryo_col_wt]]
    )
  )

message("\nWT non-expressing embryos included:")

print(
  df_wt_filtered %>%
    dplyr::count(
      EmbryoID
    ),
  n = Inf
)

df_nonexp <- compute_dispersion(
  df_wt_filtered,
  "WT Non-expressing"
) %>%
  dplyr::filter(
    is.finite(Value),
    Value > 0
  ) %>%
  dplyr::transmute(
    Group = "Non-expressing",
    EmbryoID,
    Value
  )

message(
  "\nWT Non-expressing: ",
  nrow(df_nonexp),
  " triplets from ",
  dplyr::n_distinct(df_nonexp$EmbryoID),
  " embryo(s)"
)

# ============================================================
# COMBINE DATASETS
# ============================================================

df_all <- dplyr::bind_rows(
  df_nc14,
  df_nonexp
) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = group_levels
    )
  ) %>%
  dplyr::left_join(
    group_spacing,
    by = "Group"
  )

message("\nFinal group summary:")

print(
  df_all %>%
    dplyr::group_by(
      Group
    ) %>%
    dplyr::summarise(
      n_nuclei = dplyr::n(),
      N_embryos = dplyr::n_distinct(
        EmbryoID
      ),
      median_value = stats::median(
        Value,
        na.rm = TRUE
      ),
      mean_value = mean(
        Value,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)

# ============================================================
# NUCLEUS/TRIPLET-LEVEL WILCOXON TEST
# ============================================================

nc14_vals <- df_all$Value[
  df_all$Group == "nc14"
]

nonexp_vals <- df_all$Value[
  df_all$Group == "Non-expressing"
]

p_val <- if (
  length(nc14_vals) >= 3 &&
    length(nonexp_vals) >= 3
) {
  suppressWarnings(
    stats::wilcox.test(
      nc14_vals,
      nonexp_vals,
      exact = FALSE
    )$p.value
  )
} else {
  NA_real_
}

stars <- dplyr::case_when(
  is.na(p_val)      ~ "ns",
  p_val <= 0.0001   ~ "****",
  p_val <= 0.001    ~ "***",
  p_val <= 0.01     ~ "**",
  p_val <= 0.05     ~ "*",
  TRUE              ~ "ns"
)

stats_result <- tibble::tibble(
  Comparison = "nc14_vs_Non-expressing",
  n_nc14 = length(nc14_vals),
  n_nonexpressing = length(nonexp_vals),
  p_raw = p_val,
  stars = stars
)

message("\nNucleus-level Wilcoxon result:")

print(
  stats_result,
  n = Inf
)

# ============================================================
# Y POSITIONS AND STATISTICAL BRACKET
# ============================================================

y_max_data <- max(
  df_all$Value,
  na.rm = TRUE
)

y_min_data <- min(
  df_all$Value[
    df_all$Value > 0
  ],
  na.rm = TRUE
)

y_axis_bot <- y_min_data * 0.5

annot_df <- tibble::tibble(
  bracket_id = "nc14_vs_nonexpressing",
  x1 = 1,
  x2 = 2,
  y = y_max_data * 1.35,
  y_cap = y_max_data * 1.35 / 1.035,
  y_text = y_max_data * 1.35 * 1.08,
  stars = stars
) %>%
  dplyr::filter(
    stars != "ns"
  )

bracket_path_df <- annot_df %>%
  dplyr::select(
    bracket_id,
    x1,
    x2,
    y,
    y_cap
  ) %>%
  dplyr::rowwise() %>%
  dplyr::do(
    tibble::tibble(
      bracket_id = .$bracket_id,
      x = c(
        .$x1,
        .$x1,
        .$x2,
        .$x2
      ),
      y = c(
        .$y_cap,
        .$y,
        .$y,
        .$y_cap
      ),
      point_order = 1:4
    )
  ) %>%
  dplyr::ungroup()

y_axis_top <- if (nrow(annot_df) > 0) {
  max(
    annot_df$y_text,
    na.rm = TRUE
  ) * 1.25
} else {
  y_max_data * 1.5
}

# ============================================================
# nc14 MEDIAN REFERENCE LINE
# ============================================================

nc14_median <- stats::median(
  nc14_vals,
  na.rm = TRUE
)

median_seg <- tibble::tibble(
  x_start = 0.6,
  x_end = 2.4,
  y_med = nc14_median
)

# ============================================================
# BOTTOM LABELS
# ============================================================

bottom_label_df <- df_all %>%
  dplyr::group_by(
    Group,
    xpos
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),

    N = dplyr::n_distinct(
      EmbryoID
    ),

    percent_below =
      100 *
      sum(
        Value < RADIAL_CUTOFF,
        na.rm = TRUE
      ) /
      n_total,

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(
      round(percent_below),
      "%\n",
      "n=", n_total,
      "\n",
      "N=", N
    ),

    y_label = y_axis_bot * 1.35
  )

# ============================================================
# EMBRYO MEDIANS — DISPLAY ONLY
# ============================================================

embryo_median_df <- df_all %>%
  dplyr::group_by(
    Group,
    EmbryoID,
    xpos
  ) %>%
  dplyr::summarise(
    embryo_median = stats::median(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ============================================================
# FINAL PLOT
# ============================================================

p <- ggplot2::ggplot(
  df_all,
  ggplot2::aes(
    x = xpos,
    y = Value
  )
) +

  gghalves::geom_half_violin(
    ggplot2::aes(
      group = Group
    ),
    side = "r",
    width = violin_w,
    fill = NA,
    alpha = 1,
    linewidth = 0.5,
    colour = "black",
    trim = TRUE
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      x = xpos - dot_left_shift,
      y = Value,
      colour = Group
    ),
    inherit.aes = FALSE,
    shape = 16,
    size = 0.9,
    alpha = 0.35,
    position = ggplot2::position_jitter(
      width = dot_jit_x,
      height = 0,
      seed = 1
    )
  ) +

  ggplot2::geom_boxplot(
    ggplot2::aes(
      fill = Group,
      group = Group
    ),
    width = box_w,
    outlier.shape = NA,
    linewidth = 0.5,
    colour = "black",
    alpha = 0.30
  ) +

  ggplot2::geom_point(
    data = embryo_median_df,
    ggplot2::aes(
      x = xpos,
      y = embryo_median
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.5,
    alpha = 0.85,
    stroke = 0.5,
    fill = "#686765",
    colour = "white",
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0,
      seed = 42
    )
  ) +

  ggplot2::geom_segment(
    data = median_seg,
    ggplot2::aes(
      x = x_start,
      xend = x_end,
      y = y_med,
      yend = y_med
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.55,
    colour = "black"
  ) +

  {
    if (nrow(annot_df) > 0) {
      list(
        ggplot2::geom_path(
          data = bracket_path_df,
          ggplot2::aes(
            x = x,
            y = y,
            group = bracket_id
          ),
          inherit.aes = FALSE,
          linewidth = 0.5,
          colour = "black",
          lineend = "square",
          linejoin = "mitre"
        ),

        ggplot2::geom_text(
          data = annot_df,
          ggplot2::aes(
            x = (x1 + x2) / 2,
            y = y_text,
            label = stars
          ),
          inherit.aes = FALSE,
          size = FONT_STARS,
          fontface = "bold",
          colour = "black",
          family = FONT_FAMILY
        )
      )
    }
  } +

  ggplot2::geom_text(
    data = bottom_label_df,
    ggplot2::aes(
      x = xpos,
      y = y_label,
      label = label
    ),
    inherit.aes = FALSE,
    size = FONT_N_LABEL,
    lineheight = 0.9,
    colour = "black",
    family = FONT_FAMILY
  ) +

  ggplot2::scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      "nc14",
      "Non-expressing"
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.06,
        0.08
      )
    )
  ) +

  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = SIGMA
    ),
    breaks = y_breaks,
    labels = scales::number_format(
      accuracy = 0.01
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.02,
        0
      )
    )
  ) +

  ggplot2::coord_cartesian(
    ylim = c(
      y_axis_bot,
      y_axis_top
    ),
    clip = "off"
  ) +

  ggplot2::scale_fill_manual(
    values = group_colors,
    drop = FALSE,
    name = NULL
  ) +

  ggplot2::scale_colour_manual(
    values = group_colors,
    drop = FALSE,
    name = NULL
  ) +

  ggplot2::guides(
    fill = "none",
    colour = "none"
  ) +

  ggplot2::labs(
    x = NULL,
    y = expression(
      "Hub dispersion (" *
        mu *
        "m)"
    )
  ) +

  ggplot2::theme_classic(
    base_size = 20,
    base_family = FONT_FAMILY
  ) +

  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = FONT_AXIS_TEXT,
      colour = "black"
    ),

    axis.text.y = ggplot2::element_text(
      size = FONT_AXIS_TEXT,
      colour = "black"
    ),

    axis.title.y = ggplot2::element_text(
      size = FONT_AXIS_TITLE,
      margin = ggplot2::margin(
        r = 8
      )
    ),

    axis.line = ggplot2::element_line(
      linewidth = 0.5,
      colour = "black"
    ),

    axis.ticks = ggplot2::element_line(
      linewidth = 0.4,
      colour = "black"
    ),

    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    legend.position = "none",

    plot.margin = ggplot2::margin(
      6,
      18,
      45,
      6
    )
  )

print(p)
ggplot2::ggsave(
  filename = "F:\\Thesis\\PAPER PDFS\\nc.pdf",
  plot = p,
  width = 2.5,
  height = 3.1,
  units = "in",
  device = cairo_pdf
)

