# ============================================================
# INTER-PROBE DISTANCES — NEW COMBINED DATASET
# FILTER: retain only pairwise distances >= 0.01 µm
# PLOT ORDER: Promoter-E6, Promoter-DG2, E6-DG2
# WT (SXLGFP): NON-EXPRESSING / DORSAL / VENTRAL
#
# Data source:
#   ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx
#   sheet = "Filtered_Triplets"
#
# Statistics:
# - Wilcoxon rank-sum tests at nucleus/triplet level
# - ONLY:
#     OFF vs Dorsal
#     OFF vs Ventral
# - BH correction across these 2 comparisons
#   separately within each probe-pair group
#
# Plot:
# - transparent right half-violins
# - jittered individual nuclei/triplets
# - pooled boxplots
# - embryo median dots
# - n / N labels
# - ns brackets included
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(gghalves)
  library(scales)
  library(tibble)
  library(stringr)
})

# ============================================================
# INPUT
# ============================================================

xlsx_path <- paste0(
  "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/",
  "AP_file_2026/ROI/",
  "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

sheet <- "Filtered_Triplets"
WT_GENE_ID <- "SXLGFP"

# Remove embryo-condition groups with fewer than this many triplets
LOW_N_THRESHOLD <- 8

CONTROL_COND <- "shavenbaby OFF"

ROI_ALIASES <- list(
  "shavenbaby OFF" = c(
    "nonexp",
    "non_exp",
    "non-exp",
    "nonexpressing",
    "non-expressing",
    "ne"
  ),
  "shavenbaby ON (dorsal)" = c(
    "dorsal",
    "d"
  ),
  "shavenbaby ON (ventral)" = c(
    "ventral",
    "v"
  )
)

# ============================================================
# PSEUDO-LOG SETTINGS
# ============================================================

SIGMA <- 0.001

y_breaks <- c(
  0.01,
  0.05,
  0.10,
  0.20,
  0.50,
  1.00,
  2.00
)

# ============================================================
# ORDER AND COLOURS
# ============================================================

cond_levels <- c(
  "shavenbaby OFF",
  "shavenbaby ON (dorsal)",
  "shavenbaby ON (ventral)"
)

cond_cols <- c(
  "shavenbaby OFF"          = "#686765",
  "shavenbaby ON (dorsal)"  = "#cc2829",
  "shavenbaby ON (ventral)" = "#2b3f99"
)

pair_levels <- c(
  "Promoter-E6",
  "Promoter-DG2",
  "E6-DG2"
)

group_spacing <- tibble::tibble(
  PairGroup = factor(
    pair_levels,
    levels = pair_levels
  ),
  x_center = c(
    1.00,
    2.40,
    3.80
  )
)

dodge_offset <- 0.42

cond_offsets <- c(
  "shavenbaby OFF"          = -dodge_offset,
  "shavenbaby ON (dorsal)"  = 0,
  "shavenbaby ON (ventral)" = dodge_offset
)

# ============================================================
# STYLE
# ============================================================

box_w <- 0.18
violin_w <- 0.37

dot_band_w <- box_w / 2
dot_left_shift <- dot_band_w
dot_jit_x <- dot_band_w * 0.9

FONT_FAMILY <- "Arial"

FONT_AXIS_TITLE <- 8
FONT_AXIS_TEXT <- 7
FONT_N_LABEL <- 6
FONT_STARS <- 7

# ============================================================
# HELPERS
# ============================================================

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    NA_character_
  } else {
    hit[[1]]
  }
}

find_embryo_col <- function(df) {
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
        "Could not find embryo/image ID column.\n",
        "Columns present:\n",
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  col
}

d3 <- function(x1, y1, z1, x2, y2, z2) {
  sqrt(
    (x1 - x2)^2 +
      (y1 - y2)^2 +
      (z1 - z2)^2
  )
}

pairwise_wilcox <- function(d, c1, c2) {
  x <- d$Value[d$Condition == c1]
  y <- d$Value[d$Condition == c2]

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) >= 3 && length(y) >= 3) {
    suppressWarnings(
      stats::wilcox.test(
        x,
        y,
        exact = FALSE
      )$p.value
    )
  } else {
    NA_real_
  }
}

get_xpos <- function(pair, cond) {
  xc <- group_spacing$x_center[
    as.character(group_spacing$PairGroup) == pair
  ][1]

  xc + cond_offsets[[cond]]
}

# ============================================================
# LOAD COMBINED DATA
# ============================================================

message("Reading combined Excel file:")
message(xlsx_path)

df_raw_all <- suppressMessages(
  readxl::read_excel(
    xlsx_path,
    sheet = sheet,
    guess_max = 100000
  )
)

if (!"gene_id" %in% names(df_raw_all)) {
  stop("Column 'gene_id' was not found.", call. = FALSE)
}

if (!"roi" %in% names(df_raw_all)) {
  stop("Column 'roi' was not found.", call. = FALSE)
}

message(
  "gene_id values present: ",
  paste(
    sort(unique(df_raw_all$gene_id)),
    collapse = ", "
  )
)

df_wt <- df_raw_all %>%
  dplyr::filter(
    gene_id == WT_GENE_ID
  )

if (nrow(df_wt) == 0) {
  stop(
    paste0(
      "No rows found for gene_id == '",
      WT_GENE_ID,
      "'."
    ),
    call. = FALSE
  )
}

message(
  "ROI values present for ",
  WT_GENE_ID,
  ": ",
  paste(
    sort(unique(df_wt$roi)),
    collapse = ", "
  )
)

# ============================================================
# IDENTIFY EMBRYO COLUMN
# ============================================================

embryo_col <- find_embryo_col(df_wt)

df_wt <- df_wt %>%
  dplyr::mutate(
    EmbryoID = as.character(
      .data[[embryo_col]]
    ),
    Experiment = stringr::str_extract(
      EmbryoID,
      "Experiment-[0-9]+"
    )
  )

# ============================================================
# ASSIGN CONDITION FROM ROI
# ============================================================

roi_lower <- tolower(
  trimws(
    as.character(df_wt$roi)
  )
)

condition_vector <- rep(
  NA_character_,
  nrow(df_wt)
)

for (condition_name in names(ROI_ALIASES)) {
  condition_vector[
    roi_lower %in% ROI_ALIASES[[condition_name]]
  ] <- condition_name
}

df_wt$Condition <- condition_vector

unmatched_roi <- unique(
  df_wt$roi[
    is.na(df_wt$Condition)
  ]
)

if (length(unmatched_roi) > 0) {
  message(
    "Unmatched ROI values excluded: ",
    paste(
      unmatched_roi,
      collapse = ", "
    )
  )
}

# ============================================================
# FIND DISTANCE OR COORDINATE COLUMNS
# ============================================================

dist_DG_E <- pick_col(
  df_wt,
  c(
    "Dist_CoM_DG_E",
    "dist_DG_E_um",
    "DG_E_distance_um",
    "distance_DG_E"
  )
)

dist_DG_P <- pick_col(
  df_wt,
  c(
    "Dist_CoM_DG_P",
    "dist_DG_P_um",
    "DG_P_distance_um",
    "distance_DG_P"
  )
)

dist_E_P <- pick_col(
  df_wt,
  c(
    "Dist_CoM_E_P",
    "dist_E_P_um",
    "E_P_distance_um",
    "distance_E_P"
  )
)

use_precomputed_distances <- !any(
  is.na(
    c(
      dist_DG_E,
      dist_DG_P,
      dist_E_P
    )
  )
)

if (use_precomputed_distances) {
  message("\nUsing precomputed distance columns:")

  print(
    tibble::tibble(
      PairGroup = c(
        "E6-DG2",
        "Promoter-DG2",
        "Promoter-E6"
      ),
      selected_column = c(
        dist_DG_E,
        dist_DG_P,
        dist_E_P
      )
    ),
    n = Inf
  )
} else {
  message(
    "\nComplete precomputed distance columns were not found. ",
    "Distances will be calculated from coordinates."
  )

  x_DG <- pick_col(
    df_wt,
    c(
      "x_DG_um",
      "DG_CoM_X",
      "DG_x_um",
      "x_DG",
      "DG_x"
    )
  )

  y_DG <- pick_col(
    df_wt,
    c(
      "y_DG_um",
      "DG_CoM_Y",
      "DG_y_um",
      "y_DG",
      "DG_y"
    )
  )

  z_DG <- pick_col(
    df_wt,
    c(
      "z_DG_um",
      "DG_CoM_Z",
      "DG_z_um",
      "z_DG",
      "DG_z"
    )
  )

  x_E <- pick_col(
    df_wt,
    c(
      "x_E_um",
      "E_CoM_X",
      "E_x_um",
      "x_E",
      "E_x"
    )
  )

  y_E <- pick_col(
    df_wt,
    c(
      "y_E_um",
      "E_CoM_Y",
      "E_y_um",
      "y_E",
      "E_y"
    )
  )

  z_E <- pick_col(
    df_wt,
    c(
      "z_E_um",
      "E_CoM_Z",
      "E_z_um",
      "z_E",
      "E_z"
    )
  )

  x_P <- pick_col(
    df_wt,
    c(
      "x_svb_um",
      "P_CoM_X",
      "x_P_um",
      "x_promoter_um",
      "x_svb",
      "x_promoter"
    )
  )

  y_P <- pick_col(
    df_wt,
    c(
      "y_svb_um",
      "P_CoM_Y",
      "y_P_um",
      "y_promoter_um",
      "y_svb",
      "y_promoter"
    )
  )

  z_P <- pick_col(
    df_wt,
    c(
      "z_svb_um",
      "P_CoM_Z",
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
        "Neither complete precomputed distance columns nor complete ",
        "coordinate columns were found.\nColumns present:\n",
        paste(names(df_wt), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  message("\nCoordinate columns selected:")

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
}

# ============================================================
# CREATE LONG PAIRWISE-DISTANCE TABLE
# ============================================================

if (use_precomputed_distances) {
  df_all <- df_wt %>%
    dplyr::filter(
      !is.na(Condition)
    ) %>%
    dplyr::transmute(
      Condition,
      EmbryoID,
      Experiment,

      DG2_E6 = as.numeric(
        .data[[dist_DG_E]]
      ),

      DG2_Promoter = as.numeric(
        .data[[dist_DG_P]]
      ),

      E6_Promoter = as.numeric(
        .data[[dist_E_P]]
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        DG2_E6,
        DG2_Promoter,
        E6_Promoter
      ),
      names_to = "PairRaw",
      values_to = "Value"
    )
} else {
  df_all <- df_wt %>%
    dplyr::filter(
      !is.na(Condition)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      DG2_E6 = d3(
        .data[[x_DG]],
        .data[[y_DG]],
        .data[[z_DG]],
        .data[[x_E]],
        .data[[y_E]],
        .data[[z_E]]
      ),

      DG2_Promoter = d3(
        .data[[x_DG]],
        .data[[y_DG]],
        .data[[z_DG]],
        .data[[x_P]],
        .data[[y_P]],
        .data[[z_P]]
      ),

      E6_Promoter = d3(
        .data[[x_E]],
        .data[[y_E]],
        .data[[z_E]],
        .data[[x_P]],
        .data[[y_P]],
        .data[[z_P]]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      Condition,
      EmbryoID,
      Experiment,
      DG2_E6,
      DG2_Promoter,
      E6_Promoter
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        DG2_E6,
        DG2_Promoter,
        E6_Promoter
      ),
      names_to = "PairRaw",
      values_to = "Value"
    )
}

df_all <- df_all %>%
  dplyr::mutate(
    PairGroup = dplyr::case_when(
      PairRaw == "E6_Promoter" ~ "Promoter-E6",
      PairRaw == "DG2_Promoter" ~ "Promoter-DG2",
      PairRaw == "DG2_E6" ~ "E6-DG2",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(
    !is.na(PairGroup),
    is.finite(Value),
    Value >= 0.01
  ) %>%
  dplyr::mutate(
    Condition = factor(
      Condition,
      levels = cond_levels
    ),

    PairGroup = factor(
      PairGroup,
      levels = pair_levels
    )
  ) %>%
  dplyr::left_join(
    group_spacing,
    by = "PairGroup"
  ) %>%
  dplyr::mutate(
    xpos = x_center +
      unname(
        cond_offsets[
          as.character(Condition)
        ]
      ),

    group_id = interaction(
      PairGroup,
      Condition,
      drop = TRUE
    )
  )

message(
  "\nAll pairwise distances below 0.01 µm were removed before QC, ",
  "statistics, medians and plotting."
)

message("\nDistance-value counts after the 0.01 µm distance filter and before low-n QC:")

print(
  df_all %>%
    dplyr::count(
      PairGroup,
      Condition
    ),
  n = Inf
)

# ============================================================
# LOW-N QC FILTER
# Applied by embryo-condition group
# ============================================================

# Count the number of retained triplets from one probe-pair
# so that each triplet is counted once rather than three times.
qc_source <- df_all %>%
  dplyr::filter(
    PairGroup == pair_levels[[1]]
  ) %>%
  dplyr::group_by(
    Condition,
    EmbryoID,
    Experiment
  ) %>%
  dplyr::summarise(
    n_triplets = dplyr::n(),
    .groups = "drop"
  )

low_n_pairs <- qc_source %>%
  dplyr::filter(
    n_triplets < LOW_N_THRESHOLD
  ) %>%
  dplyr::select(
    Condition,
    EmbryoID
  )

message(
  "\nExcluded embryo-condition groups with fewer than ",
  LOW_N_THRESHOLD,
  " triplets:"
)

print(
  low_n_pairs,
  n = Inf
)

df_all <- df_all %>%
  dplyr::anti_join(
    low_n_pairs,
    by = c(
      "Condition",
      "EmbryoID"
    )
  )

message("\nDistance-value counts after low-n QC:")

print(
  df_all %>%
    dplyr::count(
      PairGroup,
      Condition
    ),
  n = Inf
)

message("\nEmbryos retained per condition:")

print(
  df_all %>%
    dplyr::group_by(
      Condition
    ) %>%
    dplyr::summarise(
      N_embryos = dplyr::n_distinct(
        EmbryoID
      ),
      .groups = "drop"
    )
)

# ============================================================
# STATISTICS
# ONLY OFF vs Dorsal and OFF vs Ventral
# BH correction across these 2 tests within each PairGroup
# ns brackets retained
# ============================================================

comparisons_tbl <- tibble::tibble(
  Comparison = c(
    "OFF_vs_Dorsal",
    "OFF_vs_Ventral"
  ),

  c1 = c(
    "shavenbaby OFF",
    "shavenbaby OFF"
  ),

  c2 = c(
    "shavenbaby ON (dorsal)",
    "shavenbaby ON (ventral)"
  ),

  bracket_level = c(
    1,
    2
  )
)

stats_all <- purrr::map_dfr(
  pair_levels,
  function(pg) {
    d <- df_all %>%
      dplyr::filter(
        PairGroup == pg
      )

    comparisons_tbl %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n_c1 = sum(
          d$Condition == c1,
          na.rm = TRUE
        ),

        n_c2 = sum(
          d$Condition == c2,
          na.rm = TRUE
        ),

        p_raw = pairwise_wilcox(
          d,
          c1,
          c2
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        # BH correction across the two requested comparisons
        # within this probe-pair group
        p_adj = stats::p.adjust(
          p_raw,
          method = "BH"
        ),

        stars = dplyr::case_when(
          is.na(p_adj)      ~ "ns",
          p_adj <= 0.0001   ~ "****",
          p_adj <= 0.001    ~ "***",
          p_adj <= 0.01     ~ "**",
          p_adj <= 0.05     ~ "*",
          TRUE              ~ "ns"
        ),

        PairGroup = pg,

        x1 = purrr::map2_dbl(
          pg,
          c1,
          get_xpos
        ),

        x2 = purrr::map2_dbl(
          pg,
          c2,
          get_xpos
        )
      )
  }
)

message("\nNucleus/triplet-level Wilcoxon results:")

print(
  stats_all %>%
    dplyr::select(
      PairGroup,
      Comparison,
      n_c1,
      n_c2,
      p_raw,
      p_adj,
      stars
    ),
  n = Inf
)

# ============================================================
# Y POSITIONS AND STATISTICAL BRACKETS
# ============================================================

tops <- df_all %>%
  dplyr::group_by(
    PairGroup
  ) %>%
  dplyr::summarise(
    top = max(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

annot_df <- stats_all %>%
  dplyr::left_join(
    tops,
    by = "PairGroup"
  ) %>%
  dplyr::mutate(
    y = top *
      dplyr::case_when(
        bracket_level == 1 ~ 1.35,
        bracket_level == 2 ~ 1.75,
        TRUE ~ 1.35
      ),

    y_cap = y / 1.035,
    y_text = y * 1.08,

    bracket_id = paste(
      PairGroup,
      Comparison,
      sep = "_"
    )
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

y_axis_bot <- min(
  df_all$Value,
  na.rm = TRUE
) * 0.5

y_axis_top <- max(
  annot_df$y_text,
  na.rm = TRUE
) * 1.35

# ============================================================
# OFF MEDIAN DOTTED LINE PER PROBE PAIR
# ============================================================

off_medians <- df_all %>%
  dplyr::filter(
    Condition == CONTROL_COND
  ) %>%
  dplyr::group_by(
    PairGroup
  ) %>%
  dplyr::summarise(
    y_med = stats::median(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    x = purrr::map_dbl(
      as.character(PairGroup),
      ~ get_xpos(.x, CONTROL_COND)
    ),

    x_start = x - box_w / 2,
    x_end = x + box_w / 2
  )

# ============================================================
# n / N LABELS
# ============================================================

bottom_label_df <- df_all %>%
  dplyr::group_by(
    PairGroup,
    Condition,
    xpos
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),

    N = dplyr::n_distinct(
      EmbryoID
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(
      "n=",
      n_total,
      "\nN=",
      N
    ),

    y_label = y_axis_bot * 1.35
  )

# ============================================================
# EMBRYO MEDIAN DOTS
# ============================================================

embryo_median_df <- df_all %>%
  dplyr::group_by(
    PairGroup,
    Condition,
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
# SEPARATORS BETWEEN PAIR GROUPS
# ============================================================

seps <- c(
  mean(
    c(
      group_spacing$x_center[1],
      group_spacing$x_center[2]
    )
  ),

  mean(
    c(
      group_spacing$x_center[2],
      group_spacing$x_center[3]
    )
  )
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

  ggplot2::geom_vline(
    xintercept = seps,
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.3
  ) +

  gghalves::geom_half_violin(
    ggplot2::aes(
      group = group_id
    ),
    side = "r",
    width = violin_w,
    fill = NA,
    alpha = 1,
    linewidth = 0.3,
    colour = "black",
    trim = TRUE
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      x = xpos - dot_left_shift,
      y = Value,
      colour = Condition
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
      fill = Condition,
      group = group_id
    ),
    width = box_w,
    outlier.shape = NA,
    linewidth = 0.3,
    colour = "black",
    alpha = 0.30
  ) +

  ggplot2::geom_point(
    data = embryo_median_df,
    ggplot2::aes(
      x = xpos,
      y = embryo_median,
      fill = Condition
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 0.8,
    alpha = 0.85,
    stroke = 0.3,
    colour = "white",
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0,
      seed = 42
    )
  ) +

  ggplot2::geom_segment(
    data = off_medians,
    ggplot2::aes(
      x = x_start,
      xend = x_end,
      y = y_med,
      yend = y_med
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.4,
    colour = "black"
  ) +

  ggplot2::geom_text(
    data = bottom_label_df,
    ggplot2::aes(
      x = xpos,
      y = y_label,
      label = label
    ),
    inherit.aes = FALSE,
    size = FONT_N_LABEL / ggplot2::.pt,
    lineheight = 0.9,
    colour = "black",
    family = FONT_FAMILY
  ) +

  ggplot2::geom_path(
    data = bracket_path_df,
    ggplot2::aes(
      x = x,
      y = y,
      group = bracket_id
    ),
    inherit.aes = FALSE,
    linewidth = 0.3,
    colour = "black",
    lineend = "square",
    linejoin = "mitre"
  ) +

  ggplot2::geom_text(
    data = annot_df,
    ggplot2::aes(
      x = (x1 + x2) / 2,
      y = y_text,
      label = stars
    ),
    inherit.aes = FALSE,
    size = FONT_STARS / ggplot2::.pt,
    fontface = "bold",
    colour = "black",
    family = FONT_FAMILY
  ) +

  ggplot2::scale_x_continuous(
    breaks = group_spacing$x_center,
    labels = as.character(
      group_spacing$PairGroup
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.02,
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
    values = cond_cols,
    drop = FALSE,
    name = NULL
  ) +

  ggplot2::scale_colour_manual(
    values = cond_cols,
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
      "Inter-probe distance (" *
        mu *
        "m)"
    )
  ) +

  ggplot2::theme_classic(
    base_size = FONT_AXIS_TITLE,
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
      linewidth = 0.3,
      colour = "black"
    ),

    axis.ticks = ggplot2::element_line(
      linewidth = 0.3,
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

# ============================================================
# SAVE OUTPUTS
# ============================================================

ggplot2::ggsave(
  filename = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_new_dataset_interprobe_ge_0.01_order_PromoterE6_PromoterDG2_E6DG2.pdf"
  ),
  plot = p,
  width = 6.8,
  height = 4.4,
  units = "in",
  device = cairo_pdf
)

write.csv(
  stats_all,
  file = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_new_dataset_interprobe_ge_0.01_order_PromoterE6_PromoterDG2_E6DG2_stats.csv"
  ),
  row.names = FALSE
)

write.csv(
  qc_source,
  file = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_new_dataset_interprobe_ge_0.01_order_PromoterE6_PromoterDG2_E6DG2_QC.csv"
  ),
  row.names = FALSE
)

