# ============================================================
# WT HUB DISPERSION
# NON-EXPRESSING vs DORSAL AND VENTRAL
#
# NO MANUAL EMBRYO EXCLUSIONS
#
# Statistics:
# - nucleus/triplet-level Wilcoxon tests
# - only OFF vs Dorsal and OFF vs Ventral
# - BH correction across these 2 comparisons only
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(gghalves)
  library(tibble)
  library(scales)
})

xlsx_path <- paste0(
  "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/",
  "AP_file_2026/ROI/",
  "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

sheet <- "Filtered_Triplets"
WT_GENE_ID <- "SXLGFP"
LOW_N_THRESHOLD <- 8
CONTROL_COND <- "shavenbaby OFF"
RADIAL_CUTOFF <- 0.5

ROI_ALIASES <- list(
  "shavenbaby OFF" = c("nonexp", "non_exp", "non-exp", "nonexpressing", "non-expressing", "ne"),
  "shavenbaby ON (dorsal)" = c("dorsal", "d"),
  "shavenbaby ON (ventral)" = c("ventral", "v")
)

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

condition_spacing <- tibble::tibble(
  Condition = factor(cond_levels, levels = cond_levels),
  xpos = c(1, 2, 3)
)

SIGMA <- 0.001
y_breaks <- c(0.1, 0.2, 0.5, 1, 2)
box_w <- 0.35
violin_w <- 1.0
dot_band_w <- box_w / 2
dot_left_shift <- dot_band_w
dot_jit_x <- dot_band_w * 0.9
FONT_FAMILY <- "Arial"
FONT_AXIS_TITLE <- 8
FONT_AXIS_TEXT <- 7
FONT_N_LABEL <- 6
FONT_STARS <- 7

d3 <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

find_embryo_col <- function(df, label) {
  col <- pick_col(df, c(
    "FileName_Image", "filename_image", "FileName", "filename",
    "EmbryoID", "embryo_id", "embryo", "Embryo", "Image", "image",
    "triplet_source"
  ))
  if (is.na(col)) {
    stop(paste0(
      "Could not find embryo/image ID column in ", label,
      ".\nColumns present:\n", paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  col
}

pairwise_wilcox_nuclei <- function(d, c1, c2) {
  x <- d$Value[d$Condition == c1]
  y <- d$Value[d$Condition == c2]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
}

message("Reading Excel file: ", xlsx_path)
df_raw_all <- suppressMessages(readxl::read_excel(
  xlsx_path,
  sheet = sheet,
  guess_max = 100000
))

if (!"gene_id" %in% names(df_raw_all)) stop("Column 'gene_id' was not found.", call. = FALSE)
if (!"roi" %in% names(df_raw_all)) stop("Column 'roi' was not found.", call. = FALSE)

df_wt <- df_raw_all %>% dplyr::filter(gene_id == WT_GENE_ID)
if (nrow(df_wt) == 0) stop(paste0("No rows found for gene_id == '", WT_GENE_ID, "'."), call. = FALSE)

embryo_col <- find_embryo_col(df_wt, paste0("gene_id = ", WT_GENE_ID))
df_wt <- df_wt %>%
  dplyr::mutate(
    EmbryoID = as.character(.data[[embryo_col]]),
    Experiment = stringr::str_extract(EmbryoID, "Experiment-[0-9]+")
  )

roi_lower <- tolower(trimws(as.character(df_wt$roi)))
condition_vector <- rep(NA_character_, nrow(df_wt))
for (condition_name in names(ROI_ALIASES)) {
  condition_vector[roi_lower %in% ROI_ALIASES[[condition_name]]] <- condition_name
}
df_wt$Condition <- condition_vector

unmatched_roi <- unique(df_wt$roi[is.na(df_wt$Condition)])
if (length(unmatched_roi) > 0) {
  message("Unmatched ROI values excluded: ", paste(unmatched_roi, collapse = ", "))
}

x_DG <- pick_col(df_wt, c("x_DG_um", "DG_CoM_X", "DG_x_um", "x_DG", "DG_x"))
y_DG <- pick_col(df_wt, c("y_DG_um", "DG_CoM_Y", "DG_y_um", "y_DG", "DG_y"))
z_DG <- pick_col(df_wt, c("z_DG_um", "DG_CoM_Z", "DG_z_um", "z_DG", "DG_z"))
x_E  <- pick_col(df_wt, c("x_E_um", "E_CoM_X", "E_x_um", "x_E", "E_x"))
y_E  <- pick_col(df_wt, c("y_E_um", "E_CoM_Y", "E_y_um", "y_E", "E_y"))
z_E  <- pick_col(df_wt, c("z_E_um", "E_CoM_Z", "E_z_um", "z_E", "E_z"))
x_P  <- pick_col(df_wt, c("x_svb_um", "P_CoM_X", "x_P_um", "x_promoter_um", "x_svb", "x_promoter"))
y_P  <- pick_col(df_wt, c("y_svb_um", "P_CoM_Y", "y_P_um", "y_promoter_um", "y_svb", "y_promoter"))
z_P  <- pick_col(df_wt, c("z_svb_um", "P_CoM_Z", "z_P_um", "z_promoter_um", "z_svb", "z_promoter"))

coordinate_cols <- c(x_DG, y_DG, z_DG, x_E, y_E, z_E, x_P, y_P, z_P)
if (any(is.na(coordinate_cols))) {
  stop(paste0(
    "One or more coordinate columns could not be found.\nColumns present:\n",
    paste(names(df_wt), collapse = ", ")
  ), call. = FALSE)
}

message("Coordinate columns selected:")
print(tibble::tibble(
  coordinate = c("DG x", "DG y", "DG z", "E x", "E y", "E z", "P x", "P y", "P z"),
  selected_column = coordinate_cols
), n = Inf)

df_all <- df_wt %>%
  dplyr::filter(!is.na(Condition)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    centroid_x = (.data[[x_DG]] + .data[[x_E]] + .data[[x_P]]) / 3,
    centroid_y = (.data[[y_DG]] + .data[[y_E]] + .data[[y_P]]) / 3,
    centroid_z = (.data[[z_DG]] + .data[[z_E]] + .data[[z_P]]) / 3,
    d_DG = d3(.data[[x_DG]], .data[[y_DG]], .data[[z_DG]], centroid_x, centroid_y, centroid_z),
    d_E  = d3(.data[[x_E]], .data[[y_E]], .data[[z_E]], centroid_x, centroid_y, centroid_z),
    d_P  = d3(.data[[x_P]], .data[[y_P]], .data[[z_P]], centroid_x, centroid_y, centroid_z),
    Value = d_DG + d_E + d_P
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(Condition, EmbryoID, Experiment, Value) %>%
  dplyr::filter(is.finite(Value), Value > 0) %>%
  dplyr::mutate(Condition = factor(Condition, levels = cond_levels)) %>%
  dplyr::left_join(condition_spacing, by = "Condition") %>%
  dplyr::mutate(group_id = Condition)

message("Nuclei before low-n QC:")
print(df_all %>% dplyr::count(Condition))

qc_table <- df_all %>%
  dplyr::group_by(Condition, EmbryoID) %>%
  dplyr::summarise(
    n_triplets = dplyr::n(),
    embryo_median = stats::median(Value, na.rm = TRUE),
    .groups = "drop"
  )

low_n_pairs <- qc_table %>%
  dplyr::filter(n_triplets < LOW_N_THRESHOLD) %>%
  dplyr::select(Condition, EmbryoID)

message("Excluded embryo-condition groups with fewer than ", LOW_N_THRESHOLD, " nuclei:")
print(low_n_pairs, n = Inf)

df_all <- df_all %>%
  dplyr::anti_join(low_n_pairs, by = c("Condition", "EmbryoID"))

message("Final nuclei retained per condition:")
print(df_all %>% dplyr::count(Condition))

message("Final embryos retained per condition:")
print(df_all %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(N_embryos = dplyr::n_distinct(EmbryoID), .groups = "drop"))

embryo_median_df <- df_all %>%
  dplyr::group_by(Condition, EmbryoID, xpos) %>%
  dplyr::summarise(
    embryo_median = stats::median(Value, na.rm = TRUE),
    .groups = "drop"
  )

comparisons_tbl <- tibble::tibble(
  Comparison = c("OFF_vs_Dorsal", "OFF_vs_Ventral"),
  c1 = c("shavenbaby OFF", "shavenbaby OFF"),
  c2 = c("shavenbaby ON (dorsal)", "shavenbaby ON (ventral)"),
  x1 = c(1, 1),
  x2 = c(2, 3),
  bracket_level = c(1, 2)
)

stats_all <- comparisons_tbl %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_c1 = sum(df_all$Condition == c1, na.rm = TRUE),
    n_c2 = sum(df_all$Condition == c2, na.rm = TRUE),
    p_raw = pairwise_wilcox_nuclei(df_all, c1, c2)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_raw, method = "BH"),
    stars = dplyr::case_when(
      is.na(p_adj)    ~ "ns",
      p_adj <= 0.0001 ~ "****",
      p_adj <= 0.001  ~ "***",
      p_adj <= 0.01   ~ "**",
      p_adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  )

message("Nucleus-level Wilcoxon results:")
print(stats_all %>% dplyr::select(Comparison, n_c1, n_c2, p_raw, p_adj, stars), n = Inf)

stats_plot <- stats_all %>% dplyr::filter(stars != "ns")
y_max_data <- max(df_all$Value, na.rm = TRUE)
y_min_data <- min(df_all$Value, na.rm = TRUE)
y_axis_bot <- y_min_data * 0.5

annot_df <- stats_plot %>%
  dplyr::mutate(
    y = y_max_data * dplyr::case_when(
      bracket_level == 1 ~ 1.35,
      bracket_level == 2 ~ 1.80,
      TRUE ~ 1.35
    ),
    y_cap = y / 1.035,
    y_text = y * 1.08,
    bracket_id = Comparison
  )

bracket_path_df <- annot_df %>%
  dplyr::select(bracket_id, x1, x2, y, y_cap) %>%
  dplyr::rowwise() %>%
  dplyr::do(tibble::tibble(
    bracket_id = .$bracket_id,
    x = c(.$x1, .$x1, .$x2, .$x2),
    y = c(.$y_cap, .$y, .$y, .$y_cap),
    point_order = 1:4
  )) %>%
  dplyr::ungroup()

y_axis_top <- if (nrow(annot_df) > 0) max(annot_df$y_text, na.rm = TRUE) * 1.30 else y_max_data * 1.50

off_median <- df_all %>%
  dplyr::filter(Condition == CONTROL_COND) %>%
  dplyr::summarise(y_med = stats::median(Value, na.rm = TRUE)) %>%
  dplyr::pull(y_med)

off_seg <- tibble::tibble(x_start = 0.6, x_end = 3.4, y_med = off_median)

bottom_label_df <- df_all %>%
  dplyr::group_by(Condition, xpos) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    N = dplyr::n_distinct(EmbryoID),
    percent_below = 100 * sum(Value < RADIAL_CUTOFF, na.rm = TRUE) / n_total,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(round(percent_below), "%\n", "n=", n_total, "\n", "N=", N),
    y_label = y_axis_bot * 1.35
  )

p <- ggplot2::ggplot(df_all, ggplot2::aes(x = xpos, y = Value)) +
  gghalves::geom_half_violin(
    ggplot2::aes(group = group_id),
    side = "r", width = violin_w, fill = NA, alpha = 1,
    linewidth = 0.3, colour = "black", trim = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = xpos - dot_left_shift, y = Value, colour = Condition),
    inherit.aes = FALSE, shape = 16, size = 0.9, alpha = 0.35,
    position = ggplot2::position_jitter(width = dot_jit_x, height = 0, seed = 1)
  ) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = Condition, group = group_id),
    width = box_w, outlier.shape = NA, linewidth = 0.3,
    colour = "black", alpha = 0.30
  ) +
  ggplot2::geom_point(
    data = embryo_median_df,
    ggplot2::aes(x = xpos, y = embryo_median, fill = Condition),
    inherit.aes = FALSE, shape = 21, size = 0.8, alpha = 0.85,
    stroke = 0.3, colour = "white",
    position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 42)
  ) +
  ggplot2::geom_segment(
    data = off_seg,
    ggplot2::aes(x = x_start, xend = x_end, y = y_med, yend = y_med),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.4, colour = "black"
  ) +
  {
    if (nrow(annot_df) > 0) {
      list(
        ggplot2::geom_path(
          data = bracket_path_df,
          ggplot2::aes(x = x, y = y, group = bracket_id),
          inherit.aes = FALSE, linewidth = 0.3, colour = "black",
          lineend = "square", linejoin = "mitre"
        ),
        ggplot2::geom_text(
          data = annot_df,
          ggplot2::aes(x = (x1 + x2) / 2, y = y_text, label = stars),
          inherit.aes = FALSE, size = FONT_STARS / ggplot2::.pt,
          fontface = "bold", colour = "black", family = FONT_FAMILY
        )
      )
    }
  } +
  ggplot2::geom_text(
    data = bottom_label_df,
    ggplot2::aes(x = xpos, y = y_label, label = label),
    inherit.aes = FALSE, size = FONT_N_LABEL / ggplot2::.pt,
    lineheight = 0.9, colour = "black", family = FONT_FAMILY
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Non-exp", "Dorsal", "Ventral"),
    expand = ggplot2::expansion(mult = c(0.06, 0.08))
  ) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = SIGMA),
    breaks = y_breaks,
    labels = scales::number_format(accuracy = 0.01),
    expand = ggplot2::expansion(mult = c(0.02, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(y_axis_bot, y_axis_top), clip = "off") +
  ggplot2::scale_fill_manual(values = cond_cols, drop = FALSE, name = NULL) +
  ggplot2::scale_colour_manual(values = cond_cols, drop = FALSE, name = NULL) +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::labs(x = NULL, y = expression("Hub dispersion (" * mu * "m)")) +
  ggplot2::theme_classic(base_size = FONT_AXIS_TITLE, base_family = FONT_FAMILY) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = FONT_AXIS_TEXT, colour = "black"),
    axis.text.y = ggplot2::element_text(size = FONT_AXIS_TEXT, colour = "black"),
    axis.title.y = ggplot2::element_text(size = FONT_AXIS_TITLE, margin = ggplot2::margin(r = 8)),
    axis.line = ggplot2::element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "black"),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(6, 18, 45, 6)
  )

print(p)

ggplot2::ggsave(
  filename = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_hub_dispersion_OFF_vs_Dorsal_and_Ventral_no_manual_exclusions.pdf"
  ),
  plot = p,
  width = 2.5,
  height = 3.1,
  units = "in",
  device = cairo_pdf
)

write.csv(
  stats_all,
  file = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_hub_dispersion_OFF_vs_Dorsal_and_Ventral_no_manual_exclusions_stats.csv"
  ),
  row.names = FALSE
)

write.csv(
  qc_table,
  file = paste0(
    "F:/Thesis/PAPER PDFS/",
    "WT_hub_dispersion_no_manual_exclusions_QC_table.csv"
  ),
  row.names = FALSE
)

