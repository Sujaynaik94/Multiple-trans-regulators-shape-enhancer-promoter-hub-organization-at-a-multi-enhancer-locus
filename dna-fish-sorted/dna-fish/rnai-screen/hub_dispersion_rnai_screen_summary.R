# ============================================================
# RNAi — HUB DISPERSION PANELS + SUMMARY BARPLOT
# SINGLE VERTICAL PLOT
# ORDER: DORSAL, VENTRAL, NON-EXPRESSING, SUMMARY
# SUMMARY X-AXIS ALIGNED TO UPPER PANELS USING BLANK LacZ POSITION
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(ggtext)
  library(gghalves)
  library(tibble)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

# ====== INPUT ======
xlsx_path <- paste0(
  "Z:/Sujay/RNAi screen/DNA-FISH_RNAi_screen/",
  "AP_file_2026/ROI/",
  "ALL_GENES_ALL_ROIS_TRIPLETS_WITH_SUMMARY.xlsx"
)

sheet <- "Filtered_Triplets"
CONTROL_GENE <- "LacZ"
RADIAL_CUTOFF <- 0.5

# ====== PSEUDO-LOG ======
SIGMA <- 0.001
y_breaks <- c(0.1, 0.2, 0.5, 1, 2)

# ====== CONDITIONS ======
condition_levels <- c("Non-expressing", "Dorsal", "Ventral")

roi_to_condition <- c(
  "nonexp"  = "Non-expressing",
  "dorsal"  = "Dorsal",
  "ventral" = "Ventral"
)

# ====== GENE ORDER ======
# Fixed plotting order requested by the user.
# Only these genes are retained; CP190 and 3-hour Sxl-GFP are excluded.
plot_order <- c(
  "LacZ", "Awh", "pnr", "Ubx", "Blimp-1", "ttk",
  "lola", "Exd", "bowl", "hth", "pnt", "Med13",
  "CTCF", "BEAF-32", "Su(Hw)", "Mod(mdg4)"
)

all_gene_order <- plot_order

tf_screen <- c(
  "Awh", "pnr", "Ubx", "Blimp-1", "ttk",
  "lola", "Exd", "bowl", "hth", "pnt"
)

mediator_screen <- "Med13"

boundary_screen <- c(
  "CTCF", "BEAF-32", "Su(Hw)", "Mod(mdg4)"
)

# ====== GENE GROUPS ======
gene_groups <- tibble::tibble(
  gene = c(
    "LacZ",
    "CTCF", "Su(Hw)", "CP190", "BEAF-32", "Mod(mdg4)",
    "pnr", "Blimp-1", "pnt", "ttk", "bowl", "lola", "Awh",
    "Ubx", "Exd", "hth", "Med13"
  ),
  Group = c(
    "Control",
    rep("Boundary elements", 5),
    rep("TF + Hox/cofactors", 10),
    "Mediator"
  )
)

# ====== COLORS ======
gene_colors <- c(
  "LacZ"      = "#686765",
  "Awh"       = "#008088",
  "pnr"       = "#008088",
  "Ubx"       = "#008088",
  "Blimp-1"   = "#008088",
  "ttk"       = "#008088",
  "lola"      = "#008088",
  "Exd"       = "#008088",
  "bowl"      = "#008088",
  "hth"       = "#008088",
  "pnt"       = "#008088",
  "Med13"     = "#7B4FA3",
  "CTCF"      = "#E07B39",
  "BEAF-32"   = "#E07B39",
  "Su(Hw)"    = "#E07B39",
  "Mod(mdg4)" = "#E07B39"
)

side_bar_colors <- c(
  "shavenbaby OFF"          = "#686765",
  "shavenbaby ON (dorsal)"  = "#cc2829",
  "shavenbaby ON (ventral)" = "#2b3f99"
)

# ====== STYLE ======
box_w <- 0.35
violin_w <- 1.0
dot_band_w <- box_w / 2
dot_left_shift <- dot_band_w
dot_jit_x <- dot_band_w * 0.9

# ====== HELPERS ======
d3 <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

label_fun <- function(x) {
  sapply(x, function(g) {
    if (g == CONTROL_GENE) {
      "<span style='color:black;'>Control</span>"
    } else {
      paste0("<span style='color:black;'><i>", g, "</i></span>")
    }
  })
}

# ====== LOAD DATA ======
df_raw <- suppressMessages(
  readxl::read_excel(xlsx_path, sheet = sheet)
)

embryo_col <- pick_col(
  df_raw,
  c("FileName_Image", "filename_image", "embryo", "Embryo", "Image", "image")
)

if (is.na(embryo_col)) {
  stop(
    paste0(
      "Could not find embryo/image ID column.\nColumns present:\n",
      paste(names(df_raw), collapse = ", ")
    ),
    call. = FALSE
  )
}

normalize_roi_condition <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  x_clean <- gsub("[[:space:]_-]+", "", x_clean)

  dplyr::case_when(
    x_clean == "dorsal" ~ "Dorsal",
    x_clean == "ventral" ~ "Ventral",
    x_clean %in% c("nonexp", "nonexpressing", "off") ~ "Non-expressing",
    TRUE ~ NA_character_
  )
}

clean_gene_name <- function(x) {
  x_clean <- toupper(trimws(as.character(x)))
  x_clean <- sub("_[0-9]+$", "", x_clean)
  x_clean <- sub("_MASKS?$", "", x_clean)
  x_clean <- gsub("[[:space:]_-]+", "", x_clean)

  dplyr::recode(
    x_clean,
    "LACZ"   = "LacZ",
    "CTCF"   = "CTCF",
    "SUHW"   = "Su(Hw)",
    "SU(HW)" = "Su(Hw)",
    "MDG4"   = "Mod(mdg4)",
    "MOD(MDG4)" = "Mod(mdg4)",
    "BEAF32" = "BEAF-32",
    "CP190"  = "CP190",
    "BLIMP1" = "Blimp-1",
    "LOLAJ"  = "lola",
    "LOLA"   = "lola",
    "BOWL"   = "bowl",
    "PNR"    = "pnr",
    "TTK"    = "ttk",
    "PNT"    = "pnt",
    "UBX"    = "Ubx",
    "EXD"    = "Exd",
    "HTH"    = "hth",
    "MED13"  = "Med13",
    "AWH"    = "Awh",
    .default = x_clean
  )
}

df <- df_raw %>%
  dplyr::mutate(
    gene = clean_gene_name(gene_id),
    Condition = normalize_roi_condition(roi),
    EmbryoID = as.character(.data[[embryo_col]])
  ) %>%
  dplyr::filter(
    !is.na(Condition),
    !is.na(gene),
    gene %in% plot_order
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = plot_order),
    Condition = factor(Condition, levels = condition_levels)
  ) %>%
  dplyr::left_join(gene_groups, by = "gene")

cat("\nFixed plotting order:\n")
print(plot_order)

cat("\nRows retained after gene and ROI cleaning:", nrow(df), "\n")
cat("\nRetained data by condition and gene:\n")
print(
  df %>%
    dplyr::count(Condition, gene, name = "n") %>%
    tidyr::pivot_wider(
      names_from = gene,
      values_from = n,
      values_fill = 0
    )
)

control_counts <- df %>%
  dplyr::filter(gene == CONTROL_GENE) %>%
  dplyr::count(Condition, name = "n_control")

cat("\nLacZ control rows by condition:\n")
print(control_counts)

missing_control_conditions <- setdiff(
  condition_levels,
  as.character(control_counts$Condition[control_counts$n_control > 0])
)

if (length(missing_control_conditions) > 0) {
  stop(
    paste0(
      "LacZ control is missing after cleaning for: ",
      paste(missing_control_conditions, collapse = ", "),
      ". Raw workbook control is expected to be LACZ_MASKS."
    ),
    call. = FALSE
  )
}

# ====== COMPUTE HUB DISPERSION ======
df_area <- df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    centroid_x = (x_DG_um + x_E_um + x_svb_um) / 3,
    centroid_y = (y_DG_um + y_E_um + y_svb_um) / 3,
    centroid_z = (z_DG_um + z_E_um + z_svb_um) / 3,
    d_DG  = d3(x_DG_um,  y_DG_um,  z_DG_um,  centroid_x, centroid_y, centroid_z),
    d_E   = d3(x_E_um,   y_E_um,   z_E_um,   centroid_x, centroid_y, centroid_z),
    d_svb = d3(x_svb_um, y_svb_um, z_svb_um, centroid_x, centroid_y, centroid_z),
    Value = d_DG + d_E + d_svb
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, Group, Condition, EmbryoID, Value) %>%
  dplyr::filter(is.finite(Value), Value > 0)

# ============================================================
# CONDITION PANEL FUNCTION
# ============================================================
make_panel <- function(df_in, condition_name) {

  df_p <- df_in %>%
    dplyr::filter(
      Condition == condition_name,
      as.character(gene) %in% plot_order
    ) %>%
    dplyr::mutate(
      gene = factor(as.character(gene), levels = plot_order),
      gene = droplevels(gene),
      x_num = as.numeric(gene),
      gene_col = gene_colors[as.character(gene)]
    )

  if (!(CONTROL_GENE %in% as.character(df_p$gene))) return(NULL)

  y_max_data <- max(df_p$Value, na.rm = TRUE)
  y_min_data <- min(df_p$Value, na.rm = TRUE)

  y_axis_bot <- y_min_data * 0.5
  y_axis_top <- y_max_data * 3
  y_lab <- y_max_data * 1.6

  ctrl_vals <- df_p %>%
    dplyr::filter(as.character(gene) == CONTROL_GENE) %>%
    dplyr::pull(Value)

  stats_tbl <- purrr::map_dfr(
    setdiff(levels(df_p$gene), CONTROL_GENE),
    function(g) {
      gene_vals <- df_p %>%
        dplyr::filter(as.character(gene) == g) %>%
        dplyr::pull(Value)

      p <- if (length(gene_vals) >= 3 && length(ctrl_vals) >= 3) {
        suppressWarnings(
          stats::wilcox.test(gene_vals, ctrl_vals, exact = FALSE)$p.value
        )
      } else {
        NA_real_
      }

      tibble::tibble(gene = g, p_raw = p)
    }
  ) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_raw, method = "BH"),
      stars = dplyr::case_when(
        is.na(p_adj)   ~ "ns",
        p_adj <= 0.001 ~ "***",
        p_adj <= 0.01  ~ "**",
        p_adj <= 0.05  ~ "*",
        TRUE           ~ "ns"
      ),
      gene = factor(gene, levels = levels(df_p$gene))
    )

  ctrl_median <- stats::median(ctrl_vals, na.rm = TRUE)

  bottom_label_df <- df_p %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      N = dplyr::n_distinct(EmbryoID),
      percent_below = 100 * sum(Value < RADIAL_CUTOFF, na.rm = TRUE) / n_total,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene = factor(as.character(gene), levels = plot_order),
      label = paste0(round(percent_below), "%\nn=", n_total, "\nN=", N),
      y_label = y_axis_bot * 1.35
    )

  embryo_median_df <- df_p %>%
    dplyr::group_by(gene, EmbryoID, x_num, gene_col) %>%
    dplyr::summarise(
      embryo_median = stats::median(Value, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(df_p, ggplot2::aes(x = gene, y = Value)) +
    gghalves::geom_half_violin(
      ggplot2::aes(group = gene),
      side = "r",
      width = violin_w,
      fill = NA,
      alpha = 1,
      linewidth = 0.5,
      colour = "black",
      trim = TRUE
    ) +
    ggplot2::geom_point(
      data = df_p,
      ggplot2::aes(x = x_num - dot_left_shift, y = Value, color = gene_col),
      inherit.aes = FALSE,
      shape = 16,
      size = 0.9,
      alpha = 0.22,
      position = ggplot2::position_jitter(width = dot_jit_x, height = 0, seed = 1)
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = gene, fill = gene_col),
      outlier.shape = NA,
      linewidth = 0.5,
      width = box_w,
      colour = "black",
      alpha = 0.30
    ) +
    ggplot2::geom_point(
      data = embryo_median_df,
      ggplot2::aes(x = x_num, y = embryo_median, fill = gene_col),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.5,
      alpha = 0.85,
      stroke = 0.5,
      colour = "white",
      position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 42)
    ) +
    ggplot2::geom_hline(
      yintercept = ctrl_median,
      linetype = "dotted",
      linewidth = 0.55,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = dplyr::left_join(
        stats_tbl,
        tibble::tibble(
          gene = factor(names(gene_colors), levels = plot_order),
          gene_col = gene_colors
        ),
        by = "gene"
      ),
      ggplot2::aes(x = gene, y = y_lab, label = stars, color = gene_col),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold",
      family = "Arial"
    ) +
    ggplot2::geom_text(
      data = bottom_label_df,
      ggplot2::aes(x = gene, y = y_label, label = label),
      inherit.aes = FALSE,
      size = 3.5,
      lineheight = 0.9,
      color = "black",
      family = "Arial"
    ) +
    ggplot2::scale_x_discrete(
      labels = label_fun,
      expand = ggplot2::expansion(add = c(0.55, 0.85))
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = SIGMA),
      breaks = y_breaks,
      labels = scales::number_format(accuracy = 0.01),
      expand = ggplot2::expansion(mult = c(0.02, 0))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_axis_bot, y_axis_top),
      clip = "off"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(
      title = condition_name,
      x = NULL,
      y = expression("Hub dispersion (" * mu * "m)")
    ) +
    ggplot2::theme_classic(base_size = 20, base_family = "Arial") +
    ggplot2::theme(
      axis.text.x = ggtext::element_markdown(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        color = "black",
        size = 16
      ),
      axis.text.y = ggplot2::element_text(size = 14, colour = "black"),
      axis.title.y = ggplot2::element_text(size = 16, margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
      plot.margin = ggplot2::margin(2, 18, 10, 6)
    )
}

# ============================================================
# SUMMARY BARPLOT DATA
# ============================================================
df_disp <- df_area %>%
  dplyr::rename(Dispersion = Value)

summary_tbl <- df_disp %>%
  dplyr::group_by(gene, Group, Condition) %>%
  dplyr::summarise(
    mean_disp = mean(Dispersion, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

lacz_tbl <- summary_tbl %>%
  dplyr::filter(gene == CONTROL_GENE) %>%
  dplyr::select(Condition, lacz_mean_disp = mean_disp)

stats_summary_tbl <- df_disp %>%
  dplyr::filter(gene != CONTROL_GENE) %>%
  dplyr::group_by(Condition) %>%
  dplyr::group_modify(function(.x, .y) {

    ctrl_vals <- df_disp %>%
      dplyr::filter(
        Condition == .y$Condition,
        gene == CONTROL_GENE
      ) %>%
      dplyr::pull(Dispersion)

    .x %>%
      dplyr::distinct(gene) %>%
      dplyr::mutate(
        p_raw = purrr::map_dbl(gene, function(g) {
          gene_vals <- .x %>%
            dplyr::filter(gene == g) %>%
            dplyr::pull(Dispersion)

          if (length(ctrl_vals) >= 3 && length(gene_vals) >= 3) {
            suppressWarnings(
              stats::wilcox.test(gene_vals, ctrl_vals, exact = FALSE)$p.value
            )
          } else {
            NA_real_
          }
        }),
        p_adj = stats::p.adjust(p_raw, method = "BH"),
        signif_label = dplyr::case_when(
          is.na(p_adj)   ~ "",
          p_adj <= 0.001 ~ "***",
          p_adj <= 0.01  ~ "**",
          p_adj <= 0.05  ~ "*",
          TRUE           ~ ""
        ),
        is_sig = signif_label != ""
      )
  }) %>%
  dplyr::ungroup()

effect_tbl <- summary_tbl %>%
  dplyr::filter(gene != CONTROL_GENE) %>%
  dplyr::left_join(lacz_tbl, by = "Condition") %>%
  dplyr::mutate(
    pct_change_vs_lacz = 100 * ((mean_disp / lacz_mean_disp) - 1)
  ) %>%
  dplyr::left_join(stats_summary_tbl, by = c("gene", "Condition"))

# Keep same x-axis categories as upper panels, including blank LacZ.
all_gene_order <- all_gene_order[
  all_gene_order %in% c("LacZ", unique(df_disp$gene))
]

tf_screen <- tf_screen[tf_screen %in% all_gene_order]
mediator_screen <- mediator_screen[mediator_screen %in% all_gene_order]
boundary_screen <- boundary_screen[boundary_screen %in% all_gene_order]

gene_group_map <- tibble::tibble(
  gene = all_gene_order,
  Group = dplyr::case_when(
    all_gene_order == "LacZ" ~ "Control placeholder",
    all_gene_order %in% tf_screen ~ "TF + Hox/cofactors",
    all_gene_order %in% mediator_screen ~ "Mediator",
    all_gene_order %in% boundary_screen ~ "Boundary elements",
    TRUE ~ NA_character_
  )
)

dodge_w <- 0.78
bar_w <- 0.70

offset_map <- c(
  "shavenbaby OFF" = -dodge_w / 3,
  "shavenbaby ON (dorsal)" = 0,
  "shavenbaby ON (ventral)" = dodge_w / 3
)

plot_df <- effect_tbl %>%
  dplyr::filter(gene %in% setdiff(all_gene_order, "LacZ")) %>%
  dplyr::mutate(
    gene = factor(gene, levels = all_gene_order),
    side_condition = dplyr::recode(
      as.character(Condition),
      "Non-expressing" = "shavenbaby OFF",
      "Dorsal" = "shavenbaby ON (dorsal)",
      "Ventral" = "shavenbaby ON (ventral)"
    ),
    side_condition = factor(side_condition, levels = names(side_bar_colors)),
    x_num = as.numeric(gene),
    x_center = x_num + unname(offset_map[as.character(side_condition)])
  )

star_df <- plot_df %>%
  dplyr::filter(is_sig) %>%
  dplyr::mutate(
    star_y = dplyr::case_when(
      pct_change_vs_lacz > 0 ~ pct_change_vs_lacz + 3,
      pct_change_vs_lacz < 0 ~ pct_change_vs_lacz - 3,
      TRUE ~ 3
    ),
    star_vjust = dplyr::case_when(
      pct_change_vs_lacz > 0 ~ 0,
      pct_change_vs_lacz < 0 ~ 1,
      TRUE ~ 0
    )
  )

xpos_tbl <- tibble::tibble(
  gene = factor(all_gene_order, levels = all_gene_order),
  x_num = seq_along(all_gene_order)
) %>%
  dplyr::left_join(gene_group_map, by = "gene")

bg_tbl <- xpos_tbl %>%
  dplyr::filter(Group != "Control placeholder") %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    xmin = min(x_num) - 0.5,
    xmax = max(x_num) + 0.5,
    .groups = "drop"
  )

separator_x <- c(
  length(tf_screen) + 1.5,
  length(tf_screen) + length(mediator_screen) + 1.5
)

y_values <- c(plot_df$pct_change_vs_lacz, star_df$star_y)
y_top <- ceiling(max(y_values, na.rm = TRUE) / 10) * 10 + 10
y_bot <- ceiling(abs(min(y_values, na.rm = TRUE)) / 10) * 10 + 10

gene_label_map <- stats::setNames(
  ifelse(
    all_gene_order == "LacZ",
    "",
    paste0("<span style='color:black;'><i>", all_gene_order, "</i></span>")
  ),
  all_gene_order
)

p_summary <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = bg_tbl,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = Group
    ),
    alpha = 0.45,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_col(
    data = plot_df,
    ggplot2::aes(
      x = gene,
      y = pct_change_vs_lacz,
      fill = side_condition,
      group = side_condition
    ),
    position = ggplot2::position_dodge(width = dodge_w),
    width = bar_w,
    color = "black",
    linewidth = 0.2
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.7,
    color = "black"
  ) +
  ggplot2::geom_vline(
    xintercept = separator_x,
    linewidth = 0.8,
    color = "grey55"
  ) +
  ggplot2::geom_text(
    data = star_df,
    ggplot2::aes(
      x = x_center,
      y = star_y,
      label = signif_label,
      vjust = star_vjust
    ),
    size = 5,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "TF + Hox/cofactors" = "#008088",
      "Mediator" = "#7B4FA3",
      "Boundary elements" = "#E07B39",
      side_bar_colors
    ),
    breaks = names(side_bar_colors),
    name = NULL
  ) +
  ggplot2::scale_x_discrete(
    limits = all_gene_order,
    labels = gene_label_map,
    expand = ggplot2::expansion(add = c(0.55, 0.85))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-y_bot, y_top),
    breaks = pretty(c(-y_bot, y_top), n = 8),
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Summary: hub dispersion — % change vs LacZ control",
    x = NULL,
    y = "% change in hub dispersion vs LacZ"
  ) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    axis.text.x = ggtext::element_markdown(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 18,
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(size = 16, color = "black"),
    axis.title.y = ggplot2::element_text(size = 18, color = "black"),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 14),
    panel.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 18, 8, 6)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
  )

# ============================================================
# FINAL SINGLE VERTICAL PLOT
# ORDER: DORSAL, VENTRAL, NON-EXPRESSING, SUMMARY
# ============================================================
p_all_dorsal <- make_panel(df_area, "Dorsal")
p_all_ventral <- make_panel(df_area, "Ventral")
p_all_nonexp <- make_panel(df_area, "Non-expressing")

plot_objects <- list(
  Dorsal = p_all_dorsal,
  Ventral = p_all_ventral,
  Non_expressing = p_all_nonexp,
  Summary = p_summary
)

invalid_plots <- names(plot_objects)[
  !vapply(
    plot_objects,
    function(x) inherits(x, c("gg", "ggplot", "patchwork")),
    logical(1)
  )
]

if (length(invalid_plots) > 0) {
  stop(
    paste0(
      "These plot objects were not created correctly: ",
      paste(invalid_plots, collapse = ", "),
      ". Review the printed condition × gene counts above."
    ),
    call. = FALSE
  )
}

p_single <- patchwork::wrap_plots(
  plot_objects,
  ncol = 1,
  heights = c(1, 1, 1, 0.9)
)

print(p_single)
ggplot2::ggsave(
  filename = "F:\\Thesis\\PAPER PDFS\\rnai.pdf",
  plot = p_single,
  width = 16,
  height = 20,
  units = "in",
  device = cairo_pdf
)
