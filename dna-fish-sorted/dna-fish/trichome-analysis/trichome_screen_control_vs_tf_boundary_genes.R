####for the paper tf and boundary together+compacted
# =========================================================
# RNAi TRICHOME SCREEN
# Standalone paper plot:
#   One combined panel
#   Control -> selected TF/Hox genes -> boundary elements
#   Gene labels at 45 degrees
#   Reduced spacing between genes
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(ggtext)
  library(ggh4x)
})

# =========================================================
# PARAMETERS
# =========================================================
raw_point_size   <- 2
raw_point_alpha  <- 0.55
raw_jitter_width <- 0.015
raw_x_shift      <- 0.06

median_dot_size <- 3
median_x_shift  <- -0.04

whisker_line_width <- 1.0

gene_spacing <- 0.78

Y_MIN_DORSAL_LATERAL <- 200
Y_MIN_VENTRAL        <- 100

Y_MAX_DORSAL_LATERAL <- 550
Y_MAX_VENTRAL        <- 280

# =========================================================
# FUNCTIONS
# =========================================================
box_stats <- function(x) {
  x <- x[!is.na(x)]
  qs <- grDevices::boxplot.stats(x)$stats
  
  tibble::tibble(
    ymin   = qs[1],
    lower  = qs[2],
    median = qs[3],
    upper  = qs[4],
    ymax   = qs[5]
  )
}

# =========================================================
# LOAD DATA
# =========================================================
setwd("Z:/Sujay/RNAi screen/Trichome_RNAi_Screen")

df_DL <- readxl::read_excel(
  "Trichome_quantification_final.xlsx",
  sheet = 1
)

df_V <- readxl::read_excel(
  "Trichome_quantification_final.xlsx",
  sheet = 2
)

df_long_DL <- df_DL %>%
  dplyr::mutate(`Dorsal+Lateral` = Dorsal + Lateral) %>%
  tidyr::pivot_longer(
    cols = c(`Dorsal+Lateral`),
    names_to = "Region",
    values_to = "Trichome_Count"
  )

df_long_V <- df_V %>%
  tidyr::pivot_longer(
    cols = c(Ventral),
    names_to = "Region",
    values_to = "Trichome_Count"
  )

df_long <- dplyr::bind_rows(df_long_DL, df_long_V) %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::mutate(
    Gene = trimws(as.character(Gene)),
    Gene = dplyr::recode(
      Gene,
      "LacZ"       = "Control",
      
      "suhw"       = "Su(Hw)",
      "Suhw"       = "Su(Hw)",
      "SUHW"       = "Su(Hw)",
      "su(Hw)"     = "Su(Hw)",
      "Su(Hw)"     = "Su(Hw)",
      "SuHw"       = "Su(Hw)",
      
      "Blimp1"     = "Blimp-1",
      "BLIMP1"     = "Blimp-1",
      
      "TTK"        = "ttk",
      
      "Bowl"       = "bowl",
      "BOWL"       = "bowl",
      
      "ubx"        = "Ubx",
      "UBX"        = "Ubx",
      
      "Exd"        = "Exd",
      "exd"        = "Exd",
      "EXD"        = "Exd",
      
      "Hth"        = "hth",
      "HTH"        = "hth",
      
      "LolaJ"      = "lola",
      "LOLAJ"      = "lola",
      "Lola"       = "lola",
      "LOLA"       = "lola",
      
      "awh"        = "Awh",
      "AWH"        = "Awh",
      
      "mod(mdg4)"  = "Mod(mdg4)",
      "MDG4"       = "Mod(mdg4)",
      "Mod(mdg4)"  = "Mod(mdg4)",
      
      "BEAF32"     = "BEAF-32",
      "BEAF-32"    = "BEAF-32",
      
      .default = Gene
    ),
    Region = as.character(Region)
  )

# =========================================================
# GENE GROUPS + ORDER
# =========================================================
plot_order <- c(
  "Control",
  
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
  "Mod(mdg4)"
)

gene_groups <- tibble::tibble(
  Gene = plot_order,
  Group = c(
    "Control",
    rep("Transcription factors", 7),
    rep("Boundary elements", 4)
  )
)

df_long <- df_long %>%
  dplyr::left_join(gene_groups, by = "Gene") %>%
  dplyr::filter(Gene %in% plot_order) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = plot_order),
    Group = factor(
      Group,
      levels = c(
        "Control",
        "Transcription factors",
        "Boundary elements"
      )
    ),
    Region = factor(
      Region,
      levels = c("Dorsal+Lateral", "Ventral"),
      labels = c("Dorsal", "Ventral")
    )
  )

# =========================================================
# COLORS + LABELS
# =========================================================
group_colors <- c(
  "Control"               = "black",
  "Transcription factors" = "#008088",
  "Boundary elements"     = "#E07B39"
)

x_labels_md <- paste0(
  "<span style='color:black;'><i>",
  plot_order,
  "</i></span>"
)

x_labels_md[plot_order == "Control"] <-
  "<span style='color:black;'>Control</span>"

names(x_labels_md) <- plot_order

# =========================================================
# X POSITIONS WITH REDUCED SPACING
# =========================================================
x_position_tbl <- tibble::tibble(
  Gene = factor(plot_order, levels = plot_order),
  x_base = seq_along(plot_order) * gene_spacing
)

df_long <- df_long %>%
  dplyr::left_join(x_position_tbl, by = "Gene") %>%
  dplyr::mutate(
    x_raw = x_base + raw_x_shift,
    x_median = x_base + median_x_shift
  )

# =========================================================
# WILCOXON TESTS  [UPDATED: BH correction added within each Region]
# =========================================================
control_name <- "Control"

pval_df <- expand.grid(
  Gene   = setdiff(plot_order, control_name),
  Region = levels(df_long$Region),
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_value = {
      g <- Gene
      r <- Region
      
      x <- df_long %>%
        dplyr::filter(Gene == g, Region == r) %>%
        dplyr::pull(Trichome_Count)
      
      y <- df_long %>%
        dplyr::filter(Gene == control_name, Region == r) %>%
        dplyr::pull(Trichome_Count)
      
      if (length(x) < 2 || length(y) < 2) {
        NA_real_
      } else {
        suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
      }
    }
  ) %>%
  dplyr::ungroup() %>%
  # ---- BH correction across genes tested, within each Region ----
dplyr::group_by(Region) %>%
  dplyr::mutate(
    p_adj = stats::p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_signif = dplyr::case_when(
      is.na(p_adj)   ~ "ns",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  dplyr::left_join(x_position_tbl, by = "Gene") %>%
  dplyr::mutate(
    x_median = x_base + median_x_shift
  )

# =========================================================
# SUMMARY STATS
# =========================================================
summary_df <- df_long %>%
  dplyr::group_by(Gene, Region, Group, x_median) %>%
  dplyr::group_modify(~ box_stats(.x$Trichome_Count)) %>%
  dplyr::ungroup()

control_median_df <- summary_df %>%
  dplyr::filter(Gene == control_name) %>%
  dplyr::select(Region, control_median = median)

star_y_sub <- tibble::tibble(
  Region = factor(
    c("Dorsal", "Ventral"),
    levels = c("Dorsal", "Ventral")
  ),
  y = c(
    Y_MAX_DORSAL_LATERAL - 25,
    Y_MAX_VENTRAL - 15
  )
)

pval_plot <- pval_df %>%
  dplyr::left_join(star_y_sub, by = "Region")

# =========================================================
# PLOT
# =========================================================
p_combined <- ggplot2::ggplot() +
  
  ggplot2::geom_jitter(
    data = df_long,
    ggplot2::aes(
      x = x_raw,
      y = Trichome_Count
    ),
    inherit.aes = FALSE,
    color = "grey35",
    width = raw_jitter_width,
    shape = 1,
    stroke = 1,
    height = 0,
    size = raw_point_size,
    alpha = raw_point_alpha
  ) +
  
  ggplot2::geom_linerange(
    data = summary_df,
    ggplot2::aes(
      x = x_median,
      ymin = ymin,
      ymax = ymax,
      color = Group
    ),
    inherit.aes = FALSE,
    linewidth = whisker_line_width
  ) +
  
  ggplot2::geom_hline(
    data = control_median_df,
    ggplot2::aes(yintercept = control_median),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    data = summary_df,
    ggplot2::aes(
      x = x_median,
      y = median,
      color = Group
    ),
    inherit.aes = FALSE,
    shape = 16,
    size = median_dot_size
  ) +
  
  ggplot2::geom_text(
    data = pval_plot %>% dplyr::filter(p_signif != "ns"),
    ggplot2::aes(
      x = x_median,
      y = y,
      label = p_signif
    ),
    inherit.aes = FALSE,
    size = 7.5,
    fontface = "bold",
    color = "black"
  ) +
  
  ggplot2::geom_text(
    data = pval_plot %>% dplyr::filter(p_signif == "ns"),
    ggplot2::aes(
      x = x_median,
      y = y,
      label = p_signif
    ),
    inherit.aes = FALSE,
    size = 5.5,
    color = "grey60"
  ) +
  
  ggplot2::facet_grid(
    rows = ggplot2::vars(Region),
    switch = "y",
    scales = "free_y"
  ) +
  
  ggh4x::facetted_pos_scales(
    y = list(
      Region == "Dorsal" ~
        ggplot2::scale_y_continuous(
          limits = c(Y_MIN_DORSAL_LATERAL, Y_MAX_DORSAL_LATERAL),
          expand = ggplot2::expansion(mult = c(0, 0))
        ),
      Region == "Ventral" ~
        ggplot2::scale_y_continuous(
          limits = c(Y_MIN_VENTRAL, Y_MAX_VENTRAL),
          expand = ggplot2::expansion(mult = c(0, 0))
        )
    )
  ) +
  
  ggplot2::scale_color_manual(
    values = group_colors,
    name = "Gene class"
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = x_position_tbl$x_base,
    labels = x_labels_md[plot_order],
    expand = ggplot2::expansion(add = 0.25)
  ) +
  
  ggplot2::labs(
    title = "Transcription factors and boundary elements",
    x = NULL,
    y = "Number of Trichomes"
  ) +
  
  ggplot2::theme_classic(base_size = 24) +
  
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.y = unit(1.5, "lines"),
    
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(
      face = "plain",
      angle = 90,
      size = 22
    ),
    
    axis.text.x = ggtext::element_markdown(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 20,
      color = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 22,
      color = "black"
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 26,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 0.9,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.8,
      color = "black"
    ),
    
    axis.ticks.length = ggplot2::unit(4, "pt"),
    
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 12),
    legend.text = ggplot2::element_text(size = 18),
    
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      size = 26
    )
  )

print(p_combined)
