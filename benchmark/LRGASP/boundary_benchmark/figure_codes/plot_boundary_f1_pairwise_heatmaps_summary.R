# ============================================================
# Combined heatmaps: TSS + TES side-by-side, one shared legend
# ============================================================
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

chem_order_full <- c("ONT-dRNA", "ONT-cDNA", "ONT-CapTrap", "PacBio", "PacBio-CapTrap")
other_order <- c("IsoQuant", "bambu", "FLAIR", "NAGATA")

make_lortia_heat <- function(pw, header_text) {
  
  pw <- pw %>%
    mutate(
      Annotator_1 = if_else(Annotator_1 == "Flair", "FLAIR", Annotator_1),
      Annotator_2 = if_else(Annotator_2 == "Flair", "FLAIR", Annotator_2),
      Annotator_1 = if_else(Annotator_1 == "BAMBU", "bambu", Annotator_1),
      Annotator_2 = if_else(Annotator_2 == "BAMBU", "bambu", Annotator_2)
    )
  
  df <- pw %>%
    filter(Annotator_1 == "LoRTIA" | Annotator_2 == "LoRTIA") %>%
    mutate(
      Other = if_else(Annotator_1 == "LoRTIA", Annotator_2, Annotator_1),
      delta_F1 = if_else(
        Annotator_1 == "LoRTIA",
        `mean_diff_F1_(a1_minus_a2)`,
        -`mean_diff_F1_(a1_minus_a2)`
      ),
      p_adj = p_adj_BH,
      Header = header_text
    ) %>%
    mutate(
      Other = factor(Other, levels = other_order),
      Chemistry = factor(Chemistry, levels = chem_order_full[chem_order_full %in% unique(Chemistry)])
    )
  
  ggplot(df, aes(x = Chemistry, y = Other, fill = delta_F1)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", delta_F1)), size = 2.6) +
    geom_point(
      data = df %>% filter(!is.na(p_adj) & p_adj < 0.05),
      shape = 8, size = 2.0, colour = "black",
      position = position_nudge(y = -0.25)
    ) +
    scale_fill_gradient2(
      low  = "steelblue",
      mid  = "white",
      high = "firebrick",
      midpoint = 0,
      oob = scales::squish,
      name = expression(Delta*"F1 (LoRTIA - other)")
    ) +
    facet_grid(. ~ Header) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      text             = element_text(family = "sans"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "black", linewidth = 0.4),
      strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
      strip.text       = element_text(face = "bold", size = 9),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = element_text(size = 8),
      legend.position  = "right",
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# ---- load pairwise ----
pw_tss <- read_tsv("PCRMethods__pairwise.tsv", show_col_types = FALSE)
pw_tes <- read_tsv(file.path( "TES_metrics.tsv_pairwise.tsv"), show_col_types = FALSE)

# ---- build plots ----
p_tss <- make_lortia_heat(pw_tss, "TSS: LoRTIA vs other annotators")
p_tes <- make_lortia_heat(pw_tes, "TES: LoRTIA vs other annotators")

# IMPORTANT: use the SAME fill limits for both panels (shared legend scale)
# Compute global limit across both
get_lim <- function(pw) {
  df <- pw %>%
    filter(Annotator_1 == "LoRTIA" | Annotator_2 == "LoRTIA") %>%
    mutate(delta = if_else(Annotator_1 == "LoRTIA",
                           `mean_diff_F1_(a1_minus_a2)`,
                           -`mean_diff_F1_(a1_minus_a2)`))
  max(abs(df$delta), na.rm = TRUE)
}
lim_global <- max(get_lim(pw_tss), get_lim(pw_tes))
lim_global <- max(0.10, lim_global)

# Re-apply identical scale to both plots
p_tss <- p_tss + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick",
                                      midpoint=0, limits=c(-lim_global, lim_global),
                                      oob=scales::squish,
                                      name=expression(Delta*"F1 (LoRTIA - other)"))
p_tes <- p_tes + scale_fill_gradient2(low="steelblue", mid="white", high="firebrick",
                                      midpoint=0, limits=c(-lim_global, lim_global),
                                      oob=scales::squish,
                                      name=expression(Delta*"F1 (LoRTIA - other)"))

# ---- combine with one shared legend ----
combined <- (p_tss | p_tes) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# preview
combined

# save
ggsave("TSS_TES_pairwise_heatmaps.pdf", combined, width=9.5, height=2.8, dpi=300)