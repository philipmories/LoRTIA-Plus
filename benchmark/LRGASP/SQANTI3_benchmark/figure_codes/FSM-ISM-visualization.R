# ============================================================
# Combined SQANTI figure (A/B forests + C trade-off)
# with legends for A/B, grey header for C, centered titles for A/B
# ============================================================
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

chem_order <- c("ONT-dRNA","ONT-cDNA","ONT-CapTrap","PacBio","PacBio-CapTrap")
other_order <- c("IsoQuant","bambu","FLAIR","NAGATA")
annot_order <- c("LoRTIA","IsoQuant","bambu","FLAIR","NAGATA")

# comparator colors (forest A/B)
other_cols <- c(
  "IsoQuant" = "#009E73",
  "bambu"    = "#E69F00",
  "FLAIR"    = "#CC79A7",
  "NAGATA"   = "#999999"
)

# scatter colors
annot_cols <- c(
  "LoRTIA"   = "#0072B2",
  "IsoQuant" = "#009E73",
  "bambu"    = "#E69F00",
  "FLAIR"    = "#CC79A7",
  "NAGATA"   = "#999999"
)

chem_shapes <- c(
  "ONT-dRNA"       = 15,
  "ONT-cDNA"       = 17,
  "ONT-CapTrap"    = 16,
  "PacBio"         = 3,
  "PacBio-CapTrap" = 4
)

norm_annot <- function(x) {
  x <- ifelse(x == "Flair", "FLAIR", x)
  x <- ifelse(x == "BAMBU", "bambu", x)
  x
}

theme_nar <- function() {
  theme_bw(base_size = 9) +
    theme(
      text             = element_text(family = "sans"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "black", linewidth = 0.4),
      strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
      strip.text       = element_text(face = "bold", size = 9),
      axis.text.x      = element_text(size = 8),
      axis.text.y      = element_text(size = 8),
      plot.margin      = margin(4,4,4,4)
    )
}

# ============================================================
# 1) Forest plots (A/B) from posthoc table
# ============================================================
ph <- read_tsv("sqanti_fsmism_lortia_vs_others.tsv", show_col_types = FALSE) %>%
  mutate(
    other = norm_annot(Annotator_2),
    Chemistry = factor(Chemistry, levels = chem_order),
    other = factor(other, levels = other_order),
    sig = !is.na(p_adj_BH) & p_adj_BH < 0.05,
    se = sd_delta / sqrt(n_cells),
    tcrit = qt(0.975, df = pmax(n_cells - 1, 1)),
    ci_low = mean_delta - tcrit * se,
    ci_high = mean_delta + tcrit * se
  )

plot_forest <- function(metric_name, title_text, show_legend = TRUE) {
  df <- ph %>% filter(Metric == metric_name)
  
  p <- ggplot(df, aes(x = mean_delta, y = other)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.18, linewidth = 0.6, colour = "grey35") +
    geom_point(aes(colour = other, shape = sig), size = 2.7, stroke = 0.9) +
    scale_colour_manual(values = other_cols, name = "Comparator") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       name = "BH-adjusted\nP < 0.05") +
    facet_grid(. ~ Chemistry) +
    labs(x = expression(Delta*" (LoRTIA − other) with 95% CI"), y = NULL) +
    ggtitle(title_text) +
    theme_nar() +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7)
    ) +
    guides(
      colour = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(size = 3)),
      shape  = guide_legend(order = 2, nrow = 1, byrow = TRUE, override.aes = list(size = 3))
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}

pA <- plot_forest(
  metric_name = "Recall_FSM_ISM",
  title_text  = "Effect sizes by chemistry: active-catalog recovery (FSM+ISM)",
  show_legend = FALSE
)

pB <- plot_forest(
  metric_name = "FL_ratio_FSM_over_FSM_ISM",
  title_text  = "Effect sizes by chemistry: full-lengthness among recovered transcripts",
  show_legend = TRUE
)

# ============================================================
# 2) Trade-off scatter (C) with grey header strip
# ============================================================
recall <- read_tsv("SQANTI_Recall_FSM_ISM.tsv", show_col_types = FALSE) %>%
  mutate(
    Annotator = norm_annot(Annotator),
    Chemistry = factor(Chemistry, levels = chem_order),
    Annotator = factor(Annotator, levels = annot_order)
  ) %>%
  rename(Recovery = F1)

fl <- read_tsv("SQANTI_FL_ratio.tsv", show_col_types = FALSE) %>%
  mutate(
    Annotator = norm_annot(Annotator),
    Chemistry = factor(Chemistry, levels = chem_order),
    Annotator = factor(Annotator, levels = annot_order)
  ) %>%
  rename(FullLength = F1)

tt <- recall %>%
  inner_join(fl, by = c("Chemistry","Cell line","Annotator")) %>%
  mutate(Header = "Trade-off between recovery and full-lengthness (SQANTI FSM/ISM)")

cent <- tt %>%
  group_by(Chemistry, Annotator) %>%
  summarise(Recovery = mean(Recovery), FullLength = mean(FullLength), .groups="drop")

pC <- ggplot(tt, aes(Recovery, FullLength)) +
  geom_point(aes(colour = Annotator, shape = Chemistry), alpha = 0.35, size = 1.8) +
  geom_point(data = cent, aes(colour = Annotator, shape = Chemistry),
             alpha = 1, size = 3.2, stroke = 0.9) +
  scale_colour_manual(values = annot_cols, name = "Annotator") +
  scale_shape_manual(values = chem_shapes, name = "Chemistry") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1), expand=expansion(mult=c(0,0.02))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1), expand=expansion(mult=c(0,0.02))) +
  coord_equal() +
  facet_grid(. ~ Header) +   # grey header
  labs(x = "FSM+ISM recovery (recall)", y = "Full-lengthness: FSM/(FSM+ISM)") +
  theme_nar() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 3)),
    shape  = guide_legend(order = 2, nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 3))
  )

# ============================================================
# 3) Combine into one figure
# Layout: A/B stacked, C on the right (wider)
# ============================================================
left_col <- pA / pB + plot_layout(heights = c(1, 1))
combined <- (left_col | pC) + plot_layout(widths = c(2.1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold", size=12))

combined
ggsave("SQANTI_Figure_ForestA_ForestB_TradeoffC.pdf",
       combined, width = 13.5, height = 6.4, dpi = 300)