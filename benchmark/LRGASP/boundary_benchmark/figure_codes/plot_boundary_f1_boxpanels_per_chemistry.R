library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)

## ------------------------------------------------------------
## 0) Input fájlok
## ------------------------------------------------------------
file_tss <- "PCRMethods.tsv"
file_tes <- "PCRMethodsTES.tsv"

## ------------------------------------------------------------
## 1) Előfeldolgozás
## ------------------------------------------------------------
prep_f1 <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      Chemistry = if_else(is.na(Chemistry) | Chemistry == "NA", "PacBio", Chemistry),
      Annotator = if_else(Annotator == "Flair", "FLAIR", Annotator)
    )
  
  annotators  <- c("LoRTIA", setdiff(sort(unique(df$Annotator)), "LoRTIA"))
  chem_levels <- c("ONT-dRNA", "ONT-cDNA", "ONT-CapTrap", "PacBio", "PacBio-CapTrap")
  
  df %>%
    mutate(
      Annotator = factor(Annotator, levels = annotators),
      Chemistry = factor(Chemistry, levels = chem_levels),
      `Cell line` = as.factor(`Cell line`)
    )
}

f1_tss <- prep_f1(file_tss)
f1_tes <- prep_f1(file_tes)

## ------------------------------------------------------------
## 2) Színek / shape
## ------------------------------------------------------------
annot_cols <- c(
  "LoRTIA"   = "#0072B2",
  "IsoQuant" = "#009E73",
  "bambu"    = "#E69F00",
  "FLAIR"    = "#CC79A7",
  "NAGATA"   = "#999999"
)

cell_cols <- c(
  "H1"      = "#1b9e77",
  "H1-endo" = "#d95f02",
  "H1-DE"   = "#d95f02",
  "WTC11"   = "#7570b3"
)

cell_shapes <- c(
  "H1"      = 16,
  "H1-endo" = 17,
  "H1-DE"   = 17,
  "WTC11"   = 15
)

chem_order_common <- c("ONT-dRNA", "ONT-cDNA", "ONT-CapTrap", "PacBio", "PacBio-CapTrap")

## ------------------------------------------------------------
## 3) Szürke keretes "header" (facet nélkül, ezért bombabiztos)
## ------------------------------------------------------------
make_header <- function(label) {
  ggplot() +
    annotate("rect",
             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "grey90", colour = "black", linewidth = 0.4
    ) +
    annotate("text", x = 0, y = 0, label = label, fontface = "bold", size = 3) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(
      plot.margin = margin(0, 2, 0, 2)
    )
}

## ------------------------------------------------------------
## 4) Paneltest (box + pontok), üresre "no data"
## ------------------------------------------------------------
make_body <- function(df, chem_label, y_lab) {
  df_sub <- df %>% filter(as.character(Chemistry) == chem_label)
  
  if (nrow(df_sub) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "no data", size = 3) +
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
        theme_bw(base_size = 10) +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.4),
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.margin = margin(2, 2, 4, 2)
        )
    )
  }
  
  ymin <- min(df_sub$F1, na.rm = TRUE)
  ymax <- max(df_sub$F1, na.rm = TRUE)
  brks <- pretty(c(ymin, ymax), n = 4)
  y_low  <- min(brks)
  y_high <- max(brks)
  
  ggplot(df_sub, aes(x = Annotator, y = F1)) +
    geom_boxplot(
      aes(fill = Annotator),
      position      = position_nudge(x = -0.15),
      width         = 0.45,
      colour        = "black",
      linewidth     = 0.35,
      outlier.shape = NA,
      alpha         = 0.85
    ) +
    geom_point(
      aes(colour = `Cell line`, shape = `Cell line`),
      position = position_nudge(x = 0.25),
      size     = 2.0,
      alpha    = 0.95
    ) +
    scale_fill_manual(values = annot_cols, name = "Annotator") +
    scale_colour_manual(values = cell_cols,  name = "Cell line") +
    scale_shape_manual(values = cell_shapes, name = "Cell line") +
    scale_y_continuous(
      limits = c(y_low, y_high),
      breaks = brks,
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = y_lab) +
    theme_bw(base_size = 10) +
    theme(
      text  = element_text(family = "sans"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      panel.border       = element_rect(colour = "black", linewidth = 0.4),
      
      axis.title.y       = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 8),
      
      legend.position    = "right",
      plot.margin        = margin(2, 2, 4, 2)
    )
}

## ------------------------------------------------------------
## 5) Komplett panel (header + body)
## ------------------------------------------------------------
make_panel <- function(df, chem_label, y_lab) {
  make_header(chem_label) /
    make_body(df, chem_label, y_lab) +
    plot_layout(heights = c(0.14, 1))
}

## ------------------------------------------------------------
## 6) Panelek soronként
## ------------------------------------------------------------
panels_tss <- lapply(chem_order_common, function(ch) make_panel(f1_tss, ch, "TSS F1"))
panels_tes <- lapply(chem_order_common, function(ch) make_panel(f1_tes, ch, "TES F1"))

row_tss <- wrap_plots(panels_tss, nrow = 1)
row_tes <- wrap_plots(panels_tes, nrow = 1)

## ------------------------------------------------------------
## 7) 2 sor + közös legenda alul
## ------------------------------------------------------------
full_fig <- (row_tss / row_tes) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box      = "vertical"
  )

full_fig

## ------------------------------------------------------------
## 8) Mentés (ha kell)
## ------------------------------------------------------------
# ggsave("PCRMethods_TSS_TES_F1.pdf", full_fig, width = 14, height = 6, units = "in")
# ggsave("PCRMethods_TSS_TES_F1.png", full_fig, width = 14, height = 6, units = "in", dpi = 300)