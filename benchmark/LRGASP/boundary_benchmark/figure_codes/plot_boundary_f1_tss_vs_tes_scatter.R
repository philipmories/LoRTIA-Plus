# ============================================================
# Scatter only: TSS vs TES F1 (with grey header strip)
# ============================================================
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

tss <- read_tsv("PCRMethods.tsv", show_col_types = FALSE) %>%
  mutate(
    Chemistry = if_else(Chemistry == "NA", "PacBio", Chemistry),
    Annotator = if_else(Annotator == "Flair", "FLAIR", Annotator),
    Annotator = if_else(Annotator == "BAMBU", "bambu", Annotator)
  )

tes <- read_tsv("PCRMethodsTES.tsv", show_col_types = FALSE) %>%
  mutate(
    Chemistry = if_else(Chemistry == "NA", "PacBio", Chemistry),
    Annotator = if_else(Annotator == "Flair", "FLAIR", Annotator),
    Annotator = if_else(Annotator == "BAMBU", "bambu", Annotator)
  )

chem_order <- c("ONT-dRNA", "ONT-cDNA", "ONT-CapTrap", "PacBio", "PacBio-CapTrap")
annotators <- c("LoRTIA", "IsoQuant", "bambu", "FLAIR", "NAGATA")

chem_shapes <- c(
  "ONT-dRNA"       = 15,
  "ONT-cDNA"       = 17,
  "ONT-CapTrap"    = 16,
  "PacBio"         = 3,
  "PacBio-CapTrap" = 4
)

annot_cols <- c(
  "LoRTIA"   = "#0072B2",
  "IsoQuant" = "#009E73",
  "bambu"    = "#E69F00",
  "FLAIR"    = "#CC79A7",
  "NAGATA"   = "#999999"
)

tss <- tss %>% mutate(Chemistry = factor(Chemistry, levels = chem_order),
                      Annotator = factor(Annotator, levels = annotators))
tes <- tes %>% mutate(Chemistry = factor(Chemistry, levels = chem_order),
                      Annotator = factor(Annotator, levels = annotators))

tt <- tss %>%
  select(`Cell line`, Chemistry, Annotator, F1_TSS = F1) %>%
  inner_join(
    tes %>% select(`Cell line`, Chemistry, Annotator, F1_TES = F1),
    by = c("Cell line", "Chemistry", "Annotator")
  ) %>%
  mutate(Header = "TSS vs TES F1 across annotators and chemistries")

p_tss_tes <- ggplot(tt, aes(x = F1_TSS, y = F1_TES,
                            colour = Annotator, shape = Chemistry)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = annot_cols, name = "Annotator") +
  scale_shape_manual(values = chem_shapes, name = "Chemistry") +
  scale_x_continuous(
    limits = c(0.10, 0.60),
    breaks = seq(0.10, 0.60, 0.10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0.10, 0.80),
    breaks = seq(0.10, 0.80, 0.10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  facet_grid(. ~ Header) +
  labs(x = "TSS F1", y = "TES F1") +
  theme_bw(base_size = 9) +
  theme(
    text             = element_text(family = "sans"),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", linewidth = 0.4),
    
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 9),
    
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 8),
    
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    
    plot.margin      = margin(4, 4, 4, 4)
  )

p_tss_tes

ggsave("TSS_vs_TES_scatter.pdf", p_tss_tes, width = 7, height = 3.2, dpi = 300)