# =========================================================
# Figure X
# Catalog-external transcript discovery and support structure (NIC/NNC)
# =========================================================

library(tidyverse)
library(scales)
library(patchwork)

# -------------------------
# Input files
# -------------------------
burden_file    <- "nic_nnc_burden_by_condition.tsv"
support_file   <- "novel_junction_support_by_annotator_chemistry.tsv"
plaus_file     <- "novel_junction_plausibility_by_annotator_chemistry.tsv"
iso_file       <- "novel_isoform_canonicality_by_annotator_chemistry_category.tsv"
locus_file     <- "novel_locus_consistency_by_annotator_chemistry_category.tsv"

stop_if_missing <- function(files) {
  miss <- files[!file.exists(files)]
  if (length(miss) > 0) {
    stop(
      "Hiányzó fájl(ok):\n", paste(miss, collapse = "\n"),
      "\n\nAktuális mappa:\n", getwd(),
      "\n\nElérhető TSV-k:\n", paste(list.files(pattern = "tsv$", ignore.case = TRUE), collapse = "\n")
    )
  }
}

stop_if_missing(c(burden_file, support_file, plaus_file, iso_file, locus_file))

# -------------------------
# Ordering and labels
# -------------------------
annotator_levels <- c("LoRTIA", "BAMBU", "FLAIR", "IsoQuant", "NAGATA")
chemistry_levels <- c("ONT-CapTrap", "ONT-cDNA", "ONT-dRNA", "PacBio", "PacBio-CapTrap")

annotator_labels <- c(
  "LoRTIA" = "LoRTIA Plus",
  "BAMBU" = "bambu",
  "FLAIR" = "FLAIR",
  "IsoQuant" = "IsoQuant",
  "NAGATA" = "NAGATA"
)

chemistry_labels <- c(
  "ONT-CapTrap" = "ONT-CapTrap",
  "ONT-cDNA" = "ONT-cDNA",
  "ONT-dRNA" = "ONT-dRNA",
  "PacBio" = "PacBio cDNA",
  "PacBio-CapTrap" = "PacBio-CapTrap"
)

# -------------------------
# Theme
# -------------------------
base_theme <- theme_classic(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = "grey60", linewidth = 0.5),
    strip.text = element_text(face = "bold", colour = "black"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
    plot.title = element_text(face = "bold", hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# -------------------------
# Colors
# -------------------------
cat_cols <- c(
  "NIC" = "#F8766D",
  "NNC" = "#619CFF"
)

support_cols <- c(
  "Cross-tool" = "#4DAF4A",
  "Multi-chemistry" = "#377EB8",
  "Multi-cell-line" = "#984EA3",
  "Single-condition" = "#BDBDBD"
)

metric_cols <- c(
  "Junction canonicality" = "#1B9E77",
  "Isoform canonicality" = "#D95F02",
  ">=2 NNC isoforms / gene" = "#7570B3",
  "Independent-supported NNC junctions / gene" = "#E7298A"
)

# =========================================================
# PANEL A — NIC/NNC burden
# =========================================================
burden_df <- readr::read_tsv(burden_file, show_col_types = FALSE) %>%
  mutate(
    annotator = factor(annotator, levels = annotator_levels),
    chemistry = factor(chemistry, levels = chemistry_levels)
  ) %>%
  group_by(annotator, chemistry) %>%
  summarise(
    n_total_isoforms = sum(n_total_isoforms, na.rm = TRUE),
    n_nic = sum(n_nic, na.rm = TRUE),
    n_nnc = sum(n_nnc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    frac_nic = n_nic / n_total_isoforms,
    frac_nnc = n_nnc / n_total_isoforms
  ) %>%
  pivot_longer(
    cols = c(frac_nic, frac_nnc),
    names_to = "category",
    values_to = "fraction"
  ) %>%
  mutate(
    category = recode(category, frac_nic = "NIC", frac_nnc = "NNC"),
    category = factor(category, levels = c("NIC", "NNC"))
  )

pA <- ggplot(burden_df, aes(x = annotator, y = fraction, fill = category)) +
  geom_col(width = 0.72) +
  facet_wrap(~chemistry, nrow = 1, labeller = labeller(chemistry = chemistry_labels)) +
  scale_x_discrete(labels = annotator_labels) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Novelty burden",
    x = "Annotator",
    y = "Fraction of reconstructed isoforms"
  ) +
  base_theme

# =========================================================
# PANEL B — Novel junction support structure
# =========================================================
support_df <- readr::read_tsv(support_file, show_col_types = FALSE) %>%
  mutate(
    annotator = factor(annotator, levels = annotator_levels),
    chemistry = factor(chemistry, levels = chemistry_levels)
  ) %>%
  pivot_longer(
    cols = c(
      frac_cross_tool,
      frac_multi_chemistry_single_tool,
      frac_multi_cell_line_single_tool,
      frac_single_condition
    ),
    names_to = "support_class",
    values_to = "fraction"
  ) %>%
  mutate(
    support_class = recode(
      support_class,
      frac_cross_tool = "Cross-tool",
      frac_multi_chemistry_single_tool = "Multi-chemistry",
      frac_multi_cell_line_single_tool = "Multi-cell-line",
      frac_single_condition = "Single-condition"
    ),
    support_class = factor(
      support_class,
      levels = c("Cross-tool", "Multi-chemistry", "Multi-cell-line", "Single-condition")
    )
  )

pB <- ggplot(support_df, aes(x = annotator, y = fraction, fill = support_class)) +
  geom_col(width = 0.72) +
  facet_wrap(~chemistry, nrow = 1, labeller = labeller(chemistry = chemistry_labels)) +
  scale_x_discrete(labels = annotator_labels) +
  scale_fill_manual(values = support_cols, name = NULL) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Novel junction support structure",
    x = "Annotator",
    y = "Fraction of novel junctions"
  ) +
  base_theme

# =========================================================
# PANEL C — NNC canonicality
# =========================================================
plaus_df <- readr::read_tsv(plaus_file, show_col_types = FALSE) %>%
  filter(scope == "nnc_only") %>%
  transmute(
    annotator,
    chemistry,
    metric = "Junction canonicality",
    value = canonical_fraction
  )

iso_df <- readr::read_tsv(iso_file, show_col_types = FALSE) %>%
  filter(category == "NNC") %>%
  transmute(
    annotator,
    chemistry,
    metric = "Isoform canonicality",
    value = canonical_isoform_fraction
  )

canon_df <- bind_rows(plaus_df, iso_df) %>%
  mutate(
    annotator = factor(annotator, levels = annotator_levels),
    chemistry = factor(chemistry, levels = chemistry_levels),
    metric = factor(metric, levels = c("Junction canonicality", "Isoform canonicality"))
  )

pC <- ggplot(canon_df, aes(x = annotator, y = value, colour = metric, group = metric)) +
  geom_point(position = position_dodge(width = 0.45), size = 2.7) +
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.5) +
  facet_wrap(~chemistry, nrow = 1, labeller = labeller(chemistry = chemistry_labels)) +
  scale_x_discrete(labels = annotator_labels) +
  scale_colour_manual(values = metric_cols[c("Junction canonicality", "Isoform canonicality")], name = NULL) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Canonicality of NNC novelty",
    x = "Annotator",
    y = "Canonical fraction"
  ) +
  base_theme

# =========================================================
# PANEL D — NNC locus-level consistency
# =========================================================
locus_df <- readr::read_tsv(locus_file, show_col_types = FALSE) %>%
  filter(category == "NNC") %>%
  mutate(
    frac_independent_supported_genes = genes_with_independent_supported_novel_junction / genes_with_novelty
  ) %>%
  select(
    annotator, chemistry,
    frac_genes_with_ge_2_novel_isoforms,
    frac_independent_supported_genes
  ) %>%
  pivot_longer(
    cols = c(frac_genes_with_ge_2_novel_isoforms, frac_independent_supported_genes),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    annotator = factor(annotator, levels = annotator_levels),
    chemistry = factor(chemistry, levels = chemistry_levels),
    metric = recode(
      metric,
      frac_genes_with_ge_2_novel_isoforms = ">=2 NNC isoforms / gene",
      frac_independent_supported_genes = "Independent-supported NNC junctions / gene"
    ),
    metric = factor(
      metric,
      levels = c(">=2 NNC isoforms / gene", "Independent-supported NNC junctions / gene")
    )
  )

pD <- ggplot(locus_df, aes(x = annotator, y = value, colour = metric, group = metric)) +
  geom_point(position = position_dodge(width = 0.45), size = 2.7) +
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.5) +
  facet_wrap(~chemistry, nrow = 1, labeller = labeller(chemistry = chemistry_labels)) +
  scale_x_discrete(labels = annotator_labels) +
  scale_colour_manual(
    values = metric_cols[c(">=2 NNC isoforms / gene", "Independent-supported NNC junctions / gene")],
    name = NULL
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Locus-level consistency of NNC novelty",
    x = "Annotator",
    y = "Fraction of genes"
  ) +
  base_theme

# =========================================================
# Assemble figure
# =========================================================
p_all <- (pA / pB / pC / pD) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 13, colour = "black"),
    plot.tag.position = c(0.01, 0.99)
  )

print(p_all)

ggsave(
  filename = "FigureX_NIC_NNC_discovery.pdf",
  plot = p_all,
  width = 14,
  height = 12,
  dpi = 300
)

ggsave(
  filename = "FigureX_NIC_NNC_discovery.png",
  plot = p_all,
  width = 14,
  height = 12,
  dpi = 300
)