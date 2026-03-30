library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(grid)

# ==========================================================
# BEÁLLÍTÁSOK
# ==========================================================
input_file <- "SQANTI3_basic_stat.tsv"
outdir <- "Ablock_IsoformYield_SQANTI3"

plot_versions <- c("main_avg_cells", "detailed_cellrows")

collapse_cells_method <- "weighted" # "weighted" vagy "simple"

default_categories <- c(
  "FSM","ISM","NIC","NNC",
  "Antisense","Fusion","GenicGenomic","GenicIntron","Intergenic"
)
categories_use <- default_categories

chem_levels  <- c("ONT-dRNA","ONT-cDNA","ONT-CapTrap","PacBio","PacBio-CapTrap")
cell_levels  <- c("H1-hES","H1-DE","WTC11")
annot_levels <- c("BAMBU","FLAIR","IsoQuant","LoRTIA","NAGATA")

# ====== FELSŐ yield plot chemistry jelölése ======
chem_label_mode_main    <- "legend" # "none" | "legend" | "ticks"
chem_label_mode_cellrows <- "legend" # "none" | "legend" | "ticks"

# ====== Színek ======
chem_cols <- c(
  "ONT-dRNA"       = "#4C78A8",
  "ONT-cDNA"       = "#F58518",
  "ONT-CapTrap"    = "#E45756",
  "PacBio"         = "#72B7B2",
  "PacBio-CapTrap" = "#54A24B"
)

cat_cols <- c(
  "FSM"          = "#4C78A8",
  "ISM"          = "#F58518",
  "NIC"          = "#54A24B",
  "NNC"          = "#E45756",
  "Antisense"    = "#72B7B2",
  "Fusion"       = "#B279A2",
  "GenicGenomic" = "#9D755D",
  "GenicIntron"  = "#8C8C8C",
  "Intergenic"   = "#FF9DA6"
)

repeat_y_axis_each_panel <- TRUE
yield_free_y <- FALSE
yield_y_expand_top <- 0.04

save_plots <- TRUE
pdf_width <- 18
pdf_height <- 7
extra_height_for_cellrows <- 1.2

# -------- GLOBÁLIS FEJLÉC --------
header_text <- "Isoform yield and SQANTI3 class composition across chemistries and annotators"
header_size_pt <- 14
header_bold <- TRUE
header_fill <- "grey90"
header_border <- "grey35"
header_border_lwd <- 0.6
header_rel_height <- 0.09

# ==========================================================
# patchwork + ggh4x
# ==========================================================
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!has_patchwork) {
  install.packages("patchwork")
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
}
if (!has_patchwork) stop("Nem sikerült telepíteni a 'patchwork' csomagot.")
library(patchwork)

has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (repeat_y_axis_each_panel && !has_ggh4x) {
  install.packages("ggh4x")
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}
if (repeat_y_axis_each_panel && !has_ggh4x) {
  warning("A 'ggh4x' nem elérhető, ezért a y-tengely csak a bal szélső panelen fog látszani.")
}

# ==========================================================
# THEME + HEADER
# ==========================================================
theme_Astyle <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width  = unit(0.70, "cm"),
      legend.spacing.x  = unit(0.25, "cm"),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      panel.spacing = unit(0.7, "lines"),
      plot.margin = margin(t = 6, r = 10, b = 10, l = 10)
    )
}

add_global_header <- function(p) {
  header <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
              fill = header_fill, color = header_border, linewidth = header_border_lwd) +
    annotate(
      "text", x = 0.5, y = 0.5, label = header_text,
      fontface = if (header_bold) "bold" else "plain",
      size = header_size_pt / ggplot2::.pt
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  header / p + patchwork::plot_layout(heights = c(header_rel_height, 1))
}

save_plot <- function(filename, p, w, h) {
  fn <- file.path(outdir, filename)
  ggsave(fn, p, width = w, height = h)
  message("Mentve: ", fn)
}

# ==========================================================
# 1) BEOLVASÁS + NORMALIZÁLÁS
# ==========================================================
df0 <- read_tsv(input_file, show_col_types = FALSE)
need <- c("chemistry","cell","category","annotator","n_isoforms")
miss <- setdiff(need, names(df0))
if (length(miss) > 0) stop("Hiányzó oszlop(ok): ", paste(miss, collapse = ", "))

canon_chem <- function(x) str_replace_all(x, regex("captrap", ignore_case = TRUE), "CapTrap")
canon_annot <- function(x) {
  xlow <- tolower(x)
  dplyr::case_when(
    xlow == "bambu" ~ "BAMBU",
    xlow == "flair" ~ "FLAIR",
    xlow == "isoquant" ~ "IsoQuant",
    xlow == "lortia" ~ "LoRTIA",
    xlow == "nagata" ~ "NAGATA",
    TRUE ~ x
  )
}

df <- df0 %>%
  transmute(
    chemistry = canon_chem(as.character(chemistry)),
    cell = as.character(cell),
    category = as.character(category),
    annotator = canon_annot(as.character(annotator)),
    n_isoforms = as.numeric(n_isoforms)
  ) %>%
  filter(category %in% categories_use) %>%
  mutate(
    chemistry = factor(chemistry, levels = chem_levels),
    cell      = factor(cell, levels = cell_levels),
    annotator = factor(annotator, levels = annot_levels),
    category  = factor(category, levels = categories_use)
  )

if (save_plots) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 2) TOP: isoform yield
# ==========================================================
yield_cell <- df %>%
  group_by(annotator, chemistry, cell) %>%
  summarise(n_isoforms_total = sum(n_isoforms, na.rm = TRUE), .groups = "drop")

yield_avg <- yield_cell %>%
  group_by(annotator, chemistry) %>%
  summarise(n_isoforms_avg = mean(n_isoforms_total, na.rm = TRUE), .groups = "drop") %>%
  mutate(cell = factor("AVG", levels = c("AVG")))

# ==========================================================
# 3) BOTTOM: SQANTI3 composition
# ==========================================================
comp_cell <- df %>%
  group_by(annotator, chemistry, cell, category) %>%
  summarise(n_isoforms = sum(n_isoforms, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(annotator, chemistry, cell, category, fill = list(n_isoforms = 0))

collapse_cells_comp <- function(d, method = c("weighted","simple")) {
  method <- match.arg(method)
  if (method == "weighted") {
    d %>%
      group_by(annotator, chemistry, category) %>%
      summarise(n_isoforms = sum(n_isoforms, na.rm = TRUE), .groups = "drop") %>%
      mutate(cell = factor("AVG", levels = c("AVG")))
  } else {
    prop_cell <- d %>%
      group_by(annotator, chemistry, cell) %>%
      mutate(total = sum(n_isoforms, na.rm = TRUE),
             prop = ifelse(total > 0, n_isoforms / total, NA_real_)) %>%
      ungroup()
    
    prop_cell %>%
      group_by(annotator, chemistry, category) %>%
      summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
      group_by(annotator, chemistry) %>%
      mutate(prop = prop / sum(prop, na.rm = TRUE)) %>%
      ungroup() %>%
      transmute(
        annotator, chemistry, category,
        n_isoforms = prop,
        cell = factor("AVG", levels = c("AVG"))
      )
  }
}

comp_avg <- collapse_cells_comp(comp_cell, method = collapse_cells_method)

# ==========================================================
# 4) Facet helper
# ==========================================================
facet_wrap_annot <- function(p, free_y = FALSE) {
  if (repeat_y_axis_each_panel && has_ggh4x) {
    p + ggh4x::facet_wrap2(
      ~ annotator, nrow = 1,
      axes = "all", remove_labels = "none",
      scales = if (free_y) "free_y" else "fixed"
    )
  } else {
    p + facet_wrap(~ annotator, nrow = 1, scales = if (free_y) "free_y" else "fixed")
  }
}

facet_grid_cell_annot <- function(p, free_y = FALSE) {
  if (repeat_y_axis_each_panel && has_ggh4x) {
    p + ggh4x::facet_grid2(
      cell ~ annotator,
      axes = "all", remove_labels = "none",
      scales = if (free_y) "free_y" else "fixed"
    )
  } else {
    p + facet_grid(cell ~ annotator, scales = if (free_y) "free_y" else "fixed")
  }
}

# ==========================================================
# 5) TOP plot (yield) – chemistry legend a felső blokk ALÁ
# ==========================================================
plot_yield <- function(d,
                       facet_mode = c("wrap","cellrows"),
                       chem_label_mode = c("none","legend","ticks"),
                       ylab = "Average number of isoforms / sample") {
  facet_mode <- match.arg(facet_mode)
  chem_label_mode <- match.arg(chem_label_mode)
  
  if ("n_isoforms_avg" %in% names(d)) {
    d2 <- d %>% mutate(y = n_isoforms_avg)
  } else if ("n_isoforms_total" %in% names(d)) {
    d2 <- d %>% mutate(y = n_isoforms_total)
  } else if ("n_isoforms" %in% names(d)) {
    d2 <- d %>% mutate(y = n_isoforms)
  } else {
    stop("plot_yield(): Nem találok n_isoforms_avg / n_isoforms_total / n_isoforms oszlopot.")
  }
  
  p <- ggplot(d2, aes(x = chemistry, y = y, fill = chemistry)) +
    geom_col(width = 0.82, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = chem_cols, breaks = chem_levels, limits = chem_levels, name = "Chemistry") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, yield_y_expand_top))) +
    labs(x = NULL, y = ylab) +
    theme_Astyle() +
    theme(axis.title.x = element_blank())
  
  if (chem_label_mode == "none") {
    p <- p +
      guides(fill = "none") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none")
  } else if (chem_label_mode == "legend") {
    p <- p +
      guides(fill = guide_legend(title = "Chemistry", nrow = 1, byrow = TRUE)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom")
  } else if (chem_label_mode == "ticks") {
    p <- p +
      guides(fill = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
            legend.position = "none")
  }
  
  if (facet_mode == "wrap") {
    p <- facet_wrap_annot(p, free_y = yield_free_y)
  } else {
    p <- facet_grid_cell_annot(p, free_y = yield_free_y)
  }
  
  p
}

# ==========================================================
# 6) BOTTOM plot (SQANTI3 stacked) – legend marad legalul
# ==========================================================
plot_sqanti <- function(d, facet_mode = c("wrap","cellrows")) {
  facet_mode <- match.arg(facet_mode)
  
  p <- ggplot(d, aes(x = chemistry, y = n_isoforms, fill = category)) +
    geom_col(width = 0.82, color = "white", linewidth = 0.25, position = position_fill()) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = c(0, 0.5, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = cat_cols,
      breaks = categories_use,
      limits = categories_use,
      name = "SQANTI3 categories"
    ) +
    guides(fill = guide_legend(title = "SQANTI3 categories", nrow = 1, byrow = TRUE)) +
    labs(x = "Chemistry", y = "SQANTI3 ratio (%)") +
    theme_Astyle() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )
  
  if (facet_mode == "wrap") {
    p <- facet_wrap_annot(p, free_y = FALSE)
  } else {
    p <- facet_grid_cell_annot(p, free_y = FALSE)
  }
  
  p
}

# ==========================================================
# 7) A-blokk összerakása – NINCS guides='collect'
#    chemistry legend a felső blokk ALÁ,
#    SQANTI legend a teljes ábra ALJÁRA kerül
# ==========================================================
build_Ablock <- function(yield_plot, sqanti_plot) {
  core <- (yield_plot / sqanti_plot) + patchwork::plot_layout(heights = c(0.46, 0.54))
  add_global_header(core)
}

# ==========================================================
# 8) FUTTATÁS + MENTÉS (MAIN és DETAILED külön objektum)
# ==========================================================
if (save_plots) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Ablock_MAIN <- NULL
Ablock_DETAILED <- NULL

if ("main_avg_cells" %in% plot_versions) {
  p_top_main <- plot_yield(
    yield_avg,
    facet_mode = "wrap",
    chem_label_mode = chem_label_mode_main,
    ylab = "Average number of isoforms / sample"
  )
  p_bot_main <- plot_sqanti(comp_avg, facet_mode = "wrap")
  pA_main <- build_Ablock(p_top_main, p_bot_main)
  
  Ablock_MAIN <- pA_main
  
  if (save_plots) {
    save_plot(paste0("Ablock_MAIN_cellsAVG_", collapse_cells_method, ".pdf"),
              pA_main, pdf_width, pdf_height)
  } else {
    print(pA_main)
  }
}

if ("detailed_cellrows" %in% plot_versions) {
  p_top_det <- plot_yield(
    yield_cell,
    facet_mode = "cellrows",
    chem_label_mode = chem_label_mode_cellrows,
    ylab = "Number of isoforms / sample"
  )
  p_bot_det <- plot_sqanti(comp_cell, facet_mode = "cellrows")
  pA_det <- build_Ablock(p_top_det, p_bot_det)
  
  Ablock_DETAILED <- pA_det
  
  n_cells <- length(unique(as.character(comp_cell$cell)))
  h_cellrows <- pdf_height * n_cells + extra_height_for_cellrows
  
  if (save_plots) {
    save_plot("Ablock_DETAILED_cells_3rows.pdf", pA_det, pdf_width, h_cellrows)
  } else {
    print(pA_det)
  }
}

# kompatibilitás: pA maradjon a MAIN (ha létezik), különben a DETAILED
if (!is.null(Ablock_MAIN)) {
  pA <- Ablock_MAIN
} else if (!is.null(Ablock_DETAILED)) {
  pA <- Ablock_DETAILED
}
