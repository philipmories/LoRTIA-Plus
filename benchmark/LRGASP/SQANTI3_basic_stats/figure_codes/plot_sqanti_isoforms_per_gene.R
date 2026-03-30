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
input_file <- "isoformsgene.tsv"

plot_versions <- c("collapsed_cells", "cellrows")
avg_cells_method <- "weighted"     # "weighted" vagy "simple"

# --- STACK + LABEL ERGONÓMIA ---
reverse_stack <- FALSE             # FALSE: 1 alul, >=6 felül
label_threshold <- 0.06
label_size <- 2.8
label_fontface <- "bold"
check_overlap_labels <- FALSE

# Mean isoforms/gene felirat felül
show_mean_labels <- FALSE          # <<< kapcsold OFF-ra, ha nem kell
mean_digits <- 2


# --- Y skála headroom automatikusan ---
y_top <- if (show_mean_labels) 1.10 else 1.00  # <<< ez szedi ki a "fehér részt"
y_breaks <- c(0, 0.25, 0.5, 0.75, 1)

# Színek
col_bins <- c(
  "1"    = "#8DA0CB",
  "2-3"  = "#FC8D62",
  "4-5"  = "#66C2A5",
  ">=6"  = "#E78AC3"
)

x_text_angle <- 45
repeat_y_axis_each_panel <- TRUE

# Mentés
save_plots <- TRUE
outdir <- "IsoformsPerGene_Astyle_header"
pdf_width_single  <- 18
pdf_height_single <- 6

match_panel_size_in_cellrows <- TRUE
extra_height_for_cellrows <- 1.2

# -------- GLOBÁLIS FEJLÉC --------
header_text <- "Number of Isoforms per Gene"
header_size_pt <- 14
header_bold <- TRUE
header_fill <- "grey90"
header_border <- "grey35"
header_border_lwd <- 0.6
header_rel_height <- 0.09

# ==========================================================
# patchwork
# ==========================================================
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!has_patchwork) {
  install.packages("patchwork")
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
}
if (!has_patchwork) stop("Nem sikerült telepíteni a 'patchwork' csomagot.")
library(patchwork)

# ==========================================================
# ggh4x (axes='all')
# ==========================================================
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (repeat_y_axis_each_panel && !has_ggh4x) {
  install.packages("ggh4x")
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}
if (repeat_y_axis_each_panel && !has_ggh4x) {
  warning("A 'ggh4x' nem elérhető, ezért a y-tengely csak a bal szélső panelen fog látszani.")
}

# ==========================================================
# SEGÉDEK
# ==========================================================
canon_annot <- function(x) {
  xlow <- tolower(x)
  dplyr::case_when(
    xlow == "bambu"    ~ "BAMBU",
    xlow == "flair"    ~ "FLAIR",
    xlow == "isoquant" ~ "IsoQuant",
    xlow == "lortia"   ~ "LoRTIA",
    xlow == "nagata"   ~ "NAGATA",
    TRUE ~ x
  )
}
canon_chem <- function(x) str_replace_all(x, regex("captrap", ignore_case = TRUE), "CapTrap")

theme_Astyle <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
      axis.text.x  = element_text(angle = x_text_angle, hjust = 1, vjust = 1),
      
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold"),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.75, "cm"),
      
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
# 1) BEOLVASÁS
# ==========================================================
df0 <- read_tsv(input_file, show_col_types = FALSE)

need <- c("Annotator","Chemical","Cell","n_genes_total","isoform_bin","n_genes_in_bin",
          "frac_genes_in_bin","mean_isoforms_per_gene","median_isoforms_per_gene")
miss <- setdiff(need, names(df0))
if (length(miss) > 0) stop("Hiányzó oszlop(ok): ", paste(miss, collapse = ", "))

chem_levels  <- c("ONT-dRNA","ONT-cDNA","ONT-CapTrap","PacBio","PacBio-CapTrap")
cell_levels  <- c("H1-hES","H1-DE","WTC11")
annot_levels <- c("BAMBU","FLAIR","IsoQuant","LoRTIA","NAGATA")
bin_levels   <- c("1","2-3","4-5",">=6")

df <- df0 %>%
  transmute(
    annotator = canon_annot(as.character(Annotator)),
    chemistry = canon_chem(as.character(Chemical)),
    cell      = as.character(Cell),
    n_genes_total = as.numeric(n_genes_total),
    isoform_bin   = as.character(isoform_bin),
    n_genes_in_bin = as.numeric(n_genes_in_bin),
    frac_genes_in_bin = as.numeric(frac_genes_in_bin),
    mean_isoforms_per_gene = as.numeric(mean_isoforms_per_gene),
    median_isoforms_per_gene = as.numeric(median_isoforms_per_gene)
  ) %>%
  mutate(
    chemistry   = factor(chemistry, levels = chem_levels),
    cell        = factor(cell, levels = cell_levels),
    annotator   = factor(annotator, levels = annot_levels),
    isoform_bin = factor(isoform_bin, levels = bin_levels)
  )

# ==========================================================
# 2) Stabilizálás: complete bin + újraszámolt arány
# ==========================================================
d_unit <- df %>%
  group_by(annotator, chemistry, cell) %>%
  summarise(
    n_genes_total = max(n_genes_total, na.rm = TRUE),
    mean_isoforms_per_gene = max(mean_isoforms_per_gene, na.rm = TRUE),
    median_isoforms_per_gene = max(median_isoforms_per_gene, na.rm = TRUE),
    .groups = "drop"
  )

d_bins <- df %>%
  select(annotator, chemistry, cell, isoform_bin, n_genes_in_bin) %>%
  tidyr::complete(annotator, chemistry, cell, isoform_bin, fill = list(n_genes_in_bin = 0)) %>%
  left_join(d_unit, by = c("annotator","chemistry","cell")) %>%
  mutate(
    frac_genes_in_bin = ifelse(n_genes_total > 0, n_genes_in_bin / n_genes_total, NA_real_)
  )

# ==========================================================
# 3) Cell összevonás (AVG)
# ==========================================================
collapse_cells <- function(d, method = c("weighted","simple")) {
  method <- match.arg(method)
  
  if (method == "weighted") {
    out <- d %>%
      group_by(annotator, chemistry, isoform_bin) %>%
      summarise(
        n_genes_in_bin = sum(n_genes_in_bin, na.rm = TRUE),
        n_genes_total  = sum(n_genes_total,  na.rm = TRUE),
        mean_isoforms_per_gene =
          sum(mean_isoforms_per_gene * n_genes_total, na.rm = TRUE) / sum(n_genes_total, na.rm = TRUE),
        median_isoforms_per_gene = median(median_isoforms_per_gene, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        frac_genes_in_bin = ifelse(n_genes_total > 0, n_genes_in_bin / n_genes_total, NA_real_),
        cell = factor("AVG", levels = c("AVG"))
      )
  } else {
    out <- d %>%
      group_by(cell, annotator, chemistry, isoform_bin) %>%
      summarise(prop_cell = mean(frac_genes_in_bin, na.rm = TRUE), .groups = "drop") %>%
      group_by(annotator, chemistry, isoform_bin) %>%
      summarise(frac_genes_in_bin = mean(prop_cell, na.rm = TRUE), .groups = "drop") %>%
      group_by(annotator, chemistry) %>%
      mutate(frac_genes_in_bin = frac_genes_in_bin / sum(frac_genes_in_bin, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(
        n_genes_in_bin = frac_genes_in_bin,
        n_genes_total  = 1,
        mean_isoforms_per_gene = NA_real_,
        median_isoforms_per_gene = NA_real_,
        cell = factor("AVG", levels = c("AVG"))
      )
  }
  out
}

# ==========================================================
# 4) Plot
# ==========================================================
plot_isoforms_per_gene <- function(d, facet_mode = c("wrap","cellrows")) {
  facet_mode <- match.arg(facet_mode)
  
  d_lab <- d %>%
    mutate(
      lab = ifelse(!is.na(frac_genes_in_bin) & frac_genes_in_bin >= label_threshold,
                   percent(frac_genes_in_bin, accuracy = 1), NA_character_)
    )
  
  p <- ggplot(d, aes(x = chemistry, y = frac_genes_in_bin, fill = isoform_bin)) +
    geom_col(
      width = 0.82, color = "white", linewidth = 0.25,
      position = position_stack(reverse = reverse_stack)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = y_breaks,
      limits = c(0, y_top),                 # <<< itt nincs több fölös headroom
      expand = expansion(mult = c(0, 0.00)) # <<< és itt sincs extra "felső levegő"
    ) +
    scale_fill_manual(
      values = col_bins,
      breaks = bin_levels,
      limits = bin_levels,
      name = "Number of isoforms"
    ) +
    labs(x = "Chemistry", y = "Proportion of genes") +
    theme_Astyle()
  
  if (facet_mode == "wrap") {
    if (repeat_y_axis_each_panel && has_ggh4x) {
      p <- p + ggh4x::facet_wrap2(~ annotator, nrow = 1, axes = "all", remove_labels = "none")
    } else {
      p <- p + facet_wrap(~ annotator, nrow = 1)
    }
  } else {
    if (repeat_y_axis_each_panel && has_ggh4x) {
      p <- p + ggh4x::facet_grid2(cell ~ annotator, axes = "all", remove_labels = "none")
    } else {
      p <- p + facet_grid(cell ~ annotator)
    }
  }
  
  p <- p +
    geom_text(
      data = d_lab,
      aes(label = lab),
      position = position_stack(vjust = 0.5, reverse = reverse_stack),
      size = label_size,
      fontface = label_fontface,
      color = "black",
      lineheight = 0.9,
      check_overlap = check_overlap_labels,
      na.rm = TRUE
    )
  
  # mean isoforms/gene felirat felül (csak ha ON)
  if (show_mean_labels) {
    meanlab <- d %>%
      distinct(cell, annotator, chemistry, mean_isoforms_per_gene) %>%
      mutate(mean_txt = ifelse(is.na(mean_isoforms_per_gene), NA_character_,
                               sprintf(paste0("%.", mean_digits, "f"), mean_isoforms_per_gene)))
    
    p <- p +
      geom_text(
        data = meanlab,
        aes(x = chemistry, y = 1.03, label = mean_txt),
        inherit.aes = FALSE,
        size = 3,
        vjust = 0
      ) +
      coord_cartesian(clip = "off")
  }
  
  add_global_header(p)
}

# ==========================================================
# 5) FUTTATÁS + MENTÉS  (+ EXPORT OBJEKTUMOK: pD)
# ==========================================================
if (save_plots) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pD_main   <- NULL   # collapsed_cells (AVG) -> master D panel
pD_detail <- NULL   # cellrows (3 rows)     -> detailed PDF

if ("collapsed_cells" %in% plot_versions) {
  d_avg <- collapse_cells(d_bins, method = avg_cells_method)
  pD_main <- plot_isoforms_per_gene(d_avg, facet_mode = "wrap")
  
  if (save_plots) {
    save_plot(
      paste0("IsoformsPerGene_cellsAVG_", avg_cells_method, "_header.pdf"),
      pD_main, pdf_width_single, pdf_height_single
    )
  } else {
    print(pD_main)
  }
}

if ("cellrows" %in% plot_versions) {
  pD_detail <- plot_isoforms_per_gene(d_bins, facet_mode = "cellrows")
  
  n_cells <- length(unique(as.character(d_bins$cell)))
  h_cellrows <- if (match_panel_size_in_cellrows) {
    pdf_height_single * n_cells + extra_height_for_cellrows
  } else {
    max(pdf_height_single, 7)
  }
  
  if (save_plots) {
    save_plot(
      "IsoformsPerGene_cells_3rows_header.pdf",
      pD_detail, pdf_width_single, h_cellrows
    )
  } else {
    print(pD_detail)
  }
}

# --- EXPORT a masterhez ---
# Alapból a MAIN (AVG) menjen pD-ként:
pD <- pD_main

# Opcionális, név szerint is:
IsoformsPerGene_MAIN   <- pD_main
IsoformsPerGene_DETAIL <- pD_detail

