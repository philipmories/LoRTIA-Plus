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

default_categories <- c("FSM","ISM","NIC","NNC")
categories_use <- default_categories

plot_versions <- c("by_category", "by_category_cellrows", "collapsed_cells")

# ---- CELL-összevonás módja ----
# "cell_mean" = cellánkénti arányok egyszerű átlaga (EZ KELL NEKED)
# "weighted"  = isoform-súlyozott (összes isoform dominál)
avg_cells_method <- "cell_mean"   # <<<<<<

# % feliratok
label_threshold <- 0.06
label_size <- 3
label_alpha <- 0.75

# Figure 2 C színek
col_coding    <- "#F08080"  # Coding (felül)
col_noncoding <- "#87CEFA"  # Non-coding (alul)

repeat_y_axis_each_panel <- TRUE

# Mentés
save_plots <- TRUE
outdir <- "coding_composition_plots_Astyle_header"
pdf_width_single  <- 18
pdf_height_single <- 5

match_panel_size_in_cellrows <- TRUE
extra_height_for_cellrows <- 1.2

# -------- GLOBÁLIS FEJLÉC --------
header_text <- "Coding/non-coding composition per annotator"
header_size_pt <- 14
header_bold <- TRUE
header_fill <- "grey90"
header_border <- "grey35"
header_border_lwd <- 0.6
header_rel_height <- 0.09

# ---- faktor szintek (LoRTIA elöl) ----
annot_levels <- c("LoRTIA","Flair","IsoQuant","NAGATA","bambu")
cell_levels  <- c("H1-hES","H1-DE","WTC11")

# ==========================================================
# patchwork (fejlécsávhoz)
# ==========================================================
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!has_patchwork) {
  install.packages("patchwork")
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
}
if (!has_patchwork) stop("Nem sikerült telepíteni a 'patchwork' csomagot.")
library(patchwork)

# ==========================================================
# ggh4x (axes="all")
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
    xlow == "lortia"   ~ "LoRTIA",
    xlow == "flair"    ~ "Flair",
    xlow == "isoquant" ~ "IsoQuant",
    xlow == "nagata"   ~ "NAGATA",
    xlow == "bambu"    ~ "bambu",
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
      axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1),
      
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold"),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
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
# 1) BEOLVASÁS + TISZTÍTÁS
# ==========================================================
df0 <- read_tsv(input_file, show_col_types = FALSE)

need <- c("chemistry","cell","category","annotator","coding_rate","n_isoforms")
miss <- setdiff(need, names(df0))
if (length(miss) > 0) stop("Hiányzó oszlop(ok): ", paste(miss, collapse = ", "))

df <- df0 %>%
  transmute(
    chemistry   = canon_chem(as.character(chemistry)),
    cell        = as.character(cell),
    category    = as.character(category),
    annotator   = canon_annot(as.character(annotator)),
    coding_rate = as.numeric(coding_rate),
    n_isoforms  = as.numeric(n_isoforms)
  )

avail_cats <- sort(unique(df$category))
categories_use <- intersect(categories_use, avail_cats)
if (length(categories_use) == 0) stop("Nincs ábrázolható kategória.")
df <- df %>% filter(category %in% categories_use)

df <- df %>%
  mutate(
    category  = factor(category, levels = categories_use),
    annotator = factor(annotator, levels = annot_levels),
    cell      = factor(cell, levels = cell_levels)
  )

if (save_plots) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 2) ELŐKÉSZÍTÉS (SZÁMOLÁSOK)
#    - Non-coding alul
#    - Coding felül
#    - avg_cells_method: "cell_mean" vagy "weighted"
# ==========================================================
prep_composition <- function(df, keep_cell = FALSE, avg_method = c("cell_mean","weighted")) {
  avg_method <- match.arg(avg_method)
  
  # cell×annot×category×chem -> összevonjuk kémián belül, majd cellen belül arányt képzünk
  d_cellchem <- df %>%
    group_by(cell, annotator, category, chemistry) %>%
    summarise(
      n_isoforms_sum = sum(n_isoforms, na.rm = TRUE),
      coding_n_sum   = sum(coding_rate * n_isoforms, na.rm = TRUE),
      .groups = "drop"
    )
  
  # cell×annot×category szint (chem-ek összevonása)
  d_cell <- d_cellchem %>%
    group_by(cell, annotator, category) %>%
    summarise(
      n_isoforms_total = sum(n_isoforms_sum, na.rm = TRUE),
      coding_n_cell    = sum(coding_n_sum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      coding_prop_cell = ifelse(n_isoforms_total > 0, coding_n_cell / n_isoforms_total, NA_real_),
      noncoding_prop_cell = 1 - coding_prop_cell
    )
  
  if (keep_cell) {
    # cellrows: cellenként külön bar-ok (itt nem átlagolunk cellák között)
    d_sum <- d_cell %>%
      transmute(
        cell, annotator, category,
        coding_prop = coding_prop_cell,
        noncoding_prop = noncoding_prop_cell
      )
  } else {
    if (avg_method == "weighted") {
      # isoform-súlyozott (összeadjuk cellákon át a számlálót és nevezőt)
      d_sum <- d_cell %>%
        group_by(annotator, category) %>%
        summarise(
          coding_prop = ifelse(sum(n_isoforms_total, na.rm = TRUE) > 0,
                               sum(coding_n_cell, na.rm = TRUE) / sum(n_isoforms_total, na.rm = TRUE),
                               NA_real_),
          .groups = "drop"
        ) %>%
        mutate(noncoding_prop = 1 - coding_prop)
    } else {
      # cell-átlag: cellánkénti arányok egyszerű átlaga (minden cell = 1 szavazat)
      d_sum <- d_cell %>%
        group_by(annotator, category) %>%
        summarise(
          coding_prop = mean(coding_prop_cell, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(noncoding_prop = 1 - coding_prop)
    }
  }
  
  # hosszú: a rajzoláshoz (count lehet prop is, position_fill úgyis arányt csinál)
  d_long <- d_sum %>%
    pivot_longer(
      cols = c(coding_prop, noncoding_prop),
      names_to = "status", values_to = "prop"
    ) %>%
    mutate(
      status = recode(status, noncoding_prop = "Non-coding", coding_prop = "Coding"),
      status = factor(status, levels = c("Coding","Non-coding")),
      count  = prop,  # “pszeudo count”
      lab    = ifelse(!is.na(prop) & prop >= label_threshold, percent(prop, accuracy = 1), NA_character_)
    )
  
  # CÍMKE POZÍCIÓK: mindig a saját sáv KÖZEPÉN
  # Non-coding: 0 .. noncoding_prop  -> közép = noncoding_prop/2
  # Coding: noncoding_prop .. 1      -> közép = noncoding_prop + coding_prop/2
  d_labs <- d_sum %>%
    mutate(
      lab_noncoding = ifelse(!is.na(noncoding_prop) & noncoding_prop >= label_threshold,
                             percent(noncoding_prop, accuracy = 1), NA_character_),
      lab_coding = ifelse(!is.na(coding_prop) & coding_prop >= label_threshold,
                          percent(coding_prop, accuracy = 1), NA_character_),
      y_noncoding = noncoding_prop / 2,
      y_coding = noncoding_prop + coding_prop / 2
    )
  
  list(sum = d_sum, long = d_long, labs = d_labs)
}

# ==========================================================
# 3) PLOT
# ==========================================================
plot_composition <- function(prep, facet_mode = c("wrap","cellrows")) {
  facet_mode <- match.arg(facet_mode)
  
  d_long <- prep$long
  d_labs <- prep$labs
  
  p <- ggplot(d_long, aes(x = category, y = count, fill = status)) +
    geom_col(
      position = position_fill(reverse = FALSE),
      width = 0.82,
      color = "white",
      linewidth = 0.25
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = c("Non-coding" = col_noncoding, "Coding" = col_coding),
      breaks = c("Non-coding","Coding"),
      limits = c("Non-coding","Coding")
    ) +
    labs(
      x = "SQANTI+ category",
      y = "Proportion of isoforms"
    ) +
    theme_Astyle()
  
  # facet
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
  
  # címkék: saját sáv közepén
  p <- p +
    geom_label(
      data = d_labs,
      aes(x = category, y = y_noncoding, label = lab_noncoding),
      inherit.aes = FALSE,
      size = label_size,
      label.size = 0,
      fill = alpha("white", label_alpha),
      color = "black",
      label.padding = unit(0.08, "lines"),
      na.rm = TRUE
    ) +
    geom_label(
      data = d_labs,
      aes(x = category, y = y_coding, label = lab_coding),
      inherit.aes = FALSE,
      size = label_size,
      label.size = 0,
      fill = alpha("white", label_alpha),
      color = "black",
      label.padding = unit(0.08, "lines"),
      na.rm = TRUE
    )
  
  add_global_header(p)
}

# ==========================================================
# 4) FUTTATÁS + EXPORT OBJEKTUMOK
# ==========================================================
pC_by_category <- NULL
pC_cellrows    <- NULL
pC_avgcells    <- NULL

if ("by_category" %in% plot_versions) {
  prep <- prep_composition(df, keep_cell = FALSE, avg_method = avg_cells_method)
  pC_by_category <- plot_composition(prep, facet_mode = "wrap")
  if (save_plots) save_plot(paste0("Cstyle_by_category_", avg_cells_method, "_header.pdf"),
                            pC_by_category, pdf_width_single, pdf_height_single) else print(pC_by_category)
}

if ("by_category_cellrows" %in% plot_versions) {
  prep <- prep_composition(df, keep_cell = TRUE, avg_method = avg_cells_method)
  pC_cellrows <- plot_composition(prep, facet_mode = "cellrows")
  
  n_cells <- length(unique(as.character(df$cell)))
  h_cellrows <- if (match_panel_size_in_cellrows) pdf_height_single * n_cells + extra_height_for_cellrows else max(pdf_height_single, 7)
  
  if (save_plots) save_plot("Cstyle_by_category_cells_3rows_header.pdf",
                            pC_cellrows, pdf_width_single, h_cellrows) else print(pC_cellrows)
}

if ("collapsed_cells" %in% plot_versions) {
  prep <- prep_composition(df, keep_cell = FALSE, avg_method = avg_cells_method)
  pC_avgcells <- plot_composition(prep, facet_mode = "wrap")
  
  Coding_MAIN <- pC_avgcells  # <<< ezt fogod használni tovább
  
  if (save_plots) save_plot(paste0("Cstyle_cellsAVG_", avg_cells_method, "_header.pdf"),
                            pC_avgcells, pdf_width_single, pdf_height_single) else print(pC_avgcells)
}

# kényelmi alias (ha kell)
pC <- pC_avgcells
