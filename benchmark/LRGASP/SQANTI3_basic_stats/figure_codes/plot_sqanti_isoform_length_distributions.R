library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(scales)
library(grid)

# ==========================================================
# BEÁLLÍTÁSOK
# ==========================================================
input_file <- "sqanti_length_boxplot_stats.tsv"
outdir     <- "Bblock_IsoLen_Boxplots"

cap <- 20000
categories <- "ALL"  # pl: c("FSM","ISM","NIC","NNC")

filter_cells   <- NULL
filter_chems   <- NULL
filter_annots  <- NULL
filter_file_re <- NULL

# Melyik nézeteket mentsük:
make_plot_cell_x_chem <- TRUE
make_plot_cell_avg    <- TRUE
make_plot_cells_only  <- TRUE

cells_only_side_by_side <- TRUE

# Y skála
y_scale_mode <- "panel_max"   # "panel_max" vagy "fixed_cap"
repeat_y_axis_each_panel <- TRUE

# Mentés
save_plots <- TRUE

# Vászon / panel arány
base_width_in <- 18
panel_aspect  <- 0.78
extra_w_in    <- 1.0
extra_h_in    <- 1.7
min_height_in <- 4.6
max_height_in <- 20

# Cell×chem: csak MAGASSÁGOT növelünk, hogy ne legyen összenyomva
keep_width_in_cellxchem <- TRUE

# -------- A-STÍLUSÚ GLOBÁLIS FEJLÉC (patchwork-sáv) --------
header_text <- "Transcript length distributions of toolkit"
header_size_pt <- 14
header_bold <- TRUE
header_fill <- "grey90"
header_border <- "grey35"
header_border_lwd <- 0.6
header_rel_height <- 0.09

# Faktor szintek (a TE fájlod alapján ez a jó!)
cell_levels  <- c("H1-hES","H1-DE","WTC11")
chem_levels  <- c("ONT-dRNA","ONT-cDNA","ONT-CapTrap","PacBio","PacBio-CapTrap")
annot_levels <- c("LoRTIA","Flair","IsoQuant","NAGATA","bambu")

# ==========================================================
# MULTI (FSM/ISM/NIC/NNC egy oldalon)
# ==========================================================
make_multi_page <- TRUE
multi_categories <- c("FSM","ISM","NIC","NNC")

# "byChem"  -> structural_category ~ chemistry (4x5 panel)
# "chemALL" -> chemistry összevonva, facet_wrap kategóriánként (szebb/nagyobb panelek)
multi_mode <- "chemALL"

multi_pdf_width <- 18  # szélesség
# magasság automatikus lesz (lásd lent)

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
# ggh4x (axes='all' + independent y)
# ==========================================================
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (repeat_y_axis_each_panel && !has_ggh4x) {
  install.packages("ggh4x")
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}
if (repeat_y_axis_each_panel && !has_ggh4x) {
  warning("A 'ggh4x' nem elérhető -> y-tengely nem lesz minden panelen.")
}

# ==========================================================
# SEGÉDEK
# ==========================================================
safe_filename <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

canon_annotator <- function(x) {
  xlow <- str_to_lower(x)
  dplyr::case_when(
    xlow == "lortia"   ~ "LoRTIA",
    xlow == "flair"    ~ "Flair",
    xlow == "isoquant" ~ "IsoQuant",
    xlow == "nagata"   ~ "NAGATA",
    xlow == "bambu"    ~ "bambu",
    TRUE ~ x
  )
}

# >>> FONTOS: a TE fájlodban H1 / H1-endo / WTC11 van <<<
parse_meta_from_file <- function(file_vec) {
  f <- str_replace_all(file_vec, regex("captrap", ignore_case = TRUE), "CapTrap")
  cell      <- str_extract(f, "^(H1-hES|H1-DE|WTC11)")
  chemistry <- str_extract(f, "(ONT-dRNA|ONT-cDNA|ONT-CapTrap|PacBio-CapTrap|PacBio)")
  annotator <- str_extract(f, regex("(LoRTIA|Flair|IsoQuant|NAGATA|bambu)", ignore_case = TRUE))
  annotator <- canon_annotator(annotator)
  tibble(cell = cell, chemistry = chemistry, annotator = annotator, file_norm = f)
}

make_file_key_rm_cell <- function(file_norm, cell) {
  str_replace(file_norm, paste0("^", cell, "[._-]*"), "")
}

make_file_key_rm_cell_chem <- function(file_norm, cell, chemistry) {
  x <- make_file_key_rm_cell(file_norm, cell)
  x <- str_replace_all(x, fixed(chemistry), "")
  x <- str_replace_all(x, "___+", "_")
  x <- str_replace_all(x, "__+", "_")
  x <- str_replace_all(x, "^[_\\.-]+|[_\\.-]+$", "")
  x
}

apply_cap_cols <- function(d, cap) {
  d %>%
    mutate(
      min_p    = pmin(min, cap),
      q1_p     = pmin(q1, cap),
      median_p = pmin(median, cap),
      q3_p     = pmin(q3, cap),
      max_p    = pmin(max, cap),
      max_over_cap = max > cap
    )
}

theme_Astyle <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold"),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width  = unit(0.70, "cm"),
      legend.spacing.x  = unit(0.25, "cm"),
      
      panel.spacing = unit(0.9, "lines"),
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

calc_height_from_panels <- function(width_in, ncol, nrow) {
  panel_w <- max(0.1, (width_in - extra_w_in) / ncol)
  panel_h <- panel_w * panel_aspect
  h_plot_only <- nrow * panel_h + extra_h_in
  h_plot_only <- max(min_height_in, min(max_height_in, h_plot_only))
  h_plot_only * (1 + header_rel_height)
}

# ==========================================================
# 1) BEOLVASÁS + PARSE
# ==========================================================
df0 <- read_tsv(input_file, show_col_types = FALSE)

req <- c("file","structural_category","min","q1","median","q3","max")
miss <- setdiff(req, names(df0))
if (length(miss) > 0) stop("Hiányzó oszlop(ok): ", paste(miss, collapse = ", "))

# tisztítás: structural_category trim + egységesítés (FSM/ISM/NIC/NNC stb.)
df0 <- df0 %>%
  mutate(structural_category = toupper(str_trim(as.character(structural_category))))

meta <- parse_meta_from_file(df0$file)

df <- df0 %>%
  bind_cols(meta %>% select(cell, chemistry, annotator, file_norm)) %>%
  mutate(
    chemistry = factor(chemistry, levels = chem_levels),
    annotator = factor(annotator, levels = annot_levels),
    cell      = factor(cell, levels = cell_levels)
  )

if (!is.null(filter_file_re)) {
  df <- df %>% filter(str_detect(file, regex(filter_file_re, ignore_case = TRUE)))
}

bad_n <- sum(is.na(df$cell) | is.na(df$chemistry) | is.na(df$annotator))
if (bad_n > 0) {
  bad_examples <- df %>%
    filter(is.na(cell) | is.na(chemistry) | is.na(annotator)) %>%
    distinct(file) %>% head(10)
  warning(
    "Parse-olhatatlan sorok: ", bad_n,
    "\nPéldák:\n", paste(bad_examples$file, collapse = "\n"),
    "\n>>> Ellenőrizd, hogy a fájlnevek tényleg így kezdődnek: H1 / H1-endo / WTC11, és tartalmazzák a chemistry-t + annotatort."
  )
  df <- df %>% filter(!is.na(cell), !is.na(chemistry), !is.na(annotator))
}

if (!is.null(filter_cells))  df <- df %>% filter(as.character(cell) %in% filter_cells)
if (!is.null(filter_chems))  df <- df %>% filter(as.character(chemistry) %in% filter_chems)
if (!is.null(filter_annots)) df <- df %>% filter(as.character(annotator) %in% filter_annots)

df <- df %>% apply_cap_cols(cap)

if (save_plots) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 2) AGGREGÁLÁSOK
# ==========================================================
aggregate_over_cells <- function(dfi) {
  dfi %>%
    mutate(file_key = make_file_key_rm_cell(file_norm, as.character(cell))) %>%
    group_by(structural_category, chemistry, annotator, file_key) %>%
    summarise(
      min = mean(min, na.rm = TRUE),
      q1 = mean(q1, na.rm = TRUE),
      median = mean(median, na.rm = TRUE),
      q3 = mean(q3, na.rm = TRUE),
      max = mean(max, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cell = factor("AVG", levels = c("AVG")),
      file = paste0("AVG__", file_key),
      file_norm = file,
      chemistry = factor(as.character(chemistry), levels = chem_levels),
      annotator = factor(as.character(annotator), levels = annot_levels)
    ) %>%
    apply_cap_cols(cap)
}

aggregate_over_chems_keep_cells <- function(dfi) {
  dfi %>%
    mutate(file_key2 = make_file_key_rm_cell_chem(file_norm, as.character(cell), as.character(chemistry))) %>%
    group_by(structural_category, cell, annotator, file_key2) %>%
    summarise(
      min = mean(min, na.rm = TRUE),
      q1 = mean(q1, na.rm = TRUE),
      median = mean(median, na.rm = TRUE),
      q3 = mean(q3, na.rm = TRUE),
      max = mean(max, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      chemistry = factor("ALL", levels = c("ALL")),
      file = paste0("CHEMALL__", file_key2),
      file_norm = file,
      cell = factor(as.character(cell), levels = cell_levels),
      annotator = factor(as.character(annotator), levels = annot_levels)
    ) %>%
    apply_cap_cols(cap)
}

# cells + chems összevonása (MULTI "chemALL" módhoz)
aggregate_over_cells_and_chems <- function(dfi) {
  dfi %>%
    mutate(file_key3 = make_file_key_rm_cell_chem(file_norm, as.character(cell), as.character(chemistry))) %>%
    group_by(structural_category, annotator, file_key3) %>%
    summarise(
      min = mean(min, na.rm = TRUE),
      q1 = mean(q1, na.rm = TRUE),
      median = mean(median, na.rm = TRUE),
      q3 = mean(q3, na.rm = TRUE),
      max = mean(max, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cell = factor("AVG", levels = c("AVG")),
      chemistry = factor("ALL", levels = c("ALL")),
      file = paste0("AVG_ALL__", file_key3),
      file_norm = file,
      annotator = factor(as.character(annotator), levels = annot_levels)
    ) %>%
    apply_cap_cols(cap)
}

# ==========================================================
# 3) PLOT (A-stílusú header sávval)
# ==========================================================
plot_one_category <- function(dfi, cap,
                              facet_mode = c("cell_x_chem","cell_avg_x_chem","cells_only"),
                              add_header = TRUE) {
  facet_mode <- match.arg(facet_mode)
  
  pos <- position_dodge2(width = 0.8, preserve = "single")
  dfi <- dfi %>% mutate(group_id = paste(cell, chemistry, annotator, file, sep = "|"))
  
  blank0 <- dfi %>% distinct(cell, chemistry) %>% mutate(annotator = annot_levels[1], y0 = 0)
  
  d_cap <- dfi %>% filter(max_over_cap) %>% mutate(cap_y = cap)
  
  p <- ggplot(dfi, aes(x = annotator, fill = annotator)) +
    geom_boxplot(
      aes(group = group_id, ymin = min_p, lower = q1_p, middle = median_p, upper = q3_p, ymax = max_p),
      stat = "identity", width = 0.65, position = pos
    ) +
    geom_point(
      data = d_cap,
      aes(x = annotator, y = cap_y, group = group_id),
      inherit.aes = FALSE, position = pos
    ) +
    geom_blank(data = blank0, aes(x = annotator, y = y0), inherit.aes = FALSE) +
    labs(
      x = "Annotator",
      y = "Isoform length (bp)",
      fill = "Annotator"
    ) +
    theme_Astyle() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = breaks_pretty(n = 4)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  if (identical(y_scale_mode, "fixed_cap")) {
    p <- p + coord_cartesian(ylim = c(0, cap))
  }
  
  want_independent <- identical(y_scale_mode, "panel_max") && has_ggh4x
  axes_all <- repeat_y_axis_each_panel && has_ggh4x
  
  if (facet_mode %in% c("cell_x_chem","cell_avg_x_chem")) {
    if (has_ggh4x) {
      p <- p + ggh4x::facet_grid2(
        cell ~ chemistry, drop = TRUE,
        axes = if (axes_all) "all" else "margins",
        scales = if (want_independent) "free_y" else "fixed",
        independent = if (want_independent) "y" else "none"
      )
    } else {
      p <- p + facet_grid(cell ~ chemistry, drop = TRUE,
                          scales = if (identical(y_scale_mode,"panel_max")) "free_y" else "fixed")
    }
  }
  
  if (facet_mode == "cells_only") {
    facet_formula <- if (cells_only_side_by_side) (. ~ cell) else (cell ~ .)
    if (has_ggh4x) {
      p <- p + ggh4x::facet_grid2(
        facet_formula, drop = TRUE,
        axes = if (axes_all) "all" else "margins",
        scales = if (want_independent) "free_y" else "fixed",
        independent = if (want_independent) "y" else "none"
      )
    } else {
      p <- p + facet_grid(facet_formula, drop = TRUE,
                          scales = if (identical(y_scale_mode,"panel_max")) "free_y" else "fixed")
    }
  }
  
  if (add_header) return(add_global_header(p))
  p
}

# MULTI plot: kategóriák egy oldalon
plot_multi_byChem <- function(d_avg_multi, cap, cats, add_header = TRUE) {
  pos <- position_dodge2(width = 0.8, preserve = "single")
  
  d <- d_avg_multi %>%
    mutate(
      structural_category = factor(structural_category, levels = cats),
      group_id = paste(structural_category, chemistry, annotator, file, sep = "|")
    )
  
  d_cap <- d %>% filter(max_over_cap) %>% mutate(cap_y = cap)
  
  p <- ggplot(d, aes(x = annotator, fill = annotator)) +
    geom_boxplot(
      aes(group = group_id, ymin = min_p, lower = q1_p, middle = median_p, upper = q3_p, ymax = max_p),
      stat = "identity", width = 0.65, position = pos
    ) +
    geom_point(
      data = d_cap,
      aes(x = annotator, y = cap_y, group = group_id),
      inherit.aes = FALSE, position = pos
    ) +
    labs(x = "Annotator", y = "Isoform length (bp)", fill = "Annotator") +
    theme_Astyle() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = breaks_pretty(n = 3)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  if (identical(y_scale_mode, "fixed_cap")) {
    p <- p + coord_cartesian(ylim = c(0, cap))
  }
  
  want_independent <- identical(y_scale_mode, "panel_max") && has_ggh4x
  axes_all <- repeat_y_axis_each_panel && has_ggh4x
  
  if (has_ggh4x) {
    p <- p + ggh4x::facet_grid2(
      structural_category ~ chemistry,
      axes = if (axes_all) "all" else "margins",
      scales = if (want_independent) "free_y" else "fixed",
      independent = if (want_independent) "y" else "none"
    )
  } else {
    p <- p + facet_grid(structural_category ~ chemistry,
                        scales = if (identical(y_scale_mode,"panel_max")) "free_y" else "fixed")
  }
  
  if (add_header) add_global_header(p) else p
}

plot_multi_chemALL <- function(d_all, cap, cats, add_header = TRUE) {
  pos <- position_dodge2(width = 0.8, preserve = "single")
  
  d <- d_all %>%
    mutate(
      structural_category = factor(structural_category, levels = cats),
      group_id = paste(structural_category, annotator, file, sep = "|")
    )
  
  d_cap <- d %>% filter(max_over_cap) %>% mutate(cap_y = cap)
  
  p <- ggplot(d, aes(x = annotator, fill = annotator)) +
    geom_boxplot(
      aes(group = group_id, ymin = min_p, lower = q1_p, middle = median_p, upper = q3_p, ymax = max_p),
      stat = "identity", width = 0.65, position = pos
    ) +
    geom_point(
      data = d_cap,
      aes(x = annotator, y = cap_y, group = group_id),
      inherit.aes = FALSE, position = pos
    ) +
    labs(x = "Annotator", y = "Isoform length (bp)", fill = "Annotator") +
    theme_Astyle() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = breaks_pretty(n = 4)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    facet_wrap(~ structural_category, ncol = 2)
  
  if (identical(y_scale_mode, "fixed_cap")) {
    p <- p + coord_cartesian(ylim = c(0, cap))
  } else {
    # panel_max esetén itt is hasznos a free_y, de facet_wrap-nál sima ggplot free_y
    p <- p + facet_wrap(~ structural_category, ncol = 2, scales = "free_y")
  }
  
  if (add_header) add_global_header(p) else p
}

# ==========================================================
# 4) FUTTATÁS + MENTÉS
# ==========================================================
all_cats <- sort(unique(df$structural_category))
cats_to_plot <- if (is.null(categories) || identical(categories, "ALL")) all_cats else toupper(categories)
cats_to_plot <- intersect(cats_to_plot, all_cats)
if (length(cats_to_plot) == 0) stop("Nincs ábrázolható kategória (ellenőrizd a structural_category neveket).")

IsoLen_MAIN <- NULL

for (cat in cats_to_plot) {
  dcat <- df %>% filter(structural_category == cat)
  if (nrow(dcat) == 0) next
  
  chem_n <- dcat %>% distinct(chemistry) %>% nrow()
  
  # referencia (AVG: 1 sor)
  height_avg_ref <- calc_height_from_panels(base_width_in, ncol = chem_n, nrow = 1)
  
  # (2) AVG (cell averaged)  >>> EZ A "FŐ" nézet, itt CELL-ÁTLAG van <<<
  if (make_plot_cell_avg) {
    d_avg <- aggregate_over_cells(dcat)
    p_avg <- plot_one_category(d_avg, cap, facet_mode = "cell_avg_x_chem", add_header = TRUE)
    
    fn <- file.path(outdir, paste0("B_isoLen_", safe_filename(cat), "_cellAVG_", y_scale_mode, "_cap", cap, ".pdf"))
    ggsave(fn, plot = p_avg, width = base_width_in, height = height_avg_ref, limitsize = FALSE)
    message("Mentve: ", fn)
    
    if (is.null(IsoLen_MAIN)) {
      IsoLen_MAIN <- p_avg
    }
  }
  
  # (1) cell × chemistry  >>> itt emeljük a MAGASSÁGOT, hogy ne legyen összenyomva <<<
  if (make_plot_cell_x_chem) {
    p_cxc <- plot_one_category(dcat, cap, facet_mode = "cell_x_chem", add_header = TRUE)
    
    cell_n <- dcat %>% distinct(cell) %>% nrow()
    
    if (keep_width_in_cellxchem) {
      width_cxc  <- base_width_in
      height_cxc <- calc_height_from_panels(width_cxc, ncol = chem_n, nrow = cell_n)
    } else {
      # régi jellegű skálázás (ha mégis kéne)
      width_cxc  <- base_width_in
      height_cxc <- calc_height_from_panels(width_cxc, ncol = chem_n, nrow = cell_n)
    }
    
    fn <- file.path(outdir, paste0("B_isoLen_", safe_filename(cat), "_cellXchem_", y_scale_mode, "_cap", cap, ".pdf"))
    ggsave(fn, plot = p_cxc, width = width_cxc, height = height_cxc, limitsize = FALSE)
    message("Mentve: ", fn)
  }
  
  # (3) cellsOnly (chemistry averaged)
  if (make_plot_cells_only) {
    d_cells_only <- aggregate_over_chems_keep_cells(dcat)
    p_cells <- plot_one_category(d_cells_only, cap, facet_mode = "cells_only", add_header = TRUE)
    
    fn <- file.path(outdir, paste0("B_isoLen_", safe_filename(cat), "_cellsOnly_",
                                   if (cells_only_side_by_side) "sideBySide" else "stacked",
                                   "_", y_scale_mode, "_cap", cap, ".pdf"))
    ggsave(fn, plot = p_cells, width = base_width_in, height = height_avg_ref, limitsize = FALSE)
    message("Mentve: ", fn)
  }
}

if (is.null(IsoLen_MAIN)) warning("IsoLen_MAIN nem jött létre (nincs adat / nincs kategória).")

# ==========================================================
# 5) MULTI oldal (FSM/ISM/NIC/NNC egy PDF-re)  <<< EZ NEM LESZ ÜRES >>>
# ==========================================================
if (make_multi_page) {
  cats_req <- toupper(multi_categories)
  cats_ok  <- intersect(cats_req, all_cats)
  
  if (length(cats_ok) == 0) {
    warning("MULTI: egyetlen kért kategória sincs az adatokban. Kért: ",
            paste(cats_req, collapse = ", "),
            " | Elérhető: ", paste(all_cats, collapse = ", "))
  } else {
    
    d_multi <- df %>% filter(structural_category %in% cats_ok)
    
    if (multi_mode == "byChem") {
      # cellát átlagoljuk, chemistry marad
      d_avg_multi <- aggregate_over_cells(d_multi)
      
      # biztosítjuk az order-t
      d_avg_multi <- d_avg_multi %>% mutate(structural_category = factor(structural_category, levels = cats_ok))
      
      chem_n <- d_avg_multi %>% distinct(chemistry) %>% nrow()
      h_multi <- calc_height_from_panels(multi_pdf_width, ncol = chem_n, nrow = length(cats_ok))
      
      p_multi <- plot_multi_byChem(d_avg_multi, cap, cats_ok, add_header = TRUE)
      
      fn <- file.path(outdir, paste0("B_isoLen_MULTI_byChem_", paste(cats_ok, collapse = "_"), ".pdf"))
      ggsave(fn, plot = p_multi, width = multi_pdf_width, height = h_multi, limitsize = FALSE)
      message("Mentve: ", fn, " (", multi_pdf_width, " x ", round(h_multi, 2), " in)")
      
    } else if (multi_mode == "chemALL") {
      # cell + chemistry összevonás (nagy, publikáció-barát panelek)
      d_all <- aggregate_over_cells_and_chems(d_multi) %>%
        mutate(structural_category = factor(structural_category, levels = cats_ok))
      
      # ha mégis 0 sor (nem kéne), álljunk meg magyarázattal
      if (nrow(d_all) == 0) stop("MULTI chemALL: 0 sor lett az aggregálás után. (Valószínű parse / filter gond.)")
      
      # magasság: 2 oszlopos wrap -> sorok száma
      ncol_wrap <- 2
      nrow_wrap <- ceiling(length(cats_ok) / ncol_wrap)
      h_multi <- calc_height_from_panels(multi_pdf_width, ncol = ncol_wrap, nrow = nrow_wrap)
      
      p_multi <- plot_multi_chemALL(d_all, cap, cats_ok, add_header = TRUE)
      
      fn <- file.path(outdir, paste0("B_isoLen_MULTI_chemALL_", paste(cats_ok, collapse = "_"), ".pdf"))
      ggsave(fn, plot = p_multi, width = multi_pdf_width, height = h_multi, limitsize = FALSE)
      message("Mentve: ", fn, " (", multi_pdf_width, " x ", round(h_multi, 2), " in)")
    } else {
      warning("Ismeretlen multi_mode: ", multi_mode, " (használd: 'byChem' vagy 'chemALL').")
    }
  }
}
