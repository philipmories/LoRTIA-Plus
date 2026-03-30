library(tidyverse)
library(stringr)
library(scales)
library(patchwork)
library(grid)
library(cowplot)

# =========================================================
# KSHV COMBINED FIGURE
# A = metrics panels
# B = transcript category stacked barplot
# C = upset plots from ONE combined TSV (FP.tsv)
# =========================================================

# -------------------------
# File checks
# -------------------------
stop_if_missing <- function(files) {
  miss <- files[!file.exists(files)]
  if (length(miss) > 0) {
    stop(
      "Hiányzó fájl(ok):\n", paste(miss, collapse = "\n"),
      "\n\nAktuális mappa:\n", getwd(),
      "\n\nTSV-k a mappában:\n", paste(list.files(pattern = "tsv$", ignore.case = TRUE), collapse = "\n")
    )
  }
}

validate_metric_file <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, n_max = 5)
  req <- c("TP", "FP", "FN")
  miss <- setdiff(req, names(df))
  
  if (ncol(df) < 4) {
    stop("Túl kevés oszlop a metric fájlban: ", path)
  }
  if (length(miss) > 0) {
    stop("Hiányzó kötelező oszlop(ok) a metric fájlban ", path, ": ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

validate_category_file <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, n_max = 5)
  
  if (!("Category" %in% names(df))) {
    stop("Hiányzó 'Category' oszlop a category fájlban: ", path)
  }
  
  if (!any(c("isoform", "ref_id") %in% names(df))) {
    stop("A category fájlban nincs 'isoform' vagy 'ref_id' oszlop: ", path)
  }
  
  invisible(TRUE)
}

validate_combined_upset_file <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, n_max = 5)
  
  req <- c("Start", "End", "Strand")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop("Hiányzó kötelező oszlop(ok) a combined upset fájlban: ",
         paste(miss, collapse = ", "))
  }
  
  if (!any(c("Original_Transcript_ID", "Cluster_ID") %in% names(df))) {
    stop("Nincs 'Original_Transcript_ID' vagy 'Cluster_ID' oszlop a combined upset fájlban.")
  }
  
  dRNA_cols  <- names(df)[str_detect(names(df), "^dRNA_false_positives__.*_present$")]
  dcDNA_cols <- names(df)[str_detect(names(df), "^dcDNA_false_positives__.*_present$")]
  
  dRNA_cols  <- dRNA_cols[!str_detect(dRNA_cols, "__files_present$")]
  dcDNA_cols <- dcDNA_cols[!str_detect(dcDNA_cols, "__files_present$")]
  
  if (length(dRNA_cols) == 0) {
    stop("Nem találtam dRNA annotátoroszlopokat a combined upset fájlban.")
  }
  if (length(dcDNA_cols) == 0) {
    stop("Nem találtam dcDNA annotátoroszlopokat a combined upset fájlban.")
  }
  
  bad_dRNA <- dRNA_cols[
    !map_lgl(df[dRNA_cols], ~ all(.x %in% c(0, 1, "0", "1", NA)))
  ]
  bad_dcDNA <- dcDNA_cols[
    !map_lgl(df[dcDNA_cols], ~ all(.x %in% c(0, 1, "0", "1", NA)))
  ]
  
  if (length(bad_dRNA) > 0) {
    stop("Nem bináris dRNA annotátoroszlop(ok): ", paste(bad_dRNA, collapse = ", "))
  }
  if (length(bad_dcDNA) > 0) {
    stop("Nem bináris dcDNA annotátoroszlop(ok): ", paste(bad_dcDNA, collapse = ", "))
  }
  
  invisible(TRUE)
}

# -------------------------
# Theme elements
# -------------------------
strip_box <- theme(
  strip.background = element_rect(fill = "grey92", colour = "grey55", linewidth = 0.6),
  strip.text = element_text(
    face = "bold",
    colour = "black",
    size = 9,
    margin = margin(t = 2, b = 2)
  )
)

right_border <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4)
)

right_border_stacked <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25)
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, colour = "black"),
    panel.spacing = unit(0.65, "lines")
  )

make_grey_header_small <- function(text, text_size = 3.0) {
  ggplot() +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "grey92", colour = "grey55", linewidth = 0.35
    ) +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = text,
      fontface = "bold",
      size = text_size,
      colour = "black"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_upset_header <- function(text, text_size = 3.0) {
  ggplot() +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "grey92", colour = "grey55", linewidth = 0.35
    ) +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = text,
      fontface = "bold",
      size = text_size,
      colour = "black"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 4, 0, 4))
}

# =========================================================
# HELPERS
# =========================================================

clean_tool_name <- function(x) {
  x %>%
    str_remove("\\.gff3$") %>%
    str_remove("\\.(TSS|TES|Transcript|Transcripts)$") %>%
    str_replace("^FLAIR_.*$", "FLAIR") %>%
    str_replace("^bambu$", "BAMBU") %>%
    str_replace("^Bambu$", "BAMBU")
}

# =========================================================
# PART A - TOP METRICS PANELS
# =========================================================

read_feature_file <- function(path, feature_name, chemistry_name) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  names(df)[1] <- "tool"
  
  df %>%
    mutate(
      tool = clean_tool_name(tool),
      TP = as.numeric(TP),
      FP = as.numeric(FP),
      FN = as.numeric(FN)
    ) %>%
    transmute(
      chemistry = chemistry_name,
      feature = feature_name,
      tool,
      TP, FP, FN,
      precision = if_else(TP + FP == 0, 0, TP / (TP + FP)),
      recall    = if_else(TP + FN == 0, 0, TP / (TP + FN)),
      F1        = if_else(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    )
}

metric_files <- c(
  "dcDNA_stats_TSS.tsv",
  "dcDNA_stats_TES.tsv",
  "dcDNA_stats.tsv",
  "dRNA_stats_TSS.tsv",
  "dRNA_stats_TES.tsv",
  "dRNA_stats.tsv"
)

stop_if_missing(metric_files)
walk(metric_files, validate_metric_file)

df_all <- bind_rows(
  read_feature_file("dcDNA_stats_TSS.tsv", "TSS",         "dcDNA"),
  read_feature_file("dcDNA_stats_TES.tsv", "TES",         "dcDNA"),
  read_feature_file("dcDNA_stats.tsv",     "Transcripts", "dcDNA"),
  read_feature_file("dRNA_stats_TSS.tsv",  "TSS",         "dRNA"),
  read_feature_file("dRNA_stats_TES.tsv",  "TES",         "dRNA"),
  read_feature_file("dRNA_stats.tsv",      "Transcripts", "dRNA")
)

y_top <- 1.25
label_offset <- 0.04
label_cap <- 1.20

metric_colors <- c(
  "Precision" = "#F8766D",
  "Recall"    = "#00BA38",
  "F1"        = "#619CFF"
)

df_all_long <- df_all %>%
  pivot_longer(
    cols = c(precision, recall, F1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric, precision = "Precision", recall = "Recall", F1 = "F1"),
    metric = factor(metric, levels = c("Precision", "Recall", "F1")),
    feature = factor(feature, levels = c("TSS", "TES", "Transcripts")),
    chemistry = factor(chemistry, levels = c("dcDNA", "dRNA")),
    tool = factor(tool, levels = c("LoRTIA", "BAMBU", "FLAIR", "IsoQuant", "NAGATA")),
    label_y = pmin(value + label_offset, label_cap)
  ) %>%
  filter(!is.na(tool))

make_feature_plot <- function(feature_name, show_y_title = TRUE) {
  df_sub <- df_all_long %>% filter(feature == feature_name)
  dodge <- position_dodge2(width = 0.72, preserve = "single")
  
  ggplot(df_sub, aes(x = tool, y = value, fill = metric)) +
    geom_col(position = dodge, width = 0.56) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dotted",
      colour = "grey82",
      linewidth = 0.22
    ) +
    geom_text(
      aes(
        y = label_y,
        label = sprintf("%.3f", value),
        group = metric
      ),
      position = dodge,
      angle = 90,
      hjust = -0.05,
      vjust = 0.5,
      size = 2.55,
      colour = "black"
    ) +
    scale_fill_manual(values = metric_colors, name = NULL) +
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = "Annotator",
      y = if (show_y_title) "Score" else NULL,
      title = feature_name
    ) +
    facet_wrap(~chemistry, ncol = 1) +
    base_theme +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, colour = "black"),
      legend.position = "bottom"
    ) +
    right_border +
    strip_box
}

p_tss <- make_feature_plot("TSS", TRUE)
p_tes <- make_feature_plot("TES", FALSE)
p_tr  <- make_feature_plot("Transcripts", FALSE)

p_top <- (p_tss | p_tes | p_tr) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

# =========================================================
# PART B - CATEGORY STACKED BARPLOT
# =========================================================

clean_category_annotator_names <- function(x) {
  x %>%
    str_remove("\\.gff3_match$") %>%
    str_remove("_match$") %>%
    str_replace("^FLAIR_.*$", "FLAIR") %>%
    str_replace("^bambu$", "BAMBU") %>%
    str_replace("^Bambu$", "BAMBU")
}

recode_category <- function(x) {
  x0 <- x
  xl <- str_to_lower(x0)
  out <- x0
  
  out <- ifelse(
    str_detect(xl, "intergenic") & str_detect(xl, "nc") & str_detect(xl, "ori"),
    "ncRNA", out
  )
  
  out <- ifelse(str_detect(xl, "asrna|antisens|antisense"), "antisense", out)
  
  out <- ifelse(out == x0 & str_detect(x0, "AT"), "TES isoform", out)
  out <- ifelse(out == x0 & str_detect(x0, "5'"), "TSS isoform", out)
  out <- ifelse(out == x0 & str_detect(x0, "AT") & str_detect(x0, "5'"), "TES isoform", out)
  
  out <- ifelse(out == x0 & str_detect(xl, "spliced"), "spliced", out)
  
  out <- ifelse(
    out == x0 & (
      str_detect(xl, "polycistronic") |
        (str_detect(xl, "complex") & str_detect(xl, "canonic"))
    ),
    "polycistronic", out
  )
  
  out <- ifelse(
    out == x0 & (
      str_detect(xl, "monocistronic") |
        str_detect(xl, "canonic")
    ),
    "canonic", out
  )
  
  out <- ifelse(
    out == x0 & (
      str_detect(xl, "^nc$") |
        str_detect(xl, "non-coding") |
        str_detect(xl, "\\bnc\\b")
    ),
    "ncRNA", out
  )
  
  out <- ifelse(out == x0, "other", out)
  out
}

cat_levels <- c(
  "polycistronic",
  "canonic",
  "TSS isoform",
  "TES isoform",
  "spliced",
  "antisense",
  "ncRNA",
  "other"
)

prep_one_cat <- function(path, panel_label, default_annotator = "LoRTIA") {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  
  id_col <- case_when(
    "isoform" %in% names(df) ~ "isoform",
    "ref_id"  %in% names(df) ~ "ref_id",
    TRUE ~ NA_character_
  )
  
  if (is.na(id_col)) {
    stop("Nincs 'isoform' vagy 'ref_id' oszlop: ", path)
  }
  
  if (!("Category" %in% names(df))) {
    stop("Missing 'Category' column: ", path)
  }
  
  match_cols <- names(df)[str_detect(names(df), "\\.gff3_match$|_match$")]
  match_cols <- match_cols[!str_detect(match_cols, "pred_ids")]
  
  if (length(match_cols) == 0) {
    out <- df %>%
      transmute(
        ID = .data[[id_col]],
        Category_recode = recode_category(Category),
        Panel = panel_label,
        Annotator = default_annotator
      ) %>%
      count(Panel, Annotator, Category_recode, name = "n") %>%
      group_by(Panel, Annotator) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
  } else {
    out <- df %>%
      transmute(
        ID = .data[[id_col]],
        Category_recode = recode_category(Category),
        across(all_of(match_cols))
      ) %>%
      pivot_longer(
        cols = all_of(match_cols),
        names_to = "Annotator",
        values_to = "called"
      ) %>%
      mutate(
        called = as.numeric(called),
        Annotator = clean_category_annotator_names(Annotator)
      ) %>%
      filter(called > 0) %>%
      mutate(Panel = panel_label) %>%
      count(Panel, Annotator, Category_recode, name = "n") %>%
      group_by(Panel, Annotator) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
  }
  
  out %>%
    mutate(
      Panel = factor(Panel, levels = c("dcDNA", "dRNA")),
      Annotator = factor(Annotator, levels = c("LoRTIA", "BAMBU", "FLAIR", "IsoQuant", "NAGATA")),
      Category_recode = factor(Category_recode, levels = cat_levels)
    ) %>%
    filter(!is.na(Annotator))
}

category_files <- c("category_dcdna.tsv", "category_dRNA.tsv")
stop_if_missing(category_files)
walk(category_files, validate_category_file)

cat_long <- bind_rows(
  prep_one_cat("category_dcdna.tsv", "dcDNA", default_annotator = "LoRTIA"),
  prep_one_cat("category_dRNA.tsv",  "dRNA",  default_annotator = "LoRTIA")
)

category_colors <- c(
  "polycistronic" = "#1b9e77",
  "canonic"       = "#d95f02",
  "TSS isoform"   = "#7570b3",
  "TES isoform"   = "#e7298a",
  "spliced"       = "#66a61e",
  "antisense"     = "#e6ab02",
  "ncRNA"         = "#a6761d",
  "other"         = "#666666"
)

make_stacked_body <- function(panel_name, show_y = TRUE) {
  df_sub <- cat_long %>% filter(Panel == panel_name)
  
  ggplot(df_sub, aes(x = Annotator, y = frac, fill = Category_recode)) +
    geom_col(width = 0.58, colour = NA) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_fill_manual(
      values = category_colors,
      name = "Category",
      drop = FALSE,
      guide = guide_legend(ncol = 1, byrow = TRUE)
    ) +
    labs(
      x = NULL,
      y = if (show_y) "Isoform fraction" else NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.20, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, colour = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.5),
      legend.text = element_text(size = 7.5)
    ) +
    right_border_stacked
}

p_stacked_dc <- make_grey_header_small("dcDNA") / make_stacked_body("dcDNA", TRUE) +
  plot_layout(heights = c(0.09, 0.91))

p_stacked_dr <- make_grey_header_small("dRNA") / make_stacked_body("dRNA", FALSE) +
  plot_layout(heights = c(0.09, 0.91))

p_stacked <- (p_stacked_dc | guide_area() | p_stacked_dr) +
  plot_layout(widths = c(1, 0.36, 1), guides = "collect") &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical"
  )

# =========================================================
# PART C - UPSET PLOTS FROM ONE COMBINED TSV
# =========================================================

clean_combined_upset_annotator_names <- function(x) {
  x %>%
    str_remove("^dRNA_false_positives__") %>%
    str_remove("^dcDNA_false_positives__") %>%
    str_remove("\\.gff3_present$") %>%
    str_remove("_present$") %>%
    str_replace("^FLAIR_.*$", "FLAIR") %>%
    str_replace("^bambu$", "BAMBU") %>%
    str_replace("^Bambu$", "BAMBU")
}

prepare_upset_data_from_combined <- function(path, chemistry = c("dRNA", "dcDNA")) {
  chemistry <- match.arg(chemistry)
  
  df <- readr::read_tsv(path, show_col_types = FALSE)
  
  id_col <- if ("Original_Transcript_ID" %in% names(df)) {
    "Original_Transcript_ID"
  } else if ("Cluster_ID" %in% names(df)) {
    "Cluster_ID"
  } else {
    NA_character_
  }
  
  if (is.na(id_col)) {
    stop("Nincs használható transcript ID oszlop.")
  }
  
  if (chemistry == "dRNA") {
    annot_cols_raw <- names(df)[str_detect(names(df), "^dRNA_false_positives__.*_present$")]
  } else {
    annot_cols_raw <- names(df)[str_detect(names(df), "^dcDNA_false_positives__.*_present$")]
  }
  
  annot_cols_raw <- annot_cols_raw[!str_detect(annot_cols_raw, "__files_present$")]
  annot_cols_clean <- clean_combined_upset_annotator_names(annot_cols_raw)
  
  meta_df <- df %>%
    transmute(
      Tr_name = .data[[id_col]],
      Start   = suppressWarnings(as.numeric(Start)),
      Stop    = suppressWarnings(as.numeric(End)),
      Strand  = Strand
    )
  
  annot_df <- df[, annot_cols_raw, drop = FALSE] %>%
    mutate(across(everything(), ~ as.integer(.x)))
  
  names(annot_df) <- annot_cols_clean
  
  out <- bind_cols(meta_df, annot_df)
  
  preferred_order <- c("LoRTIA", "BAMBU", "FLAIR", "IsoQuant", "NAGATA")
  annot_cols <- c(
    intersect(preferred_order, annot_cols_clean),
    setdiff(annot_cols_clean, preferred_order)
  )
  
  out <- out %>%
    select(Tr_name, Start, Stop, Strand, all_of(annot_cols))
  
  list(df = out, annot_cols = annot_cols)
}

build_upset_components <- function(df, annot_cols) {
  set_order <- annot_cols
  set_order_rev <- rev(set_order)
  
  set_positions <- tibble(
    Annotator = set_order_rev,
    y_num = seq_along(set_order_rev)
  )
  
  set_sizes <- tibble(
    Annotator = set_order,
    Size = colSums(df[annot_cols] == 1, na.rm = TRUE)
  ) %>%
    left_join(set_positions, by = "Annotator")
  
  combo_tbl <- df %>%
    mutate(
      degree = rowSums(across(all_of(annot_cols))),
      combo_key = pmap_chr(across(all_of(annot_cols)), function(...) {
        vals <- as.integer(c(...))
        active <- annot_cols[vals == 1]
        if (length(active) == 0) "None" else paste(active, collapse = " & ")
      })
    ) %>%
    filter(degree > 0) %>%
    count(combo_key, degree, across(all_of(annot_cols)), name = "n") %>%
    filter(n > 2) %>%
    arrange(desc(n), desc(degree)) %>%
    mutate(
      combo_index = row_number(),
      label_y = n + max(n) * 0.04
    )
  
  if (nrow(combo_tbl) == 0) {
    stop("Nincs olyan intersection, amelynek mérete > 2 a kiválasztott kémiára.")
  }
  
  matrix_df <- combo_tbl %>%
    select(combo_index, all_of(annot_cols)) %>%
    pivot_longer(
      cols = all_of(annot_cols),
      names_to = "Annotator",
      values_to = "present"
    ) %>%
    left_join(set_positions, by = "Annotator")
  
  bg_df <- tidyr::expand_grid(
    combo_index = combo_tbl$combo_index,
    Annotator = set_positions$Annotator
  ) %>%
    left_join(set_positions, by = "Annotator")
  
  row_bg <- set_positions %>%
    mutate(ymin = y_num - 0.45, ymax = y_num + 0.45)
  
  list(
    set_sizes = set_sizes,
    combo_tbl = combo_tbl,
    matrix_df = matrix_df,
    bg_df = bg_df,
    row_bg = row_bg,
    set_positions = set_positions
  )
}

make_upset_body <- function(path, chemistry_name = c("dRNA", "dcDNA")) {
  chemistry_name <- match.arg(chemistry_name)
  prep <- prepare_upset_data_from_combined(path, chemistry = chemistry_name)
  comps <- build_upset_components(prep$df, prep$annot_cols)
  
  max_set <- max(comps$set_sizes$Size, na.rm = TRUE)
  y_breaks <- comps$set_positions$y_num
  y_labels <- comps$set_positions$Annotator
  y_limits <- c(0.5, max(y_breaks) + 0.5)
  
  x_max <- max(comps$combo_tbl$combo_index) + 0.5
  x_min <- 0.5
  
  # LEFT-BOTTOM: set sizes
  p_sets <- ggplot() +
    geom_rect(
      data = comps$row_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey94",
      colour = NA
    ) +
    geom_rect(
      data = comps$set_sizes,
      aes(
        xmin = 0, xmax = Size,
        ymin = y_num - 0.23, ymax = y_num + 0.23
      ),
      fill = "#F06A5D",
      colour = "#4D4D4D",
      linewidth = 0.3
    ) +
    geom_text(
      data = comps$set_sizes,
      aes(x = Size + max_set * 0.03, y = y_num, label = Size),
      hjust = 0, size = 2.5, colour = "black"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.18))
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = y_limits,
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.line.x = element_line(linewidth = 0.35, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(colour = "black", size = 8),
      axis.text.y = element_text(colour = "black", size = 9),
      panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 6, b = 2, l = 2)
    ) +
    coord_cartesian(clip = "off")
  
  # RIGHT-TOP: intersection bars
  p_int <- ggplot(comps$combo_tbl, aes(x = combo_index, y = n)) +
    geom_col(width = 0.56, fill = "#6D9EEB", colour = "#4D4D4D", linewidth = 0.3) +
    geom_text(
      aes(y = label_y, label = n),
      size = 2.4, colour = "black", angle = 90, vjust = -0.15
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = NULL, y = "Intersection size") +
    theme_classic(base_size = 9) +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.y = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.y = element_line(linewidth = 0.25, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 8),
      axis.title.y = element_text(colour = "black", size = 9),
      panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 2, b = 2, l = 2)
    )
  
  # RIGHT-BOTTOM: matrix
  p_matrix <- ggplot() +
    geom_rect(
      data = comps$row_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey94",
      colour = "white",
      linewidth = 0.2
    ) +
    geom_point(
      data = comps$bg_df,
      aes(x = combo_index, y = y_num),
      size = 2.6, colour = "grey78"
    ) +
    geom_point(
      data = comps$matrix_df %>% filter(present == 1),
      aes(x = combo_index, y = y_num),
      size = 3.0, colour = "#6D9EEB"
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = rep("", length(y_breaks)),
      limits = y_limits,
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Intersection size") +
    theme_classic(base_size = 9) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = NA),
      axis.title.y = element_text(colour = NA),
      panel.border = element_blank(),
      plot.margin = margin(t = 1, r = 2, b = 1, l = 2)
    )
  
  aligned_right <- cowplot::align_plots(
    p_int, p_matrix,
    align = "v",
    axis = "lr"
  )
  
  top_left_spacer <- ggplot() + theme_void()
  
  top_row <- cowplot::plot_grid(
    top_left_spacer,
    aligned_right[[1]],
    ncol = 2,
    rel_widths = c(0.34, 0.66)
  )
  
  bottom_row <- cowplot::plot_grid(
    p_sets,
    aligned_right[[2]],
    ncol = 2,
    rel_widths = c(0.34, 0.66),
    align = "h",
    axis = "tb"
  )
  
  cowplot::plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(0.70, 0.30),
    align = "v"
  )
}

combined_upset_file <- "FP.tsv"
stop_if_missing(c(combined_upset_file))
validate_combined_upset_file(combined_upset_file)

p_upset_dcdna <- cowplot::plot_grid(
  make_upset_header("dcDNA"),
  make_upset_body(combined_upset_file, chemistry_name = "dcDNA"),
  ncol = 1,
  rel_heights = c(0.05, 0.95),
  align = "v"
)

p_upset_drna <- cowplot::plot_grid(
  make_upset_header("dRNA"),
  make_upset_body(combined_upset_file, chemistry_name = "dRNA"),
  ncol = 1,
  rel_heights = c(0.05, 0.95),
  align = "v"
)

p_upset <- cowplot::plot_grid(
  p_upset_dcdna,
  p_upset_drna,
  ncol = 2,
  rel_widths = c(1, 1)
)

# =========================================================
# FINAL LAYOUT
# =========================================================

row_A <- wrap_elements(full = p_top)
row_B <- wrap_elements(full = p_stacked)
row_C <- wrap_elements(full = p_upset)

p_all <- row_A / row_B / row_C +
  plot_layout(heights = c(1.10, 0.70, 0.74)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, colour = "black"),
    plot.tag.position = c(0.01, 0.99)
  )

print(p_all)

ggsave(
  "Figure1_KSHV_Benchmark_Combined_v14.pdf",
  p_all,
  width = 10.8,
  height = 11.1,
  dpi = 300
)