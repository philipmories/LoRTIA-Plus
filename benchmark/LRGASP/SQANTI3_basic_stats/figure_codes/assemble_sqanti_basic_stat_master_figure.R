library(ggplot2)
library(patchwork)
library(grid)

# ==========================================================
# BEÁLLÍTÁSOK
# ==========================================================
outdir_master <- "Figure2_MASTER"
outfile_pdf   <- "Figure2_ABCD_MASTER.pdf"

pdf_width  <- 18
pdf_height <- 26   # emeld/csökkentsd igény szerint

# blokkok egymáshoz viszonyított magassága
hA <- 1.35
hB <- 1.00
hC <- 1.00
hD <- 1.05

dir.create(outdir_master, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 0) ELLENŐRZÉS: objektumok megvannak-e
# ==========================================================
req_objs <- c("Ablock_MAIN", "IsoLen_MAIN", "Coding_MAIN", "IsoformsPerGene_MAIN")
missing <- req_objs[!vapply(req_objs, exists, logical(1), inherits = TRUE)]
if (length(missing) > 0) {
  stop(
    "Hiányzó objektum(ok) a sessionből: ", paste(missing, collapse = ", "), "\n",
    "Futtasd le a megfelelő script(ek)et úgy, hogy ezek *_MAIN néven exportálódjanak."
  )
}

# ==========================================================
# 1) TRÜKK: minden blokkot 1 panelként “becsomagolunk”
#    -> így NEM lesz extra betű (E/F/G), és nem csúszik el
# ==========================================================
as_one_panel <- function(x) {
  patchwork::wrap_elements(full = x)
}

pA <- as_one_panel(get("Ablock_MAIN", inherits = TRUE))
pB <- as_one_panel(get("IsoLen_MAIN", inherits = TRUE))
pC <- as_one_panel(get("Coding_MAIN", inherits = TRUE))
pD <- as_one_panel(get("IsoformsPerGene_MAIN", inherits = TRUE))

# ==========================================================
# 2) LAYOUT: egymás alá (A, B, C, D)
# ==========================================================
fig2_master <- (pA / pB / pC / pD) +
  plot_layout(heights = c(hA, hB, hC, hD)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1),
    plot.tag.margin = margin(4, 6, 0, 0)
  )

# ==========================================================
# 3) MENTÉS
# ==========================================================
outpath <- file.path(outdir_master, outfile_pdf)
ggsave(outpath, fig2_master, width = pdf_width, height = pdf_height, limitsize = FALSE)

message("Mentve: ", outpath)
