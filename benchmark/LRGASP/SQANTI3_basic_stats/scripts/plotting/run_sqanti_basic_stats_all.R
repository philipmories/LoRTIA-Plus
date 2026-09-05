#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_sqanti_basic_stats_all.R OUTPUT_DIRECTORY")
}

output_dir <- normalizePath(args[[1]], mustWork = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Cannot determine the driver-script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
figure_code_dir <- dirname(script_path)

inputs <- c(
  SQANTI_BASIC_STAT_FILE = file.path(output_dir, "Figure_input", "SQANTI3_basic_stat.tsv"),
  SQANTI_LENGTH_STATS_FILE = file.path(output_dir, "Figure_input", "sqanti_length_boxplot_stats.tsv"),
  SQANTI_ISOFORMS_GENE_FILE = file.path(output_dir, "Figure_input", "isoformsgene.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing plot input(s):\n", paste(missing_inputs, collapse = "\n"))
}

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr", "scales", "grid",
  "patchwork", "ggh4x"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them before rerunning the workflow."
  )
}

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(output_dir)
do.call(Sys.setenv, as.list(inputs))

source(file.path(figure_code_dir, "plot_sqanti_yield_and_category_composition.R"), chdir = FALSE)
source(file.path(figure_code_dir, "plot_sqanti_isoform_length_distributions.R"), chdir = FALSE)
source(file.path(figure_code_dir, "plot_sqanti_coding_noncoding_composition.R"), chdir = FALSE)
source(file.path(figure_code_dir, "plot_sqanti_isoforms_per_gene.R"), chdir = FALSE)
source(file.path(figure_code_dir, "assemble_sqanti_basic_stat_master_figure.R"), chdir = FALSE)

message("[OK] SQANTI3 basic-statistics plots completed under: ", output_dir)
