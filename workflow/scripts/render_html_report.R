#!/usr/bin/env Rscript

# Load required libraries
library(rmarkdown)

new_dir_path <- dirname(as.character(snakemake@output[["html_report"]]))
if (!dir.exists(new_dir_path)) {
  dir.create(new_dir_path, recursive = TRUE)
}

# Render the HTML report
rmarkdown::render(
  input = paste0(snakemake@input[["snakefile_dir"]], "/scripts/generate_html_report.Rmd"),
  output_file = basename(as.character(snakemake@output[["html_report"]])),
  output_dir = dirname(as.character(snakemake@output[["html_report"]])),
  quiet = TRUE
)
