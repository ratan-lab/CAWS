#!/usr/bin/env Rscript

# Load required libraries
library(rmarkdown)

# Render the HTML report
rmarkdown::render(
  input = "workflow/scripts/generate_html_report.Rmd",
  output_file = basename(as.character(snakemake@output[["html_report"]])),
  output_dir = dirname(as.character(snakemake@output[["html_report"]])),
  quiet = TRUE
)