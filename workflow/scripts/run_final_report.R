#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

filename <- file.path(args[1], "stats", "stats.txt")
stats <- read_tsv(filename, col_names=TRUE)

# SEACR files
filename <- file.path(args[1], "stats", "seacr", "peaks_frip.txt")
frip <- read_tsv(filename, col_names=TRUE)

filename <- file.path(args[1], "stats", "seacr", "peaks_num.txt")
num <- read_tsv(filename, col_names=TRUE)

seacr <- left_join(frip, num, by=c("sampleID"="sampleID", "condition"="condition")) |> select(sampleID, condition, everything()) |> distinct()
colnames(seacr) <- c("sampleID", "condition", paste0("seacr_", colnames(seacr)[3:ncol(seacr)]))

# MACS files 
filename <- file.path(args[1], "stats", "macs3", "peaks_frip.txt")
frip <- read_tsv(filename, col_names=TRUE)

filename <- file.path(args[1], "stats", "macs3", "peaks_num.txt")
num <- read_tsv(filename, col_names=TRUE)

macs3 <- left_join(frip, num, by=c("sampleID"="sampleID", "condition"="condition")) |> select(sampleID, condition, everything()) |> distinct()
colnames(macs3) <- c("sampleID", "condition", paste0("macs3_", colnames(macs3)[3:ncol(macs3)]))

filename <- file.path(args[1], "stats", "combined_stats.txt")
left_join(stats, seacr, by=c("Sample"="sampleID")) |>
    distinct() |>
    left_join(macs3, by=c("Sample"="sampleID", "condition"="condition")) |>
    distinct() |>
    write_tsv(filename)
