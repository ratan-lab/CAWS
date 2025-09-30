#!/usr/bin/env Rscript

library(tidyverse)

# Read SEACR stats files
seacr_num <- read_tsv(as.character(snakemake@input[["seacr_num"]]))
seacr_reprod <- read_tsv(as.character(snakemake@input[["seacr_reprod"]]))
seacr_frip <- read_tsv(as.character(snakemake@input[["seacr_frip"]]))

# Read MACS3 stats files
macs3_num <- read_tsv(as.character(snakemake@input[["macs3_num"]]))
macs3_reprod <- read_tsv(as.character(snakemake@input[["macs3_reprod"]]))
macs3_frip <- read_tsv(as.character(snakemake@input[["macs3_frip"]]))

# Consolidate SEACR stats
seacr_consolidated <- seacr_num %>%
  select(sampleID, condition, peakN) %>%
  rename(seacr_peakN = peakN) %>%
  left_join(
    seacr_reprod %>% select(sampleID, peakReprodNum, peakReprodRate) %>%
    rename(seacr_peakReprodNum = peakReprodNum, seacr_peakReprodRate = peakReprodRate),
    by = "sampleID"
  ) %>%
  left_join(
    seacr_frip %>% select(sampleID, FRiP) %>%
    rename(seacr_FRiP = FRiP),
    by = "sampleID"
  )

# Consolidate MACS3 stats
macs3_consolidated <- macs3_num %>%
  select(sampleID, peakN) %>%
  rename(macs3_peakN = peakN) %>%
  left_join(
    macs3_reprod %>% select(sampleID, peakReprodNum, peakReprodRate) %>%
    rename(macs3_peakReprodNum = peakReprodNum, macs3_peakReprodRate = peakReprodRate),
    by = "sampleID"
  ) %>%
  left_join(
    macs3_frip %>% select(sampleID, FRiP) %>%
    rename(macs3_FRiP = FRiP),
    by = "sampleID"
  )

# Combine both methods
consolidated_stats <- seacr_consolidated %>%
  left_join(macs3_consolidated, by = "sampleID") %>%
  select(sampleID, condition,
         seacr_peakN, seacr_peakReprodNum, seacr_peakReprodRate, seacr_FRiP,
         macs3_peakN, macs3_peakReprodNum, macs3_peakReprodRate, macs3_FRiP)

# Write consolidated stats
consolidated_stats %>% write_tsv(as.character(snakemake@output[["consolidated"]]))