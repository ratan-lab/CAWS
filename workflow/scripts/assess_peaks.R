#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(viridis)
library(chromVAR)


samplesheet <- read_tsv(as.character(snakemake@params[["samplesheet"]]))
stats <- read_tsv(as.character(snakemake@input[["stats"]]))
folder <- as.character(snakemake@params[["subdir"]])
method <- as.character(snakemake@params[["method"]])

suffix <- NA
if (method == "SEACR") {
    suffix <- ".peaks.stringent.bed"
} else if (method == "MACS") {
    suffix <- "_peaks.narrowPeak"
}
stopifnot(is.na(suffix) == FALSE)

peakN = tibble()
peakW = tibble()
peakO = tibble()
peakF = tibble()

for (absname in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname) == 0L) next
    filename = basename(absname)
    sample = str_split(filename, suffix)[[1]][1]
    peakInfo = read.table(absname, header = FALSE, fill = TRUE)  |> mutate(width = abs(V3-V2))
    peakN = tibble(peakN=nrow(peakInfo), sampleID=sample) |> bind_rows(peakN)
    peakW = tibble(width=peakInfo$width, sampleID=sample) |> bind_rows(peakW)
}

df <- left_join(peakN, samplesheet) |> select(sampleID, condition, peakN)
df |> write_tsv(as.character(snakemake@output[["peakN"]]))

df <- left_join(peakW, samplesheet) |> select(sampleID, condition, width)

# Calculate peak width statistics
peakW_stats <- df |>
    group_by(sampleID, condition) |>
    summarise(
        median_width = median(width, na.rm = TRUE),
        mean_width = mean(width, na.rm = TRUE),
        q25_width = quantile(width, 0.25, na.rm = TRUE),
        q75_width = quantile(width, 0.75, na.rm = TRUE),
        min_width = min(width, na.rm = TRUE),
        max_width = max(width, na.rm = TRUE),
        narrow_peaks = sum(width <= 500, na.rm = TRUE),
        broad_peaks = sum(width > 2000, na.rm = TRUE),
        total_peaks = n(),
        .groups = 'drop'
    ) |>
    mutate(
        pct_narrow = round(narrow_peaks * 100.0 / total_peaks, 2),
        pct_broad = round(broad_peaks * 100.0 / total_peaks, 2)
    )

# Write peak width statistics
write_tsv(peakW_stats, as.character(snakemake@output[["peakW"]]))


for (absname1 in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname1) == 0L) next
    for (absname2 in as.character(snakemake@input[["peaks"]])) {
        if (absname1 == absname2) next
        if (file.size(absname2) == 0L) next

        filename1 = basename(absname1)
        filename2 = basename(absname2)

        sample1 = str_split(filename1, suffix)[[1]][1]
        sample2 = str_split(filename2, suffix)[[1]][1]

        peakInfo1 = read.table(absname1, header = FALSE, fill = TRUE)  
        peakInfo2 = read.table(absname2, header = FALSE, fill = TRUE)  
    
        peakInfo1.gr = GRanges(peakInfo1$V1, IRanges(start = peakInfo1$V2, end = peakInfo1$V3), strand = "*")
        peakInfo2.gr = GRanges(peakInfo2$V1, IRanges(start = peakInfo2$V2, end = peakInfo2$V3), strand = "*")
        overlap.gr = peakInfo1.gr[findOverlaps(peakInfo1.gr, peakInfo2.gr)@from]

        peakO = tibble(sample1=sample1, sample2=sample2, overlap=length(overlap.gr)) |> bind_rows(peakO)
    }
}

df <- left_join(peakO, samplesheet, by=c("sample1"="sampleID")) |>
      select(sample1, condition, sample2, overlap) |>
      dplyr::rename("condition1"="condition") |>
      left_join(samplesheet, by=c("sample2"="sampleID")) |>
      select(sample1, condition1, sample2, condition, overlap) |>
      dplyr::filter(condition1 == condition) |>
      select(sample1, condition, overlap) |>
      left_join(peakN, by=c("sample1"="sampleID")) |>
      select(sample1, condition, overlap, peakN) |>
      mutate(peakReprodRate = round(overlap * 100.0 / peakN, 2)) |>
      select(sample1, condition, peakN, overlap, peakReprodRate) |>
      dplyr::rename("sampleID"="sample1") |>
      dplyr::rename("peakReprodNum"="overlap")
df |> write_tsv(as.character(snakemake@output[["peakR"]]))

# FRiPs
for (absname in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname) == 0L) next
    filename = basename(absname)
    sample = str_split(filename, suffix)[[1]][1]

    peakRes = read.table(absname, header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(folder, "/", sample, ".sorted.qflt.bam")

    # Determine if sample is paired-end from samplesheet
    # This logic must match the Python is_missing_read2() function in Snakefile
    sample_info <- samplesheet |> filter(sampleID == sample)
    is_paired <- if (nrow(sample_info) > 0) {
        read2_value <- sample_info$read2[1]
        # Check if read2 is missing/empty (indicates single-end)
        # List of values that indicate missing read2 (matches Python version)
        missing_indicators <- c("", "na", "n/a", "null", "none", "-", "nan", "nil", "empty", "0", "false", "missing", "absent", "single", "se")
        !(is.na(read2_value) || trimws(tolower(as.character(read2_value))) %in% missing_indicators)
    } else {
        stop(paste("Sample", sample, "not found in samplesheet"))
    }

    fragment_counts = getCounts(bamFile, peak.gr, paired = is_paired, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] |> sum()
    peakF = tibble(inPeakN=inPeakN, sampleID=sample) |> bind_rows(peakF)
}

df <- left_join(peakF, samplesheet) |>
      select(sampleID, condition, inPeakN) |>
      left_join(stats, by=c("sampleID"="Sample")) |>
      select(sampleID, condition, inPeakN, MappedFragNum) |>
      mutate(FRiP = round(inPeakN * 100.0 /MappedFragNum, 2))
df |> write_tsv(as.character(snakemake@output[["peakF"]]))

# GC Content Analysis
peakG = tibble()
for (gc_file in as.character(snakemake@input[["gc_files"]])) {
    if (file.exists(gc_file) && file.size(gc_file) > 0L) {
        gc_data = read.table(gc_file, header = FALSE, col.names = c("sampleID", "median_gc"))
        peakG = bind_rows(peakG, gc_data)
    }
}

df <- left_join(peakG, samplesheet) |> select(sampleID, condition, median_gc)
df |> write_tsv(as.character(snakemake@output[["peakG"]]))


 
