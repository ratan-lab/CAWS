#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(Rsamtools)
library(viridis)


samplesheet <- read_tsv(as.character(snakemake@params[["samplesheet"]]))
folder <- as.character(snakemake@params[["subdir"]])
method <- as.character(snakemake@params[["method"]])

# peak file suffix, e.g. ".peaks.stringent.bed" or "_peaks.narrowPeak" / "_peaks.broadPeak"
suffix <- as.character(snakemake@params[["suffix"]])

peakN = tibble()
peakW = tibble()
peakF = tibble()

for (absname in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname) == 0L) next
    filename = basename(absname)
    sample = str_split(filename, fixed(suffix))[[1]][1]
    peakInfo = read.table(absname, header = FALSE, fill = TRUE) |> mutate(width = V3 - V2) |> filter(width > 0)
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


# Reproducibility: peaks of each sample that overlap a peak in any other
# replicate of the same condition (one row per sample)
peak_files <- as.character(snakemake@input[["peaks"]])
peak_files <- peak_files[file.size(peak_files) > 0L]
peak_gr <- peak_files |>
    set_names(\(f) str_remove(basename(f), fixed(suffix))) |>
    map(\(f) read.table(f, header = FALSE, fill = TRUE) |> filter(V3 > V2) |>
        with(GRanges(V1, IRanges(start = V2, end = V3), strand = "*")))
cond <- samplesheet |> select(sampleID, condition) |> deframe()

df <- tibble(sampleID = names(peak_gr)) |>
      mutate(condition = cond[sampleID],
             others = map2(sampleID, condition, \(s, c) setdiff(names(peak_gr)[cond[names(peak_gr)] == c], s))) |>
      filter(lengths(others) > 0) |>
      mutate(peakReprodNum = map2_int(sampleID, others, \(s, o)
                 sum(overlapsAny(peak_gr[[s]], reduce(do.call(c, unname(peak_gr[o])))))),
             peakN = map_int(sampleID, \(s) length(peak_gr[[s]])),
             peakReprodRate = round(peakReprodNum * 100.0 / peakN, 2)) |>
      select(sampleID, condition, peakN, peakReprodNum, peakReprodRate)
df |> write_tsv(as.character(snakemake@output[["peakR"]]))

# FRiPs
missing_indicators <- c("", "na", "n/a", "null", "none", "-", "nan", "nil", "empty", "0", "false", "missing", "absent", "single", "se")

for (absname in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname) == 0L) next
    filename = basename(absname)
    sample = str_split(filename, fixed(suffix))[[1]][1]

    peakRes = read.table(absname, header = FALSE, fill = TRUE)
    peak.gr = reduce(GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*"))
    bamFile = paste0(folder, "/", sample, ".sorted.qflt.bam")
    baiFile = paste0(bamFile, ".bai")
    if (!file.exists(baiFile)) stop(paste("BAI index not found:", baiFile))

    sample_info <- samplesheet |> filter(sampleID == sample)
    if (nrow(sample_info) == 0) stop(paste("Sample", sample, "not found in samplesheet"))
    read2_value <- sample_info$read2[1]
    is_paired <- !(is.na(read2_value) || trimws(tolower(as.character(read2_value))) %in% missing_indicators)

    flag <- if (is_paired) {
        scanBamFlag(isUnmappedQuery = FALSE, isFirstMateRead = TRUE)
    } else {
        scanBamFlag(isUnmappedQuery = FALSE)
    }
    param <- ScanBamParam(which = peak.gr, flag = flag)
    inPeakN <- sum(countBam(bamFile, param = param)$records)
    totalN <- sum(idxstatsBam(bamFile)$mapped)
    if (is_paired) totalN <- totalN / 2L
    peakF = tibble(inPeakN=inPeakN, totalN=totalN, sampleID=sample) |> bind_rows(peakF)
}

df <- left_join(peakF, samplesheet) |>
      select(sampleID, condition, inPeakN, totalN) |>
      mutate(FRiP = round(inPeakN * 100.0 / totalN, 2)) |>
      select(-totalN)
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


 
