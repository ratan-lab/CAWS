#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(viridis)
library(chromVAR)
library(ggpubr)

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
    filename = str_split(absname, "/")[[1]][3]
    sample = str_split(filename, suffix)[[1]][1]
    peakInfo = read.table(absname, header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
    peakN = tibble(peakN=nrow(peakInfo), sampleID=sample) %>% bind_rows(peakN)
    peakW = tibble(width=peakInfo$width, sampleID=sample) %>% bind_rows(peakW)
}

df <- left_join(peakN, samplesheet) %>% select(sampleID, condition, peakN)
df %>% write_tsv(as.character(snakemake@output[["peakN"]]))
figA <- df %>% 
        ggplot(aes(x = condition, y = peakN, fill = condition)) +
        geom_boxplot() +
        geom_jitter(aes(), position = position_jitter(0.15)) +
        scale_fill_viridis(discrete=TRUE, begin=0.1, end=0.55, option="magma", alpha=0.8) +
        scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
        theme_bw(base_size = 18) +
        ylab("Number of Peaks") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45))

df <- left_join(peakW, samplesheet) %>% select(sampleID, condition, width)

# Calculate peak width statistics
peakW_stats <- df %>%
    group_by(sampleID, condition) %>%
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
    ) %>%
    mutate(
        pct_narrow = round(narrow_peaks * 100.0 / total_peaks, 2),
        pct_broad = round(broad_peaks * 100.0 / total_peaks, 2)
    )

# Write peak width statistics
write_tsv(peakW_stats, sub("peaks_num.txt", "peaks_width_stats.txt", as.character(snakemake@output[["peakN"]])))

figB = df %>%
        ggplot(aes(x=sampleID, y=width, fill=condition)) +
        geom_violin() +
        scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
        scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
        scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
        theme_bw(base_size = 18) +
        ylab("Width of Peaks") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45))

# Enhanced peak width distribution plot
figB2 <- df %>%
    ggplot(aes(x = width, fill = condition)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_x_continuous(trans = "log10", breaks = c(100, 500, 1000, 5000, 10000)) +
    facet_wrap(~condition, scales = "free_y") +
    theme_bw(base_size = 14) +
    xlab("Peak Width (bp)") +
    ylab("Count") +
    ggtitle("Peak Width Distribution by Condition")

# Peak width category plot
figB3 <- peakW_stats %>%
    select(sampleID, condition, pct_narrow, pct_broad) %>%
    pivot_longer(cols = c(pct_narrow, pct_broad), names_to = "category", values_to = "percentage") %>%
    mutate(category = case_when(
        category == "pct_narrow" ~ "Narrow (≤500bp)",
        category == "pct_broad" ~ "Broad (>2000bp)"
    )) %>%
    ggplot(aes(x = sampleID, y = percentage, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.8, option = "plasma") +
    theme_bw(base_size = 14) +
    ylab("Percentage of Peaks") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(fill = "Peak Category")

for (absname1 in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname1) == 0L) next
    for (absname2 in as.character(snakemake@input[["peaks"]])) {
        if (absname1 == absname2) next
        if (file.size(absname2) == 0L) next

        filename1 = str_split(absname1, "/")[[1]][3]
        filename2 = str_split(absname2, "/")[[1]][3]

        sample1 = str_split(filename1, suffix)[[1]][1]
        sample2 = str_split(filename2, suffix)[[1]][1]

        peakInfo1 = read.table(absname1, header = FALSE, fill = TRUE)  
        peakInfo2 = read.table(absname2, header = FALSE, fill = TRUE)  
    
        peakInfo.gr = GRanges(peakInfo1$V1, IRanges(start = peakInfo1$V2, end = peakInfo1$V3), strand = "*")
        overlap.gr = peakInfo.gr
        peakInfo.gr = GRanges(peakInfo2$V1, IRanges(start = peakInfo2$V2, end = peakInfo2$V3), strand = "*")
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]

        peakO = tibble(sample1=sample1, sample2=sample2, overlap=length(overlap.gr)) %>% bind_rows(peakO)
    }
}

df <- left_join(peakO, samplesheet, by=c("sample1"="sampleID")) %>%
      select(sample1, condition, sample2, overlap) %>%
      dplyr::rename("condition1"="condition") %>%
      left_join(samplesheet, by=c("sample2"="sampleID")) %>%
      select(sample1, condition1, sample2, condition, overlap) %>%
      dplyr::filter(condition1 == condition) %>%
      select(sample1, condition, overlap) %>%
      left_join(peakN, by=c("sample1"="sampleID")) %>%
      select(sample1, condition, overlap, peakN) %>%
      mutate(peakReprodRate = round(overlap * 100.0 / peakN, 2)) %>%
      select(sample1, condition, peakN, overlap, peakReprodRate) %>%
      dplyr::rename("sampleID"="sample1") %>%
      dplyr::rename("peakReprodNum"="overlap")
df %>% write_tsv(as.character(snakemake@output[["peakR"]]))
figC <- df %>% 
        ggplot(aes(x=sampleID, y=peakReprodRate, fill=condition, label=peakReprodRate)) +
        geom_bar(stat = "identity") +
        geom_text(vjust = 0.1) +
        scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
        scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
        theme_bw(base_size = 18) +
        ylab("% of Peaks Reproduced") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45))

# FRiPs
for (absname in as.character(snakemake@input[["peaks"]])) {
    if (file.size(absname) == 0L) next
    filename = str_split(absname, "/")[[1]][3]
    sample = str_split(filename, suffix)[[1]][1]
    
    peakRes = read.table(absname, header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(folder, "/", sample, ".sorted.qflt.bam")
    fragment_counts = getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    peakF = tibble(inPeakN=inPeakN, sampleID=sample) %>% bind_rows(peakF)
}

df <- left_join(peakF, samplesheet) %>%
      select(sampleID, condition, inPeakN) %>%
      left_join(stats, by=c("sampleID"="Sample")) %>%
      select(sampleID, condition, inPeakN, MappedFragNum) %>%
      mutate(FRiP = round(inPeakN * 100.0 /MappedFragNum, 2))
df %>% write_tsv(as.character(snakemake@output[["peakF"]]))
figD <- df %>% 
        ggplot(aes(x=condition, y=FRiP, fill=condition, label=FRiP)) +
        geom_boxplot() +
        geom_jitter() +
        scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
        scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
        theme_bw(base_size = 18) +
        ylab("% of Fragments in Peaks") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45))

# Create multi-page PDF with enhanced peak width analysis
pdf(as.character(snakemake@output[["fig"]]), width=16, height=20)

# Page 1: Original 4-panel plot
ggarrange(figA, figB, figC, figD, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")

# Page 2: Enhanced peak width analysis
ggarrange(figB2, figB3, ncol = 1, nrow=2, common.legend = FALSE, heights = c(1, 1))

dev.off()

# Also create separate peak width analysis PDF
pdf(sub("peaks_fig.pdf", "peaks_width_analysis.pdf", as.character(snakemake@output[["fig"]])), width=14, height=10)
ggarrange(figB2, figB3, ncol = 1, nrow=2, common.legend = FALSE, heights = c(1, 1))
dev.off()

 
