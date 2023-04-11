#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(viridis)
library(chromVAR)
library(ggpubr)

samplesheet <- read_tsv(snakemake@params[["samplesheet"]])
stats <- read_tsv(snakemake@input[["stats"]])
folder <- snakemake@params[["subdir"]]

peakN = tibble()
peakW = tibble()
peakO = tibble()
peakF = tibble()

for (absname in snakemake@input[["peaks"]]) {
    filename = str_split(absname, "/")[[1]][2]
    sample = str_split(filename, ".peaks.stringent.bed")[[1]][1]
    peakInfo = read.table(absname, header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
    peakN = tibble(peakN=nrow(peakInfo), sampleID=sample) %>% bind_rows(peakN)
    peakW = tibble(width=peakInfo$width, sampleID=sample) %>% bind_rows(peakW)
}

df <- left_join(peakN, samplesheet) %>% select(sampleID, condition, peakN)
df %>% write_tsv(snakemake@output[["peakN"]])
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

for (absname1 in snakemake@input[["peaks"]]) {
    for (absname2 in snakemake@input[["peaks"]]) {
        if (absname1 == absname2) next

        filename1 = str_split(absname1, "/")[[1]][2]
        filename2 = str_split(absname2, "/")[[1]][2]

        sample1 = str_split(filename1, ".peaks.stringent.bed")[[1]][1]
        sample2 = str_split(filename2, ".peaks.stringent.bed")[[1]][1]

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
df %>% write_tsv(snakemake@output[["peakR"]])
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
for (absname in snakemake@input[["peaks"]]) {
    filename = str_split(absname, "/")[[1]][2]
    sample = str_split(filename, ".peaks.stringent.bed")[[1]][1]
    
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
df %>% write_tsv(snakemake@output[["peakF"]])
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

pdf(snakemake@output[["fig"]], width=14, height=14)
ggarrange(figA, figB, figC, figD, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

 
