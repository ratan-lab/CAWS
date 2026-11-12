#!/usr/bin/env Rscript

library(tidyverse)

df1 = tibble()
for (absname in snakemake@input[["aln"]]) {
    filename = strsplit(absname, "/")[[1]][2]
    samplename = strsplit(filename, ".alignments.txt")[[1]][1]
   
    data = read.table(absname, fill=NA, header=FALSE)
    res = tibble(Sample=samplename, 
                 SequencingDepth=as.numeric(data$V1[1]), 
                 MappedFragNum=as.numeric(data$V1[4]) + as.numeric(data$V1[5]),
                 UniqueMappedFragNum=as.numeric(data$V1[4]),
                 AlignmentRate=data$V1[6],
                 UniqueAlignmentRate=paste0(round(UniqueMappedFragNum*100.0/SequencingDepth,2), "%"))
    df1 = bind_rows(df1, res)
}

mt = c()
for (absname in snakemake@input[["mtn"]]) {
    num = scan(absname)
    mt = c(mt, num)
}
df1$MtNumReads = mt
df1 = df1 |> mutate(MtAlignmentRate = round(MtNumReads*100.0/SequencingDepth,2)) |> mutate(MtAlignmentRate = paste0(MtAlignmentRate, "%")) 

df2 = tibble()
for (absname in snakemake@input[["dup"]]) {
    filename = strsplit(absname, "/")[[1]][2]
    samplename = strsplit(filename, ".rmdup.txt")[[1]][1]
   
    # tab-separated: ESTIMATED_LIBRARY_SIZE is empty for single-end data
    data = read.table(absname, nrows=1, header=TRUE, sep="\t")
    res = tibble(Sample=samplename,
              DuplicationFrac=as.numeric(data$PERCENT_DUPLICATION),
              EstimatedLibrarySize=as.numeric(data$ESTIMATED_LIBRARY_SIZE)) |>
          mutate(DuplicationRate = round(DuplicationFrac * 100, 2)) |>
          mutate(DuplicationRate = paste0(DuplicationRate,"%")) |>
          mutate(UniqueFragNum = data$READ_PAIRS_EXAMINED + data$UNPAIRED_READS_EXAMINED -
                                 data$READ_PAIR_DUPLICATES - data$UNPAIRED_READ_DUPLICATES) |>
          select(-DuplicationFrac)
    df2 = bind_rows(df2, res)
}

# FastQ Screen: % of reads mapping only to each contaminant (one or multiple hits)
contaminants <- unlist(snakemake@params[["contaminants"]])
df3 <- snakemake@input[["screen"]] |>
    map_dfr(\(f) read_tsv(f, skip = 1, n_max = length(contaminants) + 1, show_col_types = FALSE) |>
        filter(Genome %in% contaminants) |>
        transmute(Sample = str_remove(basename(f), "_screen.txt$"),
                  Genome,
                  OnlyRate = paste0(round(`%One_hit_one_genome` + `%Multiple_hits_one_genome`, 2), "%"))) |>
    pivot_wider(names_from = Genome, values_from = OnlyRate, names_glue = "{Genome}OnlyRate")

result <- left_join(df1, df2, by=c("Sample"="Sample")) |>
          left_join(df3, by=c("Sample"="Sample"))

write_tsv(result, file=snakemake@output[[1]])
