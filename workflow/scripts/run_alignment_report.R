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
   
    data = read.table(absname, nrows=1, header=TRUE) 
    res = tibble(Sample=samplename, 
              DuplicationFrac=as.numeric(data$PERCENT_DUPLICATION),
              EstimatedLibrarySize=as.numeric(data$ESTIMATED_LIBRARY_SIZE)) |>
          mutate(DuplicationRate = round(DuplicationFrac * 100, 2)) |>
          mutate(DuplicationRate = paste0(DuplicationRate,"%")) |>
          mutate(UniqueFragNum=round(as.numeric(data$READ_PAIRS_EXAMINED) * (1 - DuplicationFrac))) |>
          select(-DuplicationFrac)
    df2 = bind_rows(df2, res)
}

df3 <- tibble()
for (absname in snakemake@input[["eco"]]) {
    filename = strsplit(absname, "/")[[1]][2]
    samplename = strsplit(filename, ".ecoli.txt")[[1]][1]

    data = read.table(absname, fill=NA, header=FALSE)
    res = tibble(Sample=samplename,
                 EcoliMappedFrag=as.numeric(data$V1[4]) + as.numeric(data$V1[5]),
                 EcoliAlignmentRate=data$V1[6])
    df3 = bind_rows(df3, res)
}


result <- left_join(df1, df2, by=c("Sample"="Sample")) |>
          left_join(df3, by=c("Sample"="Sample"))

write_tsv(result, file=snakemake@output[[1]])
