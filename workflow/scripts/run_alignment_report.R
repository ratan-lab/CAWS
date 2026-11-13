#!/usr/bin/env Rscript

library(tidyverse)

sample_of <- \(f, ext) str_remove(basename(f), fixed(ext))

# Alignment statistics from the raw (merged) BAM: total, mapped, unique, unaligned fragments
df1 <- snakemake@input[["aln"]] |>
    map_dfr(\(f) read_tsv(f, col_names = c("key", "value"), show_col_types = FALSE) |>
        deframe() |> as.list() |> as_tibble() |>
        mutate(Sample = sample_of(f, ".alnstats.txt"))) |>
    transmute(Sample,
              SequencingDepth = total,
              MappedFragNum = mapped,
              UniqueMappedFragNum = unique,
              AlignmentRate = paste0(round(mapped * 100.0 / total, 2), "%"),
              UniqueAlignmentRate = paste0(round(unique * 100.0 / total, 2), "%"),
              UnalignedFragNum = unaligned)

mt <- snakemake@input[["mtn"]] |>
    map_dfr(\(f) tibble(Sample = sample_of(f, ".mtdna.txt"), MtNumReads = scan(f, quiet = TRUE)))
df1 <- df1 |>
    left_join(mt, by = "Sample") |>
    mutate(MtAlignmentRate = paste0(round(MtNumReads * 100.0 / SequencingDepth, 2), "%"),
           .after = UniqueAlignmentRate) |>
    relocate(MtNumReads, .before = MtAlignmentRate)

# ENCODE library complexity; PCR bottleneck graded by PBC1 (ENCODE thresholds)
complexity <- snakemake@input[["complexity"]] |>
    map_dfr(\(f) read_tsv(f, col_types = "ddd", na = "NA") |> mutate(Sample = sample_of(f, ".complexity.txt"))) |>
    mutate(PCRBottleneck = case_when(is.na(PBC1) ~ NA_character_,
                                     PBC1 < 0.5 ~ "severe",
                                     PBC1 < 0.8 ~ "moderate",
                                     PBC1 < 0.9 ~ "mild",
                                     TRUE ~ "none"))

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

# whether duplicates were removed for each sample (dedup setting / samplesheet column)
dedup_samples <- unlist(snakemake@params[["dedup_samples"]])
df2 <- df2 |>
    mutate(Deduplicated = if_else(Sample %in% dedup_samples, "yes", "no"), .after = DuplicationRate)

result <- df1 |>
    select(-UnalignedFragNum) |>
    left_join(df2, by = "Sample") |>
    left_join(complexity, by = "Sample") |>
    left_join(df1 |> select(Sample, UnalignedFragNum), by = "Sample")

# spike-in fragments (deduplicated like the host alignments), when spikein_genome is set
spikein_files <- unlist(snakemake@input[["spikein"]])
if (length(spikein_files) > 0) {
    spikein <- spikein_files |>
        map_dfr(\(f) read_tsv(f, col_types = "d") |> mutate(Sample = sample_of(f, ".spikein.txt"))) |>
        transmute(Sample, SpikeinFragNum = spikein_fragments)
    result <- result |>
        left_join(spikein, by = "Sample") |>
        mutate(SpikeinRate = paste0(round(SpikeinFragNum * 100.0 / SequencingDepth, 3), "%"))
}

result <- result |>
    left_join(df3, by = "Sample")

write_tsv(result, file=snakemake@output[[1]])
