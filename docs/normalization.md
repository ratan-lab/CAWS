# Normalization

CAWS scales the coverage tracks (`tdf/`, `bigwig/`) and the region signal profiles in one of two ways, set by `normalization`. Peak calling, FRiP, consensus peaks and the other statistics use the unscaled data and are the same in both modes.

**CPM (default).** Coverage is divided by the number of filtered fragments of the sample (excluding `mt_chrom`) and multiplied by one million. This corrects for sequencing depth. It assumes that every library contains about the same total amount of signal, so it shows *where* signal is, relative to the rest of the library. It cannot show a change in the *total* amount of signal. Libraries are amplified and sequenced to similar depths, so a sample in which the target is twice as abundant everywhere looks the same as the untreated sample after CPM.

**Spike-in.** CUT&Tag reactions contain carry-over *E. coli* DNA from the production of pA-Tn5. Its amount depends on the enzyme, not on the sample, so it can serve as a fixed reference. With `normalization: "spikein"`, coverage is multiplied by 10,000 / the number of spike-in fragments of the sample (fragments per 10,000 spike-in fragments, as in the CUT&Tag protocol of Zheng, Ahmad and Henikoff). This keeps global changes. For example, if a G4-stabilizing ligand doubles BG4 signal at all G4s, spike-in scaled tracks show the doubling, while CPM hides most of it and shows an apparent loss of background.

**When to use spike-in normalization:**
- Treatments or genotypes expected to change the total amount of the target, e.g. G4 ligands with BG4, inhibitors of a histone methyltransferase, or knockdown of the target protein.
- Only when the samples to be compared had the same number of cells and used the **same batch of pA-Tn5**. Carry-over differs between batches, and different cell numbers change the ratio of human to *E. coli* DNA for technical reasons.
- Only with enough spike-in fragments. Check `SpikeinFragNum` in `stats/stats.txt` first. Fewer than about 1,000 fragments give noisy scale factors, and CAWS writes a warning to `logs/count_spikein/` for such samples.

Use CPM for questions about *where* the target binds (peak landscapes, region enrichment) and for histone marks without expected global changes. With spike-in normalization, keep in mind that strong competition for pA-Tn5 could also change how efficiently the carry-over DNA is tagmented. An exogenous spike-in, such as a fixed amount of *Drosophila* chromatin added to each reaction, avoids this. Its genome can be given in `contaminant_genomes` and selected as `spikein_genome` in the same way.

**Configuration:**

```json
"contaminant_genomes": {"Ecoli": "/path/to/ecoli_index"},
"spikein_genome": "Ecoli",
"normalization": "spikein"
```

Setting `spikein_genome` with `"normalization": "cpm"` counts spike-in fragments for every sample (columns `SpikeinFragNum` and `SpikeinRate` in `stats/stats.txt`) without changing the tracks. This shows whether spike-in normalization would be usable before switching it on.

**How spike-in fragments are counted.** The reads of a sample that did not align to the host genome (pairs with neither mate aligned, for paired-end data) are aligned to the spike-in genome with the same Bowtie2 settings as the host. Fragments with MAPQ ≥ `minquality` are counted (paired-end: concordant pairs). Duplicates are removed from the count when the sample is deduplicated, so the spike-in is counted the same way as the host signal it scales. Because only reads that did not align to the host are used, sequences shared by the two genomes are not counted as spike-in. BAM input must contain the unaligned reads, otherwise no spike-in fragments are found, and spike-in normalization stops with an error for that sample.

**Outputs:**
- `bigwig/{sample}.bw` and `tdf/{sample}.tdf`: normalized coverage tracks (the scale factor used is in `logs/bam_to_bed_pe/` or `logs/bam_to_bed_se/`)
- `alignments/{sample}.spikein.txt`: spike-in fragment count
- `stats/stats.txt`: `SpikeinFragNum` and `SpikeinRate`
