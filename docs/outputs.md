# Outputs

## Output structure

```
output_directory/
├── qc/                  # FastQC reports; fastq_screen/ contaminant screens
├── trimmed/             # Adapter-trimmed reads
├── alignments/          # alignment statistics, complexity, spike-in counts, fragment lengths, mtDNA counts, unfiltered bigWigs; units/: Bowtie2 logs per run
├── dedupalignments/     # BAMs used downstream (deduplicated where dedup applies), Picard metrics, filtered BAMs
├── bedalignments/       # Raw-count bedgraphs (SEACR input) and 500 bp binned counts
├── peaks/
│   ├── seacr/           # SEACR peaks
│   ├── macs3/           # MACS3 peaks
│   └── consensus/       # pooled peak calls and consensus peaks per condition (seacr/, macs3/)
├── stats/               # stats.txt, samples.tsv, peaks_consolidated.txt, consensus_peaks.txt, region_enrichment.txt
│   ├── enrichment/      # per region set: profile plot, heatmap, matrix and the three enrichment measures
│   ├── seacr/           # SEACR metrics, peak GC, peak heatmaps
│   └── macs3/           # MACS3 metrics, peak GC
├── tdf/                 # IGV TDF coverage tracks, normalized (one per sample), igv_session.xml
├── bigwig/              # the same normalized coverage tracks as bigWig (one per sample)
├── reports/             # HTML report, MultiQC report, pipeline_info/ (config, samplesheet, software versions)
├── annotations/         # tss.bed; regions/: reference points and merged regions per region set
└── logs/                # Log files
```

**Note on alignments**: every sample has its BAM in `dedupalignments/`: deduplicated where `dedup` applies, otherwise an unchanged copy of the alignments. The raw `alignments/{sample}.sorted.bam` (and index) are temporary and removed once all steps that use them have finished; run with `--notemp` to keep them.

**Note on IGV sessions**: open `tdf/igv_session.xml` in IGV (File → Open Session) to load the TDF tracks of all samples (experimental samples grouped by condition, then controls), the consensus peaks and the per-sample peak calls. Paths are relative, so keep the output directory layout when copying it.

**Note on consensus peaks**: `peaks/consensus/{seacr,macs3}/{condition}.consensus.bed` lists, for each consensus peak, the number and names of supporting replicates. Conditions with a single replicate have no consensus peaks.

**Note on coverage tracks**: the TDF and bigWig tracks hold the same coverage of the filtered fragments, scaled to counts per million filtered fragments (CPM, excluding fragments on `mt_chrom` from the total) or, with `normalization: "spikein"`, to fragments per 10,000 spike-in fragments. See [Normalization](normalization.md). The raw-count bedgraphs in `bedalignments/` are kept for SEACR. The unfiltered `alignments/{sample}_raw.bw` tracks are not normalized and are used for the peak heatmaps only.

## Interactive HTML report

**Report:** `reports/cutntag_analysis_report.html`  

### Dashboard Pages  
- **Sample Details:** per-sample peak counts, FRiP and peak GC for SEACR and MACS3, plus TSS profile score (if a GTF is provided), and consensus peak counts per condition  
- **Alignment:** alignment statistics table (including whether each sample was deduplicated) and rates (alignment, duplication, mtDNA, contaminants)  
- **Fragment Length:** fragment size distributions (paired-end samples)  
- **Coverage Similarity:** correlation of binned coverage between all samples  
- **SEACR / MACS3:** peak counts, FRiP, width distributions and width categories by condition  
- **SEACR Peaks / MACS3 Peaks:** pairwise peak overlap heatmaps (shared / total merged peak intervals)  
- **Region Enrichment:** enrichment over TSSs and user region sets, with IgG ratios  

**Interactive Features:**  
- Plotly-based zoom/pan/hover  
- Sortable, searchable tables with export buttons  

A separate **MultiQC** report (`reports/multiqc_report.html`) collects the FastQC, Trim Galore, Bowtie2, Picard and FastQ Screen results in the standard MultiQC format.
