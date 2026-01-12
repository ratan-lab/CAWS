---
title: 'CAWS: A comprehensive Snakemake workflow for CUT&Tag analysis'
tags:
  - Python
  - R
  - Snakemake
  - bioinformatics
  - epigenomics
  - CUT&Tag
  - chromatin profiling
  - peak calling
authors:
  - name: Akrosh Ratan
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
affiliations:
 - name: Department of Public Health Sciences, University of Virginia, Charlottesville, VA, USA
   index: 1
   ror:
date: 12 January 2026
bibliography: paper.bib
---

# Summary

CUT&Tag (Cleavage Under Targets and Tagmentation) is an epigenomic profiling method that maps protein-DNA interactions and histone modifications with high sensitivity and low background [@Kaya-Okur:2019]. Unlike traditional chromatin immunoprecipitation sequencing (ChIP-seq), CUT&Tag requires fewer cells, produces lower background signal, and does not require crosslinking or sonication. However, analyzing CUT&Tag data requires specialized computational workflows that account for its unique characteristics, including fragment size distributions, spike-in normalization, and optimized peak calling strategies.

**CAWS (CUT&Tag Analysis Workflow System)** is a comprehensive, automated Snakemake [@Molder:2021] pipeline designed specifically for CUT&Tag data analysis. The workflow takes raw paired-end or single-end FASTQ files as input and produces publication-ready results including quality control metrics, peak calls, coverage tracks, and interactive visualizations. CAWS integrates industry-standard bioinformatics tools (FastQC [@Andrews:2010], Trim Galore, Bowtie2 [@Langmead:2012], MACS3 [@Zhang:2008], SEACR [@Meers:2019]) into a unified framework that ensures reproducibility through Conda-based environment management and provides comprehensive quality assessment tailored to CUT&Tag experiments.

# Statement of Need

The rapid adoption of CUT&Tag has created a demand for standardized, reproducible analysis workflows. While general-purpose ChIP-seq pipelines exist (e.g., nf-core/cutandrun [@Ewels:2020]), they often require significant customization for CUT&Tag data and may not implement CUT&Tag-specific quality control metrics. Researchers analyzing CUT&Tag data face several challenges:

1. **Lack of CUT&Tag-specific QC**: Generic ChIP-seq pipelines do not assess CUT&Tag-specific quality metrics such as E. coli spike-in alignment rates for normalization [@Meers:2019].

2. **Inconsistent peak calling strategies**: CUT&Tag data benefits from specialized peak callers like SEACR [@Meers:2019], which uses sparse enrichment calling optimized for low-background data, but many existing workflows only implement MACS2/3.

3. **Mixed sequencing type support**: Modern experiments often combine single-end and paired-end sequencing data, but most pipelines assume homogeneous sequencing types across all samples.

4. **Reproducibility barriers**: Ad hoc analysis scripts are difficult to share, reproduce, and adapt across different computing environments.

CAWS addresses these gaps by providing:

- **CUT&Tag-optimized QC**: Specialized quality metrics including fragment size distributions, TSS enrichment scores, and E. coli spike-in normalization support.

- **Dual peak calling**: Implements both SEACR (optimized for CUT&Tag) and MACS3 (broadly compatible) to provide robust peak identification.

- **Flexible input handling**: Native support for paired-end, single-end, and mixed datasets within a single analysis, with automatic detection and appropriate processing for each sample type.

- **Reproducible execution**: Snakemake 9.x compatibility with modern executor plugins for HPC environments, Conda-based dependency management, and comprehensive documentation.

- **Interactive reporting**: Automated generation of HTML reports with Plotly-based visualizations for exploratory data analysis.

CAWS has been designed for research laboratories performing epigenomic profiling studies, particularly those investigating chromatin modifications, transcription factor binding, and chromatin-associated proteins. The pipeline reduces the barrier to rigorous CUT&Tag analysis and ensures consistent, reproducible results across different computing environments.

# Key Features

**Comprehensive workflow**: CAWS implements a complete analysis pipeline from raw FASTQ files to publication-ready outputs, including:

- Quality control and adapter trimming (FastQC, Trim Galore)
- Read alignment with quality filtering (Bowtie2)
- E. coli spike-in quantification for normalization
- PCR duplicate removal (optional)
- Peak calling using both SEACR and MACS3
- Fraction of Reads in Peaks (FRiP) calculation
- TSS enrichment analysis
- Signal heatmaps at peaks
- Interactive HTML reports with comprehensive QC metrics

**Dual-rule architecture for mixed datasets**: CAWS implements a sophisticated dual-rule architecture that automatically detects sequencing type from sample metadata and applies appropriate processing rules. Single-end and paired-end samples are processed with different SAM filtering flags and peak calling parameters, ensuring optimal results for each data type. This design allows researchers to analyze heterogeneous datasets (e.g., public data combined with in-house experiments) without manual intervention.

**HPC-ready execution**: The pipeline includes pre-configured SLURM profiles for high-performance computing environments, with resource requirements optimized for each rule. Snakemake 9.x executor plugin support ensures compatibility with modern cluster schedulers.

**Extensibility**: The modular Snakemake design allows researchers to easily customize processing steps, add new tools, or modify parameters for specific experimental designs.

# Comparison to Related Software

CAWS differs from existing pipelines in several key aspects:

- **nf-core/cutandrun** [@Ewels:2020]: A Nextflow-based pipeline for CUT&RUN/CUT&Tag. While comprehensive, it does not natively support mixed single-end/paired-end datasets and uses Nextflow rather than Snakemake, which may be less familiar to some computational biologists.

- **General ChIP-seq pipelines** (ENCODE ChIP-seq pipeline, ChIPseeker): Designed primarily for ChIP-seq data and lack CUT&Tag-specific optimizations such as SEACR peak calling and fragment size-based quality metrics.

- **SEACR** [@Meers:2019]: A standalone peak caller optimized for CUT&Tag but requires manual upstream processing (alignment, filtering, bedGraph generation). CAWS integrates SEACR into a complete automated workflow.

CAWS combines the CUT&Tag-specific features of specialized tools with the automation and reproducibility of a complete workflow system, while adding unique capabilities for mixed sequencing type support.

# Acknowledgements

We acknowledge contributions from the computational biology community who have developed the underlying tools integrated into CAWS. Development of CAWS was supported by [funding information to be added].

# References
