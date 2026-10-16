# nf-core/rnadnavar: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarizes results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Quality Control](#quality-control) - Raw read QC and trimming, FastQC and fastp reports for input data assessment
- [Alignment](#alignment) - Read alignment to reference genome (BAM/CRAM)
- [Preprocessing](#preprocessing) - GATK preprocessing steps following best practices
- [Variant Calling](#variant-calling) - Somatic variant detection, VCF files from multiple callers (Mutect2, Strelka2, SAGE).
- [Annotation](#annotation) - Variant annotation with VEP
- [Consensus](#consensus) - Consensus variant calling across multiple callersin MAF format
- [Filtering](#filtering) - Variant filtering and quality control in MAF format
- [Realignment](#realignment) - Optional RNA-specific realignment
- [MultiQC](#multiqc) - Aggregate report describing results and QC
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Quality Control

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastp.html`: fastp report in HTML format.
  - `*.fastp.json`: fastp report in JSON format.
  - `*.fastp.log`: fastp log file.
  - `*.fastp.fastq.gz`: Trimmed FastQ files (if trimming is enabled).

</details>

[fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast, all-in-one preprocessing for FastQ files. It is developed in C++ with multithreading supported to afford high performance. fastp is used in this pipeline for quality control and adapter trimming when `--trim_fastq` is specified.

## Alignment

### BWA-MEM / BWA-MEM2 / STAR

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `*.bam`: Aligned reads in BAM format.
  - `*.bai`: BAM index files.
  - `*.cram`: Aligned reads in CRAM format (if `--save_output_as_bam` is false).
  - `*.crai`: CRAM index files.

</details>

The pipeline supports multiple aligners:

- [BWA-MEM](https://github.com/lh3/bwa) for DNA alignment
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) for faster DNA alignment
- [STAR](https://github.com/alexdobin/STAR) for RNA alignment
- [DRAGMAP](https://github.com/Illumina/DRAGMAP) for high-accuracy alignment

Alignment files are saved in BAM or CRAM format depending on the `--save_output_as_bam` parameter.

### Alignment QC

<details markdown="1">
<summary>Output files</summary>

- `alignment/samtools_stats/`
  - `*.stats`: Samtools stats output with alignment statistics.
- `alignment/mosdepth/`
  - `*.mosdepth.global.dist.txt`: Global coverage distribution.
  - `*.mosdepth.summary.txt`: Coverage summary statistics.

</details>

Quality control metrics for aligned reads are generated using [samtools stats](http://www.htslib.org/doc/samtools-stats.html) and [mosdepth](https://github.com/brentp/mosdepth).

## Preprocessing

### GATK Preprocessing

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/`
  - `markduplicates/`
    - `*.bam`: BAM files with duplicates marked.
    - `*.bai`: BAM index files.
    - `*.metrics`: Duplicate marking metrics.
  - `splitncigarreads/` (RNA samples only)
    - `*.bam`: BAM files with split N CIGAR reads.
    - `*.bai`: BAM index files.
  - `baserecalibrator/`
    - `*.table`: Base quality score recalibration tables.
  - `applybqsr/`
    - `*.bam`: BAM files with recalibrated base quality scores.
    - `*.bai`: BAM index files.

</details>

The pipeline follows [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) for preprocessing:

1. **MarkDuplicates**: Identifies and marks duplicate reads
2. **SplitNCigarReads**: Splits reads with N CIGAR operations (RNA-seq specific)
3. **BaseRecalibrator**: Generates base quality score recalibration table
4. **ApplyBQSR**: Applies base quality score recalibration

## Variant Calling

### Mutect2

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/mutect2/`
  - `*.vcf.gz`: Variant calls in VCF format.
  - `*.vcf.gz.tbi`: VCF index files.
  - `*.stats`: Mutect2 statistics.
  - `*.f1r2.tar.gz`: F1R2 counts for orientation bias filtering.

</details>

[Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) is GATK's somatic variant caller for detecting SNVs and indels in tumor samples.

### Strelka2

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/strelka/`
  - `*.somatic.snvs.vcf.gz`: Somatic SNV calls.
  - `*.somatic.indels.vcf.gz`: Somatic indel calls.
  - `*.vcf.gz.tbi`: VCF index files.

</details>

[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs.

### SAGE

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/sage/`
  - `*.vcf.gz`: Variant calls in VCF format.
  - `*.vcf.gz.tbi`: VCF index files.

</details>

[SAGE](https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md) is a precise and highly sensitive somatic SNV, MNV and INDEL caller.

## Annotation

### VEP Annotation

<details markdown="1">
<summary>Output files</summary>

- `annotation/vep/`
  - `*.vcf.gz`: Annotated variants in VCF format.
  - `*.vcf.gz.tbi`: VCF index files.
  - `*.html`: VEP summary report.

</details>

[Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) determines the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.

### VCF to MAF Conversion

<details markdown="1">
<summary>Output files</summary>

- `annotation/vcf2maf/`
  - `*.maf`: Variants in Mutation Annotation Format (MAF).

</details>

Variants are converted from VCF to [MAF format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) using [vcf2maf](https://github.com/mskcc/vcf2maf) for downstream analysis and visualization.

## Consensus

### Consensus Calling

<details markdown="1">
<summary>Output files</summary>

- `consensus/`
  - `*.consensus.vcf`: Consensus variants across all callers in VCF format.
  - `*.consensus.maf`: Consensus variants in MAF format.
  - `*.consensus_*.vcf`: Caller-specific consensus variants.
  - `*.consensus_*.maf`: Caller-specific consensus variants in MAF format.
  - `*.pdf`: Consensus analysis plots and visualizations.

</details>

The consensus module combines variant calls from multiple callers (Mutect2, Strelka2, SAGE) to improve variant calling accuracy and reduce false positives.

## Filtering

### MAF Filtering

<details markdown="1">
<summary>Output files</summary>

- `filtering/`
  - `*.filtered.maf`: Filtered variants in MAF format.

</details>

Variants are filtered based on various criteria including:

- Population frequency (gnomAD) (controlled by config of `maf_filtering` module)
- Whitelist/blacklist regions (controlled by parameters `whitelist` and `blacklist`)
- Quality metrics
- Custom filtering parameters (controlled by config of `maf_filtering` module)

### RNA-specific Filtering

<details markdown="1">
<summary>Output files</summary>

- `filtering/rna/`
  - `*.rna_filtered.maf`: RNA-specific filtered variants in MAF format.

</details>

Additional filtering steps specific to RNA-seq data to account for RNA editing, splicing artifacts, and other RNA-specific noise.

## Realignment

### RNA Realignment

<details markdown="1">
<summary>Output files</summary>

- `realignment/`
  - `*.realigned.bam`: Realigned BAM files for variant regions.
  - `*.realigned.bai`: BAM index files.
  - `*.realigned.maf`: Variants from realigned regions.

</details>

Optional realignment step using [HISAT2](https://daehwankimlab.github.io/hisat2/) for RNA samples in regions where variants were detected. This helps improve variant calling accuracy in RNA-seq data by addressing alignment artifacts.

## Normalization

### VT Normalization

<details markdown="1">
<summary>Output files</summary>

- `normalization/`
  - `*.normalized.vcf.gz`: Normalized variants in VCF format.
  - `*.stats`: Normalization statistics.

</details>

Variants are normalized using [vt](https://genome.sph.umich.edu/wiki/Vt) decompose and normalize to ensure consistent representation of variants across different callers.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`: Report describing the pipeline run.
  - `execution_timeline_*.html`: Timeline of the pipeline execution.
  - `execution_trace_*.txt`: Trace file with detailed execution information.
  - `pipeline_dag_*.html`: DAG visualization of the pipeline.
  - `software_versions.yml`: Software versions used in the pipeline.
  - `params_*.json`: Parameters used for the pipeline run.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
