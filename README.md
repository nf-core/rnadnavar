<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-rnadnavar_logo_dark.png">
    <img alt="nf-core/rnadnavar" src="docs/images/nf-core-rnadnavar_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/rnadnavar)
[![GitHub Actions CI Status](https://github.com/nf-core/rnadnavar/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/rnadnavar/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnadnavar/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/rnadnavar/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnadnavar/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rnadnavar)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnadnavar-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnadnavar)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

The **nf-core/rnadnavar** is a bioinformatics best-practice
analysis pipeline for RNA somatic mutation detection
able to perform in parallel.

Initially designed for cancer research, the pipeline
uses different variant calling algorithms and applies a
consensus approach. A final filtering stage, should
provide a set of annotated somatic variants.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across
multiple compute infrastructures in a very portable
manner. It uses Docker/Singularity containers making
installation trivial and results highly reproducible.
The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this
pipeline uses one container per process which makes it
much easier to maintain and update software
dependencies. Where possible, these processes have been
submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

Depending on the options and samples provided, the
pipeline will run different tasks. This is controlled
mainly through `--step` and `--tools` parameters. This
is a summary of the possible tasks to run with the pipeline:

- Quality control and trimming (enabled by
  `--trim_fastq` and runs [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and
  [`fastp`](https://github.com/OpenGene/fastp))
- Map Reads to Reference (BWA-mem, BWA-mem2, dragmap
  and/or STAR)
- GATK preprocessing for DNA and RNA bulk sequencing
  samples (`GATK MarkDuplicates`, `GATK SplitNCigarReads`,`GATK
BaseRecalibrator` and `GATK ApplyBQSR`)
- Summarise alignment statistics (`samtools stats`, `mosdepth`)
- Variant calling (enabled with `--tools`)
  - `Mutect2`
  - `Strelka2`
  - `SAGE`
- Annotation with `VEP` (enabled with `--tools` adding
  `vep`)
- Normalisation of VCFs with VT (enabled with `--tools`
  adding `normalisation`)
- Transformation of VCF to MAF and consensus of variant
  calling results (enabled with `--tools` adding
  `consensus`)
- Filtering of MAF files applying optional gnomad,
  whitelisting and blacklisting (enabled with `--tools`
  adding `filtering`)
- Realignment step where reads from regions where a variant
  was found will be extracted and re-processed, only for
  RNA due to higher levels of background noise (enabled
  with `--tools` adding `realignment`).
- Filtering of MAF files specific for RNA (enabled with
  `--tools` adding `rna_filtering`)

<p align="center">
    <img title="RNADNAVAR Workflow"
src="docs/images/rnadnavar_schemav3.png">
</p>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Key requirements for the analysis

- RNA tumor samples must be associated with a DNA normal sample using the same patient ID. DNA tumour sample is optional.
- Sample types are specified with status: 0 (normal DNA), 1 (tumor DNA), 2 (tumor RNA)
- Reference files are required for the analysis.

### Quick Start

The typical command for running the pipeline is:

```bash
nextflow run nf-core/rnadnavar \
   -r <VERSION> \
   -profile <PROFILE> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --genome <GENOME> \
   --tools <TOOLS>
```

### Pipeline Steps

The pipeline can start from different entry points but `tools` need to be specified to run. If you are planning to run the pipeline by stages please read [usage](docs/usage.md#starting-from-different-steps).

### Input requirements

You will need to create a samplesheet with information about the samples you want to analyze. The samplesheet must be a comma-separated file with a header row.

Essential columns are:

- `patient` - Unique patient identifier
- `status` - Sample type (0/1/2)
- `sample` - Unique sample identifier
- And one of these options:
  - if variant calling is going to be performed of these three options should be an input (_note that this is included in the `realignment` step by default_):
    - `fastq_1`, `fastq_2` - Paths to paired-end FASTQ files (for raw data) (this also needs `lane` column)
    - `bam`,`bai` - Paths to BAM files and indices (for pre-aligned data)
    - `cram`,`crai` - Paths to CRAM files and indices (for pre-aligned data)
  - if the starting step is from `consensus` and no `realignment` will be performed (add `--skip_tools 'realignment'` to your command or params file) then the input are VCF/MAF:
    - `caller` - Caller used to generate that MAF (to annotate it during the consensus approach)
    - `maf` - Paths to MAF file
  - if `realignment` is desired then please add `normal_id` and take a look to [usage](docs/usage.md#starting-from-different-steps).

Examples CSVs:

Starting from `mapping`:

```csv
patient,sample,status,lane,fastq_1,fastq_2
PT1,DNA_NORMAL_SAMPLE,0,LX,DNA_NORMAL_SAMPLE_L002_R1_001.fastq.gz,DNA_NORMAL_SAMPLE1_L002_R2_001.fastq.gz
PT1,DNA_TUMOUR_SAMPLE,1,LX,DNA_TUMOUR_SAMPLE_L002_R1_001.fastq.gz,DNA_TUMOUR_SAMPLE_L002_R2_001.fastq.gz
PT1,RNA_TUMOUR_SAMPLE,2,LX,RNA_TUMOUR_SAMPLE_L002_R1_001.fastq.gz,RNA_TUMOUR_SAMPLE_L002_R2_001.fastq.gz
```

Starting from `consensus` no realignment:

```csv
patient,sample,status,caller,maf
PT1,RNA_TUMOUR_SAMPLE_1,2,mutect,DNA_TUMOUR_SAMPLE_1.mutect.maf
PT1,RNA_TUMOUR_SAMPLE_1,2,strelka,RNA_TUMOUR_SAMPLE_1.strelka.maf
PT1,RNA_TUMOUR_SAMPLE_2,2,mutect,DNA_TUMOUR_SAMPLE_2.mutect.maf
PT1,RNA_TUMOUR_SAMPLE_2,2,strelka,RNA_TUMOUR_SAMPLE_2.strelka.maf
```

### Parameters

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rnadnavar/usage) and the [parameter documentation](https://nf-co.re/rnadnavar/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnadnavar/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rnadnavar/output).

## Credits

The nf-core/rnadnavar was originally written by Raquel
Manzano Garcia at Cancer Research UK Cambridge Institute
with the continuous support of [Maxime U Garcia](https://github.com/maxulysse). The
workflow is based on
[RNA-MuTect](https://github.com/broadinstitute/RNA_MUTECT_1.0-1) which was
originally published by [Yizhak, _et al_ 2019 (Science)](https://www.science.org/doi/10.1126/science.aaw0726)

We thank the following people for their assistance in the development of this pipeline:
TBC

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnadnavar` channel](https://nfcore.slack.com/channels/rnadnavar) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/rnadnavar for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
