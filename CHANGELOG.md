# nf-core/rnadnavar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [2026-06-23]

### `Added`

- Added support for optional Mutect2 force-calling inputs (`mutect2_alleles` and `mutect2_alleles_tbi`)
- Added `conf/empty.config` to support strict-syntax-safe config loading

### `Fixed`

- Fixed multiple strict-syntax compilation issues across local workflows and subworkflows
- Fixed samplesheet parsing and moved tumour/normal composition validation upstream of channel creation
- Fixed interval-preparation strict-syntax issues, including name collisions and duration handling
- Fixed local wrappers and callers to match updated nf-core module input/output signatures and version emits
- Fixed reference-channel consumption bugs affecting FASTA/FAI/DICT propagation in preprocessing, realignment and variant-calling paths
- Fixed realignment-specific logic in RNA filtering and downstream consensus handling
- Fixed MultiQC input/config channel handling and mixed software-version aggregation
- Fixed consensus input ordering to improve deterministic behaviour when resuming runs
- Fixed non-deterministic nf-test outputs and updated snapshots accordingly
- Fixed `includeConfig` handling in `nextflow.config` for newer Nextflow strict syntax
- Fixed VEP cache initialisation and updated unzip-dependent wiring

### `Changed`

- Updated the nf-core template and a broad set of nf-core modules and subworkflows
- Harmonised FASTA/FAI/DICT channel shapes across local subworkflows
- Updated local realignment and HISAT2 call wiring for newer module interfaces
- Updated local samtools callers (`view`, `convert`, `merge`, `faidx`) to current nf-core module contracts
- Updated local GATK and Picard integrations, including restored patches for `picard/filtersamreads` and `gatk4/splitncigarreads`
- Replaced deprecated tabix usage with current htslib-based handling where appropriate
- Updated nf-test plugin configuration and test snapshots

### `Removed`

- Removed obsolete tabix modules and deprecated module usage paths
- Removed leftover config and workflow parameters no longer used after the template/module refresh (for example `hook_url`)

---

## v1.0dev - [2025-08-12]

### `Added`

- Added full pipeline nf-test (`test_full.nf.test`)
- Added test snapshots
- Restored rnadnavar logo in pipeline output with color support

### `Fixed`

- Fixed help text and documentation URLs
- Fixed pipeline-specific parameter descriptions in schema
- Fixed workflow path references and documentation links
- Fixed consensus module caller value to avoid warnings
- Fixed VEP cache handling and empty RNA edits processing
- Fixed issue with run_consensus.R plotting

### `Changed`

- Cleaned up code formatting and style across configuration files
- Updated test configurations and `.nftignore` files
- Improved subworkflow organization and consistency (typos, remove unused code blocks, spacing, comments, readability)
- Massive speed up to run_consensus.R
- Migrated consensus module from custom Docker container to Seqera Wave containers

### `Removed`

- Cleaned up unused VCFlib and VT variant processing modules
- Removed obsolete module configurations and test files

---

## v1.0dev - [2025-08-01]

### `Added`

- Added VCF simple test data for annotation testing
- Added comprehensive nf-test framework integration
- Added GitHub Actions for automated testing with nf-test
- Added support for multiple template versions (3.2.0, 3.2.1, 3.3.1)

### `Fixed`

- **Major Fix**: Resolved ConcurrentModificationException error from Java processes
- Fixed local module implementations and configurations
- Fixed subworkflow naming and structure issues
- Fixed indentation and formatting issues across multiple files
- Fixed redundancies in workflow logic
- Fixed SAGE variant caller integration and configuration
- Fixed local utilities and helper functions
- Fixed nf-core subworkflow integrations

### `Changed`

- **Massive Module Update**: Updated all nf-core modules to latest versions including:
  - bcftools/sort and bcftools/stats
  - BWA, BWA-MEM2, and DRAGMAP aligners
  - GATK4 modules (applybqsr, baserecalibrator)
  - FastQC, FastP, and other QC tools
  - VEP annotation modules
  - All associated test files and metadata
- Updated VEP cache version parameter from string '110' to integer 110
- Cleaned up main workflow structure and organization
- Cleaned up local subworkflows and removed redundant code
- Cleaned up configuration files and parameter definitions
- Updated template to nf-core/tools version 3.3.1
- Reorganized and cleaned up local modules
- Updated GitHub workflows and CI/CD configurations
- Added new configuration parameters: `help_full`, `show_hidden`, `trace_report_suffix`, `modules_testdata_base_path`
- Improved parameter organization and formatting in nextflow.config

### `Dependencies`

- Updated all nf-core modules to their latest stable versions
- Updated container images and Conda environments for all processes
- **Updated Nextflow minimum version requirement from >=23.04.0 to >=24.04.2**
- Updated nf-core template to version 3.3.1

### `Deprecated`

- Removed redundant workflow components
- Cleaned up deprecated configuration parameters
- Removed obsolete test configurations

---

## Initial Release

Initial release of nf-core/rnadnavar, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- RNA and DNA integrated analysis pipeline for somatic mutation detection
- Support for multiple variant callers (Mutect2, Strelka2, SAGE)
- Comprehensive preprocessing with GATK4 best practices
- VEP annotation and filtering capabilities
- Consensus variant calling approach
- RNA-specific filtering and realignment steps
- MultiQC reporting and quality control
- Support for both BWA-MEM and STAR alignment
