# nf-core/rnadnavar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
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

### `Changed`

- **Updated nf-schema plugin** from v2.4.2 to v2.5.0 for logo color support
- Cleaned up code formatting and style across configuration files
- Updated test configurations and `.nftignore` files
- Improved subworkflow organization and consistency (typos, remove unused code blocks, spacing, comments, readability)


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
