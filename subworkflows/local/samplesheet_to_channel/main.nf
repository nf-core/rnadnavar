//
// Function to parse samplesheet to appropriate channel structure
//

workflow  SAMPLESHEET_TO_CHANNEL{

    take:
    ch_from_samplesheet

    main:
    ch_from_samplesheet.dump(tag:"ch_from_samplesheet")
    input_sample = ch_from_samplesheet
        .map{ meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf ->
            // generate patient_sample key to group lanes together
            [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf] ]
        }
        .tap{ ch_with_patient_sample } // save the channel
        .groupTuple() //group by patient_sample to get all lanes
        .map { patient_sample, ch_items ->
            // get number of lanes per sample
            [ patient_sample, ch_items.size() ]
        }
        .combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
        .map { patient_sample, num_lanes, ch_items ->

            (meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf) = ch_items
            if (meta.lane && fastq_2) {
                meta           = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN         = params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell   = flowcellLaneFromFastq(fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${meta.sample}.${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
                if (meta.status >= 2) { // STAR does not need '@RG'
                    read_group  = "ID:${flowcell}.${meta.sample}.${meta.lane} ${CN}PU:${meta.lane} SM:${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
                }
                meta           = meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'fastq', size: 1]

                if (params.step == 'mapping') return [ meta, [ fastq_1, fastq_2 ] ]
                else {
                    error("Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }
            // start for realignment or will do realignment later starting after pre-processing
            } else if ((maf || vcf) && (params.step=="realignment" || (params.tools && params.tools.split(',').contains("realignment")))){
                if (meta.lane == null) meta.lane = "LX"
                meta            = meta + [id: "${meta.sample}-${meta.lane}-realign".toString()]
                def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}_realign\\t${CN}PU:${meta.lane}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
                if (meta.status >= 2) { // STAR does not need '@RG'
                    read_group  = "ID:${meta.sample}_${meta.lane}_realign ${CN}PU:${meta.lane} SM:${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
                }
                if (meta.status >= 2 || meta.status==0){ // these are the files that will go through realignment
                    if (cram)  return [ meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'cram', size: 1], cram, crai, maf ]
                    else if (bam) return [ meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1], bam, bai, maf ]
                    else {
                        error("Combination error")}
                } else if (meta.status >= 1){  // either DNA (status=1) or RNA (status=2)
                    if (meta.normal_id == null){
                        error("When tool 'realigment' enabled, `normal_id` should be added to the csv for maf files.")
                        // Need to get the normal id to create the tumour_vs_normal id
                    } else {
                        meta = meta + [id: meta.sample + "_vs_" + meta.normal_id, data_type: 'maf', variantcaller: variantcaller ?: '']
                        return [ meta, maf ]
                    }
                }
            // start from BAM
            } else if (meta.lane && bam) {
                if (params.step != 'mapping' && !bai) {
                    error("BAM index (bai) should be provided.")
                }
                meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
                if (meta.status >= 2) { // STAR does not need '@RG'
                    read_group  = "ID:${meta.sample}_${meta.lane} ${CN}PU:${meta.lane} SM:${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
                }
                meta            = meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

                if (params.step != 'annotate') return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // recalibration
            } else if (table && cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            } else if (params.step == "recalibrate" && !table) {

                error("Step is `$params.step` but samplesheet has no recalibration table. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")

            }else if (table && bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
                else {
                    error("Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration or variant_calling
            } else if (cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
            } else if (bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']
                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is 2 `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // annotation
            } else if (vcf) {
                meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

                if (params.step == 'annotate' ) return [ meta - meta.subMap('lane'), vcf ]
                else if (params.step == 'norm') {
                    if (meta.status == 0){ // TODO: more specific checks on this is needed
                        error("Samplesheet contains vcf files with status 0, vcfs should only be for tumours (1|2).")
                    }
                    else if (meta.normal_id == null){
                        error("When step 'norm', `normal_id` should be added to the csv for all samples.")
                        // Need to get the normal id to create the tumour_vs_normal id
                    } else {
                        meta = meta + [id: meta.sample + "_vs_" + meta.normal_id, data_type: 'vcf', variantcaller: variantcaller ?: '']
                        return [ meta - meta.subMap('lane'), vcf ]
                    }
                }
                else {
                    error("Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }
            } else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }
    input_sample.dump(tag:"input_sample")

    if (params.step != 'annotate' && params.tools && !params.build_only_index) {
        // Two checks for ensuring that the pipeline stops with a meaningful error message if
        // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
        // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
        input_sample.filter{ it[0].status == 1 | it[0].status == 2 }.ifEmpty{ // In this case, the sample-sheet contains no tumor-samples
            if (!params.build_only_index) {
                def tools_tumor = ['mutect2', 'sage', 'strelka']
                def tools_tumor_asked = []
                tools_tumor.each{ tool ->
                    if (params.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
                }
                if (!tools_tumor_asked.isEmpty()) {
                    error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
                }
            }
        }
        input_sample.filter{ it[0].status == 0 }.ifEmpty{ // In this case, the sample-sheet contains no normal/germline-samples
            def tools_requiring_normal_samples = ['mutect', 'sage', 'strelka']
            def requested_tools_requiring_normal_samples = []
            tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
                if (params.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
            }
            if (!requested_tools_requiring_normal_samples.isEmpty()) {
                error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample at the moment: ' + requested_tools_requiring_normal_samples.join(", "))
            }
        }
    }

    // Fails when wrongfull extension for intervals file
    if (params.wes && !params.step == 'annotate') {
        if (params.intervals && !params.intervals.endsWith("bed"))  error("Target file specified with `--intervals` must be in BED format for targeted data")
        else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
    } else if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("list")) error("Intervals file must end with .bed, .list, or .interval_list")

    if (params.step == 'mapping' && params.aligner.contains("dragmap") && !(params.skip_tools && params.skip_tools.split(',').contains("baserecalibrator"))) {
        log.warn("DragMap was specified as aligner. Base recalibration is not contained in --skip_tools. It is recommended to skip baserecalibration when using DragMap\nhttps://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode")
    }


    // Warns when missing files or params for mutect2
    if (params.tools && params.tools.split(',').contains('mutect2')) {
        if (!params.pon) {
            log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
        }
        if (!params.germline_resource) {
            log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
        }
        if (params.pon && params.pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")) {
            log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
        }
    }

    // Fails when missing resources for baserecalibrator
    if (!params.dbsnp && !params.known_indels) {
        if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))) {
            error("Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command.")
        }
    }

    // Fails when missing tools for variant_calling or annotate
    if ((params.step == 'variant_calling' || params.step == 'annotate') && !params.tools) {
        error("Please specify at least one tool when using `--step ${params.step}`.\nhttps://nf-co.re/rnadnavar/parameters#tools")
    }


    if ((params.download_cache) && (params.vep_cache)) {
        error("Please specify either `--download_cache` or `--vep_cache`.\nhttps://nf-co.re/rnadnavar/usage#how-to-customise-snpeff-and-vep-annotation")
    }

    emit:
    input_sample
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}
