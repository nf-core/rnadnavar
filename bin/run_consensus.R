#!/usr/bin/env Rscript
# Date: Sun 20 Sep 2020
# Author: Raquel Manzano - @RaqManzano
# Script: Find overlaps between vcf - specific filename convention will be expected to extract the variant caller name
# Expected filename pattern: <sample>_<caller>.maf or <sample>.<caller>.maf (e.g. S001_mutect2.maf or S001.mutect2.maf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
options(warn=-1)

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stringr))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(datatable.fread.input.cmd.message=FALSE)
# Collect arguments
script_args <- commandArgs(TRUE)
# default arguments
if(length(script_args) < 1) {
    script_args <- c("--help")
}

# Help section ----
if("--help" %in% script_args) {
    cat("The consensusR script:

        Arguments:
        --id=sample_id              - character, id for the sample (default: sample)
        --input_dir=input_directory - character, Directory containing input files
        --out_prefix=output_prefix  - character, prefix for outputs
        --thr=thr                   - integer, the number of callers that support a mutation to be called consensus (default: 2)
        --interval=coord_interval   - character, chrom:start-end to analyse
        --help                      - print this text
        Note: Please ensure caller name is as follows: sample_id.caller.maf
        Example:
        ./run_consensus.R --id=sample01 --input_dir=inputs/  --out_prefix=consensus --thr=2 \n\n")

    q(save="no")
}

# Dealing with missing params ----
if(!grepl(pattern = "--input_dir=", x = paste(script_args, collapse = ""))) {
    stop("Missing --input_dir name!")
}
if(!grepl(pattern = "--id=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--id=sample")
}
if(!grepl(pattern = "--out_prefix=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--out_prefix=consensus")
}
if(!grepl(pattern = "--thr=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--thr=2")
}
if(!grepl(pattern = "--interval=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--interval=0")
}

## Parse arguments (we expect the form --arg=value) ----
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsDF       <- as.data.frame(do.call("rbind", parseArgs(script_args)))
argsDF       <- aggregate(.~V1, data=argsDF, paste, collapse=",")
argsL        <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

# Input variables ----
sampleid    <- argsL$id. # just used for the title of the plot at the end
input_files <- list.files(strsplit(argsL[["input_dir"]], split=",")[[1]], pattern="\\.maf$", full.names=TRUE)

# Extract caller names from file names - specific name convention is expected ----
get_caller <- function(filename) {
    caller_name <- sub(".*(\\.|_)(.*?)\\.maf", "\\2", filename)
    if (caller_name == filename) caller_name <- "unknown"
    return(caller_name)
}
callers <- sapply(input_files, get_caller)
names(input_files) <- callers

# check if input is a vcf or maf ----
is.vcf <- grepl(x = input_files[1], pattern = ".vcf$|.vcf.gz$", perl = T)

# init output files ----
pdf.out <- paste0(argsL$out_prefix, ".pdf")
if (is.vcf){
    vcf.out <- paste0(argsL$out_prefix, ".vcf")
} else{
    vcf.out <- paste0(argsL$out_prefix, ".maf")
}

message("- Parsing headers")
if (is.vcf){
    # First, get contigs from one of the input files
    contigs_meta <- fread(cmd=paste0("awk '/^##contig/ {print; next} !/^##contig/ {exit}' ", input_files[1]), sep = NULL, header = F)
    contigs_meta <- paste0(contigs_meta$V1, collapse = "\n")
} else {
    # this is not a vcf so no contigs
    contigs_meta <- ""
}
# First, collect vcf headers for the outputs ----
callers_meta <- list()
for ( c in callers){
    vcf_meta <- fread(cmd=paste0("awk '/^#CHROM/ || /^Hugo_Symbol/ {exit} {print}' ", input_files[c]), sep = NULL, header = F)
    if (is.vcf){
        vcf_header <- fread(cmd=paste0("zgrep -m 1 '#CHROM' ", input_files[c]), sep = NULL, header = F)
        vcf_header <- strsplit(vcf_header$V1, "\t")[[1]]
        callers_meta[[c]] <-  list(meta=paste0(vcf_meta$V1, collapse = "\n"),
                                header=vcf_header)
    } else{
        maf_header <- fread(cmd=paste0("zgrep -m 1 'Hugo_Symbol' ", input_files[c]), sep = NULL, header = F, showProgress = F)
        maf_header <- strsplit(maf_header$V1, "\t")[[1]]
        callers_meta[[c]] <-  list(meta=paste0(vcf_meta$V1, collapse = "\n"),
                                    header=maf_header)
    }
}
# Second, read and convert the input files to GenomicRanges for easy manipulation.
message("- Converting to genomic ranges")
mutsGR <- list()
muts <-  list()

## if there are intervals parse now
if (argsL$interval != '0') {

    message( ' - Subset input by interval: ', argsL$interval)
    interval_chr   = strsplit(argsL$interval, ':')[[1]][1]
    interval_start = as.numeric(strsplit(strsplit(argsL$interval, ':')[[1]][2], '-')[[1]][1])
    interval_end   = as.numeric(strsplit(strsplit(argsL$interval, ':')[[1]][2], '-')[[1]][2])

}

for(c in callers[1:length(callers)]){
    v <- input_files[c]
    message("  - ", v)
    if (!is.vcf){ # is a maf
        fread_cmd =  paste0("zgrep -v '#' ", v)
        if (argsL$interval != '0') {
            fread_cmd <- paste0(fread_cmd, "| grep -Ew 'Hugo_Symbol|", interval_chr, "'")
        }
        tmp <- fread(cmd=fread_cmd, header=TRUE)
        # tmp <- tmp[Chromosome=="chr21"]
        if (argsL$interval != '0') {
            tmp <- tmp[Chromosome == interval_chr & Start_Position >= interval_start & Start_Position < interval_end]
        }
        tmp$`#CHROM` <- tmp$Chromosome
        tmp$POS      <- tmp$Start_Position
        tmp$REF      <- tmp$Reference_Allele
        tmp$ALT      <- tmp$Tumor_Seq_Allele2
    } else { # is a vcf
        fread_cmd =  paste0("zgrep -v '##' ", v)
        if (argsL$interval != '0') {
            fread_cmd <- paste0(fread_cmd, "| grep -w ", interval_chr)
        }
        tmp <- fread(cmd=fread_cmd)
    }
    if (nrow(tmp) > 0) { # if not empty
        tmp$Caller    <- c
        tmp$mut       <- paste0(tmp$REF, ">", tmp$ALT)
        tmp$DNAchange <- paste0(tmp$`#CHROM`, ":g.", tmp$POS, tmp$REF, ">", tmp$ALT)
        tmp$start     <- tmp$POS
        # get end position to create a range
        tmp$end <- ifelse(nchar(tmp$REF) > nchar(tmp$ALT),
                        tmp$POS + nchar(tmp$REF) - 1,
                        ifelse(nchar(tmp$REF) < nchar(tmp$ALT), tmp$POS + nchar(tmp$ALT) - 1,
                                ifelse(nchar(tmp$REF) == nchar(tmp$ALT) & nchar(tmp$ALT) > 1, tmp$POS + nchar(tmp$REF),
                                        tmp$POS)
                        )
        )
        if (nrow(tmp[is.na(end)])>0){
            message("   - [WARNING] Removing ",nrow(tmp[is.na(end)]), " spurious calls with no alt. There should be 0.")
        }
        tmp <- tmp[!is.na(end)]
        muts[[c]]   <- tmp
        mutsGR[[c]] <- GenomicRanges::makeGRangesFromDataFrame(df = tmp,
                                                        ignore.strand = TRUE,
                                                        start.field = "start",
                                                        end.field = "end",
                                                        seqnames.field = "#CHROM",
                                                        keep.extra.columns = T)
    } else { # if empty add empty genomic range
        mutsGR[[c]] <- GenomicRanges::GRanges()  # empty
    }
    rm(tmp); gc(verbose = FALSE)
}

# Extract variants where their coordinates overlap between inputs ----
message("- Finding overlaps")
# Third, we find overlaps
# Split each GRange by chromosome to speed things up
mutsGR_by_chr <- lapply(mutsGR, function(gr) {
    if (!methods::is(gr, "GRanges") || length(gr) == 0) return(list())
    split(gr, seqnames(gr))
})

rm(mutsGR); gc(verbose = FALSE)
callers.names.GR <- names(mutsGR_by_chr)

overlapping.vars.list <- list()
for (c1 in callers.names.GR){
    gr1_chr <- mutsGR_by_chr[[c1]]
    if (length(gr1_chr) == 0) next  # if empty jump to next
    for (c2 in callers.names.GR){
        if (c1 == c2) next  # skip self-comparison
        gr2_chr <- mutsGR_by_chr[[c2]]
        if (length(gr2_chr) == 0) next  # if (length(gr1_chr) == 0) next
        # Only consider chromosomes both callers share
        chroms <- intersect(names(gr1_chr), names(gr2_chr))
        if (length(chroms) == 0) next
        for (chrom in chroms) {
            gr1 <- gr1_chr[[chrom]]
            gr2 <- gr2_chr[[chrom]]
            if (length(gr1) == 0 || length(gr2) == 0) next
            # The gap between 2 adjacent ranges is 0.
            hits <- GenomicRanges::findOverlaps(query = gr1, subject = gr2, maxgap = 0)
            n_hits <- length(hits)
            message(sprintf("[%s vs %s | %s] %d overlaps found", c1, c2, chrom, n_hits))
            if (n_hits == 0) next  # No overlaps found
            dnachange.hits <- muts[[c1]][Chromosome==chrom][queryHits(hits)]$DNAchange
            filt.hits      <- muts[[c1]][Chromosome==chrom][queryHits(hits)]$FILTER
            group_name <- paste0(c1, "_vs_", c2, "_", chrom)
            message("group name: ", group_name)
            # due to normalization we might find the same variant with different filters - these come from homopolymer regions
            overlapping.vars.list[[group_name]] <- unique(data.frame(DNAchange = dnachange.hits, caller = c1, FILTER = filt.hits))
        }
    }
}
overlapping.vars <- data.table::rbindlist(overlapping.vars.list, fill = TRUE)
rm(overlapping.vars.list, mutsGR_by_chr); gc(verbose = FALSE)
# Finally, extract the set of variants that will be the consensus set
overlapping.vars <- overlapping.vars[!duplicated(overlapping.vars),]
# Convert to data.table and summarise by DNAchange
setDT(overlapping.vars)
overlapping.vars <- overlapping.vars[, .(
    caller = paste(caller, collapse = "|"),
    FILTER = paste(FILTER, collapse = "|")
    ),
    by = DNAchange
]


# Overlaps that are adjacent, and only SNVs are removed if not DNP
overlapping.variants.count        <- stringr::str_count(string = overlapping.vars$caller, pattern =  stringr::fixed("|")) + 1
names(overlapping.variants.count) <- overlapping.vars$DNAchange
overlapping.variants.count.snvs   <- overlapping.variants.count[ grepl(pattern = "[0-9](A|C|G|T)>(A|C|G|T)$", x = names(overlapping.variants.count))]
overlapping.variants.count.indels <- overlapping.variants.count[!grepl(pattern = "[0-9](A|C|G|T)>(A|C|G|T)$", x = names(overlapping.variants.count))]

# only snvs with exact match will be in the consensus list
con.vars.ths.snv <- names(overlapping.variants.count.snvs[overlapping.variants.count.snvs >= argsL$thr])
# we let all indels pass as they have shown overlap
con.vars.ths.indel <- names(overlapping.variants.count.indels)

message("- There are ", prettyNum(length(con.vars.ths.snv), big.mark = ','),   " SNVs that are consensus")
message("- There are ", prettyNum(length(con.vars.ths.indel), big.mark = ','), " indels that are consensus")

con.vars.ths <- c(con.vars.ths.snv, con.vars.ths.indel)

# The next steps are for the output
# To keep the information from the consensus we extract the callers that called each mutation and its correspondent filters.
consensus_map <- overlapping.vars[, .(
                                        callers_all = paste(caller, collapse = "|"),
                                        filters_all = paste(FILTER, collapse = "|")
                                        ),
                                    by = DNAchange
    ]
# Mark which DNAchanges are consensus
consensus_map[, is_consensus := DNAchange %in% con.vars.ths]
setkey(consensus_map, DNAchange)

for (c in callers){
    message("- Annotating calls from ", c)
    if (!is.null(muts[[c]])){
        muts_dt <- as.data.table(muts[[c]])
        setkey(muts_dt, DNAchange)
        # Join consensus information
        muts_dt <- consensus_map[muts_dt, on = "DNAchange"]

        # Replace missing with self information
        muts_dt[is.na(callers_all), callers_all := Caller]
        muts_dt[is.na(filters_all), filters_all := FILTER]

        # Rename for consistency with rest of script
        setnames(muts_dt,
                old = c("callers_all", "filters_all"),
                new = c("callers", "filters"))

        muts[[c]] <- as.data.frame(muts_dt)
        rm(muts_dt); gc(FALSE)
    }
}

all.muts <- do.call(rbind.fill, muts)

# Remove duplication if consensus input came annotated with more than one caller
all.consensus.muts <- all.muts[stringi::stri_detect_regex(all.muts$Caller, "consensus", case_insensitive=TRUE),]
all.consensus.muts <- all.consensus.muts[!duplicated(all.consensus.muts$DNAchange),]
# 1) remove the consensus variants that might be duplicated 2) put back the deduplicated consensus variants
all.muts <-  all.muts[!stringi::stri_detect_regex(all.muts$Caller, "consensus", case_insensitive=TRUE),]
all.muts <- rbind(all.muts, all.consensus.muts)

message("- Preparing output")
## Prepare output
simplified.filter <- sapply(all.muts$filters, FUN = function(x){
    filt.val <- strsplit(x=x, split = "|", fixed = T)[[1]]
    filt.t   <- table(filt.val == "PASS")
    ifelse(prop.table(filt.t)["FALSE"] > 0.5, "FAIL", "PASS")
})
simplified.filter <- ifelse(is.na(simplified.filter), "PASS", simplified.filter)
all.muts$FILTER_consensus <- simplified.filter
all.muts$INFO_consensus   <-  paste0("callers=", all.muts$callers, ";filters=", all.muts$filters, ";consensus_filter=", all.muts$FILTER_consensus)

## write consensus
# I want to keep the info without duplicating the mutations
all.muts$isconsensus <- grepl(pattern = "|", x = all.muts$callers, fixed = T)

# WRITE OUTPUTS

meta_consensus <- paste0('##INFO=<ID=callers,Number=1,Type=String,Description="Variant callers that called this mutation, separated by |">\n',
                        '##INFO=<ID=filters,Number=1,Type=String,Description="Filters provided by each variant caller, separated by |">\n',
                        '##INFO=<ID=consensus_filter,Number=1,Type=String,Description="PASS if 50% or more of the callers give the mutation, otherwise FAIL.">')
extra.cols <- c()
for ( c in callers){
    if (is.vcf){
        updated_meta <-  paste(callers_meta[[c]]$meta,
                        meta_consensus,
                        sep="\n")
    vcf.out.caller <- paste0(argsL$out_prefix, "_", c,".vcf")
    } else{
        if ("Caller" %in% callers_meta[[c]]$header){
            extra.cols <- c()
        } else{
            extra.cols <- c("Caller", "callers", "filters", "FILTER_consensus", "isconsensus")
        }
        callers_meta[[c]]$header <- c(callers_meta[[c]]$header, extra.cols)
        updated_meta <- "#version 2.4" ## is a maf file
        vcf.out.caller <- paste0(argsL$out_prefix, "_", c, ".maf")
    }
    write(x = updated_meta, file = vcf.out.caller, ncolumns = 1, append = F)
    # if empty maf/vcf just write empty file with columns
    if(nrow(all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header])==0){
        empty.df <- data.frame(matrix(ncol = length(callers_meta[[c]]$header), nrow = 0))
        fwrite(x = empty.df,
            file = vcf.out.caller,
            append = T,
            sep = "\t",
            col.names = T)
        next
    }
    extra.cols <- c()
    if (is.vcf){
        fields_to_write <- all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header]
        fields_to_write$INFO <- paste(fields_to_write$INFO, all.muts[all.muts$Caller==c,]$INFO_consensus, sep=";")
        fwrite(x = fields_to_write,
            file = vcf.out.caller,
            append = T,
            sep = "\t",
            col.names = T)
    } else{
        colnames(all.muts)[colnames(all.muts)=="FILTER_consensus"] <- "FILTER_consensus"
        if (all(all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header] == callers_meta[[c]]$header)){
            to_write <- setNames(data.table(matrix(nrow = 0, ncol = length(callers_meta[[c]]$header))), callers_meta[[c]]$header)
        } else{
            to_write <- all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header]
        }
        fwrite(x = to_write,
            file = vcf.out.caller,
            append = T,
            sep = "\t",
            col.names = T)
        }
    message(" - Output in: ", vcf.out.caller)
}

# Final VCF consensus
if (is.vcf){
    meta <- paste0("##fileformat=VCFv4.2\n##source=Consensus", length(callers), "Callers (", paste0(callers, collapse = ","), ")\n", contigs_meta,
                    '##FILTER=<ID=PASS,Description="All filters passed">\n##FILTER=<ID=FAIL,Description="More than half the callers did not give a PASS">\n')
    # we need the meta contigs and the INFO
    meta <- paste0(meta, meta_consensus)
} else{
    meta <- "#version 2.4"
}

to.vcf <- all.muts[all.muts$isconsensus==T,]

if (is.vcf){
    col.out <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    to.vcf$ID <- to.vcf$DNAchange
    to.vcf$QUAL <- "."
    to.vcf$INFO <- to.vcf$INFO_consensus
    to.vcf$FORMAT <- "."
    to.vcf$FILTER <- to.vcf$FILTER_consensus
} else{
    col.out <- callers_meta[[c]]$header  # any caller header is fine
}

to.vcf <- to.vcf[,col.out][!duplicated(to.vcf),]
message("- Total variants ", prettyNum(nrow(to.vcf), big.mark = ","))
message("- Variants in consensus ", prettyNum(nrow(all.muts[(!duplicated(all.muts$DNAchange) & all.muts$isconsensus==T),]), big.mark = ","))

write(x = meta, file = vcf.out, ncolumns = 1)
fwrite(x = to.vcf, file = vcf.out, append = T, sep = "\t", col.names = T)
message("- Output in: ", vcf.out)

if (nrow(to.vcf) > 1e6){
    message("Too many variants (+1M!), plotting is skipped.")
} else {

    ## PLOTTING
    variants_list <- list()
    variants_list_pass <- list()

    for (c in callers){
        variants_list[[c]]      <- all.muts[all.muts$Caller==c,]$DNAchange
        variants_list_pass[[c]] <- all.muts[all.muts$Caller==c & all.muts$FILTER_consensus=="PASS",]$DNAchange
    }

    nonempty_variants      <- sum(lengths(variants_list)) > 0
    nonempty_variants_pass <- sum(lengths(variants_list_pass)) > 0

    if (nonempty_variants && nonempty_variants_pass){
        m <- make_comb_mat(variants_list)
        comb_order <- order(comb_size(m), decreasing = T)
        u  <- grid.grabExpr(draw(UpSet(m = m, comb_order = comb_order, column_title="All variants"), newpage = FALSE))


        m2 <- make_comb_mat(variants_list_pass)
        comb_order2 <- order(comb_size(m2), decreasing = T)
        g <- ggplot(all.muts, aes(Caller, fill=isconsensus)) +
            geom_bar() +
            coord_flip() +
            scale_fill_manual(values = c(`TRUE`='#247671', `FALSE`='#92C2B5')) +
            geom_text_repel(stat='count', aes(label=prettyNum(..count.., big.mark = ","))) +
            ggtitle(subtitle = "PASS=All filters passed (note that '.' will be considered FAIL)", label = "") +
            facet_grid(.~FILTER_consensus, scales="free")  + theme(title = element_text(color="grey40"))
        u2 <- grid.grabExpr(draw(UpSet(m = m2, comb_order = comb_order2, column_title="PASS variants"), newpage = FALSE))
        # 8.27 x 11.69
        pdf(pdf.out, width = 8.3, height = 11.7, paper = "A4")
        plot <- ggarrange(g, u, u2,
                            labels = c("A", "B", "C"),
                            ncol = 1, nrow = 3)
        print(annotate_figure(plot, top = text_grob(paste("Consensus summary for", sampleid),
                                                face = "bold", size = 10, family="Courier")))

        dev.off()
        message(" - Output in: ", pdf.out)
    }

    # check whether the unwanted file exists and remove it
    file.exists("Rplots.pdf")
    file.remove("Rplots.pdf")
}
