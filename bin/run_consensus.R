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
suppressPackageStartupMessages(library(parallel))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(datatable.fread.input.cmd.message=FALSE)
# Collect arguments
script_args <- commandArgs(TRUE)
# default arguments
if(length(script_args) < 1) {
    script_args <- c("--help")
}


# Help section
if("--help" %in% script_args) {
    cat("The consensusR script:

        Arguments:
        --input_dir=input_directory - character, Directory containing input files
        --out_prefix=output_prefix  - character, prefix for outputs
        --thr=thr                   - integer, the number of callers that support a mutation to be called consensus (default: 2)
        --cpu=cpu                   - integer, number of threads/cores to use in the parallelisation
        --interval=coord_interval   - character, chrom:start-end to analyse
        --help                      - print this text
        Note: Please ensure caller name is as follows: sample_id.caller.maf
        Example:
        ./consensusR.R --id=sample01 --input_dir=inputs/  --out_prefix=consensus --cpu=1 --thr=2 \n\n")

    q(save="no")
}
if(!grepl(pattern = "--input_dir=", x = paste(script_args, collapse = ""))) {
    stop("Missing --input_dir name!")
}
if(!grepl(pattern = "--out_prefix=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--out_prefix=consensus")
}
if(!grepl(pattern = "--thr=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--thr=2")
}
if(!grepl(pattern = "--cpu=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--cpu=1")
}
if(!grepl(pattern = "--interval=", x = paste(script_args, collapse = ""))) {
    script_args <- c(script_args, "--interval=0")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsDF       <- as.data.frame(do.call("rbind", parseArgs(script_args)))
argsDF       <- aggregate(.~V1, data=argsDF, paste, collapse=",")
argsL        <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

# Input variables
sampleid <- argsL$id
input_files <- list.files(strsplit(argsL[["input_dir"]], split=",")[[1]], pattern="\\.maf$", full.names=TRUE)

# Extract caller names from file names - specific name convention is expected
get_caller <- function(filename) {
    caller_name <- sub(".*(\\.|_)(.*?)\\.maf", "\\2", filename)
    if (caller_name == filename) caller_name <- "unknown"
    return(caller_name)

}
callers <- sapply(input_files, get_caller)

names(input_files) <- callers

is.vcf <- grepl(x = input_files[1], pattern = ".vcf$|.vcf.gz$", perl = T)

# output files
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
    contigs_meta <- ""
}
# Collecting vcf headers for the outputs
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

# Second, we read and convert the input files to GenomicRanges for easy manipulation.
message("- Converting to genomic ranges")
mutsGR <- list()
muts <-  list()
if (argsL$interval != '0') {

    message( ' - Subset input by interval: ', argsL$interval)
    interval_chr   = strsplit(argsL$interval, ':')[[1]][1]
    interval_start = as.numeric(strsplit(strsplit(argsL$interval, ':')[[1]][2], '-')[[1]][1])
    interval_end   = as.numeric(strsplit(strsplit(argsL$interval, ':')[[1]][2], '-')[[1]][2])

}

for(c in callers[1:length(callers)]){
    v <- input_files[c]
    message("  - ", v)
    if (!is.vcf){
        fread_cmd =  paste0("zgrep -v '#' ", v)
        if (argsL$interval != '0') {
            fread_cmd <- paste0(fread_cmd, "| grep -Ew 'Hugo_Symbol|", interval_chr, "'")
        }
        # message("CMD:", fread_cmd)
        tmp <- fread(cmd=fread_cmd, header=TRUE)
        if (argsL$interval != '0') {
            tmp <- tmp[Chromosome == interval_chr & Start_Position >= interval_start & Start_Position < interval_end]
        }
        tmp$`#CHROM` <- tmp$Chromosome
        tmp$POS      <- tmp$Start_Position
        tmp$REF      <- tmp$Reference_Allele
        tmp$ALT      <- tmp$Tumor_Seq_Allele2
    } else {
        fread_cmd =  paste0("zgrep -v '##' ", v)
        if (argsL$interval != '0') {
            fread_cmd <- paste0(fread_cmd, "| grep -w ", interval_chr)
        }
        tmp <- fread(cmd=fread_cmd)
    }
    if (nrow(tmp) > 0) {
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
        message("   - [WARNING] Removing ",nrow(tmp[is.na(end)]), " spurious calls with no alt. There should be 0.")
        tmp <- tmp[!is.na(end)]
        muts[[c]]   <- tmp
        mutsGR[[c]] <- GenomicRanges::makeGRangesFromDataFrame(df = tmp,
                                                        ignore.strand = TRUE,
                                                        start.field = "start",
                                                        end.field = "end",
                                                        seqnames.field = "#CHROM",
                                                        keep.extra.columns = T)
    } else {
        mutsGR[[c]] <- tmp  # empty
    }
}

# callers <- names(mutsGR) # All callers might not be present
message("- Finding overlaps")
# Third, we find overlaps
overlapping.vars <- data.frame(DNAchange=character(), caller=character(), FILTER=character())
for (c1 in names(mutsGR)){
    if (nrow(data.frame(mutsGR[[c1]]))==0) next
    for (c2 in names(mutsGR)){
        if (nrow(data.frame(mutsGR[[c2]]))==0) next
        if (c1!=c2){
            group_name <- paste0(c1, "vs", c2)
            # The gap between 2 adjacent ranges is 0.
            hits <- GenomicRanges::findOverlaps(query = mutsGR[[c1]], subject = mutsGR[[c2]], maxgap = 0)
            dnachange.hits <- muts[[c1]][queryHits(hits)]$DNAchange
            filt.hits      <- muts[[c1]][queryHits(hits)]$FILTER
            if (length(dnachange.hits) > 0) {
                # due to normalization we might find the same variant with different filters - these come from homopolymer regions
                overlapping.vars <- rbind(overlapping.vars, unique(data.frame(DNAchange = dnachange.hits, caller = c1, FILTER = filt.hits)))
            }
        }
    }
}

# Finally, extract the set of variants that will be the consensus set
overlapping.vars <- overlapping.vars[!duplicated(overlapping.vars),]
overlapping.vars <- as.data.table(overlapping.vars)[, .(caller = paste(caller, collapse = "|"), FILTER = paste(FILTER, collapse = "|")), by = DNAchange]
overlapping.vars <- as.data.frame(overlapping.vars)

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
what.caller.called <- function(row, consensus, variants){
    variant <- row["DNAchange"]
    if (variant %in% consensus){
        var.callers <- variants[variants$DNAchange==variant,]$caller
        var.callers <- paste(var.callers, collapse = "|")
        filters <- variants[variants$DNAchange==variant,]$FILTER
        filters <- paste(sub(pattern = ";",
                            replacement = ",",
                            x =filters),
                        collapse = "|")
        if (var.callers == ""){
        var.callers <- row["Caller"]
        filters <- row["FILTER"]
        }
        list(callers=var.callers, filters=filters)
    } else {
        list(callers=row["Caller"], filters=row['FILTER'])
    }
}

cpu <- ifelse(is.na(argsL$cpu)| is.null(argsL$cpu), 1, as.integer(argsL$cpu))
cl <- makeCluster(cpu)

for (c in callers){
    message("- Annotating calls from ", c)
    if (!is.null(muts[[c]])){
        values <-  parApply(cl=cl, X = muts[[c]], MARGIN = 1, FUN = what.caller.called, consensus=con.vars.ths, variants=overlapping.vars)
        muts[[c]] <- cbind( muts[[c]], as.data.frame(do.call(rbind, values)))
        muts[[c]]$callers <- unlist(muts[[c]]$callers )
        muts[[c]]$filters <- unlist(muts[[c]]$filters )
    }
}

stopCluster(cl)

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
    extra_cols <- c()
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


## PLOTTING
variants_list <- list()
variants_list_pass <- list()

for (c in callers){
    variants_list[[c]]      <- all.muts[all.muts$Caller==c,]$DNAchange
    variants_list_pass[[c]] <- all.muts[all.muts$Caller==c & all.muts$FILTER_consensus=="PASS",]$DNAchange
}

if (length(unlist(variants_list)) > 0 & length(unlist(variants_list_pass)) > 0){
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
    annotate_figure(plot, top = text_grob(paste("Consensus summary for", sampleid),
                                            face = "bold", size = 14, family="Courier"))

    dev.off()
    message(" - Output in: ", pdf.out)

}

# check whether the unwanted file exists and remove it
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")