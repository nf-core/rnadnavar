#!/usr/bin/env Rscript
# Date: Sun 20 Sep 2020
# Author: Raquel Manzano - CRUK CI Caldas lab
# Script: Find overlaps between vcf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
options(warn=-1)

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggrepel))


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
      --id=sampleID              - character, sample id to put in the plots
      --input=input_file.maf     - character, VCf/MAF file, can be specfied more than once for several inputs
      --caller=caller_name       - character, caller that was used to generate input file - has to be in the SAME order as input
      --out_prefix=output_prefix - character, preffix for outputs
      --thr=thr                  - integer, the number of callers that support a mutation to be called consensus (default: 2)
      --help                     - print this text
 
      Example:
      ./consensusR.R --id=sample01 --input=caller1.vcf --input=caller2.vcf --callers C1 --callers C2  --out_prefix=consensus \n\n")
  
  q(save="no")
}
if(!grepl(pattern = "--input=", x = paste(script_args, collapse = ""))) {
  stop("Missing --input file!")
}
if(!grepl(pattern = "--caller=", x = paste(script_args, collapse = ""))) {
  stop("Missing --caller value! Please put it in the same order as caller.")
}
if(!grepl(pattern = "--out_prefix=", x = paste(script_args, collapse = ""))) {
  script_args <- c(script_args, "--out_prefix=consensus")
}
if(!grepl(pattern = "---thr=", x = paste(script_args, collapse = ""))) {
  script_args <- c(script_args, "--thr=2")
}


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(script_args)))
argsDF <- aggregate(.~V1, data=argsDF, paste, collapse=",")
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

# Input variables
sampleid <- argsL$id
vcfs <- strsplit(argsL[["input"]], split=",")[[1]]
callers <- strsplit(argsL[["caller"]], split=",")[[1]]
names(vcfs) <- callers

# output files
pdf.out <- paste0(argsL$out_prefix, ".pdf")
if (grepl(x = vcfs[1], pattern = ".vcf", fixed = T)){
  vcf.out <- paste0(argsL$out_prefix, ".vcf")
} else{
  vcf.out <- paste0(argsL$out_prefix, ".maf")
}

message("- Parsing headers")
# First, get contigs from one of the VCFs
contigs_meta <- fread(paste0("zgrep '##contig' ", vcfs[1]), sep = NULL, header = F)
contigs_meta <- paste0(contigs_meta$V1, collapse = "\n")
# Collexting vcf headers for the outputs
callers_meta <- list()
for ( c in callers){
  vcf_meta <- fread(paste0("zgrep -E '##|#version' ", vcfs[c]), sep = NULL, header = F)
  if (grepl(pattern = ".vcf", x = vcfs[c] )){
    vcf_header <- fread(paste0("zgrep '#CHROM' ", vcfs[c]), sep = NULL, header = F)
    callers_meta[[c]] <-  list(meta=paste0(vcf_meta$V1, collapse = "\n"),
                             header=strsplit(vcf_header$V1, "\t")[[1]])
  } else{
    maf_header <- fread(paste0("zgrep 'Hugo_Symbol' ", vcfs[c]), sep = NULL, header = F)
    callers_meta[[c]] <-  list(meta=paste0(vcf_meta$V1, collapse = "\n"),
                               header=strsplit(maf_header$V1, "\t")[[1]])
  }
  
}


# Second, we read and convert the VCFs in GenomicRanges for easy manipulation.
message("- Converting to genomic ranges")
mutsGR <- list() 
muts <-  list()
for(c in callers[1:length(callers)]){
  v <- vcfs[c]
  tmp <- fread(paste0("zgrep -v '##' ", v))
  if (!grepl(x = vcfs[1], pattern = ".vcf", fixed = T)){
    tmp$`#CHROM` <- tmp$Chromosome
    tmp$POS <- tmp$Start_Position
    tmp$REF <- tmp$Reference_Allele
    tmp$ALT <- tmp$Tumor_Seq_Allele2
  } 
  
  tmp$Caller <- c
  tmp$mut <- paste0(tmp$REF, ">", tmp$ALT)
  tmp$DNAchange <- paste0(tmp$`#CHROM`, ":g.", tmp$POS, tmp$REF, ">", tmp$ALT)
  tmp$start <- tmp$POS
  # get end position to create a range
  tmp$end <- ifelse(nchar(tmp$REF) > nchar(tmp$ALT),  
                     tmp$POS + nchar(tmp$REF)-1,
                     ifelse(nchar(tmp$REF) < nchar(tmp$ALT), tmp$POS + nchar(tmp$ALT)-1,
                            ifelse(nchar(tmp$REF) == nchar(tmp$ALT) & nchar(tmp$ALT) > 1, tmp$POS + nchar(tmp$REF),
                                   tmp$POS)
                     )
  )
  muts[[c]] <- tmp
  mutsGR[[c]] <- GenomicRanges::makeGRangesFromDataFrame(df = tmp, 
                                                    ignore.strand=TRUE, 
                                                    start.field = "start", 
                                                    end.field = "end", 
                                                    seqnames.field = "#CHROM",
                                                    keep.extra.columns = T)
}

message("- Finding overlaps")
# Third, we find overlaps
overlaps <- list()
overlapping.vars <- data.frame(DNAchange=character(), caller=character(), FILTER=character())
for (c1 in callers){
  for (c2 in callers){
    if (c1!=c2){
    group_name <- paste0(c1, "vs", c2)
    # The gap between 2 adjacent ranges is 0.
    hits <- GenomicRanges::findOverlaps(query = mutsGR[[c1]], subject = mutsGR[[c2]], maxgap = 0)
    overlaps[[group_name]] <- hits
    m.hits <- muts[[c1]][queryHits(hits)]$DNAchange
    filt <- muts[[c1]][queryHits(hits)]$FILTER
    # due to normalization we might find the same variant with different filters - these come from homopolymer regions
    overlapping.vars <- rbind(overlapping.vars, data.frame(DNAchange=m.hits, caller=c1, FILTER=filt))
    }
  }
}

# Finally, extract the set of variants that will be the consensus set
overlapping.vars <- overlapping.vars[!duplicated(overlapping.vars),]
overlapping.vars <- aggregate( FILTER  ~  DNAchange + caller, overlapping.vars, paste, collapse=";" )

# Set of variants that are called in at least 2 callers
con.vars <- overlapping.vars$DNAchange
#tail(sort(table(con.vars)[table(con.vars)>=2]))
con.vars.ths <- names(table(con.vars)[table(con.vars)>=2])

message("- There are ", prettyNum(length(con.vars.ths), big.mark = ','), " variants that are consensus")


# The next steps are for the output
# To keep the information from the consensus we extract the callers that called each mutation and its correspondent filters.
what.caller.called <- function(row, consensus, variants){
  variant <- row["DNAchange"]
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
}
# Parallel step
cl <- makeCluster(detectCores())
for (c in callers){
  message("- Annotating calls from ", c)
  values <-  parApply(cl = cl, X = muts[[c]], MARGIN = 1, FUN = what.caller.called, consensus=con.vars.ths, variants=overlapping.vars)
  muts[[c]] <- cbind( muts[[c]], as.data.frame(do.call(rbind, values)))
  muts[[c]]$callers <- unlist(muts[[c]]$callers )
  muts[[c]]$filters <- unlist(muts[[c]]$filters )
}
stopCluster(cl)

all.muts <- do.call(rbind.fill, muts)

message("- Preparing output")
## Prepare output
simplified.filter <- sapply(all.muts$filters, FUN = function(x){
  filt.val <- strsplit(x=x, split = "|", fixed = T)[[1]]
  filt.t <- table(filt.val == "PASS")
  ifelse(prop.table(filt.t)["FALSE"] > 0.5, "FAIL", "PASS")
})
simplified.filter <- ifelse(is.na(simplified.filter), "PASS", simplified.filter)
all.muts$FILTER_consensus <- simplified.filter
all.muts$INFO_consensus <-  paste0("callers=", all.muts$callers, ";filters=", all.muts$filters, ";consensus_filter=", all.muts$FILTER_consensus)

## write consensus
# I want to keep the ingo without duplicating the mutations
all.muts$isconsensus <- grepl(pattern = "|", x = all.muts$callers, fixed = T)

# WRITE OUTPUTS

meta_consensus <- paste0('##INFO=<ID=callers,Number=1,Type=String,Description="Variant callers that called this mutation, separated by |">\n',
'##INFO=<ID=filters,Number=1,Type=String,Description="Filters provided by each variant caller, separated by |">\n',
'##INFO=<ID=consensus_filter,Number=1,Type=String,Description="PASS if 50% or more of the callers give the mutation, otherwise FAIL.">')
extra.cols <- c()
for ( c in callers){
  if (grepl(x = vcfs[c], pattern = ".vcf", fixed = T)){
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
  extra_cols <- c()
  if (grepl(x = vcfs[c], pattern = ".vcf", fixed = T)){
    fwrite(x = all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header], 
           file = vcf.out.caller, 
           append = T, 
           sep = "\t", 
           col.names = T)
  } else{
    colnames(all.muts)[colnames(all.muts)=="FILTER_consensus"] <- "FILTER_consensus" 
    fwrite(x = all.muts[all.muts$Caller==c,][,callers_meta[[c]]$header], 
           file = vcf.out.caller, 
           append = T, 
           sep = "\t", 
           col.names = T)
  }
  message(" - Output in: ", vcf.out.caller)
}


# Final VCF consensus
if (grepl(x = vcfs[1], pattern = ".vcf", fixed = T)){
  meta <- paste0("##fileformat=VCFv4.2\n##source=Consensus", length(callers), "Callers (", paste0(callers, collapse = ","), ")\n", contigs_meta,
                 '##FILTER=<ID=PASS,Description="All filters passed">\n##FILTER=<ID=FAIL,Description="More than half the callers did not give a PASS">\n')
  # we need the meta contigs and the INFO
  meta <- paste0(meta, 
                 meta_consensus)
} else{
  meta <- "#version 2.4"
}

# to.vcf <- all.muts[all.muts$isconsensus==T,]
if (grepl(x = vcfs[1], pattern = ".vcf", fixed = T)){
  col.out <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  to.vcf$ID <- to.vcf$DNAchange
  to.vcf$QUAL <- "."
  to.vcf$INFO <- to.vcf$INFO_consensus 
  to.vcf$FORMAT <- "."
  to.vcf$FILTER <- to.vcf$FILTER_consensus
} else{
  col.out <- callers_meta[[c]]$header
}
to.vcf <- all.muts[,col.out]
message("- Total variants ", prettyNum(nrow(to.vcf), big.mark = ","))
message("- Variants in consensus ", prettyNum(nrow(all.muts[(!duplicated(all.muts$DNAchange) & all.muts$isconsensus==T),]), big.mark = ","))


write(x = meta, file = vcf.out, ncolumns = 1)
fwrite(x = to.vcf, file = vcf.out, append = T, sep = "\t", col.names = T)
message(" - Output in: ", vcf.out)



## PLOTTING
variants_list <- list()
variants_list_pass <- list()

for (c in callers){
  variants_list[[c]]      <- all.muts[all.muts$Caller==c,]$DNAchange
  variants_list_pass[[c]] <- all.muts[all.muts$Caller==c & all.muts$FILTER_consensus=="PASS",]$DNAchange
}

m <- make_comb_mat(variants_list)
comb_order <- order(comb_size(m), decreasing = T)
m2 <- make_comb_mat(variants_list_pass)
comb_order2 <- order(comb_size(m2), decreasing = T)
g <- ggplot(all.muts, aes(Caller, fill=isconsensus)) + 
  geom_bar() + 
  coord_flip() + 
  scale_fill_manual(values = c(`TRUE`='#247671', `FALSE`='#92C2B5')) +
  geom_text_repel(stat='count', aes(label=prettyNum(..count.., big.mark = ","))) +
  ggtitle(subtitle = "PASS=All filters passed (note that '.' will be considered FAIL)", label = "") +
  facet_grid(.~FILTER_consensus, scales="free")  + theme(title = element_text(color="grey40"))
u  <- grid.grabExpr(draw(UpSet(m = m, comb_order = comb_order, column_title="All variants"), newpage = FALSE)) 
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

# check whether the unwanted file exists and remove it
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")
