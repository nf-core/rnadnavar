#!/usr/bin/env python3
"""
Author: @RaqManzano
Script: Filter variants from a MAF file producing another MAF file with the new filters added.
"""
import argparse
import pandas as pd
import subprocess
import pysam


def argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i", "--input", help="MAF(s) file input (if more than 1 consensus will be annotated)", nargs="+", required=True
    )
    parser.add_argument("-o", "--output", help="MAF file output", default="RaVeX.maf")
    parser.add_argument(
        "-g", "--gnomad_thr", help="Gnomad threshold for variants (must be annotated in MAF)", default=0.0001
    )
    parser.add_argument("--whitelist", help="BED file with variants to keep (CHROM POS REF ALT)")
    parser.add_argument("--blacklist", help="BED file with regions to remove (CHROM START END)")
    parser.add_argument("--filters", help="Other filters to be considered as PASS", default=["PASS"], nargs="+")
    parser.add_argument("--ref", help="FASTA reference file to extract context")
    parser.add_argument("--vc_priority", help="The order of priority will mark which caller annotation will be kept. Only one annotation is kept per sample and per caller.",
                        default=["mutect2", "sage", "strelka", "consensus"], nargs="+")
    return parser.parse_args()


def read_whitelist_bed(bed_file):
    """
    Columns should be chrom, start, end, ref and alt to read BED file with variants for whitelisting
    """
    if not bed_file:
        variants = []
    else:
        bed = pd.read_csv(bed_file, sep="\t", comment="#", header=None)
        if len(bed.columns) == 5:
            colnames = ["#CHROM", "POS", "END", "REF", "ALT"]
        else:
            colnames = ["#CHROM", "POS", "REF", "ALT"]
        bed.columns = colnames
        try:
            bed["DNAchange"] = bed["#CHROM"].map(str) + ":g." + bed["POS"].map(str) + bed["REF"] + ">" + bed["ALT"]
        except TypeError:
            print(
                "[ERROR] BED file for whitelist should contain CHROM, START, END, REF and ALT columns with no headers."
            )
            print(f"Please check your file: {bed_file}")
        variants = bed["DNAchange"].tolist()
    return variants


def read_blacklist_bed(bed_file):
    """
    Columns should be chrom, start, end to read BED file with regions that will be considered in blacklisting
    """

    bed = pd.read_csv(bed_file, sep="\t", comment="#", header=None)
    assert (
        len(bed.columns) >= 3
    ), "[ERROR] BED file for blacklist should at least contain CHROM, START, END columns with no headers."
    bed.rename(columns={"0": "#CHROM", "1": "POS", "2": "END"}, inplace=True)
    return bed


def read_maf(maf_file):
    if type(maf_file) == type([]):
        maf_list = []
        for m in maf_file:
            maf_list += [pd.read_csv(m, sep="\t", comment="#", low_memory=False)]
        maf = pd.concat(maf_list)
    else:
        maf = pd.read_csv(maf_file, sep="\t", comment="#", low_memory=False)
    if "DNAchange" not in maf.columns:
        maf["DNAchange"] = (
            maf["Chromosome"].map(str)
            + ":g."
            + maf["Start_Position"].map(str)
            + maf["Reference_Allele"]
            + ">"
            + maf["Tumor_Seq_Allele2"]
        )
    return maf


def noncoding(maf, noncoding):
    """
    Classifies as noncoding TRUE if first consequence is in noncoding list
    """
    try:
        maf["noncoding"] = maf["Consequence"].apply(
            lambda consequence: consequence.split("&")[0].split(",")[0] in noncoding
        )
    except AttributeError:
        maf["noncoding"] = "NA"
    return maf


def remove_ig_and_pseudo(maf):
    """
    Add IG and pseudogene filters
    """
    # fill na values
    maf[["BIOTYPE", "SYMBOL"]] = maf[["BIOTYPE", "SYMBOL"]].fillna(value="")
    maf["ig_pseudo"] = maf["BIOTYPE"].str.contains(
        "IG_C_gene|IG_D_gene|IG_J_gene|IG_V_gene|TR_C_gene|TR_J_gene|TR_V_gene|pseudogene"
    )
    return maf


def filter_homopolymer(ref_context, alt, hp_length=6):
    """
    Checks if the variant is in a homopolymer context
    """
    homopolymer = False
    if len(alt) > 1:
        alt = alt[1:]  # if insertion the check will be done
    elif alt == "-":  # if deletion no homopolymer
        return False
    elif ref_context == None:
        return None
    # we calculatate the actual length of the sequence we need to check
    # context is as follows: flank_before:REF:flank_after (flanks are same length)
    try:
        length_to_consider = int((len(ref_context) - 1) / 2)
    except TypeError:
        return None
    # we substitute the ref for the alt to check the context of the variant
    ref_context = list(ref_context)
    ref_context[length_to_consider] = alt
    ref_context = "".join(ref_context)
    # only the contact that would overlap with the expected hp_length
    # (i.e. some context does not need to be consider as it is too far)
    context_to_consider = ref_context[length_to_consider - hp_length + 1 : length_to_consider + hp_length]
    # count the bases that are equal. If count >= hp_length it will be consider a homopolymer
    for idx, base in enumerate(context_to_consider):
        context_window = context_to_consider[idx : idx + hp_length]
        if len(context_window) < hp_length:
            break  # what is left of the sequence to check is too short
        elif len(set(context_window)) == 1:
            homopolymer = True  # bingo
            break
    return homopolymer


def add_context(chrom, pos, ref, genome, flank=10):
    if pos < 10:
        flank = pos - 1
    try:
        context = genome.fetch(chrom, pos - 1 - flank, pos + flank).upper()
    except ValueError:
        print(f"[WGN] This variant has NaN in their coordinates (did liftover failed?)")
        return None
    except KeyError:
        print(
            f"[WGN] These coordinates {chrom}:{pos} are not present in the genome (are you using a different version?)"
        )
        return None
    if ref != "-":  # if it is a deletion we cannot check
        try:
            assert ref[0] == context[flank]  # check that the REF matches the context we just extracted
        except AssertionError:
            print(
                f"[WGN] This ref base ({chrom}:{int(pos)} {ref[0]} != {context[flank]}) does not correspond to its context {context[flank]}."
            )
    return context


def remove_homopolymers(maf, ref):
    """
    Check for variants in homopolymer regions (a sequence of 6 consecutive identical bases)
    """
    # read genome to get context
    genome = pysam.FastaFile(ref)
    # Add context
    maf["CONTEXT"] = maf.apply(
        lambda row: add_context(str(row["Chromosome"]), row["Start_Position"], row["Reference_Allele"], genome), axis=1
    )
    # add homopolymer True/False
    maf["homopolymer"] = maf.apply(lambda row: filter_homopolymer(row["CONTEXT"], row["Tumor_Seq_Allele2"]), axis=1)
    return maf


def remove_muts_in_range(df, blacklist):
    """
    If overlap with a blacklist region will return a dataframe with the information
    """
    df["blacklist"] = False
    df["blk_reason"] = ""
    reason = ""
    for idx, row in blacklist.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]
        if len(row) >= 3:  # optional
            reason = row[3]
        df_bkl = df[df["Chromosome"] == chrom][df["Start_Position"].between(start, end, inclusive="both")].index
        df.loc[df_bkl, "blacklist"] = True
        df.loc[df_bkl, "blk_reason"] = reason
    return df

def merge_rows(group):
    """Merge info when mut is repeated due to consensus
    """
    # Take the first row as the base
    merged_row = group.iloc[0].copy()
    # Merge the 'callers' column
    merged_row['callers'] = group['callers'].iloc[0] + "|" + group['Caller'].iloc[1] if len(group) > 1 else group['callers'].iloc[0]
    # Merge the 'filters' column
    merged_row['filters'] = group['filters'].iloc[0] + "|" +  group['FILTER_consensus'].iloc[1] if len(group) > 1 else group['callers'].iloc[0]
    # Set the Tumor_Sample_Barcode_consensus column
    merged_row['Tumor_Sample_Barcode_consensus'] = group['Tumor_Sample_Barcode'].iloc[1] if len(group) > 1 else None

    return merged_row


def filtering(maf, gnomad_thr, whitelist, blacklist, filters):
    """
    Adds filters for gnomad, blacklisting, variant calling filters. Adds muts to whitelist if match.
    """
    if "PASS" not in filters:
        filters += ["PASS"]  # a PASS is always allowed
    if whitelist:
        maf["whitelist"] = maf["DNAchange"].isin(whitelist)  # whitelist
    if not blacklist.empty:
        maf = remove_muts_in_range(df=maf, blacklist=blacklist)  # blacklist
    maf["ingnomAD"] = maf["MAX_AF"] >= gnomad_thr  # gnomad

    return maf


def add_ravex_filters(
    maf, filters, noncoding=False, homopolymer=False, ig_pseudo=False, min_alt_reads=2, blacklist=False, whitelist=False
):
    maf["RaVeX_FILTER"] = "PASS"
    maf["Existing_variation"] = maf["Existing_variation"].fillna("")
    maf["SOMATIC"] = maf["SOMATIC"].fillna("")
    if "FILTER" not in maf.columns:
        maf["FILTER"] = "PASS"  # Consider pass when no FILTER column
    maf["FILTER"] = maf["FILTER"].replace(".", "PASS")  # no filter ('.') will be treated as PASS
    ignore_consensus = False
    if "isconsensus" not in maf.columns:
        maf["isconsensus"] = True  # By default true when not known
        ignore_consensus = True
    for idx, row in maf.iterrows():
        ravex_filter = []
        if row["t_alt_count"] <= min_alt_reads:
            ravex_filter += ["min_alt_reads"]
        if not blacklist.empty:
            if row["blacklist"]:
                ravex_filter += ["blacklist"]
        if not noncoding:
            if row["noncoding"]:
                ravex_filter += ["noncoding"]
        if not homopolymer:
            if row["homopolymer"]:
                ravex_filter += ["homopolymer"]
        if not ig_pseudo:
            if row["ig_pseudo"]:
                ravex_filter += ["ig_pseudo"]
        if not row["isconsensus"]:  # of there is consensus we take the FILTER from the consensus
            if not row["FILTER"] in filters:
                ravex_filter += ["vc_filter"]
            if not ignore_consensus:
                ravex_filter += ["not_consensus"]
        elif not ignore_consensus:
            if not row["FILTER_consensus"] in filters:
                ravex_filter += ["vc_filter"]
        else:
            if not row["FILTER"] in filters:
                ravex_filter += ["vc_filter"]
        if whitelist:
            if not ravex_filter or row["whitelist"]:
                ravex_filter = ["PASS"]
        else:
            if not ravex_filter:
                ravex_filter = ["PASS"]
        ravex_filter = ";".join(ravex_filter)
        maf.at[idx, "RaVeX_FILTER"] = ravex_filter
    return maf


def deduplicate_maf(variants, vc_priority):
    deduped = []
    for caller in vc_priority:
        deduped.append(variants[variants["Caller"] == caller])
    # we deduplicate per DNAchange and Tumor_Sample_Barcode
    deduplicated = pd.concat(deduped).drop_duplicates(subset=["DNAchange", "Tumor_Sample_Barcode"], keep="first")
    # when two same mutations are still present, we assume this is due to "consensus" from different data type (DNA/RNA)
    grouped = deduplicated.groupby('DNAchange')
    # merge the info from consensus
    consensus_merged = grouped.apply(merge_rows).reset_index(drop=True)
    return consensus_merged

def write_maf(maf_df, mafin_file, mafout_file, vc_priority):
    """Write output"""
    header_lines = subprocess.getoutput(f"zgrep -Eh '#|Hugo_Symbol' {mafin_file} 2>/dev/null")
    if "Caller" in maf_df.columns:
        print("Removing duplicated variants from maf (only one entry from a caller will be kept)\nNote: `consensus` info will be merged")
        # Separate the multiallelic variants
        multiallelic_variants = maf_df[maf_df["FILTER"].str.contains("multiallelic", case=False, na=False)]
        other_variants = maf_df[~maf_df["FILTER"].str.contains("multiallelic", case=False, na=False)]

        # Combine variants with Caller from multiallelic
        multiallelic_variants["Caller"] = multiallelic_variants["Caller"] + "_multiallelic"
        maf_to_dedup = pd.concat([other_variants, multiallelic_variants])

        # Deduplicate variants
        maf_to_write = deduplicate_maf(maf_to_dedup, vc_priority + list(multiallelic_variants["Caller"].unique()))
    else:
        maf_to_write = maf_df

    # Write the header and deduplicated MAF to the output file
    with open(mafout_file, "w") as mafout:
        mafout.write(header_lines)
    maf_to_write.to_csv(mafout_file, mode="a", index=False, header=True, sep="\t")
    print(f"Done! See '{mafout_file}'.")


def main():
    noncoding_list = [
        "intron_variant",
        "intergenic_variant",
        "non_coding_transcript_variant",
        "non_coding_transcript_exon_variant",
        "mature_miRNA_variant",
        "regulatory_region_variant",
        "IGR",
        "INTRON",
        "RNA",
    ]
    args = argparser()
    maf = read_maf(args.input)
    whitelist = False
    blacklist = False
    if args.whitelist:
        whitelist = read_whitelist_bed(args.whitelist)
    if args.blacklist:
        blacklist = read_blacklist_bed(args.blacklist)
    else:
        blacklist = pd.DataFrame()
    maf = filtering(maf=maf, gnomad_thr=args.gnomad_thr, whitelist=whitelist, blacklist=blacklist, filters=args.filters)
    # tag noncoding
    maf = noncoding(maf=maf, noncoding=noncoding_list)
    # tag IG and pseudo
    maf = remove_ig_and_pseudo(maf=maf)
    # tag homopolymers
    maf = remove_homopolymers(maf=maf, ref=args.ref)
    # tag consensus
    maf = add_ravex_filters(maf=maf, filters=args.filters, blacklist=blacklist, whitelist=whitelist)
    if not args.output:
        args.output = args.input.replace(".maf", "filtered.maf")
    write_maf(maf_df=maf, mafin_file=args.input, mafout_file=args.output,vc_priority=args.vc_priority)


if __name__ == "__main__":
    main()
