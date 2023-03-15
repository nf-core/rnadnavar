#!/usr/bin/env python3
"""
Date: 13 Feb 2023
Author: @RaqManzano
Script: Filter variants from a MAF file producing another MAF file with the new filters added.
"""
import argparse
import pandas as pd
import subprocess
import pysam
import re


def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", help="MAF(s) file input (if more than 1 consensus will be annotated)",
                        nargs="+", required=True)
    parser.add_argument("-o", "--output", help="MAF file output", default="RaVeX.maf")
    parser.add_argument("-g", "--gnomad_thr",
                        help="Gnomad threshold for variants (must be annotated in MAF)",
                        default=0.0001)
    parser.add_argument("--whitelist", help="BED file with variants to keep (CHROM POS REF ALT)")
    parser.add_argument("--blacklist", help="BED file with regions to remove (CHROM START END)")
    parser.add_argument("--filters", help="Other filters to be considered as PASS", default=["PASS"], nargs="+")
    parser.add_argument("--ref", help="FASTA reference file to extract context")
    return parser.parse_args()


def read_whitelist_bed(bed_file):
    """Columns should be chrom, start, end, ref and alt"""
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
            bed["DNAchange"] = bed['#CHROM'] + ":g." + bed['POS'].map(str) + bed['REF'] + ">" + bed['ALT']
        except TypeError:
            print(
                "[ERROR] BED file for whitelist should contain CHROM, START, END, REF and ALT columns with no headers.")
            print(f"Please check your file: {bed_file}")
        variants = bed["DNAchange"].tolist()
    return variants


def read_blacklist_bed(bed_file):
    """Columns should be chrom, start, end"""

    bed = pd.read_csv(bed_file, sep="\t", comment="#", header=None)
    assert len(
        bed.columns) >= 3, "[ERROR] BED file for blacklist should at least contain CHROM, START, END columns with no headers."
    bed.rename(columns={"0": "#CHROM", "1": "POS", "2": "END"}, inplace=True)
    return bed


def read_maf(maf_file):
    if type(maf_file) == type([]):
        maf_list = []
        for m in maf_file:
            maf_list += [pd.read_csv(m, sep="\t", comment="#")]
        maf = pd.concat(maf_list)
    else:
        maf = pd.read_csv(maf_file, sep="\t", comment="#")
    if "DNAchange" not in maf.columns:
        maf["DNAchange"] = maf['Chromosome'] + ":g." + maf['Start_Position'].map(str) + \
                           maf['Reference_Allele'] + ">" + maf['Tumor_Seq_Allele2']
    return maf


def noncoding(maf, noncoding):
    maf["noncoding"] = maf["Consequence"].apply(
        lambda consequence: consequence.split('&')[0].split(',')[0] in noncoding)
    # change filter accordingly
    # maf.loc[maf["noncoding"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["noncoding"])
    return maf


def remove_ig_and_psedo(maf, ig_list):
    """Add IG and pseudogene filters"""
    # fill na values
    maf[['BIOTYPE', 'SYMBOL']] = maf[['BIOTYPE', 'SYMBOL']].fillna(value="")
    # mask variants that were annotated as IG gene or pseudogene
    # mask_ig = (maf["SYMBOL"].isin(ig_list))
    # mask_ps = (maf["BIOTYPE"].str.contains('pseudogene'))
    maf["ig_pseudo"] = (maf["SYMBOL"].isin(ig_list)) | (maf["BIOTYPE"].str.contains('pseudogene'))
    # change filter accordingly
    # maf.loc[mask_ig, 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["ig_gene"])
    # maf.loc[mask_ps, 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["pseudogene"])
    return maf


def filter_homopolymer(ref_context, alt, hp_length=6):
    """Checks if the variant is in a homopolymer context"""
    homopolymer = False
    if len(alt) > 1:
        alt = alt[1:]  # if insertion the check will be done
    elif alt == "-":  # if deletion no homopolymer
        return False
    # we calculatate the actual length of the sequence we need to check
    # context is as follows: flank_before:REF:flank_after (flanks are same length)
    length_to_consider = int((len(ref_context) - 1) / 2)
    # we substitute the ref for the alt to check the context of the variant
    ref_context = list(ref_context)
    ref_context[length_to_consider] = alt
    ref_context = "".join(ref_context)
    # only the contact that would overlap with the expected hp_length
    # (i.e. some context does not need to be consider as it is too far)
    context_to_consider = ref_context[length_to_consider - hp_length + 1:length_to_consider + hp_length]
    # count the bases that are equal. If count >= hp_length it will be consider a homopolymer
    for idx, base in enumerate(context_to_consider):
        context_window = context_to_consider[idx:idx + hp_length]
        if len(context_window) < hp_length:
            break  # what is left of the sequence to check is too short
        elif len(set(context_window)) == 1:
            homopolymer = True  # bingo
            break
    return homopolymer


def add_context(chrom, pos, ref, genome, flank=10):
    if pos < 10:
        flank = pos - 1
    context = genome.fetch(chrom, pos - 1 - flank, pos + flank).upper()
    if ref != "-":  # if it is a deletion we cannot check
        assert ref[0] == context[flank]  # check that the REF matches the context we just extracted
    return context


def add_to_filter(current_filter, new_filter):
    if new_filter:
        if current_filter == "PASS":
            current_filter = new_filter[:]
        else:
            current_filter += f';{new_filter}'
    return current_filter


def remove_homopolymers(maf, ref):
    """
    Check for variants in homopolymer regions (a sequence of 6 consecutive identical bases)
    """
    # read genome to get context
    genome = pysam.FastaFile(ref)
    # Add context
    maf["CONTEXT"] = maf.apply(lambda row: add_context(row["Chromosome"],
                                                       row["Start_Position"],
                                                       row["Reference_Allele"],
                                                       genome),
                               axis=1)
    # add homopolymer True/False
    maf["homopolymer"] = maf.apply(lambda row: filter_homopolymer(row["CONTEXT"],
                                                                  row["Tumor_Seq_Allele2"]),
                                   axis=1)
    # change filter accordingly
    # maf.loc[maf["homopolymer"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["homopolymer"])
    return maf


def remove_muts_in_range(df, blacklist):
    df["blacklist"] = False
    df["blk_reason"] = ""
    reason = ""
    for idx, row in blacklist.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]
        if len(row) >= 3:  # optional
            reason = row[3]
        df_bkl = df[df["Chromosome"] == chrom][df["Start_Position"].between(start, end, inclusive=True)].index
        df.loc[df_bkl, "blacklist"] = True
        df.loc[df_bkl, "blk_reason"] = reason
    return df


# def adjust_filters():
#     """Adjust filter for variants that have been called more than once. Some variants might FAIL/PASS
#     when when low coverage. We leverage the information of the cohort to review if this was a fair filter"""

def filtering(maf, gnomad_thr, whitelist, blacklist, filters):
    if "PASS" not in filters:
        filters += ["PASS"]
    af_columns = [x for x in maf.columns if x.endswith("_AF") and x != "t_AF" and x != "n_AF"]
    maf["whitelist"] = maf['DNAchange'].isin(whitelist)
    maf = remove_muts_in_range(df=maf, blacklist=blacklist)
    maf["ingnomAD"] = maf["MAX_AF"] >= gnomad_thr

    return maf


def read_vcf(path):
    """TODO: add source"""
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def match_ref_alt_2maf_format(pos_ref_alt):
    pos, ref, alt = pos_ref_alt
    sol = str(pos) + ref + ">" + alt
    if len(ref) < len(alt) and len(ref) == 1:  # INS
        sol = str(pos) + "->" + alt[1:]
    elif len(ref) > len(alt) and len(alt) == 1:  # DEL
        sol = str(pos + 1) + ref[1:] + ">-"
    return sol


def add_ravex_filters(maf, filters, noncoding=False, homopolymer=False, ig_pseudo=False, min_alt_reads=3,
                      consensus=True):
    maf["RaVeX_FILTER"] = "PASS"
    maf["Existing_variation"] = maf["Existing_variation"].fillna("")
    maf["SOMATIC"] = maf["SOMATIC"].fillna("")
    for idx, row in maf.iterrows():
        ravex_filter = []
        if row["t_alt_count"] <= min_alt_reads:
            ravex_filter += ["min_alt_reads"]
        if row["ingnomAD"]:
            ravex_filter += ["gnomad"]
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
            ravex_filter += ["not_consensus"]
        else:
            if not row["FILTER_consensus"] in filters:
                ravex_filter += ["vc_filter"]
        if not ravex_filter or row["whitelist"]:
            ravex_filter = ["PASS"]
        ravex_filter = ";".join(ravex_filter)
        maf.at[idx, "RaVeX_FILTER"] = ravex_filter
    return maf


def dedup_maf(dnachange, maf):
    """If more than one caller, then we need to dedup the entries. We select according to the caller we trust more:
    mutect2, strelka, sage, others for SNVs, strelka, mutect2, sage for indels
    """
    maf_variant = maf[maf["DNAchange"] == dnachange]

    if maf_variant.iloc[0]["Variant_Type"] == "SNP":
        maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("mutect")]
        if maf_variant_dedup.empty:
            maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("strelka")]
            if maf_variant_dedup.empty:
                maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("sage")]
                if maf_variant_dedup.empty:
                    maf_variant_dedup = maf_variant.drop_duplicates(subset='DNAchange', keep="first")
    else:
        maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("strelka")]
        if maf_variant_dedup.empty:
            maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("mutect")]
            if maf_variant_dedup.empty:
                maf_variant_dedup = maf_variant[maf_variant["Caller"].str.contains("sage")]
                if maf_variant_dedup.empty:
                    maf_variant_dedup = maf_variant.drop_duplicates(subset='DNAchange', keep="first")
    return maf_variant_dedup


def write_maf(maf_df, mafin_file, mafout_file, vc_priority=["mutect2", "strelka", 'sage', 'consensus']):
    """Write output"""
    header_lines = subprocess.getoutput(f"zgrep -Eh '#|Hugo_Symbol' {mafin_file} 2>/dev/null")
    print("Removing duplicated variants from maf (only one entry from a caller will be kept)")
    maf_dedup = []
    for caller in vc_priority:
        maf_dedup += [maf_df[maf_df["Caller"]==caller]]
    maf_dedup = pd.concat(maf_dedup).drop_duplicates(subset='DNAchange', keep="first")
    with open(mafout_file, "w") as mafout:
        mafout.write(header_lines)
    maf_dedup.to_csv(mafout_file, mode="a", index=False, header=True, sep="\t")
    print(f"Done! See '{mafout_file}'.")


def main():
    noncoding_list = ['intron_variant', 'intergenic_variant',
                      'non_coding_transcript_variant', 'non_coding_transcript_exon_variant',
                      'mature_miRNA_variant', 'regulatory_region_variant',
                      'IGR', 'INTRON', 'RNA']
    # we need to add COSMIC	CENTERS	CONTEXT	DBVS	NCALLERS
    args = argparser()
    maf = read_maf(args.input)
    maf = filtering(maf=maf, gnomad_thr=args.gnomad_thr,
                    whitelist=read_whitelist_bed(args.whitelist),
                    blacklist=read_blacklist_bed(args.blacklist),
                    filters=args.filters)
    # tag noncoding
    maf = noncoding(maf=maf, noncoding=noncoding_list)
    # tag IG and pseudo
    maf = remove_ig_and_psedo(maf=maf,
                              ig_list=['IGH', 'IGL', 'IGK', 'TRBV', 'TRAJ', 'pseudo'])
    # tag homopolymers
    maf = remove_homopolymers(maf=maf, ref=args.ref)
    # tag consensus
    maf = add_ravex_filters(maf=maf, filters=args.filters)
    if not args.output:
        args.output = args.input.replace(".maf", "filtered.maf")
    write_maf(maf_df=maf, mafin_file=args.input, mafout_file=args.output)


if __name__ == '__main__':
    main()
