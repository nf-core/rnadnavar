#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Filters MAF file with filters from variant calling, gnomad and whitelist/blacklist.
"""
import argparse
import pandas as pd
import subprocess


def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", help="MAF file input", required=True)
    parser.add_argument("-o", "--output", help="MAF file output")
    parser.add_argument("-g", "--gnomad_thr",
                        help="Gnomad threshold for variants (must be annotated in MAF)",
                        default=0.0001 )
    parser.add_argument("--whitelist", help="BED file with variants to keep (CHROM POS REF ALT)")
    parser.add_argument("--blacklist", help="BED file with variants to remove (CHROM POS REF ALT)")
    parser.add_argument("--filters", help="Other filters to be considered as PASS", default=["PASS"], nargs="+")
    return parser.parse_args()


def read_bed(bed_file):
    if not bed_file:
        variants = []
    else:
        bed = pd.read_csv(bed_file, sep="\t", comment="#", header=["#CHROM", "POS", "REF", "ALT"])
        bed["DNAchange"] = bed['#CHROM'] + ":g." + bed['POS'].map(str) + bed['REF'] + ">" + bed['ALT']
        variants = bed["DNAchange"].tolist()
    return variants


def read_maf(maf_file):
    maf = pd.read_csv(maf_file, sep="\t", comment="#")
    maf["DNAchange"] = maf['Chromosome'] + ":g." + maf['Start_Position'].map(str) + maf['Reference_Allele'] + ">" + maf[
        'Tumor_Seq_Allele2']
    return maf


def filtering(maf, gnomad_thr, whitelist, blacklist, filters):
    if "PASS" not in filters:
        filters += ["PASS"]
    af_columns = [x for x in maf.columns if x.endswith("_AF") and x != "t_AF" and x != "n_AF"]
    maf["KEEP"] = maf['DNAchange'].isin(whitelist)
    for af_col in af_columns:
        maf[af_col] = maf[af_col].astype(float)
        maf[af_col] = maf[af_col].fillna(0)
        maf_filtered = maf[((maf[af_col] < gnomad_thr) & (maf["KEEP"] == False))]
    maf_filtered = maf_filtered[((~maf_filtered['DNAchange'].isin(blacklist)) & (maf["KEEP"] == False))]
    maf_filtered = maf_filtered[(maf_filtered["FILTER"].isin(filters) & (maf["KEEP"] == False))]
    # print(maf_filtered[["t_alt_count", "n_alt_count", "n_AF", "t_ref_count"]])
    # print(maf_filtered["t_alt_count"].astype(int))
    maf_filtered = maf_filtered[
                                ((maf_filtered["t_alt_count"].astype(float) >= 3) &
                                 (maf_filtered["n_alt_count"].astype(float) == 0) &
                                 (maf_filtered["n_AF"].astype(float) <= 0.01) &
                                 (maf_filtered["t_ref_count"].astype(float) >= 3))
                                ]
    return maf_filtered


def write_maf(maf_df, mafin_file, mafout_file):
    header_lines = subprocess.getoutput(f"zgrep -Eh '#|Hugo_Symbol' {mafin_file} 2>/dev/null")
    with open(mafout_file, "w") as mafout:
        mafout.write(header_lines)
    maf_df.to_csv(mafout_file, mode="a", index=False, header=True, sep="\t")


def main():
    args = argparser()
    maf = read_maf(args.input)
    filtered_maf = filtering(maf=maf, gnomad_thr=args.gnomad_thr,
                             whitelist=read_bed(args.whitelist), blacklist=read_bed(args.blacklist),
                             filters=args.filters)
    if not args.output:
        args.output = args.input.replace(".maf", "filtered.maf")
    write_maf(maf_df=filtered_maf, mafin_file=args.input, mafout_file=args.output)


if __name__ == '__main__':
    main()
