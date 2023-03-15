#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Filters MAF file with filters from variant calling, gnomad and whitelist/blacklist.
TODO: investigate clustered_events in RNA
TODO: adjust filters on confident variants (known variants)
TODO: remove dubious variants from masked regions
TODO: check I am retaining as much info as possible from the vcf
TODO: add filters as options
"""
import argparse
import pandas as pd
import subprocess
import pysam


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
    parser.add_argument("--ref", help="FASTA reference file to extract context")
    return parser.parse_args()


def read_bed(bed_file):
    if not bed_file:
        variants = []
    else:
        bed = pd.read_csv(bed_file, sep="\t", comment="#", header=None)
        if len(bed.columns) == 5:
            colnames = ["#CHROM", "POS", "END", "REF", "ALT"]
        else:
            colnames = ["#CHROM", "POS", "REF", "ALT"]
        bed.columns = colnames
        bed["DNAchange"] = bed['#CHROM'] + ":g." + bed['POS'].map(str) + bed['REF'] + ">" + bed['ALT']
        variants = bed["DNAchange"].tolist()
    return variants


def read_maf(maf_file):
    maf = pd.read_csv(maf_file, sep="\t", comment="#")
    maf["DNAchange"] = maf['Chromosome'] + ":g." + maf['Start_Position'].map(str) + maf['Reference_Allele'] + ">" + maf[
        'Tumor_Seq_Allele2']
    return maf

def noncoding(maf, noncoding):
    maf["noncoding"] = maf["Consequence"].apply(
        lambda consequence: consequence.split('&')[0].split(',')[0] in noncoding)
    # change filter accordingly
    maf.loc[maf["noncoding"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["noncoding"])
    return maf


def remove_ig_and_psedo(maf, ig_list):
    """Add IG and pseudogene filters"""
    # fill na values
    maf[['BIOTYPE', 'SYMBOL']] = maf[['BIOTYPE', 'SYMBOL']].fillna(value="")
    # mask variants that were annotated as IG gene or pseudogene
    mask_ig = (maf["SYMBOL"].isin(ig_list))
    mask_ps = (maf["BIOTYPE"].str.contains('pseudogene'))
    # change filter accordingly
    maf.loc[mask_ig, 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["ig_gene"])
    maf.loc[mask_ps, 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["pseudogene"])
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
    maf.loc[maf["homopolymer"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["homopolymer"])
    return maf


def filtering(maf, gnomad_thr, whitelist, blacklist, filters):
    if "PASS" not in filters:
        filters += ["PASS"]
    af_columns = [x for x in maf.columns if x.endswith("_AF") and x != "t_AF" and x != "n_AF"]
    maf["whitelist"] = maf['DNAchange'].isin(whitelist)
    for af_col in af_columns:
        maf[af_col] = maf[af_col].astype(float)
        maf[af_col] = maf[af_col].fillna(0)
        maf_filtered = maf[((maf[af_col] < gnomad_thr) | (maf["whitelist"] == True))]
    maf_filtered = maf_filtered[((~maf_filtered['DNAchange'].isin(blacklist)) | (maf["whitelist"] == True))]
    maf_filtered = maf_filtered[(maf_filtered["FILTER"].isin(filters) | (maf["whitelist"] == True))]
    maf_filtered["KEEP"] = ((maf_filtered["t_alt_count"].astype(float) >= 3) &
                           (maf_filtered["n_alt_count"].astype(float) == 0) &
                           (maf_filtered["n_AF"].astype(float) <= 0.01) &
                           (maf_filtered["t_ref_count"].astype(float) >= 3))

    return maf_filtered


def write_maf(maf_df, mafin_file, mafout_file):
    header_lines = subprocess.getoutput(f"zgrep -Eh '#|Hugo_Symbol' {mafin_file} 2>/dev/null")
    with open(mafout_file, "w") as mafout:
        mafout.write(header_lines)
    maf_df.to_csv(mafout_file, mode="a", index=False, header=True, sep="\t")


def main():
    noncoding_list = ['intron_variant', 'intergenic_variant',
                      'non_coding_transcript_variant', 'non_coding_transcript_exon_variant',
                      'mature_miRNA_variant',
                      'IGR', 'INTRON', 'RNA']
    args = argparser()
    maf = read_maf(args.input)
    filtered_maf = filtering(maf=maf, gnomad_thr=args.gnomad_thr,
                             whitelist=read_bed(args.whitelist), blacklist=read_bed(args.blacklist),
                             filters=args.filters)
    # tag noncoding
    filtered_maf = noncoding(maf=filtered_maf, noncoding=noncoding_list)
    # tag IG and pseudo
    filtered_maf = remove_ig_and_psedo(maf=filtered_maf,
                                ig_list=['IGH', 'IGL', 'IGK', 'TRBV', 'TRAJ', 'pseudo'])
    # tag homopolymers
    filtered_maf = remove_homopolymers(maf=filtered_maf, ref=args.ref)
    if not args.output:
        args.output = args.input.replace(".maf", "filtered.maf")
    write_maf(maf_df=filtered_maf, mafin_file=args.input, mafout_file=args.output)


if __name__ == '__main__':
    main()
