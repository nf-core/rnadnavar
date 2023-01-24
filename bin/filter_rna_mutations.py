#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Filters MAF file with realignment (could it run with the consensus?), noncoding(?), homopolymers and RNA editing database
"""
import argparse
import pandas as pd
# import pysam

pd.options.mode.chained_assignment = None  # default='warn'


def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--maf", help="Input MAF file")
    parser.add_argument("--maf_2pass", help="Input MAF file after a second pass")
    parser.add_argument("--whitelist", help="Input MAF file after a second pass")
    parser.add_argument("--output", help="output file name")
    parser.add_argument("--out_suffix", help="Suffix for intermidiate output", default="withRNAfilters")
    parser.add_argument('--ref',
                        help='FASTA - out-of-date?')
    parser.add_argument('--darned',
                        default='')
    parser.add_argument('--radar',
                        default='')
    parser.add_argument('--redi',
                        default='')
    parser.add_argument('--nat',
                        default='')
    return parser.parse_args()


def realignment(maf1, maf2):
    """Get variants that are in both"""
    maf1["DNAchange"] = maf1["Chromosome"] + ":g." + \
                        maf1["Start_Position"].map(str) + \
                        maf1["Reference_Allele"] + ">" + maf1["Tumor_Seq_Allele2"]
    maf2["DNAchange"] = maf2["Chromosome"] + ":g." + \
                        maf2["Start_Position"].map(str) + \
                        maf2["Reference_Allele"] + ">" + maf2["Tumor_Seq_Allele2"]
    maf1["realignment"] = maf1["DNAchange"].isin(maf2["DNAchange"])
    maf2["realignment"] = maf2["DNAchange"].isin(maf1["DNAchange"])

    maf_intersect = maf1[maf1["realignment"] == True]
    maf1 = maf1[~maf1["DNAchange"].isin(maf_intersect["DNAchange"])]
    maf2 = maf2[~maf2["DNAchange"].isin(maf_intersect["DNAchange"])]

    return maf1, maf2, maf_intersect


def add_to_filter(current_filter, new_filter):
    if new_filter:
        if current_filter == "PASS":
            current_filter = new_filter[:]
        else:
            current_filter += f';{new_filter}'
    return current_filter


def update_filter(current_filter, new_filter):
    if current_filter != "PASS":
        new_filter = current_filter + ";" + new_filter
    return new_filter


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


def noncoding(maf, noncoding):
    maf["noncoding"] = maf["Consequence"].apply(
        lambda consequence: consequence.split('&')[0].split(',')[0] in noncoding)
    # change filter accordingly
    maf.loc[maf["noncoding"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=["noncoding"])
    return maf


def remove_rnaediting(maf, database, db_name):
    # bed format expected
    rnadb = pd.read_csv(database, sep="\t", names=["chr", "start", "end", "ref", "alt"], header=None)
    rnadb["DNAchange"] = rnadb["chr"] + ":g." + \
                         rnadb["start"].map(str) + \
                         rnadb["ref"] + ">" + rnadb["alt"]
    maf["rnaediting"] = maf["DNAchange"].isin(rnadb["DNAchange"])
    # change filter accordingly
    maf.loc[maf["rnaediting"], 'FILTER'] = maf["FILTER"].apply(add_to_filter, args=[db_name])
    return maf


def write_output(args, results, output, out_suffix):
    maf_out = args.maf.replace(".maf", f".{out_suffix}.maf").replace(".gz", "").split("/")[-1]
    maf2_out = args.maf_2pass.replace(".maf", f".{out_suffix}.maf").replace(".gz", "").split("/")[-1]
    maf12_out = output[:]
    pd.concat([results[0], results[2]]).sort_values(["Chromosome",
                                                     "Start_Position"]).to_csv(maf_out,
                                                                               sep="\t",
                                                                               index=False,
                                                                               header=True)
    pd.concat([results[1], results[2]]).sort_values(["Chromosome",
                                                     "Start_Position"]).to_csv(maf2_out,
                                                                               sep="\t",
                                                                               index=False,
                                                                               header=True)
    results[2].sort_values(["Chromosome", "Start_Position"]).to_csv(maf12_out,
                                                                    sep="\t",
                                                                    index=False,
                                                                    header=True)
    print(f"See:\n-{maf_out}\n-{maf2_out}\n-{maf12_out}")


def main():
    args = argparser()
    # should I do non-coding?
    # realignment
    # noncoding_list = ['intron_variant', 'intergenic_variant',
    #                   'non_coding_transcript_variant', 'non_coding_transcript_exon_variant',
    #                   'mature_miRNA_variant',
    #                   'IGR', 'INTRON', 'RNA']
    calls_1pass = pd.read_csv(args.maf, sep="\t", comment="#")
    # second maf from a realignment is optional
    if args.maf_2pass and args.maf != args.maf_2pass:
        calls_2pass = pd.read_csv(args.maf_2pass, sep="\t", comment="#")
        calls1, calls2, calls12 = realignment(calls_1pass, calls_2pass)

        calls = [calls1, calls2, calls12]
    else:
        calls = [calls_1pass]

    results = {}
    for idx, calls in enumerate(calls):
        if calls.empty:
            results[idx] = calls
            continue
        # # tag noncoding
        # noncoding(maf=calls, noncoding=noncoding_list)
        # # tag IG and pseudo
        # calls = remove_ig_and_psedo(maf=calls,
        #                             ig_list=['IGH', 'IGL', 'IGK', 'TRBV', 'TRAJ', 'pseudo'])
        # # tag homopolymers
        # calls = remove_homopolymers(maf=calls, ref=args.ref)
        # tag rnaediting events
        for rnadb, namedb in zip([args.darned, args.radar, args.redi, args.nat],
                                 ['darned', 'radar', 'redi', 'natEd']):
            calls = remove_rnaediting(maf=calls, database=rnadb, db_name=namedb)
        results[idx] = calls
    # write maf files
    write_output(args, results, args.output, args.out_suffix)


if __name__ == '__main__':
    main()
