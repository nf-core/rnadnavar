#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Filters MAF file with realignment (could it run with the consensus?), noncoding(?), homopolymers and RNA editing database
"""
print("Note that importing CApy, might throw some unnecessary messages.")
import argparse
import pandas as pd
import warnings
# CApy import throws FutureWarnings and unnecessary messages
warnings.simplefilter(action='ignore', category=FutureWarning)
from capy import mut # type: ignore
from liftover import ChainFile
print("-- done importing")


def argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--maf", help="Input MAF file")
    parser.add_argument("--maf_realign", help="Input MAF file after a second pass")
    parser.add_argument("--pon2", help="Path to optional second binary pon")
    parser.add_argument("--pon", help="Path to binary pon")
    parser.add_argument("--whitelist", help="BED file with variants to keep (CHROM POS REF ALT)")
    parser.add_argument("--output", help="output file name")
    parser.add_argument("--out_suffix", help="Suffix for intermediate output", default="withRNAfilters")
    parser.add_argument("--ref", help="FASTA - e.g. HG38")
    parser.add_argument("--refname", help="e.g. hg38", default="hg38")
    parser.add_argument("--ref2", help="FASTA - e.g. hg19 (hs37d5)")
    parser.add_argument("--refname2", help="e.g. HG19", default="hg19")
    parser.add_argument("--rnaedits", help="BED file(s) with known RNA editing events separated by space", nargs="+", default=[])
    parser.add_argument("--thr", default=-2.8)
    parser.add_argument("--chain", help="Chain file")

    return parser.parse_args()


def realignment(maf1, maf2):
    """
    Get variants that intersect both mafs without taking into account potential consensus variants
    """
    # Create the 'DNAchange' column
    for maf in [maf1, maf2]:
        maf["DNAchange"] = (
            maf["Chromosome"]
            + ":g."
            + maf["Start_Position"].astype(str)
            + maf["Reference_Allele"]
            + ">"
            + maf["Tumor_Seq_Allele2"]
        )
    # Filter out consensus 'DNAchange' values from both mafs (consensus is unrelated to realignment)
    maf1_non_consensus = maf1[~maf1["Caller"].str.contains("consensus", case=False)]
    maf2_non_consensus = maf2[~maf2["Caller"].str.contains("consensus", case=False)]
    # Find intersection of non-consensus 'DNAchange' values
    is_intersect = maf1_non_consensus["DNAchange"].isin(maf2_non_consensus["DNAchange"])
    intersect_changes = maf1_non_consensus["DNAchange"][is_intersect]
    # Update 'realignment' columns based on the intersection
    maf1["realignment"] = maf1["DNAchange"].isin(intersect_changes)
    maf2["realignment"] = maf2["DNAchange"].isin(intersect_changes)
    # Subset to intersected variants excluding any consensus variants
    maf1_intersect = maf1[maf1["realignment"]]
    maf2_intersect = maf2[maf2["realignment"]]
    maf_intersect = pd.concat([maf1_intersect, maf2_intersect]).drop_duplicates(subset="DNAchange")
    return maf1, maf2, maf_intersect


def add_filters(maf, rnaeditingsites, realignment, whitelist):
    """
    Check for RNA editing sites in the MAF table
    """
    print("- Annotating RNA filters")
    if not rnaeditingsites:
        rnaeditingsites = pd.DataFrame()
    if not rnaeditingsites.empty:
        # if mut is in a editing position T>C; A>G; G>A; C>T
        maf = maf.assign(mut=maf["Reference_Allele"] + ">" + maf["Tumor_Seq_Allele2"])
        maf = maf.assign(chr_start=maf['Chromosome'] + maf['Start_Position'].astype(str))
        rnaeditingsites = rnaeditingsites.assign(chr_start = rnaeditingsites['chr'] + rnaeditingsites['start'].astype(str))
        rnaedits = pd.DataFrame({"rnaediting": maf['chr_start'].isin(rnaeditingsites['chr_start']) & maf['mut'].str.contains("^T>C$|^C>T$|^A>G$|^G>A$")})
        maf.drop("chr_start", axis=1,inplace=True)
        maf = pd.merge(maf, rnaedits, left_index=True, right_index=True)
    else:
        maf["rnaediting"] = False
    # change filter accordingly
    pon_cols = [x for x in maf.columns if "pon_thr" in x]
    for idx, row in maf.iterrows():
        ravex_filter = row["RaVeX_FILTER"].split(";")
        if ravex_filter == ["PASS"]:
            ravex_filter = []
        if realignment:
            if not row["realignment"]:
                ravex_filter += ["realignment"]
        if row["rnaediting"]:
            ravex_filter += ["rnaediting"]
        for pon_col in pon_cols:
            if row[pon_col]:
                ravex_filter += ["rna_pon" + pon_col[-5:]]
        if whitelist:
            if row["whitelist"]:
                ravex_filter = ["PASS"]
        if not ravex_filter:
            ravex_filter = ["PASS"]
        ravex_filter = ";".join(ravex_filter)
        maf.at[idx, "RaVeX_FILTER"] = ravex_filter
    # column cleanup
    maf.drop(["coordinates_hg19"], axis=1, inplace=True)
    return maf


def extract_coords(df):
    """
    Extract coordinates from maf and stores them in a bed format file
    """
    coords = df["Chromosome"] + ":" + df["Start_Position"].astype(str) + "-" + df["Start_Position"].astype(str)
    return coords


def run_capy(M, pon, ref, thr, chroms, suffix="_hg38", refname2="hg19"):
    """
    Runs CApy with RNA PoN generated with tokenizer (https://github.com/getzlab/aggregate_tokens_files_TOOL)
    """
    print("- Running CApy" + suffix)
    chrom = "Chromosome"
    start = "Start_Position"
    if refname2:
        chrom += "_" + refname2
        start += "_" + refname2
    # Remove alt contigs
    M = M[M[chrom].isin(chroms)]
    # Perform string replacements and type casting
    M = M.assign(
        chr=M[chrom].str.replace("chr", "").replace({"X": "23", "Y": "24"}),
        pos=M[start],
        n_alt=M["t_alt_count"],
        n_ref=M["t_ref_count"]
    )
    M = M[M["chr"] != "M"]
    M = M.astype({"chr": "int32", "pos": "int32"})
    M.reset_index(inplace=True, drop=True)
    # Get PoN scores
    pon_scores = mut.filter_mutations_against_token_PoN(M=M, ponfile=pon, ref=ref)
    assert len(pon_scores) == M.shape[0], "PoN scores are not same length as nrow"
    # Create a DataFrame with the necessary columns including PoN scores and threshold
    pon_scores_df = pd.DataFrame({
        "pon_score" + suffix: pon_scores,
        "pon_thr" + suffix: [score >= thr for score in pon_scores]
    })
    pon_scores_df.reset_index(inplace=True, drop=True)
    assert pon_scores_df.index.equals(M.index)
    # Ensure the threshold column is of boolean type
    pon_scores_df["pon_thr" + suffix] = pon_scores_df["pon_thr" + suffix].astype(bool)
    # Merge the original DataFrame with the PoN scores DataFrame
    M = M.merge(pon_scores_df, left_index=True, right_index=True)

    return M


def add_coords2_with_liftover(M, chain_file,ref1="hg38",ref2="hg19"):
    """
    Liftover ref1 coordinates from maf to ref2
    """
    print("- Liftovering coordinates as second PoN was provided with chain file")
    converter = ChainFile(chain_file, ref1, ref2)
    na_nr = M[M["Chromosome"].isnull()].shape[0]
    if na_nr > 0:
        print(f"[WGN] Removing {na_nr} variants where Chromosome is NA")
        M = M[
            ~M["Chromosome"].isna()
        ].reindex()  # remove positions where coordinates are not present (maybe liftover went wrong)

    starting_size = M.shape[0]
    M["coordinates_" + ref2] = M.apply(lambda x: converter[x["Chromosome"]][x["Start_Position"]], axis=1)
    # Replace with tuple of None's when position has been deleted in new reference genome
    tmp = pd.DataFrame(M["coordinates_"+ref2].tolist(), index=M.index).apply(
        lambda ds: ds.map(lambda x: x if x != None else (None, None, None))
    )
    tmp[["Chromosome_" + ref2, "Start_Position_" + ref2, "STRAND_" + ref2]] = pd.DataFrame(tmp[0].tolist(), index=M.index)
    liftover_size = tmp.shape[0]
    assert starting_size == liftover_size, "Liftovered positions are not the same number as in original MAF"
    M = pd.concat([M, tmp], axis=1).drop(0, axis=1)

    return M


def write_output(args, results, output, out_suffix):
    """Write filtered maf(s). If realignment was provided then there will be three files:
    one for each input file with the new filters annotated and a third one with the merged results.
    """
    if len(results.keys()) > 1:
        maf_out = args.maf.replace(".maf", f".{out_suffix}.maf").replace(".gz", "").split("/")[-1]
        maf2_out = args.maf_realign.replace(".maf", f".{out_suffix}.maf").replace(".gz", "").split("/")[-1]
        maf12_out = output[:]
        pd.concat([results[0], results[2]]).sort_values(["Chromosome", "Start_Position"]).to_csv(
            maf_out, sep="\t", index=False, header=True
        )
        pd.concat([results[1], results[2]]).sort_values(["Chromosome", "Start_Position"]).to_csv(
            maf2_out, sep="\t", index=False, header=True
        )
        results[2].sort_values(["Chromosome", "Start_Position"]).to_csv(maf12_out, sep="\t", index=False, header=True)
        print(f"See:\n- {maf_out}\n- {maf2_out}\n- {maf12_out}")
    else:
        results[0].sort_values(["Chromosome", "Start_Position"]).to_csv(output, sep="\t", index=False, header=True)
        print(f"See: {output}")

def check_rnaediting(rnaedits):
    print("- Obtaining RNA editing sites from bed(s)")
    rnadbs = []
    for rnadb_file in rnaedits:
        print(f" - Reading {rnadb_file}")
        rnadb = pd.read_csv(rnadb_file, sep="\s+", names=["chr", "start", "end", "ref", "alt"], header=None, low_memory=False)
        rnadbs += [rnadb.assign(DNAchange=rnadb["chr"] + ":g." + rnadb["start"].map(str) + rnadb["ref"] + ">" + rnadb["alt"])]
    rnadbs_concat = pd.concat(rnadbs)
    return rnadbs_concat

def main():
    args = argparser()
    print("- Running script")
    if args.rnaedits:
        if "," in args.rnaedits[0]:
            args.rnaedits = args.rnaedits[0].strip().split(",")
    # chromosomes
    chroms = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]]

    # realignment
    try:
        calls_1pass = pd.read_csv(args.maf, sep="\t", comment="#", low_memory=False)
    except pd.errors.EmptyDataError:
        print("[Error] No columns to parse from file, is your input empty?")
        exit()
    # If REALIGNMENT provide intersect
    if args.maf_realign and args.maf != args.maf_realign:
        didrealignment = True
        calls_2pass = pd.read_csv(args.maf_realign, sep="\t", comment="#",low_memory=False)
        calls1, calls2, calls12 = realignment(calls_1pass, calls_2pass)
        calls = [calls1, calls2, calls12]
    else:
        didrealignment = False
        calls = [calls_1pass]
    results = {}
    for idx, calls in enumerate(calls):
        print(f"> Annotating {idx+1}/3")
        if calls.empty:
            results[idx] = calls
            continue
        if "Chromosome"+args.refname2 in calls.columns:
            calls.drop(["Chromosome"+args.refname2, "Start_Position"+args.refname2], axis=1, inplace=True)
        # RNA panel of normals
        if args.pon2 is not None and args.chain is not None and args.ref2 is not None:
            calls = add_coords2_with_liftover(calls, chain_file=args.chain,ref1=args.refname,ref2=args.refname2)
            calls = run_capy(M=calls, pon=args.pon2, ref=args.ref2, thr=args.thr, suffix='_'+args.refname2, chroms=chroms, refname2=args.refname2)
        if args.pon:
            calls = run_capy(M=calls, pon=args.pon,  ref=args.ref,  thr=args.thr, suffix='_'+args.refname, chroms=chroms, refname2=False)
        # Annotate known RNA editing
        rnadbs = pd.DataFrame()
        if args.rnaedits:
            rnadbs = check_rnaediting(args.rnaedits)
        calls = add_filters(maf=calls, rnaeditingsites=rnadbs, realignment=didrealignment, whitelist=args.whitelist)
        results[idx] = calls
    # write maf files
    write_output(args, results, args.output, args.out_suffix)


if __name__ == "__main__":
    main()
