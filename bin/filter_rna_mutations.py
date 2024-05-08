#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Filters MAF file with realignment (could it run with the consensus?), noncoding(?), homopolymers and RNA editing database
"""
import argparse
import pandas as pd
from liftover import ChainFile
from capy import mut

pd.options.mode.chained_assignment = None  # default='warn'


def argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--maf", help="Input MAF file")
    parser.add_argument("--maf_realign", help="Input MAF file after a second pass")
    parser.add_argument("--pon2", help="Path to optional second binary pon")
    parser.add_argument("--pon", help="Path to binary pon")
    parser.add_argument("--whitelist", help="Input MAF file after a second pass")
    parser.add_argument("--output", help="output file name")
    parser.add_argument("--out_suffix", help="Suffix for intermediate output", default="withRNAfilters")
    parser.add_argument("--ref", help="FASTA - e.g. HG38")
    parser.add_argument("--refname", help="e.g. hg38", default="hg38")
    parser.add_argument("--ref2", help="FASTA - e.g. hg19 (hs37d5)")
    parser.add_argument("--refname2", help="e.g. HG19", default="hg19")
    parser.add_argument("--rnaedits", help="BED file(s) with known RNA editing events separated by space", nargs="+")
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


def add_rnaediting_sites(maf, rnaeditingsites, realignment):
    """
    Check for RNA editing sites in the MAF table
    """
    # bed format expected
    maf["rnaediting"] = maf["DNAchange"].isin(rnaeditingsites["DNAchange"])
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
        if not ravex_filter or row["whitelist"]:
            ravex_filter = ["PASS"]
        ravex_filter = ";".join(ravex_filter)
        maf.at[idx, "RaVeX_FILTER"] = ravex_filter
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
        chrom += refname2
        start += refname2
    M = M[M[chrom].isin(chroms)]  # remove alt contigs
    M["chr"] = M[chrom].str.replace("chr", "").replace("X", "23").replace("Y", "24")
    M["pos"] = M[start]
    M = M[M["chr"] != "M"]
    M["n_alt"] = M["t_alt_count"]
    M["n_ref"] = M["t_ref_count"]
    M = M.astype({"chr": "int32", "pos": "int32"})
    pon_scores = mut.filter_mutations_against_token_PoN(M=M, ponfile=pon, ref=ref)
    assert len(pon_scores) == M.shape[0], "PoN scores are not same length as nrow"
    M["pon_score" + suffix] = pon_scores
    M["pon_thr" + suffix] = M["pon_score" + suffix] >= thr
    return M


def add_coords2_with_liftover(M, chain_file,ref1="hg38",ref2="hg19"):
    """
    Liftover ref1 coordinates from maf to ref2
    """
    converter = ChainFile(chain_file, ref1, ref2)
    na_nr = M[M["Chromosome"].isnull()].shape[0]
    if na_nr > 0:
        print(f"[WGN] Removing {na_nr} variants where Chromosome is NA")
        M = M[
            ~M["Chromosome"].isna()
        ].reindex()  # remove positions where coordinates are not present (maybe liftover went wrong)

    starting_size = M.shape[0]
    M["coordinates_"+ref2] = M.apply(lambda x: converter[x["Chromosome"]][x["Start_Position"]], axis=1)
    # Replace with tuple of None's when position has been deleted in new reference genome
    tmp = pd.DataFrame(M["coordinates19"].tolist(), index=M.index).apply(
        lambda ds: ds.map(lambda x: x if x != None else (None, None, None))
    )
    tmp[["Chromosome_"+ref2, "Start_Position_"+ref2, "STRAND_"+ref2]] = pd.DataFrame(tmp[0].tolist(), index=M.index)
    liftover_size = tmp.shape[0]
    assert starting_size == liftover_size, "Liftovered positions are not the same number as in original MAF"
    M = pd.concat([M, tmp], axis=1).drop(0, axis=1)
    return M


def write_output(args, results, output, out_suffix):
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
        print(f"See:\n-{maf_out}\n-{maf2_out}\n-{maf12_out}")
    else:
        results[0].sort_values(["Chromosome", "Start_Position"]).to_csv(output, sep="\t", index=False, header=True)
        print(f"See: {output}")


def main():
    args = argparser()
    if args.rnaedits:
        if "," in args.rnaedits[0]:
            args.rnaedits = args.rnaedits[0].strip().split(",")
    # chromosomes
    chroms = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]]

    # realignment
    calls_1pass = pd.read_csv(args.maf, sep="\t", comment="#")
    # If REALIGNMENT provide intersect
    if args.maf_realign and args.maf != args.maf_realign:
        didrealignment = True
        calls_2pass = pd.read_csv(args.maf_realign, sep="\t", comment="#")
        calls1, calls2, calls12 = realignment(calls_1pass, calls_2pass)
        calls = [calls1, calls2, calls12]
    else:
        didrealignment = False
        calls = [calls_1pass]
    # Annotate known RNA editing
    results = {}
    for idx, calls in enumerate(calls):
        if calls.empty:
            results[idx] = calls
            continue
        if "Chromosome"+args.refname2 in calls.columns:
            calls.drop(["Chromosome"+args.refname2, "Start_Position"+args.refname2], axis=1, inplace=True)
        # RNA panel of normals
        if args.pon2 and args.chain and args.ref2:
            calls = add_coords2_with_liftover(calls, chain_file=args.chain,ref1=args.refname,ref2=args.refname2)
            calls = run_capy(M=calls, pon=args.pon2, ref=args.ref2, thr=args.thr, suffix='_'+args.refname2, chroms=chroms)
        if args.pon:
            calls = run_capy(M=calls, pon=args.pon, ref=args.ref, thr=args.thr, suffix='_'+args.refname, chroms=chroms)
        # Annotate known RNA editing
        if args.rnaedits:
            rnadbs = []
            for rnadb_file in args.rnaedits:
                rnadb = pd.read_csv(rnadb_file, sep="\t", names=["chr", "start", "end", "ref", "alt"], header=None)
                rnadb["DNAchange"] = rnadb["chr"] + ":g." + rnadb["start"].map(str) + rnadb["ref"] + ">" + rnadb["alt"]
                rnadbs += [rnadb]
            rnadbs = pd.concat(rnadbs)
            calls = add_rnaediting_sites(maf=calls, rnaeditingsites=rnadbs, realignment=didrealignment)
        results[idx] = calls
    # write maf files
    write_output(args, results, args.output, args.out_suffix)


if __name__ == "__main__":
    main()
