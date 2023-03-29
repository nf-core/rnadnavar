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
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--maf", help="Input MAF file")
    parser.add_argument("--maf_2pass", help="Input MAF file after a second pass")
    parser.add_argument('--pon19', help='Path to binary pon',
                        default='/rds/project/rds-70r4qFasPsQ/PBCP/Work/rm889/projects/rna_mutectv2/inputs/GTeX/77759/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WES.panel-of-normals.c1/PoN_GTEx')
    parser.add_argument('--pon', help='Path to binary pon',
                        default='/rds/project/rds-70r4qFasPsQ/PBCP/Work/rm889/projects/rna_mutectv2/inputs/GTeX/77759/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WES.panel-of-normals.c1/PoN_GTEx')
    parser.add_argument("--whitelist", help="Input MAF file after a second pass")
    parser.add_argument("--output", help="output file name")
    parser.add_argument("--out_suffix", help="Suffix for intermidiate output", default="withRNAfilters")
    parser.add_argument('--ref',
                        help='FASTA - HG38 (GDC)')
    parser.add_argument('--ref19',
                        help='FASTA - HG19 (hs37d5)')
    parser.add_argument('--rnaedits',
                        help="BED file(s) with known RNA editing events separated by space",
                        nargs="+")
    parser.add_argument("--thr", default=-2.8)
    parser.add_argument('--chain', help='Chain file',
                        default='/rds/project/rds-70r4qFasPsQ/PBCP/Work/rm889/resources/hg38ToHg19.over.chain.gz')

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


def add_rnaediting_sites(maf, rnaeditingsites, realignment):
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
                ravex_filter += ["rnaediting"]
        if row["rnaediting"]:
            ravex_filter += ["rnaediting"]
        for pon_col in pon_cols:
            if row[pon_col]:
                ravex_filter += ["rna_pon" + pon_col[-5:]]
        if not ravex_filter:
            ravex_filter = ["PASS"]
        ravex_filter = ";".join(ravex_filter)
        maf.at[idx, "RaVeX_FILTER"] = ravex_filter
    return maf


def extract_coords(df):
    """Extract coordinates from maf and stores them in a bed format file"""
    coords = df["Chromosome"] + ":" + df["Start_Position"].astype(str) + "-" + df["Start_Position"].astype(str)
    return coords


def run_capy(M, pon, ref, thr, chroms, suffix="_hg38"):
    """Runs capy with hg38 PoN"""
    print('- Running CApy' + suffix)
    chrom = 'Chromosome'
    start = 'Start_Position'
    if "hg19" in suffix:
        chrom += "19"
        start += "19"
    M = M[M[chrom].isin(chroms)]  # remove alt contigs
    M['chr'] = M[chrom].str.replace('chr', '').replace('X', '23').replace('Y', '24')
    M['pos'] = M[start]
    M = M[M['chr'] != 'M']
    M['n_alt'] = M["t_alt_count"]
    M['n_ref'] = M["t_ref_count"]
    M = M.astype({'chr': 'int32', 'pos': 'int32'})
    M["pon_score" + suffix] = mut.filter_mutations_against_token_PoN(M=M,
                                                                     ponfile=pon,
                                                                     ref=ref)
    M['pon_thr' + suffix] = M['pon_score' + suffix] >= thr
    return M


def add_hg19_coords_with_liftover(M, chain_file):
    converter = ChainFile(chain_file, 'hg38', 'hg19')
    starting_size = M.shape[0]
    M["coordinates19"] = M.apply(lambda x: converter[x["Chromosome"]][x["Start_Position"]], axis=1)
    tmp = pd.DataFrame(M['coordinates19'].tolist(), index=M.index)
    tmp[['Chromosome19', 'Start_Position19', "STRAND19"]] = pd.DataFrame(tmp[0].tolist(), index=M.index)
    liftover_size = tmp.shape[0]
    assert starting_size == liftover_size, "Liftovered positions are not the same numer as in original MAF"
    M = pd.concat([M, tmp], axis=1).drop(0, axis=1)
    return M


def write_output(args, results, output, out_suffix):
    if len(results.keys()) > 1:
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
    else:
        results[0].sort_values(["Chromosome", "Start_Position"]).to_csv(output,
                                                                        sep="\t",
                                                                        index=False,
                                                                        header=True)
        print(f"See: {output}")


def main():
    args = argparser()
    if "," in args.rnaedits[0]:
        args.rnaedits = args.rnaedits[0].strip().split(",")
    # chromosomes
    chroms = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']]

    # realignment
    calls_1pass = pd.read_csv(args.maf, sep="\t", comment="#")
    # If REALIGNMENT provide intersect
    if args.maf_2pass and args.maf != args.maf_2pass:
        didrealignment = True
        calls_2pass = pd.read_csv(args.maf_2pass, sep="\t", comment="#")
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
        # RNA panel of normals
        if args.pon19 and args.chain and args.ref19:
            calls = add_hg19_coords_with_liftover(calls, chain_file=args.chain)
            calls = run_capy(M=calls, pon=args.pon19, ref=args.ref19, thr=args.thr, suffix="_hg19", chroms=chroms)
        calls = run_capy(M=calls, pon=args.pon, ref=args.ref, thr=args.thr, suffix="_hg38", chroms=chroms)
        # Annotate known RNA editing
        rnadbs = []
        for rnadb_file in args.rnaedits:
            rnadb = pd.read_csv(rnadb_file, sep="\t", names=["chr", "start", "end", "ref", "alt"], header=None)
            rnadb["DNAchange"] = rnadb["chr"] + ":g." + \
                                 rnadb["start"].map(str) + \
                                 rnadb["ref"] + ">" + rnadb["alt"]
            rnadbs += [rnadb]
        rnadbs = pd.concat(rnadbs)
        calls = add_rnaediting_sites(maf=calls, rnaeditingsites=rnadbs, realignment=didrealignment)
        results[idx] = calls
    # write maf files
    write_output(args, results, args.output, args.out_suffix)


if __name__ == '__main__':
    main()
