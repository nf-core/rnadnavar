#!/usr/bin/env python
"""
Author: Raquel Manzano - @RaqManzano
Script: Convert MAF to BED format keeping ref and alt info
"""
import argparse
import pandas as pd


def argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-maf", "--mafin", help="MAF input file", required=True)
    parser.add_argument("-bed", "--bedout", help="BED input file", required=True)
    parser.add_argument(
        "--extra", help="Extra columns to keep (space separated list)", nargs="+", required=False, default=[]
    )
    return parser.parse_args()


def maf2bed(maf_file, bed_file, extra):
    maf = pd.read_csv(maf_file, sep="\t", comment="#")
    bed = maf[["Chromosome", "Start_Position", "End_Position"] + extra]
    bed.to_csv(bed_file, sep="\t", index=False, header=False)


def main():
    args = argparser()
    maf2bed(maf_file=args.mafin, bed_file=args.bedout, extra=args.extra)


if __name__ == "__main__":
    main()
