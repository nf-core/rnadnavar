# !/usr/bin/env python3
"""
Date: 08 Oct 2022
Author: Raquel Manzano - Caldas lab
Script: Add info drom run.consensus.py output to a maf file
"""
import argparse
import pandas as pd
import subprocess


def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--maf")
    parser.add_argument("--txt")
    parser.add_argument("--status", choices=[1, 2])
    parser.add_argument("--status_name", choices=["DNA", "RNA"])
    return parser.parse_args()


def match_and_add(maf, txt, status_name):
    maf_out = maf.replace(".maf", ".consensus.maf")
    head_maf = subprocess.getoutput(f"zgrep '#' {maf}")
    maf = pd.read_csv(maf, sep="\t", comment="#")
    txt = pd.read_csv(txt, sep="\t", comment="#")
    txt.drop(labels="FILTER", inplace=True, axis=1)
    merged = pd.merge(maf, txt, on=['DNAchange'], how='left')
    merged["consensus_callers"] = merged["consensus_callers"].replace({"consensus": f"in{status_name}"}, regex=True)
    with open(maf_out, "w") as out:
        out.write(head_maf + "\n")
    merged.fillna("").to_csv(maf_out, sep="\t", mode="a", index=False, header=True)


def main():
    args = argparser()
    match_and_add(maf=args.maf, txt=args.txt, status_name=args.status_name)


if __name__ == '__main__':
    main()
