# !/usr/bin/env python3
"""
Author: Raquel Manzano - @RaqManzano
Script: Apply panel of normals filtering for RNA-seq data
"""
import argparse
import pandas as pd
from capy import mut


def argparser():
    parser = argparse.ArgumentParser(description='Apply panel of normals filtering for RNA-seq data')
    parser.add_argument("--maf")
    parser.add_argument("--pon")
    parser.add_argument("--thr", default=-2.8)
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--ref', help='Path to binary pon',
                        default='/rds/project/rds-70r4qFasPsQ/PBCP/Work/rm889/resources/hg19/hs37d5.fa.gz')

    return parser.parse_args()


def run_capy(maf_file, pon, ref, output, thr):
    """Runs capy with hg38 PoN"""
    print('- Running CApy HG38')
    M = pd.read_csv(maf_file, sep='\t')
    M = M[~M['Chromosome'].str.contains('_')]  # remove alt contigs
    M['chr'] = M['Chromosome'].str.replace('chr', '').replace('X', '23').replace('Y', '24')
    M['pos'] = M['Start_Position']
    M = M[M['chr'] != 'M']
    M['n_alt'] = M["alt_count"]
    M['n_ref'] = M["ref_count"]
    M["pon_score"] = mut.filter_mutations_against_token_PoN(M=M,
                                                            ponfile=pon,
                                                            ref=ref)
    M['pon_thr'] = M['pon_score_hg38'] >= thr
    M.to_csv(output, sep='\t', index=False)
    print(f'Done! See {output}')


def main():
    args = argparser()
    run_capy()


if __name__ == '__main__':
    main()
