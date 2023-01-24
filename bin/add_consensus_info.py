# !/usr/bin/env python3
"""
Date: 04 Oct 2022
Author: Raquel Manzano - Caldas lab
Script: Add consensus info to maf
"""
import argparse


def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--vcfs", help="vcfs that were used for the consensus")
    parser.add_argument("--callers", help="Name of the caller per vcf")
    parser.add_argument("--maf_prefiltering", help="MAF file to add info before filtering")
    parser.add_argument("--maf_posfiltering", help="MAF file to add info after filtering")
    return parser.parse_args()


def extract_info_from_vcf(vcf_file):
    vcf_dict = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if not line.startswith("#"):
                values = line.split("\t")
                var_id = values[2]
                var_info = values[7]
                info_values = var_info.split(";")
                info_dict = {"overlap_type": "", "overlap_id": "", "consensus_count": "", "consensus_callers": ""}
                for val in info_values:
                    val_name, val_value = val.split("=")
                    if val_name in info_dict.keys():
                        info_dict[val_name] = val_value
                vcf_dict[var_id] = info_dict

    return vcf_dict


def main():
    args = argparser()
    function()


if __name__ == '__main__':
    main()
