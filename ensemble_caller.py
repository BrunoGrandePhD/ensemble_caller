#!/usr/bin/env python2

"""
ensemble_caller.py
==================
This script implements the recommendation described and supported
in Ewing at al. (doi: 10.1038/nmeth.3407), namely to perform ensemble
single nucleotide variant (SNV) calling with a voting scheme.
Specifically, given sets of SNV calls from various algorithms,
this script will match up corresponding variants and output SNV calls
that are supported by a majority of algorithms.

See `python ensemble_caller.py --help` for a list of arguments.

Dependencies
------------
- python>=2.7
- pyvcf>=0.6.7

Known Issues
------------
- None
"""

import argparse
import vcf as pyvcf

__version__ = 0.1
__desc__ = "Perform ensemble SNV calling based on multiple algorithms"


def main():
    """Perform ensemble SNV calling."""
    # Parse command-line arguments
    args = parse_args()
    # Sort checking
    if not args.skip_sort_check:
        pass


def parse_args(args=None):
    """Parse command-line arguments and prepare them.

    Arguments:
        args (list): List of command-line arguments (for testing)

    Returns:
        Dictionary (dict) of argument-value pairs
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description=__desc__)
    parser.add_argument("vcf_files", nargs="+", type=argparse.FileType("r"), help="VCF files from multiple algorithms")
    parser.add_argument("--skip_sort_check", "-s", action="store_true", help="Skip sort check")
    args_parsed = parser.parse_args(args)
    # Convert VCF files into pyvcf objects
    args_parsed.vcf_files = [pyvcf.Reader(f) for f in args_parsed.vcf_files]
    return args_parsed


def reset_vcf_files(vcf_files):
    """Reset the VCF objects to the beginning.

    Arguments:
        vcf_files (list): List of `vcf.Reader` objects

    Returns:
        None
    """
    for vcf in vcf_files:
        vcf.reader.seek(0)


def check_sort(vcf_files):
    pass


if __name__ == '__main__':
    main()
