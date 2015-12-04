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
import vcf

__version__ = 0.1
__desc__ = "Perform ensemble SNV calling based on multiple algorithms"


def main():
    pass


def parse_args(args=None):
    """Parse command-line arguments.

    Arguments:
        args (list): List of command-line arguments (for testing)

    Returns:
        Dictionary (dict) of argument-value pairs
    """
    parser = argparse.ArgumentParser(description=__desc__)
    parser.add_argument("vcf_files", n="+", help="VCF files from multiple algorithms")
    args_parsed = parser.parse_args(args)
    return args_parsed


if __name__ == '__main__':
    main()
