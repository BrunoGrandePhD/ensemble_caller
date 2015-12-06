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

import logging
import sys
import argparse
import os
from collections import OrderedDict
import re
import vcf as pyvcf
import vcf.utils as pyvcf_utils

__version__ = 0.1
__desc__ = "Perform ensemble SNV calling based on multiple algorithms"


def main():
    """Main function."""
    # Parse command-line arguments
    args = parse_args()
    # Sort checking
    if not args.skip_sort_check:
        are_sorted(args.vcf_files)
        reset_vcf_files(args.vcf_files)
    # Create walk_together generator
    vcf_wt = create_vcf_walktogether(*args.vcf_files)
    # Extract source (method) names
    names = args.names or extract_names(args.vcf_files)
    # Iterate over records
    for record in vcf_wt:
        print record


def setup_logging():
    """Setup basic config for logging."""
    log_format = '%(asctime)s - %(levelname)s (%(module)s.%(funcName)s):  %(message)s'
    date_format = '%Y/%m/%d %H:%M:%S'  # 2010/12/12 13:46:36
    logging.basicConfig(format=log_format, level=logging.INFO, datefmt=date_format,
                        stream=sys.stderr)


def parse_args(args=None):
    """Parse command-line arguments, validates them and prepares them.

    Arguments:
        args: List of command-line arguments (for testing)

    Returns:
        Dictionary of argument-value pairs
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description=__desc__)
    parser.add_argument("vcf_files", nargs="+", help="VCF files from multiple algorithms")
    parser.add_argument("--skip_sort_check", "-s", action="store_true", help="Skip sort check")
    parser.add_argument("--names", "-n", help="Method names (comma-separated) for VCF files "
                        "(same order and length)")
    args = parser.parse_args(args)
    # Confirm that all VCF files exist
    if not all(os.path.exists(vcf) for vcf in args.vcf_files):
        raise IOError("One of the specified VCF files doesn't exist.")
    # Convert VCF files into pyvcf objects
    args.vcf_files = [pyvcf.Reader(open(f)) for f in args.vcf_files]
    # Split names on the comma
    if args.names:
        args.names = args.names.split(",")
    # Check that names is the same length as the number of VCF files
    if args.names:
        assert len(args.names) == len(args.vcf_files), "Number of names and VCF files don't match."
    return args


def reset_vcf_files(vcf_files):
    """Reset the VCF objects to the beginning.

    Arguments:
        vcf_files: List of `vcf.Reader` objects
    """
    for vcf in vcf_files:
        vcf._reader.seek(0)  # Reset the underlying file object to zero
        vcf.reader = (line.strip() for line in vcf._reader if line.strip())  # Recreate generator
        vcf._parse_metainfo()  # Reparse the meta-information to skip ahead to the VCF records


class NotSortedException(Exception):
    """Input VCF file isn't sorted"""


def are_sorted(vcf_files):
    """Check if the VCF files are sorted. If not, issue a warning
    to the user that data might be lost during VCF traversal.

    Arguments:
        vcf_files: List of `vcf.Reader` objects

    Returns:
        A boolean indicating whether the following is true:
        - The chromosome order is the same in all VCF files
        - The positions within each chromosome are numerically
          sorted
    """
    chrom_lists = []
    for vcf in vcf_files:
        chrom_lists.append(parse_order(vcf))
    is_same_order = compare_orders(chrom_lists)
    if not is_same_order:
        logging.warn("The chromosomes across VCF files aren't all sorted "
                     "lexicographically or numerically. Please ensure that "
                     "they are in the same order. Use --skip_sort_check to "
                     "suppress this warning.")
    return is_same_order


def parse_order(vcf_file):
    """Parse a VCF file to obtain a list of chromosomes as ordered
    in the VCF file. Raises a NotSortedException when two consecutive
    records from the same chromosome aren't in order. Also raises the
    same exception if there is more than one block for a given chromosome.

    Arguments:
        vcf_file: A `vcf.Reader` object

    Returns:
        A list of chromosomes as ordered in the VCF file, with
        consecutive duplicates removed.
    """
    chroms_seen = OrderedDict()
    chrom_prev = None
    pos_prev = None
    for record in vcf_file:
        chrom, pos = record.CHROM, record.POS
        if chrom != chrom_prev:
            if chrom in chroms_seen:
                raise NotSortedException("More than one consecutive block "
                                         "for chromosome {}".format(chrom))
            else:
                chroms_seen[chrom] = True
                chrom_prev = chrom
                pos_prev = 0
        else:
            if pos < pos_prev:
                raise NotSortedException("Unsorted positions within chromosome "
                                         "{} around position {}".format(chrom, pos))
            pos_prev = pos
    return chroms_seen.keys()


def compare_orders(lists):
    """Compare the order of different lists while allowing for
    missing values. Raises a NotSortedException if the elements
    aren't in the same order across the lists.

    Checks whether all lists are sorted either lexicographically or
    numerically in increasing order.

    Arguments:
        lists: List of lists consisting of objects that
            can be ordered

    Returns:
        A boolean indicated whether the lists follow the same order
            while allowing for missing values
    """

    def check_lexico(l):
        """Check if given list l is sorted lexicographically"""
        return all(l[i] <= l[i+1] for i in range(len(l)-1))

    def check_num(l):
        """Check if given list l is sorted numerically"""
        l_split = [re.findall(r"([A-Za-z._]+)(\d+)|([A-Za-z._]+)|(\d+)", str(s)) for s in l]
        for i in range(len(l_split)):
            for j in range(len(l_split[i])):
                l_split[i][j] = [x for x in l_split[i][j] if x != ""]
                for k in range(len(l_split[i][j])):
                    try:
                        l_split[i][j][k] = int(l_split[i][j][k])
                    except ValueError:
                        pass
        return all(l_split[i] <= l_split[i+1] for i in range(len(l_split)-1))

    is_sorted_lexico = all(check_lexico(l) for l in lists)
    is_sorted_num = all(check_num(l) for l in lists)
    return is_sorted_lexico or is_sorted_num


def create_vcf_walktogether(vcf_files):
    """Create a VCF walk-together generator.

    Arguments:
        vcf_files: A list of `vcf.Reader` objects

    Returns:
        A `vcf.utils.walk_together` generator
    """
    def vcf_record_sort_key(r):
        return r.CHROM, r.POS, r.REF, r.ALT
    return pyvcf_utils.walk_together(*vcf_files, vcf_record_sort_key=vcf_record_sort_key)


def extract_names(vcf_files):
    """Extract source (method) names from the meta-information.
    If the source name is absent, a name is given based on its order.

    Arguments:
        vcf_files: List of `vcf.Reader` objects

    Returns:
        List of source (method) names
    """
    names = []
    for i, vcf in enumerate(vcf_files, 1):
        # Attempt to extract it from VCF file
        try:
            name = vcf.metadata["source"][0]
        except KeyError:
            name = ""
        # If not in VCF file or empty, set it generically
        if name == "":
            name = "method_{}".format(str(i))
        names.append(name)
    return names


if __name__ == '__main__':
    main()
