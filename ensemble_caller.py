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
import os

from collections import OrderedDict
from collections import deque
import vcf as pyvcf
import vcf.utils as pyvcf_utils

__version__ = 0.1
__desc__ = "Perform ensemble SNV calling based on multiple algorithms"


def main():
    """Perform ensemble SNV calling."""
    # Parse command-line arguments
    args = parse_args()
    # Sort checking
    if not args.skip_sort_check:
        are_sorted(args.vcf_files)
        reset_vcf_files(args.vcf_files)
    # Create walk_together generator
    vcf_wt = create_vcf_walktogether(*args.vcf_files)
    # Extract source (method) names
    if args.names:
        names = args.names
    else:
        names = extract_names(args.vcf_files)
    # Iterate over records
    for record in vcf_wt:
        print record


def parse_args(args=None):
    """Parse command-line arguments, validates them and prepares them.

    Arguments:
        args (list): List of command-line arguments (for testing)

    Returns:
        Dictionary (dict) of argument-value pairs
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
        vcf_files (list): List of `vcf.Reader` objects
    """
    for vcf in vcf_files:
        vcf._reader.seek(0)  # Reset the underlying file object to zero
        vcf.reader = (line.strip() for line in vcf._reader if line.strip())  # Recreate generator
        vcf._parse_metainfo()  # Reparse the meta-information to skip ahead to the VCF records


class NotSortedException(Exception):
    """Input VCF file isn't sorted"""


def are_sorted(vcf_files):
    """Check if the VCF files are sorted.

    Arguments:
        vcf_files (list): List of `vcf.Reader` objects

    Returns:
        A boolean indicating whether the following is true:
        - The chromosome order is the same in all VCF files
        - The positions within each chromosome are numerically
          sorted
    """
    chrom_lists = []
    for vcf in vcf_files:
        chrom_lists.append(parse_order(vcf))
    return compare_orders(chrom_lists)


def parse_order(vcf_file):
    """Parse a VCF file to obtain a list of chromosomes as ordered
    in the VCF file. Raises a NotSortedException when two consecutive
    records from the same chromosome aren't in order. Also raises the
    same exception if there is more than one block for a given chromosome.

    Arguments:
        vcf_file (`vcf.Reader`): A `vcf.Reader` object

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

    Arguments:
        lists (list): List of lists consisting of objects that
            can be ordered

    Returns:
        A boolean indicated whether the lists follow the same order
            while allowing for missing values
    """
    seen = set()
    deques = [deque(l) for l in lists]
    while sum(len(d) for d in deques) > 0:
        heads = [d[0] for d in deques if len(d) > 0]
        head = sorted(heads)[0]
        if head in seen:
            raise NotSortedException("Chromosomes aren't in the same order"
                                     "across VCF files.")
        else:
            seen.add(head)
            for d in deques:
                if len(d) > 0 and d[0] == head:
                    d.popleft()
    return True


def create_vcf_walktogether(vcf_files):
    """Create a VCF walk-together generator.

    Arguments:
        vcf_files (list): A list of `vcf.Reader` objects

    Returns:
        A `vcf.utils.walk_together` generator
    """
    def vcf_record_sort_key(r):
        return r.CHROM, r.POS, r.REF, r.ALT
    return pyvcf_utils.walk_together(*vcf_files, vcf_record_sort_key=vcf_record_sort_key)


def extract_names(vcf_files):
    """Extract source (method) names from the meta-information.

    Arguments:
        vcf_files (list): List of `vcf.Reader` objects

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
