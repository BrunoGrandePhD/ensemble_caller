#!/usr/bin/env python2

"""
test_ensemble_caller.py
=======================
This script contains all unit tests for ensemble_calling.py.
"""

import os
import pytest
import ensemble_caller
import vcf as pyvcf

VCF1 = os.path.join(os.path.dirname(__file__), 'strelka.vcf')
VCF2 = os.path.join(os.path.dirname(__file__), 'museq.vcf')
VCF_SORTED1 = os.path.join(os.path.dirname(__file__), 'strelka.sort.vcf')
VCF_SORTED2 = os.path.join(os.path.dirname(__file__), 'museq.sort.vcf')


def test_parse_args():
    """Test parse_args function."""
    # Test 1: no VCF files
    test_1 = []
    with pytest.raises(SystemExit):
        ensemble_caller.parse_args(test_1)
    # Test 2: one VCF file
    test_2 = [VCF1]
    args_2 = ensemble_caller.parse_args(test_2)
    assert isinstance(args_2.vcf_files, list)
    assert len(args_2.vcf_files) == 1
    assert all(isinstance(f, pyvcf.Reader) for f in args_2.vcf_files)
    assert args_2.skip_sort_check is False
    # Test 3: two VCF files
    test_3 = [VCF1, VCF2]
    args_3 = ensemble_caller.parse_args(test_3)
    assert isinstance(args_3.vcf_files, list)
    assert len(args_3.vcf_files) == 2
    assert all(isinstance(f, pyvcf.Reader) for f in args_3.vcf_files)
    # Test 4: Skip sort check
    test_4 = ["--skip_sort_check", VCF1]
    args_4 = ensemble_caller.parse_args(test_4)
    assert args_4.skip_sort_check is True


def test_reset_vcf_files():
    """Test reset_vcf_files function"""
    vcf_file = pyvcf.Reader(open(VCF1))
    record = next(vcf_file)
    assert record.CHROM == "1" and record.POS == 2317235
    for i in range(5):
        record = next(vcf_file)
    assert not (record.CHROM == "1" and record.POS == 2317235)
    ensemble_caller.reset_vcf_files([vcf_file])
    record = next(vcf_file)
    assert record.CHROM == "1" and record.POS == 2317235


def test_check_sort():
    """Test check_sort function"""
    # Compare unsorted files
    vcf_file1 = pyvcf.Reader(open(VCF1))
    vcf_file2 = pyvcf.Reader(open(VCF2))
    is_sorted = ensemble_caller.check_sort([vcf_file1, vcf_file2])
    assert is_sorted is False
    # Compare sorted files
    vcf_file_sorted1 = pyvcf.Reader(open(VCF_SORTED1))
    vcf_file_sorted2 = pyvcf.Reader(open(VCF_SORTED2))
    is_sorted = ensemble_caller.check_sort([vcf_file_sorted1, vcf_file_sorted2])
    assert is_sorted is True
