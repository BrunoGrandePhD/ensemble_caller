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
VCF_UNSORTED1 = os.path.join(os.path.dirname(__file__), 'strelka.unsorted1.vcf')
VCF_UNSORTED2 = os.path.join(os.path.dirname(__file__), 'strelka.unsorted2.vcf')


@pytest.fixture
def setup_vcf_files():
    vcf_file1 = pyvcf.Reader(open(VCF1))
    vcf_file2 = pyvcf.Reader(open(VCF2))
    return vcf_file1, vcf_file2


@pytest.fixture
def setup_unsorted_vcf_files():
    vcf_file_unsorted1 = pyvcf.Reader(open(VCF_UNSORTED1))
    vcf_file_unsorted2 = pyvcf.Reader(open(VCF_UNSORTED2))
    return vcf_file_unsorted1, vcf_file_unsorted2


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


def test_reset_vcf_files(setup_vcf_files):
    """Test reset_vcf_files function"""
    vcf_file, _ = setup_vcf_files
    record = next(vcf_file)
    assert record.CHROM == "1" and record.POS == 2317235
    for i in range(5):
        record = next(vcf_file)
    assert not (record.CHROM == "1" and record.POS == 2317235)
    ensemble_caller.reset_vcf_files([vcf_file])
    record = next(vcf_file)
    assert record.CHROM == "1" and record.POS == 2317235


def test_are_sorted(setup_vcf_files, setup_unsorted_vcf_files):
    """Test are_sorted function"""
    # Compare sorted files
    vcf_file1, vcf_file2 = setup_vcf_files
    are_sorted = ensemble_caller.are_sorted([vcf_file1, vcf_file2])
    assert are_sorted is True
    # Compare sorted and unsorted files
    vcf_file_unsorted1, vcf_file_unsorted2 = setup_unsorted_vcf_files
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.are_sorted([vcf_file1, vcf_file_unsorted1])
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.are_sorted([vcf_file1, vcf_file_unsorted2])


def test_parse_order(setup_vcf_files, setup_unsorted_vcf_files):
    """Test parse_order function"""
    # Test 1: Obtain chromosome order from sorted file
    vcf_file, _ = setup_vcf_files
    expect_1 = ["1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X"]
    chrom_order = ensemble_caller.parse_order(vcf_file)
    assert chrom_order == expect_1
    # Test 2: Obtain chromosome order from unsorted file (2 chrom blocks)
    # Expectation: NotSortedException
    vcf_file_unsorted1, vcf_file_unsorted2 = setup_unsorted_vcf_files
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.parse_order(vcf_file_unsorted1)
    # Test 3: Obtain chromsome order from unsorted file (permuted pos)
    # Expectation: NotSortedException
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.parse_order(vcf_file_unsorted2)


def test_compare_orders():
    """Test compare_orders function"""
    # Test 1: Compare identical lists
    test_1 = [[1, 2, 3], [1, 2, 3]]
    assert ensemble_caller.compare_orders(test_1) is True
    # Test 2: Compare ordered lists with a missing value (middle)
    test_2 = [[1, 2, 3], [1, 3]]
    assert ensemble_caller.compare_orders(test_2) is True
    # Test 3: Compare ordered lists with a missing value (end)
    test_3 = [[1, 2, 3], [1, 2]]
    assert ensemble_caller.compare_orders(test_3) is True
    # Test 4: Compare unordered lists
    test_4 = [[1, 2, 3], [1, 3, 2]]
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.compare_orders(test_4)
