#!/usr/bin/env python2

"""
test_ensemble_caller.py
=======================
This script contains all unit tests for ensemble_calling.py.
"""

import os
import sys
import pytest
import ensemble_caller
import vcf as pyvcf

TEST_DIR = os.path.dirname(__file__)
VCF1 = os.path.join(TEST_DIR, 'strelka.vcf')
VCF2 = os.path.join(TEST_DIR, 'museq.vcf')
VCF_UNSORTED1 = os.path.join(TEST_DIR, 'strelka.unsorted1.vcf')
VCF_UNSORTED2 = os.path.join(TEST_DIR, 'strelka.unsorted2.vcf')
VCF_UNSORTED3 = os.path.join(TEST_DIR, 'museq.unsorted.vcf')
VCF_NOSOURCE = os.path.join(TEST_DIR, 'strelka.nosource.vcf')

ensemble_caller.setup_logging()


def open_vcf_files(*vcf_files):
    vcf_files_opened = tuple([pyvcf.Reader(open(vcf)) for vcf in vcf_files])
    if len(vcf_files_opened) == 1:
        return vcf_files_opened[0]
    else:
        return vcf_files_opened


@pytest.fixture
def setup_vcf_files():
    return open_vcf_files(VCF1, VCF2)


@pytest.fixture
def setup_unsorted_vcf_files():
    return open_vcf_files(VCF_UNSORTED1, VCF_UNSORTED2, VCF_UNSORTED3)


@pytest.fixture
def setup_nosource_vcf_file():
    return open_vcf_files(VCF_NOSOURCE)


def test_parse_args():
    """Test parse_args function."""
    # Test 1: no VCF files
    test_1 = []
    with pytest.raises(SystemExit):
        ensemble_caller.parse_args(test_1)
    # Test 2: one VCF file
    # Expectation: all default values
    test_2 = [VCF1]
    args_2 = ensemble_caller.parse_args(test_2)
    assert isinstance(args_2.vcf_files, list)
    assert len(args_2.vcf_files) == 1
    assert all(isinstance(f, pyvcf.Reader) for f in args_2.vcf_files)
    assert args_2.skip_sort_check is False
    assert isinstance(args_2.output_vcf, pyvcf.Writer)
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
    # Test 5: attempt with non-existent file
    test_5 = [VCF1, VCF2, "oups.vcf"]
    with pytest.raises(IOError):
        ensemble_caller.parse_args(test_5)
    # Test 6: Manually specify method names
    test_6 = ["--names", "strelka,mutationseq", VCF1, VCF2]
    args_6 = ensemble_caller.parse_args(test_6)
    assert args_6.names == ["strelka", "mutationseq"]
    # Test 7: Manually specify method names, but wrong number
    test_7 = ["--names", "strelka,mutationseq,mutect", VCF1, VCF2]
    with pytest.raises(AssertionError):
        ensemble_caller.parse_args(test_7)
    # Test 8: specify output_vcf in existing directory
    output_vcf_path1 = os.path.join(TEST_DIR, "output.vcf")
    test_8 = ["--output_vcf", output_vcf_path1, VCF1, VCF2]
    args_8 = ensemble_caller.parse_args(test_8)
    assert isinstance(args_8.output_vcf, pyvcf.Writer)
    os.remove(output_vcf_path1)
    # Test 9: specify output_vcf in non-existing directory
    # Expectation: non-existing directory is created
    output_vcf_path2 = os.path.join(TEST_DIR, "output_dir", "output.vcf")
    test_9 = ["--output_vcf", output_vcf_path2, VCF1, VCF2]
    ensemble_caller.parse_args(test_9)
    assert os.path.exists(os.path.join(TEST_DIR, "output_dir"))
    os.remove(output_vcf_path2)
    os.rmdir(os.path.dirname(output_vcf_path2))


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
    vcf_file_unsorted1, vcf_file_unsorted2, vcf_file_unsorted3 = setup_unsorted_vcf_files
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.are_sorted([vcf_file1, vcf_file_unsorted1])
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.are_sorted([vcf_file1, vcf_file_unsorted2])
    assert ensemble_caller.are_sorted([vcf_file1, vcf_file_unsorted3]) is False


def test_parse_order(setup_vcf_files, setup_unsorted_vcf_files):
    """Test parse_order function"""
    # Test 1: Obtain chromosome order from sorted file
    vcf_file, _ = setup_vcf_files
    expect = ["1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
              "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X"]
    chrom_order = ensemble_caller.parse_order(vcf_file)
    assert chrom_order == expect
    # Test 2: Obtain chromosome order from unsorted file (2 chrom blocks)
    # Expectation: NotSortedException
    vcf_file_unsorted1, vcf_file_unsorted2, vcf_file_unsorted3 = setup_unsorted_vcf_files
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.parse_order(vcf_file_unsorted1)
    # Test 3: Obtain chromsome order from unsorted file (permuted pos)
    # Expectation: NotSortedException
    with pytest.raises(ensemble_caller.NotSortedException):
        ensemble_caller.parse_order(vcf_file_unsorted2)
    # Test 4: Obtain chromsome order from weirdly sorted file, but still
    # in contiguous blocks
    # Expectation: No error
    assert len(ensemble_caller.parse_order(vcf_file_unsorted3)) > 0


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
    assert ensemble_caller.compare_orders(test_4) is False
    # Test 5: Compare unordered lists
    test_5 = [[3, 2, 1], [1, 2]]
    assert ensemble_caller.compare_orders(test_5) is False
    # Test 6: Compare lexicographically sorted lists
    test_6 = [["chr1", "chr10", "chr2"], ["chr1", "chr10", "chr3"]]
    assert ensemble_caller.compare_orders(test_6) is True
    # Test 7: Compare numerically sorted lists (of strings)
    test_7 = [["chr1", "chr2", "chr10"], ["chr1", "chr3", "chr10"]]
    assert ensemble_caller.compare_orders(test_7) is True
    # Test 8: Compare inconsistently sorted lists
    test_8 = [["chr1", "chr2", "chr10"], ["chr1", "chr10", "chr3"]]
    assert ensemble_caller.compare_orders(test_8) is False


def test_create_vcf_walktogether(setup_vcf_files):
    """Test create_vcf_walktogether function."""
    vcf_files = list(setup_vcf_files)
    vcf_wt = ensemble_caller.create_vcf_walktogether(vcf_files)
    length = len([r for r in vcf_wt])
    assert length == 410


def test_extract_names(setup_vcf_files, setup_nosource_vcf_file):
    """Test extract_names function."""
    # Test 1: VCF files with source
    vcf_files = list(setup_vcf_files)
    expect_1 = ["strelka", "mutationSeq_4.3.6"]
    names_1 = ensemble_caller.extract_names(vcf_files)
    assert names_1 == expect_1
    # Test 2: VCF file without source
    vcf_file_nosource = setup_nosource_vcf_file
    expect_2 = ["strelka", "method_2"]
    names_2 = ensemble_caller.extract_names([vcf_files[0], vcf_file_nosource])
    assert names_2 == expect_2
