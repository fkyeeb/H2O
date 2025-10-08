from h2o import bp2pie,cli
import argparse
import os
import shutil
import pytest
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the bp2pie command. See if arguments are parsed correctly."""

    original_main = bp2pie.main

    captured_args = None

    def mock_main(args):
        nonlocal captured_args
        captured_args = args
    
    monkeypatch.setattr(bp2pie, 'main', mock_main)

    test_args = ["h2o", "bp2pie", "-f", "tests/test_data/bp_output.txt", "-s", "tests/test_data/species_2.tre", "-od", "tests/test_data/other_output", "-p"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.bp_output_file == "tests/test_data/bp_output.txt"
    assert captured_args.summary_tree_file == "tests/test_data/species_2.tre"
    assert captured_args.output_directory == "tests/test_data/other_output"
    assert captured_args.pie_option == True

    test_args = ["h2o", "bp2pie", "-f", "tests/test_data/bp_output.txt"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.bp_output_file == "tests/test_data/bp_output.txt"
    assert captured_args.output_directory == None
    assert captured_args.pie_option == False
    assert captured_args.summary_tree_file == None

    monkeypatch.setattr(bp2pie, 'main', original_main)

def test_bp2pie():
    """Test the bp2pie function"""
    args = argparse.Namespace(
        bp_output_file="tests/test_data/bp_output.txt",
        output_directory="tests/test_data/other_output",
        summary_tree_file="tests/test_data/species_2.tre",
        pie_option=True
    )
    bp2pie.main(args)

    with open("tests/test_data/other_output/bp_data.tsv","r") as f:
        lines = f.readlines()
        assert lines[0] == "node_number\tconflict\tconcord\tunsupported\n"
        assert lines[5] == "5\t1\t1\t2\n"
    with open("tests/test_data/other_output/bp_output.tre","r") as f:
        lines = f.readlines()
        assert lines[2] == "(o,((e,f)2,((a,b),(c,d)2)));\n"
    with open("tests/test_data/other_output/bp_summary_tree_numbered.tre","r") as f:
        assert f.readline() == "(o:1.0,(((c:1.0,d:1.0)5:0.0,(a:1.0,b:1.0)4:1.0)3:1.0,(e:1.0,f:1.0)2:1.0)1:1.0):0.0;\n"

def test_bp2pie_2():
    args = argparse.Namespace(
        bp_output_file="tests/test_data/bp_output.txt",
        output_directory="tests/test_data/other_output",
        summary_tree_file=None,
        pie_option=False
    )
    bp2pie.main(args)
    with open("tests/test_data/other_output/bp_summary_tree_numbered.tre","r") as f:
        assert f.readline() == "(o:0.0,((e:0.0,f:0.0)2:0.0,((a:0.0,b:0.0)4:0.0,(c:0.0,d:0.0)5:0.0)3:0.0)1:0.0)0:0.0;\n"
    with open("tests/test_data/other_output/gokstad_pie.tre","r") as f:
        assert f.readline() == "(o:0.0,((e:0.0,f:0.0)[&pie=0,2]:0.0,((a:0.0,b:0.0)[&pie=2,2]:0.0,(c:0.0,d:0.0)[&pie=1,1]:0.0)[&pie=4,0]:0.0)[&pie=1,3]:0.0)0:0.0;\n"
    
    shutil.rmtree("tests/test_data/other_output/")