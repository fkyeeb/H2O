from h2o import phyparts2pie,cli
import argparse
import shutil
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the phyparts2pie command. See if arguments are parsed correctly."""

    original_main = phyparts2pie.main

    captured_args = None

    def mock_main(args):
        nonlocal captured_args
        captured_args = args
    
    monkeypatch.setattr(phyparts2pie, 'main', mock_main)

    test_args = ["h2o", "phyparts2pie", "-k", "tests/test_data/phyparts.node.key", "-f", "tests/test_data/phyparts.hist", "-s", "tests/test_data/phyparts.tre", "-n", "7", "-od", "tests/test_data/other_output"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.phyparts_node_key_file == "tests/test_data/phyparts.node.key"
    assert captured_args.phyparts_hist_file == "tests/test_data/phyparts.hist"
    assert captured_args.consensus_tree_file == "tests/test_data/phyparts.tre"
    assert captured_args.total_tree_number == "7"
    assert captured_args.output_directory == "tests/test_data/other_output"

    test_args = ["h2o", "phyparts2pie", "-k", "tests/test_data/phyparts.node.key", "-f", "tests/test_data/phyparts.hist", "-s", "tests/test_data/phyparts.tre", "-n", "7"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)
    assert captured_args.output_directory == None

    monkeypatch.setattr(phyparts2pie, 'main', original_main)

def test_phyparts2pie():
    """Test the phyparts2pie function"""
    args = argparse.Namespace(
        phyparts_node_key_file="tests/test_data/phyparts.node.key",
        phyparts_hist_file="tests/test_data/phyparts.hist",
        consensus_tree_file="tests/test_data/phyparts.tre",
        total_tree_number="7",
        output_directory="tests/test_data/other_output"
    )
    phyparts2pie.main(args)

    with open("tests/test_data/other_output/phyparts_summary.tsv","r") as f:
        lines = f.readlines()
        assert lines[2] == "1	6.0	1.0	0.0	0\n"
    with open("tests/test_data/other_output/gokstad_pie.tre","r") as f:
        assert f.readline() == "((((A:0.0,B:0.0)[&pie=4.0,2.0,0.0,1]:0.0,C:0.0)[&pie=6.0,1.0,0.0,0]:0.0,D:0.0)[&pie=2.0,2.0,0.0,3]:0.0,E:0.0):0.0;\n"