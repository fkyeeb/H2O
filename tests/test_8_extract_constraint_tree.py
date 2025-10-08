from h2o import extract_constraint_tree,cli
import argparse
import shutil
import os
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the extract_constraint_tree command. See if arguments are parsed correctly."""

    original_main = extract_constraint_tree.main
    
    captured_args = None

    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    monkeypatch.setattr(extract_constraint_tree, 'main', mock_main)

    test_args = ["h2o", "constraint", "-s", "tests/test_data/species.tre", "-n", "141,132,115,103,77,101,87,2", "-od", "tests/test_data/other_output"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.summary_tree_file == "tests/test_data/species.tre"
    assert captured_args.nodes == "141,132,115,103,77,101,87,2"
    assert captured_args.output_directory == "tests/test_data/other_output"

    test_args = ["h2o", "constraint", "-s", "tests/test_data/species.tre", "-t", "tests/test_data/tips.txt", "-od", "tests/test_data/other_output"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.tips_file == "tests/test_data/tips.txt"

    monkeypatch.setattr(extract_constraint_tree, 'main', original_main)

def test_extract_constraint_tree_node():
    """Test the extract_constraint_tree function"""
    args = argparse.Namespace(
        summary_tree_file="bp_summary_tree_numbered.tre",
        nodes="4,5,o,2",
        output_directory=None,
        tips_file=None
    )
    os.chdir("tests/test_data/other_output/")

    extract_constraint_tree.main(args)

    os.chdir("../../..")

    with open("tests/test_data/other_output/constraint_tree.tre","r") as f:
        tree = f.readline()
        assert "o" in tree
        assert ("a" in tree) or ("b" in tree)
        assert ("c" in tree) or ("d" in tree)
        assert ("e" in tree) or ("f" in tree)

def test_extract_constraint_tree_tip():
    """Test the extract_constraint_tree function"""
    args = argparse.Namespace(
        summary_tree_file="tests/test_data/other_output/bp_summary_tree_numbered.tre",
        nodes=None,
        tips_file="tests/test_data/tips2keep.txt",
        output_directory="tests/test_data/other_output"
    )
    extract_constraint_tree.main(args)

    with open("tests/test_data/other_output/constraint_tree.tre","r") as f:
        tree = f.readline()
        assert "e" in tree
        assert "b" in tree
        assert "d" in tree
        assert "o" not in tree
    
    shutil.rmtree("tests/test_data/other_output/")
