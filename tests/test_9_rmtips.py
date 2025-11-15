from h2o import rmtips,cli
import argparse
import shutil
import os
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the rmtips command. See if arguments are parsed correctly."""

    original_main = rmtips.main
    
    captured_args = None

    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    monkeypatch.setattr(rmtips, 'main', mock_main)

    test_args = ["h2o", "rmtips", "-t", "tests/test_data/homolog_trees", "-rm", "file1.txt", "-od", "tests/test_data/pruned_homologs", "-e", ".tre","-m","10"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.tree_dir == "tests/test_data/homolog_trees"
    assert captured_args.tips2remove == "file1.txt"
    assert captured_args.output_directory == "tests/test_data/pruned_homologs"
    assert captured_args.tree_file_ending == ".tre"
    assert captured_args.minimum_taxa == "10"

    test_args = ["h2o", "rmtips", "-t", "tests/test_data/homolog_trees", "-sv", "file1.txt", "-od", "tests/test_data/pruned_homologs", "-e", ".tre"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.tips2save == "file1.txt"
    assert captured_args.tips2remove == None

    monkeypatch.setattr(rmtips, 'main', original_main)

def test_rmtips():
    """Test the rmtips function"""
    args = argparse.Namespace(
        tree_dir="tests/test_data/homolog_trees",
        tips2remove="tests/test_data/tips2rm.txt",
        output_directory="tests/test_data/pruned_homologs/",
        tips2save=None,
        tree_file_ending=".tre",
        minimum_taxa=None
    )

    rmtips.main(args)
    with open("tests/test_data/pruned_homologs/dup_loss.tre", "r") as file:
        assert file.read() == "(o:0.0,(((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0,((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0):0.0):0.0;\n"
    shutil.rmtree("tests/test_data/pruned_homologs/")

def test_svtips():
    """Test the rmtips function"""
    args = argparse.Namespace(
        tree_dir="tests/test_data/homolog_trees",
        tips2remove=None,
        output_directory="tests/test_data/pruned_homologs/",
        tips2save="tests/test_data/tips2sv.txt",
        tree_file_ending=".tre",
        minimum_taxa="7"
    )

    rmtips.main(args)
    with open("tests/test_data/pruned_homologs/dup_loss.tre", "r") as file:
        assert file.read() == "(o:0.0,(((a:0.0,b:0.0):0.0,c:0.0):0.0,((a:0.0,b:0.0):0.0,c:0.0):0.0):0.0):0.0;\n"
    assert not os.path.exists("tests/test_data/pruned_homologs/dup_difficult.tre")
    shutil.rmtree("tests/test_data/pruned_homologs/")