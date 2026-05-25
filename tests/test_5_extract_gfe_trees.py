from h2o import extract_gfe_trees,cli
import argparse
import os
import shutil
import pytest
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the extract_gfe_trees command. See if arguments are parsed correctly."""

    original_main = extract_gfe_trees.main

    captured_args = None
    
    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    monkeypatch.setattr(extract_gfe_trees, 'main', mock_main)

    test_args = ["h2o", "extract_gfe_trees", "-t", "tests/test_data/processed_trees", "-n", "3", "-od", "tests/test_data/extract_gfe_trees_output","-d", "tests/test_data/duplication_counts"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.processed_tree_dir == "tests/test_data/processed_trees"
    assert captured_args.gfe_nodes == "3"
    assert captured_args.output_directory == "tests/test_data/extract_gfe_trees_output"
    assert captured_args.duplication_counts_dir == "tests/test_data/duplication_counts"

    monkeypatch.setattr(extract_gfe_trees, 'main', original_main)


def test_extract_gfe_trees():
    """Test the extract_gfe_trees function"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        gfe_nodes="3",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/"
    )
    extract_gfe_trees.main(args)

    assert os.path.exists('tests/test_data/other_output/ASTRAL_in_unpruned_gfe_n3.tre')
    assert os.path.exists('tests/test_data/other_output/ASTRAL_in_pruned_gfe_n3.tre')

    os.remove('tests/test_data/other_output/cat_unpruned_gfe_trees.sh')
    os.remove('tests/test_data/other_output/cat_pruned_gfe_trees.sh')
    os.remove('tests/test_data/other_output/ASTRAL_in_unpruned_gfe_n3.tre')
    os.remove('tests/test_data/other_output/ASTRAL_in_pruned_gfe_n3.tre')

def test_gfe_node_error(capsys):
    """Test the extract_gfe_trees function with gfe node error"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        gfe_nodes="3.4",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/"
    )
    with pytest.raises(SystemExit) as exc_info:
        extract_gfe_trees.main(args)
    assert exc_info.value.code == 2

    captured = capsys.readouterr()
    assert "Error: gfe node numbers must be integers." in captured.out

def test_gfe_node_out_of_range(capsys):
    """Test the extract_gfe_trees function with gfe node out of range"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        gfe_nodes="8",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/"
    )
    with pytest.raises(SystemExit) as exc_info:
        extract_gfe_trees.main(args)
    assert exc_info.value.code == 2

    captured = capsys.readouterr()
    assert "Error: gfe node numbers provided is out of range for the consensus tree." in captured.out

def test_no_pruned_tree_folder(capsys):
    """Test the extract_gfe_trees function with no ortholog tree folder"""

    shutil.rmtree('tests/test_data/processed_trees/pruned')

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        gfe_nodes="3",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/"
    )
    extract_gfe_trees.main(args)

    captured = capsys.readouterr()
    assert "pruned processed tree folder not found." in captured.out

def test_no_ortholog_tree_folder(capsys):
    """Test the extract_gfe_trees function with no ortholog tree folder"""

    shutil.rmtree('tests/test_data/processed_trees/unpruned')

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        gfe_nodes="3",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/"
    )
    with pytest.raises(SystemExit) as exc_info:
        extract_gfe_trees.main(args)
    assert exc_info.value.code == 2

    captured = capsys.readouterr()
    assert "Error: No processed ortholog tree folder found." in captured.out

    shutil.rmtree('tests/test_data/processed_trees')
    shutil.rmtree('tests/test_data/other_output')