import sys
import os
import pytest
import shutil

from h2o import cli

def test_parse_arguments_infer_ortho(monkeypatch):
    """Test the parse_arguments function for the infer_ortho command. See if arguments are parsed correctly."""

    test_args = ["h2o", "infer_ortho", "-d", "tests/test_data/homolog_trees", "-of", "tests/test_data/outgroups.txt", "-t", ".tre", "-m", "3", "-od", "tests/test_data/ortholog_trees"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    args = cli.parse_arguments()
    assert args.homolog_tree_dir == "tests/test_data/homolog_trees"
    assert args.outgroup_file == "tests/test_data/outgroups.txt"
    assert args.tree_file_ending == ".tre"
    assert args.min_ingroup_taxa == 3
    assert args.output_directory == "tests/test_data/ortholog_trees"

def test_no_subcommand(monkeypatch, capsys):
    """Test errors for no subcommand."""

    test_args = ["h2o"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    # Test that the main function prints the expected message and exits
    with pytest.raises(SystemExit) as exc_info:
        cli.main()
    
    # Assert the exit code
    assert exc_info.value.code == 1
    
    # Capture the printed output
    captured = capsys.readouterr()
    assert "No subcommand provided. Use -h for help." in captured.out


def test_main(monkeypatch):
    """Test the main function for the infer_ortho command. See if infer_orthology.py is being ran correctly."""

    test_args = ["h2o", "infer_ortho", "-d", "tests/test_data/homolog_trees", "-o", "o", "-t", ".tre", "-od", "tests/test_data/ortholog_trees"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/ortholog_trees/dup_rooted.tre')
    assert os.path.exists('tests/test_data/ortholog_trees/unpruned/dup_rooted_pruned.tre')
    assert os.path.exists('tests/test_data/ortholog_trees/pruned/dup_rooted_pruned.tre')

    # remove all the output files created by unit test
    shutil.rmtree('tests/test_data/ortholog_trees')

def test_only_pruning(monkeypatch):
    """Test pruning options"""

    test_args = ["h2o", "infer_ortho", "-d", "tests/test_data/homolog_trees", "-o", "o", "-t", ".tre", "-od", "tests/test_data/ortholog_trees", "-p"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/ortholog_trees/pruned/dup_rooted_pruned.tre')
    assert not os.path.exists('tests/test_data/ortholog_trees/unpruned/')

    shutil.rmtree('tests/test_data/ortholog_trees')

def test_only_no_pruning(monkeypatch):
    """Test pruning options"""

    test_args = ["h2o", "infer_ortho", "-d", "tests/test_data/homolog_trees", "-o", "o", "-t", ".tre", "-od", "tests/test_data/ortholog_trees", "-np"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/ortholog_trees/unpruned/dup_rooted_pruned.tre')
    assert not os.path.exists('tests/test_data/ortholog_trees/pruned/')
    
    shutil.rmtree('tests/test_data/ortholog_trees')