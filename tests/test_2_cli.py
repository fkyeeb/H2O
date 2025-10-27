import sys
import os
import pytest
import shutil

from h2o import cli

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


def test_infer_ortho_main(monkeypatch):
    """Test the main function for the infer_ortho command. See if infer_orthology.py is being ran correctly."""

    test_args = ["h2o", "infer_ortho", "-t", "tests/test_data/homolog_trees", "-o", "o", "-e", ".tre", "-od", "tests/test_data/processed_trees"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/processed_trees/dup_rooted.tre')
    assert os.path.exists('tests/test_data/processed_trees/unpruned/dup_rooted_processed.tre')
    assert os.path.exists('tests/test_data/processed_trees/pruned/dup_rooted_processed.tre')

    # remove all the output files created by unit test
    shutil.rmtree('tests/test_data/processed_trees')

def test_infer_ortho_only_pruning(monkeypatch):
    """Test pruning options"""

    test_args = ["h2o", "infer_ortho", "-t", "tests/test_data/homolog_trees", "-o", "o", "-e", ".tre", "-od", "tests/test_data/processed_trees", "-p"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/processed_trees/pruned/dup_rooted_processed.tre')
    assert not os.path.exists('tests/test_data/processed_trees/unpruned/')

    shutil.rmtree('tests/test_data/processed_trees')

def test_infer_ortho_only_no_pruning(monkeypatch):
    """Test pruning options"""

    test_args = ["h2o", "infer_ortho", "-t", "tests/test_data/homolog_trees", "-o", "o", "-e", ".tre", "-np","-od", "tests/test_data/processed_trees"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    cli.main()

    assert os.path.exists('tests/test_data/processed_trees/unpruned/dup_rooted_processed.tre')
    assert not os.path.exists('tests/test_data/processed_trees/pruned/')
    
    shutil.rmtree('tests/test_data/processed_trees')