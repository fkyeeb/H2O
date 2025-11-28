import argparse
import os
import sys
import shutil
from h2o import map_duplications,cli

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the map_dupl command. See if arguments are parsed correctly."""

    original_main = map_duplications.main

    captured_args = None
    
    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    monkeypatch.setattr(map_duplications, 'main', mock_main)

    test_args = ["h2o", "map_dupl", "-t", "tests/test_data/processed_trees", "-s", "tests/test_data/species.tre", "-od", "tests/test_data/map_dupl_output"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.processed_tree_dir == "tests/test_data/processed_trees"
    assert captured_args.species_tree_file == "tests/test_data/species.tre"
    assert captured_args.output_directory == "tests/test_data/map_dupl_output"

    monkeypatch.setattr(map_duplications, 'main', original_main)

def test_map_dupl():
    """Test the map_dupl function with advance"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        output_directory="tests/test_data/map_dupl_output",
        species_tree_file="tests/test_data/species.tre"
    )
    map_duplications.main(args)

    with open('tests/test_data/map_dupl_output/duplication_counts.tsv', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "dup_bl" in line:
                assert line == 'dup_bl\t0\t0\t0\t3\t0\t0\t0\t\n'
                break

    with open('tests/test_data/map_dupl_output/consensus_tree_numbered.tre', 'r') as f:
        assert f.read() == '(o:0.0,((e:0.0,f:0.0)2:0.0,((a:0.0,b:0.0)4:0.0,(c:0.0,d:0.0)5:0.0)3:0.0)1:0.0)0:0.0;\n(o:0.0,((e:0.0,f:0.0)0:0.0,((a:0.0,b:0.0)0:0.0,(c:0.0,d:0.0)1:0.0)4:0.0)2:0.0)0:0.0;\n'

    shutil.rmtree('tests/test_data/map_dupl_output')

def test_map_dupl_no_output_directory():
    """Test the map_dupl function with no output directory"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        species_tree_file="tests/test_data/species.tre",
        output_directory="tests/test_data/other_output/"
    )
    map_duplications.main(args)

    assert os.path.exists('tests/test_data/other_output/duplication_counts.tsv')

    # shutil.rmtree('tests/test_data/processed_trees')
    # shutil.rmtree('tests/test_data/other_output')
