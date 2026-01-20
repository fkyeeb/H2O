from h2o import map_gene_loss,cli,infer_orthology,map_duplications
import argparse
import shutil
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the map_gene_loss command. See if arguments are parsed correctly."""

    original_main = map_gene_loss.main

    captured_args = None

    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    monkeypatch.setattr(map_gene_loss, 'main', mock_main)

    test_args = ["h2o", "gene_loss", "-t", "tests/test_data/processed_trees", "-n", "3", "-od", "tests/test_data/gene_loss_output","-d","tests/test_data/duplication_counts"]
    monkeypatch.setattr(sys, "argv", test_args)

    args = cli.parse_arguments()
    args.func(args)

    assert captured_args.processed_tree_dir == "tests/test_data/processed_trees"
    assert captured_args.wgd_nodes == "3"
    assert captured_args.output_directory == "tests/test_data/gene_loss_output"
    assert captured_args.duplication_counts_dir == "tests/test_data/duplication_counts"

    monkeypatch.setattr(map_gene_loss, 'main', original_main)

def test_map_gene_loss():
    """Test the map_gene_loss function"""

    # reruning previous analyses
    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees",
        outgroup_list="o",
        tree_file_ending=".tre",
        min_dupl_tip_overlap=None,
        min_dupl_percentage_overlap=None,
        output_directory="tests/test_data/processed_trees",
        no_pruning=False,
        just_pruning=False,
        id2sp_file=None
    )
    infer_orthology.main(args)

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        species_tree_file="tests/test_data/species.tre",
        output_directory="tests/test_data/other_output/",
        id2sp_file=None
    )
    map_duplications.main(args)

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        output_directory="tests/test_data/other_output/",
        duplication_counts_dir="tests/test_data/other_output/",
        wgd_nodes="3,1"
    )
    map_gene_loss.main(args)

    with open('tests/test_data/other_output/gene_loss_counts.tsv', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "dup_loss" in line:
                assert line == 'dup_loss\t1\t1\t2\t\n'
                break
    
    shutil.rmtree('tests/test_data/other_output')
    shutil.rmtree('tests/test_data/processed_trees')