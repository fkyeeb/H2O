'''
Entry point of the package, command line interface
'''

import argparse
import sys
from importlib.metadata import version, PackageNotFoundError

from h2o import (
    infer_orthology,
    map_duplications,
    extract_wgd_trees,
    map_gene_loss,
    bp2pie,
    extract_constraint_tree
)

def infer_ortho_main(args):
    """
    direct to infer_orthology
    """
    infer_orthology.main(args)

def map_duplications_main(args):
    """
    direct to map_duplications
    """
    map_duplications.main(args)

def extract_wgd_trees_main(args):
    """
    direct to extract_wgd_trees
    """
    extract_wgd_trees.main(args)

def gene_loss_main(args):
    """
    direct to gene_loss
    """
    map_gene_loss.main(args)

def bp2pie_main(args):
    """
    direct to bp2pie
    """
    bp2pie.main(args)

def extract_constraint_tree_main(args):
    """
    direct to extract_constraint_tree
    """
    extract_constraint_tree.main(args)

def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments
    """
    parser = argparse.ArgumentParser(prog="h2o", description="H2O CLI")
    subparsers = parser.add_subparsers(title='Subcommands for H2O', help="Subcommands", dest="subcommand")

    # Subcommand: h2o infer_ortho
    parser_infer_ortho = subparsers.add_parser("infer_ortho", help="Infer orthology and label gene duplications")
    parser_infer_ortho.add_argument("-t", "--homolog_tree_dir", type=str, help="Folder containing homolog trees", required=True)
    
    # Create mutually exclusive group for outgroup options
    outgroup_group = parser_infer_ortho.add_mutually_exclusive_group(required=True)
    outgroup_group.add_argument("-o", "--outgroup_list", type=str, help="List of outgroup taxa, separated by commas, no spaces")
    outgroup_group.add_argument("-of", "--outgroup_file", type=str, help="File containing the outgroup taxa, each line is a taxon")
    
    parser_infer_ortho.add_argument("-e", "--tree_file_ending", type=str, help="File ending of the homolog trees", required=True)
    parser_infer_ortho.add_argument("-m", "--min_dupl_tip_overlap", type=int, help="Minimum number of tip overlap between two child clades to be considered as a duplication node, default is 2")
    parser_infer_ortho.add_argument("-mp", "--min_dupl_percentage_overlap", type=float, help="Minimum percentage overlap between two child clades to be considered as a duplication node, default is 0.1")
    parser_infer_ortho.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_infer_ortho.add_argument("-p", "--just_pruning", action="store_true", help="Only produce pruned ortholog trees")
    parser_infer_ortho.add_argument("-np", "--no_pruning", action="store_true", help="Only produce unpruned ortholog trees")
    parser_infer_ortho.set_defaults(func=infer_ortho_main)

    # Subcommand: h2o map_dupl
    parser_map_duplications = subparsers.add_parser("map_dupl", help="Map gene duplications onto the consensus tree")
    parser_map_duplications.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_map_duplications.add_argument("-od", "--output_directory", type=str, help="Output directory, default is creating an `other_output/` directory in the currnet directory if not exist")
    parser_map_duplications.add_argument("-s", "--species_tree_file", type=str, help="Species tree file", required=True)
    parser_map_duplications.set_defaults(func=map_duplications_main)

    # Subcommand: h2o extract_wgd_trees
    parser_extract_wgd_trees = subparsers.add_parser("extract_wgd_trees", help="Extract and concatenate homolog trees that shows gene duplications at WGD events")
    parser_extract_wgd_trees.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_extract_wgd_trees.add_argument("-n", "--wgd_nodes", type=str, help="List of WGD node numbers, separated by commas, no spaces", required=True)
    parser_extract_wgd_trees.add_argument("-d", "--duplication_counts_dir", type=str, help="Duplication counts directory, default is other_output/ directory in the current directory")
    parser_extract_wgd_trees.add_argument("-od", "--output_directory", type=str, help="Output directory,default is other_output/ directory in the current directory, create if not exist")
    parser_extract_wgd_trees.set_defaults(func=extract_wgd_trees_main)

    # Subcommand: h2o gene_loss
    parser_gene_loss = subparsers.add_parser("gene_loss", help="Map gene copy losses after given WGD events")
    parser_gene_loss.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_gene_loss.add_argument("-n", "--wgd_nodes", type=str, help="List of WGD node numbers, separated by commas, no spaces", required=True)
    parser_gene_loss.add_argument("-d", "--duplication_counts_dir", type=str, help="Duplication counts directory, default is other_output directory in the current directory")
    parser_gene_loss.add_argument("-od", "--output_directory", type=str, help="Output directory, default is creating the other_output/ directory in the current directory, create if not exist")
    parser_gene_loss.set_defaults(func=gene_loss_main)

    # Subcommand: h2o bp2pie
    parser_bp2pie = subparsers.add_parser("bp2pie", help="Extract bp conflict result for R or gokstad plotting")
    parser_bp2pie.add_argument("-f", "--bp_output_file", type=str, help="bp output file, `-tv` has to be flagged when running bp, if multiple, separate by commas, no spaces", required=True)
    parser_bp2pie.add_argument("-s", "--consensus_tree_file", type=str, help="Consensus tree file, provide if branch length different from from the tree used to run bp")
    parser_bp2pie.add_argument("-od", "--output_directory", type=str, help="Output directory, default is the current directory")
    parser_bp2pie.add_argument("-p", "--pie_option", action="store_true", help="Flag to include unsupported counts in the gokstad pie tree")
    parser_bp2pie.add_argument("-n", "--run_name", type=str, help="Name of the run, to be added to output file name, default is empty string")
    parser_bp2pie.set_defaults(func=bp2pie_main)

    # Subcommand: h2o constraint
    parser_extract_constraint_tree = subparsers.add_parser("constraint", help="Extract constraint tree")
    parser_extract_constraint_tree.add_argument("-s", "--consensus_tree_file", type=str, help="Consensus tree file", required=True)
    parser_extract_constraint_tree.add_argument("-od", "--output_directory", type=str, help="Output directory, default is the current directory")

    # Create mutually exclusive group for nodes and tips options
    nodes_tips_group = parser_extract_constraint_tree.add_mutually_exclusive_group(required=True)
    nodes_tips_group.add_argument("-n", "--nodes", type=str, help="List of nodes to keep, node numbers separated by commas, no spaces")
    nodes_tips_group.add_argument("-t", "--tips_file", type=str, help="File containing the tips to keep, each line is a tip")
    
    parser_extract_constraint_tree.set_defaults(func=extract_constraint_tree_main)

    # version
    try:
        pkg_version = version("H2O")  # ← name from pyproject.toml [project.name]
    except PackageNotFoundError:
        pkg_version = "unknown"

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {pkg_version}")

    args = parser.parse_args()
    return args

def main():
    """
    Main function to parse arguments and call the appropriate function
    """
    print("   ___                          ")
    print("  (   )                         ")
    print("   | | .-.     .--.      .--.   ")
    print("   | |/   \\   ;  _  \\   /    \\  ")
    print("   |  .-. .  (___)` |  |  .-. ; ")
    print("   | |  | |       ' '  | |  | | ")
    print("   | |  | |      / /   | |  | | ")
    print("   | |  | |     / /    | |  | | ")
    print("   | |  | |    / /     | '  | | ")
    print("   | |  | |   / '____  '  `-' / ")
    print("  (___)(___) (_______)  `.__.'  ")
    print("\nThe ancient-WGD-aware homolog to ortholog trees command-line tool, for plant phylogenomics.\n")

    args = parse_arguments()
    if args.subcommand is None:
        print("No subcommand provided. Use -h for help.")
        sys.exit(1)
    
    args.func(args)

if __name__ == "__main__":
    main()