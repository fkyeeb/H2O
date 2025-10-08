'''
Entry point of the package, command line interface
'''

import argparse
import sys

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
    parser_infer_ortho = subparsers.add_parser("infer_ortho", help="Orthology and Gene duplication Inference")
    parser_infer_ortho.add_argument("-t", "--homolog_tree_dir", type=str, help="Folder containing homolog trees", required=True)
    
    # Create mutually exclusive group for outgroup options
    outgroup_group = parser_infer_ortho.add_mutually_exclusive_group(required=True)
    outgroup_group.add_argument("-o", "--outgroup_list", type=str, help="List of outgroup taxa, separated by commas, no spaces")
    outgroup_group.add_argument("-of", "--outgroup_file", type=str, help="File containing the outgroup taxa, each line is a taxon")
    
    parser_infer_ortho.add_argument("-e", "--tree_file_ending", type=str, help="File ending of the homolog trees", required=True)
    parser_infer_ortho.add_argument("-m", "--min_ingroup_taxa", type=int, help="Minimum number of ingroup taxa, default is 3")
    parser_infer_ortho.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_infer_ortho.add_argument("-p", "--just_pruning", action="store_true", help="Only produce pruned ortholog trees")
    parser_infer_ortho.add_argument("-np", "--no_pruning", action="store_true", help="Only produce unpruned ortholog trees")
    parser_infer_ortho.set_defaults(func=infer_ortho_main)

    # Subcommand: h2o map_dupl
    parser_map_duplications = subparsers.add_parser("map_dupl", help="Map gene duplications")
    parser_map_duplications.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_map_duplications.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_map_duplications.add_argument("-s", "--species_tree_file", type=str, help="Species tree file", required=True)
    parser_map_duplications.set_defaults(func=map_duplications_main)

    # Subcommand: h2o extract_wgd_trees
    parser_extract_wgd_trees = subparsers.add_parser("extract_wgd_trees", help="Extract and concatenate homolog trees that shows gene duplications at WGD events")
    parser_extract_wgd_trees.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_extract_wgd_trees.add_argument("-n", "--wgd_nodes", type=str, help="List of WGD node numbers, separated by commas, no spaces", required=True)
    parser_extract_wgd_trees.add_argument("-d", "--duplication_counts_dir", type=str, help="Duplication counts directory")
    parser_extract_wgd_trees.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_extract_wgd_trees.set_defaults(func=extract_wgd_trees_main)

    # Subcommand: h2o gene_loss
    parser_gene_loss = subparsers.add_parser("gene_loss", help="map gene losses after given WGD events")
    parser_gene_loss.add_argument("-t", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_gene_loss.add_argument("-n", "--wgd_nodes", type=str, help="List of WGD node numbers, separated by commas, no spaces", required=True)
    parser_gene_loss.add_argument("-d", "--duplication_counts_dir", type=str, help="Duplication counts directory")
    parser_gene_loss.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_gene_loss.set_defaults(func=gene_loss_main)

    # Subcommand: h2o bp2pie
    parser_bp2pie = subparsers.add_parser("bp2pie", help="extract bp conflict result for ploting")
    parser_bp2pie.add_argument("-f", "--bp_output_file", type=str, help="bp output file, `-tv` has to be flagged when running bp", required=True)
    parser_bp2pie.add_argument("-s", "--summary_tree_file", type=str, help="Summary tree file, provide if branch length different from bp tree")
    parser_bp2pie.add_argument("-od", "--output_directory", type=str, help="Output directory, default is the current directory")
    parser_bp2pie.add_argument("-p", "--pie_option", action="store_true", help="Flag to include unsupported counts in the gokstad pie tree")
    parser_bp2pie.set_defaults(func=bp2pie_main)

    # Subcommand: h2o constraint
    parser_extract_constraint_tree = subparsers.add_parser("constraint", help="extract constraint tree")
    parser_extract_constraint_tree.add_argument("-s", "--summary_tree_file", type=str, help="Summary tree file", required=True)
    parser_extract_constraint_tree.add_argument("-od", "--output_directory", type=str, help="Output directory, default is the current directory")

    # Create mutually exclusive group for nodes and tips options
    nodes_tips_group = parser_extract_constraint_tree.add_mutually_exclusive_group(required=True)
    nodes_tips_group.add_argument("-n", "--nodes", type=str, help="List of nodes to keep, node numbers separated by commas, no spaces")
    nodes_tips_group.add_argument("-t", "--tips_file", type=str, help="File containing the tips to keep, each line is a tip")
    
    parser_extract_constraint_tree.set_defaults(func=extract_constraint_tree_main)

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
    print("\n the welcome message goes here")

    args = parse_arguments()
    if args.subcommand is None:
        print("No subcommand provided. Use -h for help.")
        sys.exit(1)
    
    args.func(args)

if __name__ == "__main__":
    main()