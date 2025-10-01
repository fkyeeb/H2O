'''
Entry point of the package, command line interface
'''

import argparse
import sys

from h2o import infer_orthology, map_duplications, extract_wgd_trees

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

def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments
    """
    parser = argparse.ArgumentParser(prog="h2o", description="H2O CLI")
    subparsers = parser.add_subparsers(title='Subcommands for H2O', help="Subcommands", dest="subcommand")

    # Subcommand: h2o infer_ortho
    parser_infer_ortho = subparsers.add_parser("infer_ortho", help="Orthology and Gene duplication Inference")
    parser_infer_ortho.add_argument("-d", "--homolog_tree_dir", type=str, help="Folder containing homolog trees", required=True)
    
    # Create mutually exclusive group for outgroup options
    outgroup_group = parser_infer_ortho.add_mutually_exclusive_group(required=True)
    outgroup_group.add_argument("-o", "--outgroup_list", type=str, help="List of outgroup taxa, separated by commas, no spaces")
    outgroup_group.add_argument("-of", "--outgroup_file", type=str, help="File containing the outgroup taxa, each line is a taxon")
    
    parser_infer_ortho.add_argument("-t", "--tree_file_ending", type=str, help="File ending of the homolog trees", required=True)
    parser_infer_ortho.add_argument("-m", "--min_ingroup_taxa", type=int, help="Minimum number of ingroup taxa, default is 3")
    parser_infer_ortho.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_infer_ortho.add_argument("-p", "--just_pruning", action="store_true", help="Only produce pruned ortholog trees")
    parser_infer_ortho.add_argument("-np", "--no_pruning", action="store_true", help="Only produce unpruned ortholog trees")
    parser_infer_ortho.set_defaults(func=infer_ortho_main)

    # Subcommand: h2o map_dupl
    parser_map_duplications = subparsers.add_parser("map_dupl", help="Map gene duplications")
    parser_map_duplications.add_argument("-d", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_map_duplications.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_map_duplications.add_argument("-t", "--species_tree_file", type=str, help="Species tree file", required=True)
    parser_map_duplications.set_defaults(func=map_duplications_main)

    # Subcommand: h2o extract_wgd_trees
    parser_extract_wgd_trees = subparsers.add_parser("extract_wgd_trees", help="Extract and concatenate homolog trees that shows gene duplications at WGD events")
    parser_extract_wgd_trees.add_argument("-d", "--processed_tree_dir", type=str, help="Folder containing processed trees", required=True)
    parser_extract_wgd_trees.add_argument("-n", "--wgd_nodes", type=str, help="List of WGD node numbers, separated by commas, no spaces", required=True)
    parser_extract_wgd_trees.add_argument("-f", "--duplication_counts_file", type=str, help="Duplication counts file")
    parser_extract_wgd_trees.add_argument("-od", "--output_directory", type=str, help="Output directory")
    parser_extract_wgd_trees.set_defaults(func=extract_wgd_trees_main)

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