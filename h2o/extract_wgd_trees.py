"""
This is the script for the subcommand extract_wgd_trees
it extracts and concatenates homolog trees that shows gene duplications at WGD events
"""
from h2o.utils import (
    check_folder,
    run_shell_command
)
import os
import sys
from pathlib import Path

def cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, pruning):
    if os.path.exists(processed_tree_folder +  pruning + "/"):
        with open(output_directory + "cat_" + pruning + "_wgd_trees.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            for i, tree in enumerate(wgd_trees):
                if i == 0:
                    f.write("cat " + processed_tree_folder + pruning + "/" + tree + "_ortho[0-9].tre > " + output_directory + "ASTRAL_in_" + pruning + "_wgd_trees.tre\n")
                else:
                    f.write("cat " + processed_tree_folder + pruning + "/" + tree + "_ortho[0-9].tre >> " + output_directory + "ASTRAL_in_" + pruning + "_wgd_trees.tre\n")
        run_shell_command("bash " + output_directory + "cat_" + pruning + "_wgd_trees.sh")  
        return True
    else:
        print(pruning + " processed tree folder not found.\n")
        return False

def main(args):
    """
    Main function to extract WGD trees
    """
    
    wgd_nodes = args.wgd_nodes.split(",")
    try:
        wgd_nodes = [int(node) for node in wgd_nodes]
    except ValueError:
        print("Error: WGD node numbers must be integers.")
        sys.exit(2)

    processed_tree_folder = check_folder(args.processed_tree_dir,error_if_not_exists=True)
    
    default_dup_dir = str(Path(processed_tree_folder).parent) + "/other_output/"
    duplication_counts_file_directory = check_folder(args.duplication_counts_file_directory,default_folder=default_dup_dir,error_if_not_exists=True)

    if os.path.exists(duplication_counts_file_directory + "duplication_counts.tsv"):
        duplication_counts_file = duplication_counts_file_directory + "duplication_counts.tsv"
    if not os.path.exists(duplication_counts_file):
        print("Error: No duplication counts file found.")
        sys.exit(2)
    
    output_directory = check_folder(args.output_directory,default_folder=duplication_counts_file_directory,create_if_not_exists=True)
    
    wgd_trees = []
    with open(duplication_counts_file, "r") as f:
        f.readline()
        for line in f:
            splt = line.strip().split("\t")
            tree_name = splt.pop(0)

            # Check if any wgd_node index is out of range for the current line
            if any(node >= len(splt) or node < 0 for node in wgd_nodes):
                print("Error: WGD node numbers provided is out of range for the summary tree.")
                sys.exit(2)

            if any(int(splt[node]) > 0 for node in wgd_nodes):
                wgd_trees.append(tree_name)
    
    unpruned = cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, pruning="unpruned")
    
    pruned = cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, pruning="pruned")

    if not unpruned and not pruned:
        print("Error: No processed ortholog tree folder found.")
        sys.exit(2)