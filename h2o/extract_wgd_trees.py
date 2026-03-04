"""
This is the script for the subcommand extract_wgd_trees
it extracts and concatenates homolog trees that shows gene duplications at WGD events
"""
from h2o.utils import (
    check_path,
    run_shell_command,
    transform_elapsed_time
)
import os
import sys
import time

def cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, pruning, tree_file_ending):
    if os.path.exists(processed_tree_folder +  pruning + "/"):
        with open(output_directory + "cat_" + pruning + "_wgd_trees.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            for i, tree in enumerate(wgd_trees):
                if i == 0:
                    f.write("cat " + processed_tree_folder + pruning + "/" + tree + "*ortho*.tre > " + output_directory + "ASTRAL_in_" + pruning + "_wgd" + tree_file_ending + ".tre\n")
                else:
                    f.write("cat " + processed_tree_folder + pruning + "/" + tree + "*ortho*.tre >> " + output_directory + "ASTRAL_in_" + pruning + "_wgd" + tree_file_ending + ".tre\n")
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
    wgd_nodes.sort()
    tree_file_ending = "_n" + str(wgd_nodes[0])

    processed_tree_folder = check_path(args.processed_tree_dir,error_if_not_exists=True)
    
    default_dup_dir = "other_output/"
    duplication_counts_dir = check_path(args.duplication_counts_dir,default_path=default_dup_dir,error_if_not_exists=True)

    duplication_counts_file = check_path(duplication_counts_dir + "duplication_counts.tsv",is_folder=False,error_if_not_exists=True)
    
    output_directory = check_path(args.output_directory,default_path=duplication_counts_dir,create_if_not_exists=True)
    
    print("------------------------------------------------------------\n")
    print(time.ctime() + "\n")
    start_time = time.time()

    wgd_trees = []
    with open(duplication_counts_file, "r") as f:
        f.readline()
        for line in f:
            splt = line.strip().split("\t")
            tree_name = splt.pop(0)

            # Check if any wgd_node index is out of range for the current line
            if any(node >= len(splt) or node < 0 for node in wgd_nodes):
                print("Error: WGD node numbers provided is out of range for the consensus tree.")
                sys.exit(2)

            if any(int(splt[node]) > 0 for node in wgd_nodes):
                wgd_trees.append(tree_name)
    
    unpruned = cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, "unpruned",tree_file_ending)
    
    pruned = cat_n_run_bash(wgd_trees, processed_tree_folder, output_directory, "pruned",tree_file_ending)

    if not unpruned and not pruned:
        print("Error: No processed ortholog tree folder found.")
        sys.exit(2)
    
    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"Done with extracting WGD trees. Total time elapsed: {elapsed}")

    print("\n------------------------------------------------------------\n\n")