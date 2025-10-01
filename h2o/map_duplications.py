"""
This is the script for the subcommand map_dupl
it maps gene duplications to species tree
and output the results in a tsv file
"""

from h2o import tree_reader as t
from h2o.utils import (
    check_folder,
    precompute_leaf_names_number_nodes,
    transform_elapsed_time
)
import os
import sys
import time
from pathlib import Path

def get_deepest_dup_parent(node):
    """
    This function takes a node and returns its deepest duplication parent.
    If node has no parent, return node
    """
    current_node = node

    while current_node.parent != None:
        if current_node.parent.label != "D":
            return current_node
        else:
            current_node = current_node.parent

    return current_node

def number_all_nodes(tree,output_folder):
    """
    Give each node a number
    Writes a new summary tree with node numbers
    And also writes the bipartitions of the tree

    :param [Node] tree: root node of the summary tree
    :param str output_folder: output directory path
    :return [dict]: dictionary with node numbers as keys and to store dupl counts as values later
    :return [dict]: dictionary with node labels as keys and tips [list] as values
    :return [file]: tsv file handle for continued writing
    """

    tf = open(output_folder + "summary_tree_numbered" + ".tre", 'w')
    tsv_file = open(output_folder + "duplication_counts" + ".tsv", 'w')
    tsv_file.write("tree\t")
    
    node_numbers = {}
    
    # Pre-compute leaf names to avoid repeated lvsnms() calls
    leaf_cache = precompute_leaf_names_number_nodes(tree,label=True)
    
    for node_label in leaf_cache:
        node_numbers[node_label] = 0
        tsv_file.write(node_label + "\t")
    
    tf.write(tree.get_newick_repr(showbl = True) + ";\n")
    tf.close()
    tsv_file.write("n/a\n")
    
    node_numbers["not_found"] = 0

    return node_numbers,leaf_cache,tsv_file


def map_dup(dup_tree,node_numbers,sp_tree_leaf_cache,tsv_file,tree_file):
    """
    Map duplications from dup_tree to species tree
    Write a line in the duplication counts file with the number of duplications at each node for this tree
    """

    previous_count = node_numbers.copy()
    
    # Pre-compute leaf names to avoid repeated lvsnms() calls
    leaf_cache = precompute_leaf_names_number_nodes(dup_tree)
    all_dup_tips = set(leaf_cache["0"])
    all_sp_tree_tips = set(sp_tree_leaf_cache["0"])

    for node in dup_tree.iternodes():
        if node.label == "D":
            deepest_dup_parent = get_deepest_dup_parent(node)
            dup_tips = set(leaf_cache[deepest_dup_parent.cache_label])
            other_dup_tips = all_dup_tips - dup_tips

            # Optimize bipartition matching with early termination
            best_match = None
            best_size = float('inf')
            
            for node_label in sp_tree_leaf_cache:
                node_ingroup_tips = sp_tree_leaf_cache[node_label]
                node_outgroup_tips = all_sp_tree_tips - node_ingroup_tips
                if dup_tips <= node_ingroup_tips and other_dup_tips <= node_outgroup_tips:
                    bipartition_size = len(node_ingroup_tips)
                    if bipartition_size < best_size:
                        best_size = bipartition_size
                        best_match = node_label

            if best_match is None:
                node_numbers["not_found"] += 1
                continue
                
            node_numbers[best_match] += 1
    
    # Write to file using provided file handle
    tsv_file.write(tree_file.split("_rooted_processed.tre")[0] + "\t")
    for node in node_numbers:
        tsv_file.write(str(node_numbers[node]-previous_count[node]) + "\t")
    tsv_file.write("\n")

    return node_numbers

def main(args):
    """
    Main function to map duplications
    """
    species_tree_file = args.species_tree_file
    if not os.path.exists(species_tree_file):
        print("Error: The species tree file " + species_tree_file + " does not exist.")
        sys.exit(2)

    processed_tree_folder = check_folder(args.processed_tree_dir) + "unpruned/"
    processed_tree_folder = check_folder(processed_tree_folder,error_if_not_exists=True)

    default_output_folder = str(Path(processed_tree_folder).parent.parent) + "/other_output/"
    output_folder = check_folder(args.output_directory,default_folder=default_output_folder,create_if_not_exists=True)

    with open(species_tree_file,"r") as f:
        sp_tree =  t.read_tree_string(f.readline().strip())

    print("------------------------------------------------------------\n\n")
    print(time.ctime() + "\n")
    start_time = time.time()

        # Initialize with consolidated file operations
    node_numbers,sp_tree_leaf_cache,tsv_file = number_all_nodes(sp_tree,output_folder)
    
    for tree_file in os.listdir(processed_tree_folder):
        if tree_file.endswith("rooted_processed.tre"):
            with open(processed_tree_folder + tree_file,"r") as f:
                dup_tree = t.read_tree_string(f.readline().strip())
            node_numbers = map_dup(dup_tree,node_numbers,sp_tree_leaf_cache,tsv_file,tree_file)

    tsv_file.close()

    # Writes numbers of gene duplications at each node in the numbered species tree
    with open(output_folder + "summary_tree_numbered" + ".tre", 'a') as file:
        for node in sp_tree.iternodes():
            if not node.istip and node.label != "":
                node.label = str(node_numbers[node.label])
        file.write(sp_tree.get_newick_repr(showbl = True) + ";\n")

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"Done with duplication mapping. Total time elapsed: {elapsed}")

    print("\n------------------------------------------------------------\n\n")