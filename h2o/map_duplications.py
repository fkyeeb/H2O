"""
This is the script for the subcommand map_dupl
it maps gene duplications to species tree
and output the results in a tsv and a tree file
"""

from h2o import tree_reader as t
from h2o.utils import (
    check_path,
    precompute_leaf_names_number_nodes,
    transform_elapsed_time,
    get_deepest_dup_parent,
    get_deepest_non_dup_parent,
    get_tips_in_ascending_order
)
import os
import time
import copy
import sys

def number_all_nodes(tree,output_folder):
    """
    Give each node a number
    Writes a new consensus tree with node numbers
    And also writes the bipartitions of the tree

    :param [Node] tree: root node of the consensus tree
    :param str output_folder: output directory path
    :return [dict]: dictionary with node numbers as keys and to store dupl counts as values later
    :return [dict]: dictionary with node labels as keys and tips [list] as values
    :return [file]: tsv file handle for continued writing
    """

    tf = open(output_folder + "consensus_tree_numbered.tre", 'w')
    tsv_file = open(output_folder + "duplication_counts.tsv", 'w')
    tsv_c_file = open(output_folder + "duplication_counts_ils_corrected.tsv", 'w')
    tsv_file.write("tree\t")
    tsv_c_file.write("tree\t")
    
    dup_node_counts = {}
    dup_node_counts_corrected = {}
    
    # Pre-compute leaf names to avoid repeated lvsnms() calls
    leaf_cache, parent2children = precompute_leaf_names_number_nodes(tree,label=True,return_set=True,return_children=True)
    
    for node_label in leaf_cache:
        dup_node_counts[node_label] = 0
        dup_node_counts_corrected[node_label] = 0
        tsv_file.write(node_label + "\t")
        tsv_c_file.write(node_label + "\t")
    
    tf.write(tree.get_newick_repr(showbl = True) + ";\n")
    tf.close()
    tsv_file.write("n/a\n")
    tsv_c_file.write("n/a\n")
    
    dup_node_counts["not_found"] = 0
    dup_node_counts_corrected["not_found"] = 0

    return dup_node_counts,dup_node_counts_corrected,leaf_cache,parent2children,tsv_file,tsv_c_file


def map_dup(dup_tree,dup_node_counts,dup_node_counts_corrected,sp_tree_leaf_cache,tsv_file,tsv_c_file,tree_file,sp_tree_parent2children,id2sp):
    """
    Map duplications from dup_tree to species tree
    Write a line in the duplication counts file with the number of duplications at each node for this tree
    """

    previous_count = dup_node_counts.copy()
    previous_count_corrected = dup_node_counts_corrected.copy()
    ks_pairs = []
    ks_pairs_ils = []
    
    # Pre-compute leaf names to avoid repeated lvsnms() calls
    leaf_cache = precompute_leaf_names_number_nodes(dup_tree,id2sp=id2sp)
    all_sp_tree_tips = set(sp_tree_leaf_cache["0"])

    for node in dup_tree.iternodes():
        if node.label == "D":
            deepest_dup_parent = get_deepest_dup_parent(node)
            deepest_non_dup_parent = get_deepest_non_dup_parent(deepest_dup_parent)
            dup_tips = set(leaf_cache[deepest_dup_parent.cache_label])
            deepest_non_dup_parent_tips = set(leaf_cache[deepest_non_dup_parent.cache_label])
            other_dup_tips = deepest_non_dup_parent_tips - dup_tips

            # Optimize bipartition matching with early termination
            best_match = None
            best_size = float('inf')

            # finding the best match on the species tree
            for node_label in sp_tree_leaf_cache:
                node_ingroup_tips = sp_tree_leaf_cache[node_label]
                node_outgroup_tips = all_sp_tree_tips - node_ingroup_tips
                if dup_tips <= node_ingroup_tips and other_dup_tips <= node_outgroup_tips:
                    bipartition_size = len(node_ingroup_tips)
                    if bipartition_size < best_size:
                        best_size = bipartition_size
                        best_match = node_label
            if best_match is None:
                dup_node_counts["not_found"] += 1
                dup_node_counts_corrected["not_found"] += 1
                continue
            dup_node_counts[best_match] += 1

            # counting tips for each duplication event
            child_tips = []
            for child in node.children:
                if child.istip:
                    tips = set([child.label])
                    child_tips.append(tips)
                else:
                    tips = set(leaf_cache[child.cache_label]) 
                    child_tips.append(tips)
            
            # ILS correction
            best_match_children =  sp_tree_parent2children[best_match]
            dup_children_tips = []
            best_match_corrected = best_match
            for child in deepest_dup_parent.children:
                if child.istip:
                    dup_children_tips.append(set([child.label]))
                else:
                    dup_children_tips.append(set(leaf_cache[child.cache_label]))
            sp_children_tips = []
            for index in range(2):
                sp_child = best_match_children[index]
                if type(sp_child) is int: # a node not a tip
                    sp_child = str(sp_child)
                    sp_children_tips.append(sp_tree_leaf_cache[sp_child])
                else: # a tip
                    sp_children_tips.append(set([sp_child]))
            for index in range(2):
                if len(sp_children_tips[index-1]) > 1:
                    sp_child_tips = sp_children_tips[index]
                    overlap_size = 1
                    if len(sp_child_tips) <= 2:
                        overlap_size = 0
                        # more restrictive for 1 and 2 tip clades
                    len1 = len(sp_child_tips & dup_children_tips[0])
                    len2 = len(sp_child_tips & dup_children_tips[1])
                    if len1 <= overlap_size or len2 <= overlap_size:
                        best_match_corrected = str(best_match_children[index-1])
                        break
            dup_node_counts_corrected[best_match_corrected] += 1

            # recording ks pairs for duplicated genes
            if id2sp:
                ks_pairs_temporary = {}
                cluster = tree_file.split("_rooted_processed_id.tre")[0]
                cluster = cluster.split(".")[0]
                if node == deepest_dup_parent:
                    overlap = dup_children_tips[0] & dup_children_tips[1]
                    for index in range(2):
                        child = deepest_dup_parent.children[index]
                        for node in child.iternodes():
                            if node.label in overlap:
                                if node.label not in ks_pairs_temporary:
                                    ks_pairs_temporary[node.label] = [[],[]]
                                ks_pairs_temporary[node.label][index].append(node.cache_label)
                for pairs in ks_pairs_temporary.values():
                    for gene1 in pairs[0]:
                        for gene2 in pairs[1]:
                            combine = gene1 + "," + gene2
                            ks_pairs.append([int(best_match),cluster,combine])
                            ks_pairs_ils.append([int(best_match_corrected),cluster,combine])
        
    # Write to file using provided file handle
    if id2sp:
        tree_file_ending = "_rooted_processed_id.tre"
    else:
        tree_file_ending = "_rooted_processed.tre"
    tsv_file.write(tree_file.split(tree_file_ending)[0] + "\t")
    tsv_c_file.write(tree_file.split(tree_file_ending)[0] + "\t")
    for node in dup_node_counts:
        tsv_file.write(str(dup_node_counts[node]-previous_count[node]) + "\t")
        tsv_c_file.write(str(dup_node_counts_corrected[node]-previous_count_corrected[node]) + "\t")
    tsv_file.write("\n")
    tsv_c_file.write("\n")

    return dup_node_counts,dup_node_counts_corrected,ks_pairs,ks_pairs_ils

def main(args):
    """
    Main function to map duplications
    """
    species_tree_file = check_path(args.species_tree_file,is_folder=False,error_if_not_exists=True)

    processed_tree_folder = check_path(args.processed_tree_dir) + "unpruned/"
    processed_tree_folder = check_path(processed_tree_folder,error_if_not_exists=True)

    default_output_folder = "other_output/"
    output_folder = check_path(args.output_directory,default_path=default_output_folder,create_if_not_exists=True)

    if args.id2sp_file:
        id2sp = check_path(args.id2sp_file,error_if_not_exists=True,is_folder=False)
    else:
        id2sp = None

    with open(species_tree_file,"r") as f:
        sp_tree = t.read_tree_string(f.readline().strip())

        print("------------------------------------------------------------\n")
        print(time.ctime() + "\n")
        start_time = time.time()

        # Initialize with consolidated file operations
        (dup_node_counts, dup_node_counts_corrected,  
        sp_tree_leaf_cache, sp_tree_parent2children, tsv_file, 
        tsv_c_file) = number_all_nodes(sp_tree, output_folder)
        
        ks_pairs = []
        ks_pairs_ils_corrected = []
        if id2sp:
            tree_file_ending = "_rooted_processed_id.tre"
        else:
            tree_file_ending = "_rooted_processed.tre"
        for tree_file in os.listdir(processed_tree_folder):
            if tree_file.endswith(tree_file_ending):
                with open(processed_tree_folder + tree_file, "r") as f:
                    dup_tree = t.read_tree_string(f.readline().strip())
                    (dup_node_counts, dup_node_counts_corrected, 
                    ks_pairs_extend, ks_pairs_ils_extend) = map_dup(
                    dup_tree, dup_node_counts, dup_node_counts_corrected, 
                    sp_tree_leaf_cache, tsv_file, tsv_c_file, tree_file, 
                    sp_tree_parent2children,id2sp)
                    ks_pairs.extend(ks_pairs_extend)
                    ks_pairs_ils_corrected.extend(ks_pairs_ils_extend)

        tsv_file.close()
        tsv_c_file.close()

        # Write ks pairs to file
        if id2sp:
            ks_pairs.sort(key=lambda x: (x[0], x[1]))
            with open(output_folder + "ks_pairs.tsv", 'w') as file:
                file.write("node\tcluster\tks_pair\n")
                for data in ks_pairs:
                    data[0] = str(data[0])
                    file.write("\t".join(data) + "\n")

            ks_pairs_ils_corrected.sort(key=lambda x: (x[0], x[1]))
            with open(output_folder + "ks_pairs_ils_corrected.tsv", 'w') as file:
                file.write("node\tcluster\tks_pair\n")
                for data in ks_pairs_ils_corrected:
                    data[0] = str(data[0])
                    file.write("\t".join(data) + "\n")

        # Writes numbers of gene duplications at each node in the numbered species tree
        with open(output_folder + "consensus_tree_numbered.tre", 'a') as file:
            sp_tree_copy = copy.deepcopy(sp_tree)
            for node in sp_tree.iternodes():
                if not node.istip and node.label != "":
                    node.label = str(dup_node_counts[node.label])
            file.write(sp_tree.get_newick_repr(showbl=True) + ";\n")
            for node in sp_tree_copy.iternodes():
                if not node.istip and node.label != "":
                    node.label = str(dup_node_counts_corrected[node.label])
            file.write(sp_tree_copy.get_newick_repr(showbl=True) + ";\n")

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"\nDone with duplication mapping. Total time elapsed: {elapsed}\n")
    print(f"Output files write to {output_folder}\n")

    print("------------------------------------------------------------\n\n")