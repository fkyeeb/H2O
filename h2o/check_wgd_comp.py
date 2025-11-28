"""
check structure of WGD and some near nodes
"""

from h2o import tree_reader as t
from h2o.utils import (
    check_path,
    precompute_leaf_names_number_nodes,
    transform_elapsed_time,
    get_deepest_dup_parent,
    get_deepest_non_dup_parent
)
import sys
import time
from collections import Counter

def make_counter_results_string(counter):
    """
    Return the message of the counts in the counter object

    :param [Counter] counter: Counter object with counts
    """
    message = "  "
    for count in sorted(counter.keys()):
        message += f"{counter[count]} trees has {count} gene duplications;  "
    return message[:-3] + "\n"

def find_duplication_node(tree,leaf_cache,node2find,tips_dict):
    """
    Find the duplication node with node2find label

    :param [Node] tree: root node of the tree
    :param [dict] leaf_cache: dictionary with node numbers as keys and set(tips) as values
    :param node2find: node.cache_label of the node to find
    :param tips_dict: dictionary with node labels as keys and list of sets of tips as values
    :return [Node]: the duplication node found, None if not found
    """
    nodes2return = []

    for node in tree.iternodes():
        if node.label == "D":
            deepest_dup_parent = get_deepest_dup_parent(node)
            deepest_non_dup_parent = get_deepest_non_dup_parent(deepest_dup_parent)
            dup_tips = set(leaf_cache[deepest_dup_parent.cache_label])
            deepest_non_dup_parent_tips = set(leaf_cache[deepest_non_dup_parent.cache_label])
            other_dup_tips = deepest_non_dup_parent_tips - dup_tips

            tips = tips_dict[node2find][0]
            other_tips = tips_dict[node2find][1]
            if dup_tips <= tips and other_dup_tips <= other_tips:
                nodes2return.append(node)
    
    if len(nodes2return) == 0:
        return None
    else:
        return nodes2return

def get_node_of_wgd_tips(node,leaf_cache,wgd_tip_sets):
    """
    This function takes tree and return a set of tips under the node corresponding to wgd node
    """
    current_node = node

    while current_node != None:
        current_tips = set(leaf_cache[current_node.cache_label])
        deepest_non_dup_parent = get_deepest_non_dup_parent(current_node)
        deepest_non_dup_parent_tips = set(leaf_cache[deepest_non_dup_parent.cache_label])
        current_other_tips = deepest_non_dup_parent_tips - current_tips

        if current_tips <= wgd_tip_sets[0] and current_other_tips <= wgd_tip_sets[1]:
            if current_node == node:
                return True
            else:
                return False
        
        current_node = current_node.parent
    
    return None

def main(args):
    """
    Main function to map gene losses
    """

    processed_tree_folder = check_path(args.processed_tree_dir) + "unpruned/"
    processed_tree_folder = check_path(processed_tree_folder,error_if_not_exists=True)

    default_output_folder = "other_output/"
    output_folder = check_path(args.output_directory,default_path=default_output_folder,create_if_not_exists=True)
    
    wgd_node = args.wgd_node
    nodes2tree_names = {wgd_node:[]}
    if args.connected_nodes is None:
        connected_nodes = None
    else:
        connected_nodes = args.connected_nodes.split(",")
        for index in range(len(connected_nodes)):
            nodes2tree_names[connected_nodes[index]] = []
        try:
            for node in nodes2tree_names:
                int(node)
        except ValueError:
            print("Error: Node numbers must be integers.")
            sys.exit(2)

    default_dup_dir = "other_output/"
    duplication_counts_dir = check_path(args.duplication_counts_dir,default_path=default_dup_dir,error_if_not_exists=True)

    duplication_counts_file = check_path(duplication_counts_dir + "duplication_counts.tsv",is_folder=False,error_if_not_exists=True)
    numbered_tree_file = check_path(duplication_counts_dir + "consensus_tree_numbered.tre",is_folder=False,error_if_not_exists=True)

    print("------------------------------------------------------------\n")
    print(time.ctime() + "\n")
    start_time = time.time()
    
    with open(numbered_tree_file,"r") as f:
        numbered_tree = t.read_tree_string(f.readline().strip())
    
    sp_tree_leaf_cache = precompute_leaf_names_number_nodes(numbered_tree,use_label=True,return_set=True)
    all_sp_tree_tips = sp_tree_leaf_cache["0"]
    tips_dict = {}
    duplication_counts_data = {}
    for node in nodes2tree_names:
        duplication_counts_data[node] = []
        tips = sp_tree_leaf_cache[node]
        tips_dict[node] = [tips,all_sp_tree_tips - tips]

    with open(duplication_counts_file,"r") as f:
        f.readline()
        for line in f:
            splt = line.strip().split("\t")
            tree_name = splt.pop(0)

            if any(int(node) >= len(splt) or int(node) < 0 for node in nodes2tree_names):
                print("Error: Node numbers provided is out of range for the consensus tree.")
                sys.exit(2)
            for node in nodes2tree_names:
                count = int(splt[int(node)])
                if count > 0:
                    nodes2tree_names[node].append(tree_name)
                    duplication_counts_data[node].append(count)
    
    # printing overall summary
    print(f"Duplication in WGD node {wgd_node} is found in {len(nodes2tree_names[wgd_node])} trees.")
    counter_message = make_counter_results_string(Counter(duplication_counts_data[wgd_node]))
    print(counter_message)

    for node in connected_nodes:
        print(f"Duplication in connected node {node} is found in {len(nodes2tree_names[node])} trees.")
        overlap_trees = set(nodes2tree_names[wgd_node]) & set(nodes2tree_names[node])
        print(f"  Overlapping trees with WGD node: {len(overlap_trees)}")
        counter_message = make_counter_results_string(Counter(duplication_counts_data[node]))
        print(counter_message)

    ### wgd node
    wgd_retain_node_counts = {}
    for node in nodes2tree_names:
        if node != wgd_node:
            wgd_retain_node_counts[node] = 0
    for tree_name in nodes2tree_names[wgd_node]:
        with open(processed_tree_folder + tree_name + "_rooted_processed.tre","r") as f:
            tree = t.read_tree_string(f.readline().strip())
        
        leaf_cache = precompute_leaf_names_number_nodes(tree)
        wgd_node_in_this_tree = find_duplication_node(tree,leaf_cache,wgd_node,tips_dict)
        if not wgd_node_in_this_tree:
            print(f"Warning: WGD node {wgd_node} not found in tree {tree_name}.")
            continue
        else:
            for node in wgd_node_in_this_tree:
                for child in node.children:
                    for node in connected_nodes:
                        if set(leaf_cache[child.cache_label]) <= tips_dict[node][0]:
                            wgd_retain_node_counts[node] += 1
                            break
    print(f"\nAmong {len(nodes2tree_names[wgd_node])} trees that has gene duplication at WGD node {wgd_node}:")
    for node in wgd_retain_node_counts:
        print(f"  Number of times where only connected node {node} is retained in one duplicated gene copy: {wgd_retain_node_counts[node]}")
    
    # compare within connected nodes
    same_w_wgd_counts = [0] * len(connected_nodes)
    for index in range(len(connected_nodes)):
        node_label = connected_nodes[index]
        for tree_name in nodes2tree_names[node_label]:
            with open(processed_tree_folder + tree_name + "_rooted_processed.tre","r") as f:
                tree = t.read_tree_string(f.readline().strip())
            
            leaf_cache = precompute_leaf_names_number_nodes(tree)
            this_node_in_this_tree = find_duplication_node(tree,leaf_cache,node_label,tips_dict)
            if not this_node_in_this_tree:
                print(f"Warning: connected node {node_label} not found in tree {tree_name}.")
                continue
            else:
                for node in this_node_in_this_tree:
                    deepest_dup_parent = get_deepest_dup_parent(node)
                    if deepest_dup_parent.parent is None:
                        continue
                    bool = get_node_of_wgd_tips(deepest_dup_parent,leaf_cache,tips_dict[wgd_node])
                    if bool is None:
                        print(f"Warning: WGD node {wgd_node} not found in tree {tree_name} when checking connected node {node_label}.")
                        continue
                    if bool is True:
                        same_w_wgd_counts[index] += 1    
        print(f"\nAmong {len(nodes2tree_names[node_label])} trees that has gene duplication at connected node {node_label}:")
        print(f"  Number of times where no tips from the wgd node is present in the tree other than those of node {node_label}: {same_w_wgd_counts[index]}")

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"\n\nDone with WGD composition check. Total time elapsed: {elapsed}")

    print("\n------------------------------------------------------------\n\n")