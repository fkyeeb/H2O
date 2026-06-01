"""
This is the script for the subcommand infer_ortho
It roots the homolog trees
infers gene duplications
and outputs the ortholog trees
"""

from h2o import tree_reader as t
from h2o.node import Node
from h2o.utils import (
    get_sister,
    check_path,
    precompute_leaf_names_number_nodes,
    transform_elapsed_time
)

import copy
import os
import sys
import time
from collections import Counter

def check_taxa_counts(filename,outgroup,all_tips):
    """
    This function checks if the outgroups are in the tree
    and if the outgroup + ingroup both duplicated. 
    If so, skip the tree for now.
    Also skip if the tree have less than 3 ingroup taxa.

    :param filename: name of the tree file
    :param outgroup: list of outgroups
    :param all_tips: list of all tips in the tree
    :return: True if the outgroup is not duplicated, False otherwise
    :rtype: bool
    """
    # record outgroup number
    outgroup_count = Counter(tip for tip in all_tips if tip in outgroup)
    ingroup_count = set(tip for tip in all_tips if tip not in outgroup)

    # check if there are at least 3 ingroup taxa
    if len(ingroup_count) < 3:
        print("Less than 3 ingroup taxa, " + filename + " is skipped.\n")
        return False

    # check if outgroup exists in the tree
    if len(outgroup_count) == 0:
        print("None of the outgroups is in the tree, " + filename + " is skipped.\n")
        return False
    # elif len(outgroup_count) == 1:
    #     print("Only one outgroup is in the tree, " + filename + " is skipped.\n")
    #     return False
    
    # check if outgroup is duplicated
    limit = len(outgroup) * 2 // 3
    duplicated_outgroup = sum(1 for count in outgroup_count.values() if count > 1)  
    if duplicated_outgroup !=0 and duplicated_outgroup >= limit:
        print("Outgroup is duplicated, " + filename + " is skipped.\n")
        return False
    return True

def locate_outgroup(file_name,tree,outgroup,leaf_cache):
    """
    This function takes a tree and a list of outgroups and returns the node of the clade.

    :param file_name: name of the tree file
    :param tree: root node of the tree
    :param outgroup: list of outgroups
    :param leaf_cache: precomputed leaf names cache
    :return: node of the clade and a boolean indicating if the tree needs to be rerooted
    :rtype: Node, bool
    """
    # make sure to only save the biggest outgroup clade
    to_root = True
    for node in tree.iternodes(order="postorder"):
        if not node.istip:
            child_tips = leaf_cache[node.cache_label]
        else:
            child_tips = [node.label]
        if node.parent != None:
            other_tip_counts = Counter(leaf_cache["0"]) - Counter(child_tips)
            other_tips = set(other_tip_counts.elements())
            child_tips = set(child_tips)
            if child_tips <= outgroup:
                if other_tips.isdisjoint(outgroup):
                    if node.parent.parent == None:
                        # tree is already rooted at the outgroup, no need to reroot
                        to_root = False
                        return tree, to_root
                    # found monophyletic outgroup, root normally
                    return node,to_root
            elif child_tips & outgroup == set():
                if other_tips <= outgroup:
                    if node.parent.parent == None:
                        to_root = False
                        # tree is already rooted at the outgroup, no need to reroot
                        return tree,to_root
                    else:
                        # if the tree is rooted within the outgroup, cut off the outgroup
                        node.parent = None
                        to_root = False
                        return node,to_root

    return "Outgroup is polyphyletic, " + file_name + " is skipped.\n",to_root

def correct_parent(node):
    """
    After copying a node, this function makes sure that the children of the copied node has the correct parent.
    """
    for child in node.children:
        if child.parent != node:
            child.parent = node

def root_tree(file_name,tree,outgroup,leaf_cache):
    """
    This function takes a tree and a list of outgroups and returns the tree rooted at the outgroup.

    :param tree: root node of the tree to be rooted
    :param outgroup: one outgroup for now
    :param leaf_cache: precomputed leaf names cache
    """
    new_root = Node()
    dummy = new_root

    outgroup_node,to_root = locate_outgroup(file_name,tree,outgroup,leaf_cache)

    if isinstance(outgroup_node,str):
        # when outgroup is polyphyletic, print the error message
        print(outgroup_node)
        return None
    elif not to_root:
        # when outgroup is already rooted
        # outgroup_node should already be the root
        return outgroup_node
    current = outgroup_node

    # set up root
    dummy.add_child(copy.copy(current))
    dummy.add_child(Node())
    outgroup_length = current.length/2
    for child in dummy.children:
        child.length = outgroup_length
    
    # add internal nodes
    while current.parent.parent != None:
        dummy = dummy.children[1]
        if dummy.parent.parent != None:
            dummy.length = current.parent.length
        sis = get_sister(current)
        dummy.add_child(copy.copy(sis))
        current = current.parent
        dummy.add_child(Node())
    
    # add the last node
    sis = get_sister(current)
    dummy.remove_child(dummy.children[1])  # remove the last dummy node
    dummy.add_child(copy.copy(sis))

    # parent is not reassigned for the untouched nodes, correct that
    [correct_parent(node) for node in new_root.iternodes() if not node.istip]
    
    return new_root

def prune_selected_taxa(node,selected_taxa):
    """
    This function takes a duplication node and prunes the taxa that are not duplicated.
    """
    node2prune = [n for n in node.iternodes() if n.istip and n.label in selected_taxa]

    # make sure no tip is skipped
    [n.prune() for n in node2prune]

def check_prune_dup(node,min_dupl_tip_overlap,min_dupl_percentage_overlap,bool,leaf_cache,keep_ss_dup):
    """
    This function checks if the node is a duplication node.
    And also saves any taxa that are not duplicated.
    """
    child1_node = node.children[0]
    child2_node = node.children[1]
    if child1_node.istip:
        child1 = set([child1_node.label])
    else:
        child1 = leaf_cache[child1_node.cache_label]

    if child2_node.istip:
        child2 = set([child2_node.label])
    else:
        child2 = leaf_cache[child2_node.cache_label]

    common = child1 & child2
    if len(child1) <= len(child2):
        child2prune = child1_node
        percentage_common = len(common) / len(child2)
    else:
        child2prune = child2_node
        percentage_common = len(common) / len(child1)

    if len(common) > 0:
        # prune overlap taxa if they are less than the minimum duplication tip overlap
        # or the percentage overlap to the bigger clade is less than the minimum duplication percentage overlap
        if len(common) < min_dupl_tip_overlap or percentage_common <= min_dupl_percentage_overlap:
            if keep_ss_dup:
                if common == child1 and common == child2:
                    return
            prune_selected_taxa(child2prune,common)
            return
        # recognize duplication node if the overlap is big enough
        elif len(common) >= min_dupl_tip_overlap and percentage_common > min_dupl_percentage_overlap:
            node.duplication = True
            node.label = "D"
            if bool:
                node.missing_dup = child1 ^ child2
                prune_selected_taxa(node,node.missing_dup)

def label_duplication_node(tree,min_dupl_tip_overlap,min_dupl_percentage_overlap,bool,leaf_cache,keep_ss_dup):
    """
    This function takes a node and checks if it is a duplication node if it is not a tip.
    It will prune the duplications that are smaller than the specified minimum clade size.
    It also saves any taxa that are not duplicated in the duplication node.
    """
    for node in tree.iternodes():
        if not node.istip:
            if len(node.children) == 2:
                # pruning happens at the same time as the duplication check
                check_prune_dup(node,min_dupl_tip_overlap,min_dupl_percentage_overlap,bool,leaf_cache,keep_ss_dup)
            else:
                print("Warning: Node with " + str(len(node.children)) + " children found.\n")
        
def get_orthologs(root):
    """
    This function splits all duplication nodes and into separate ortholog trees
    if ortholog clade size is at least 3, it will be outputted
    Else, it will be pruned off
    returns a list of ortholog trees
    """
    ortho_trees = []
    for node in root.iternodes(order="postorder"):
        if node.duplication:
            # split off the smaller clade
            clade1_size = len(node.children[0].lvsnms())
            clade2_size = len(node.children[1].lvsnms())
            if clade1_size <= clade2_size:
                if clade1_size >= 3:
                    ortho_trees.append(node.children[0])
                node.children[0].prune()
            else:
                if clade2_size >= 3:
                    ortho_trees.append(node.children[1])
                node.children[1].prune()
    if len(root.children) == 1:
        root = root.children[0]
    ortho_trees.append(root)

    return ortho_trees

def tree2id_tree(tree):
    id_tree = copy.deepcopy(tree)
    for n in id_tree.iternodes():
        if n.istip:
            n.label = n.cache_label
    return id_tree

def prune_or_not(tree,min_dupl_tip_overlap,min_dupl_percentage_overlap,output_directory,tree_name,bool,leaf_cache,id2sp,keep_ss_dup):

    label_duplication_node(tree,min_dupl_tip_overlap,min_dupl_percentage_overlap,bool,leaf_cache,keep_ss_dup)
    
    # remove node when there is only one child left
    if len(tree.children) == 1:
        tree = tree.children[0]
        tree.parent = None
    for n in tree.iternodes(order="postorder"):
        if len(n.children) == 1:
            p = n.parent
            if p != None:
                n.children[0].length += n.length
                p.add_child(n.children[0])
                p.remove_child(n)
            else:
                tree = n.children[0]
                tree.parent = None

    with open(output_directory + tree_name + "_rooted_processed.tre","w") as f:
        f.write(tree.get_newick_repr(showbl=True) + ";\n")
    
    if id2sp:
        id_tree = tree2id_tree(tree)

        with open(output_directory + tree_name + "_rooted_processed_id.tre","w") as f:
            f.write(id_tree.get_newick_repr(showbl=True) + ";\n")

    # get orthologs
    ortho_trees = get_orthologs(tree)

    l = len(ortho_trees)

    # Batch file operations for better performance
    ortho_output_contents = []
    for i in range(l-1,-1,-1):
        filename = output_directory + tree_name + "_ortho" + str(l-i) + ".tre"
        content = ortho_trees[i].get_newick_repr(showbl=True) + ";\n"
        ortho_output_contents.append((filename, content))
    
    # Write all files at once
    for filename, content in ortho_output_contents:
        with open(filename, "w") as f:
            f.write(content)
    
    if id2sp:
        ortho_id_output_contents = []
        for i in range(l-1,-1,-1):
            filename = output_directory + tree_name + "_ortho" + str(l-i) + "_id.tre"
            id_ortho_tree = tree2id_tree(ortho_trees[i])
            content = id_ortho_tree.get_newick_repr(showbl=True) + ";\n"
            ortho_id_output_contents.append((filename, content))
        
        for filename, content in ortho_id_output_contents:
            with open(filename, "w") as f:
                f.write(content)

def process_trees(tree,outgroup_list,tree_name,output_directory,min_dupl_tip_overlap,min_dupl_percentage_overlap,pruning,id2sp,keep_ss_dup):
    # root tree

    # if tree is not rooted, root it arbitrarily first
    if len(tree.children) == 3:
        half_height = tree.children[0].length/2
        tree.children[0].length = half_height
        child2 = tree.children[1]
        child3 = tree.children[2]
        tree.add_child(Node())
        tree.remove_child(child2)
        tree.remove_child(child3)
        tree.children[1].length = half_height
        tree.children[1].add_child(child2)
        tree.children[1].add_child(child3)

    # Precompute leaf names for efficiency
    leaf_cache = precompute_leaf_names_number_nodes(tree,return_set=False,id2sp=id2sp)

    if check_taxa_counts(tree_name,outgroup_list,leaf_cache["0"]):
        rooted_tree = root_tree(tree_name,tree,outgroup_list,leaf_cache)
    else:
        # if less than 2 ourgroups or if outgroup is duplicated, skip the tree
        rooted_tree = None

    if rooted_tree:
        with open(output_directory + tree_name + "_rooted.tre","w") as f:
            f.write(rooted_tree.get_newick_repr(showbl=True) + ";\n")

        # Precompute leaf names for the rooted tree
        rooted_leaf_cache = precompute_leaf_names_number_nodes(rooted_tree,return_set=True)

        # label duplication node and prune missing taxa
        rooted_tree_2prune = copy.deepcopy(rooted_tree)

        if pruning != False:
            prune_or_not(rooted_tree_2prune,min_dupl_tip_overlap,min_dupl_percentage_overlap,output_directory + "pruned/",tree_name,True,rooted_leaf_cache,id2sp,keep_ss_dup)
        if pruning != True:
            prune_or_not(rooted_tree,min_dupl_tip_overlap,min_dupl_percentage_overlap,output_directory + "unpruned/",tree_name,False,rooted_leaf_cache,id2sp,keep_ss_dup)


def main(args):
    # read tree folder and tree file ending
    tree_folder = check_path(args.homolog_tree_dir,error_if_not_exists=True)

    tree_file_ending = args.tree_file_ending
    keep_ss_dup = args.single_sample_duplications

    # read outgroup list
    if args.outgroup_list:
        outgroup_list = args.outgroup_list.split(",")
    elif args.outgroup_file:
        if not os.path.exists(args.outgroup_file):
            print("Error: The file " + args.outgroup_file + " does not exist.")
            sys.exit(2)
        with open(args.outgroup_file,"r") as f:
            outgroup_list = f.read().splitlines()
    outgroup_list = set(outgroup_list)

    # read min dupl overlap
    if args.min_dupl_tip_overlap:
        min_dupl_tip_overlap = args.min_dupl_tip_overlap
    else:
        min_dupl_tip_overlap = 2
    if args.min_dupl_percentage_overlap:
        min_dupl_percentage_overlap = args.min_dupl_percentage_overlap
    else:
        min_dupl_percentage_overlap = 0.1

    # check pruning options
    if args.just_pruning:
        pruning = True
    elif args.no_pruning:
        pruning = False
    else: # do both pruning and no pruning
        pruning = "both"
    
    # id2sp file
    if args.id2sp_file:
        id2sp = check_path(args.id2sp_file,error_if_not_exists=True,is_folder=False)
    else:
        id2sp = None

    # read output folder
    default_output_folder =  "processed_trees/"
    output_directory = check_path(args.output_directory,default_path=default_output_folder,create_if_not_exists=True)

    if pruning != False:
        if not os.path.exists(output_directory + "pruned/"):
            os.makedirs(output_directory + "pruned/")
    if pruning != True:
        if not os.path.exists(output_directory + "unpruned/"):
            os.makedirs(output_directory + "unpruned/")

    # run the script
    print("------------------------------------------------------------\n")
    print(time.ctime() + "\n")
    print("Processing homolog trees\n")
    start_time = time.time()

    # Filter tree files more efficiently using list comprehension
    tree_files = [f for f in os.listdir(tree_folder) if f.endswith(tree_file_ending)]
    
    all_tree_count = len(tree_files)
    for tree_file in tree_files:
        print("Processing tree: " + tree_file + "\n")
        with open(tree_folder + tree_file,"r") as f:
            tree = t.read_tree_string(f.readline().strip())
        tree_name = tree_file.split(tree_file_ending)[0]
        process_trees(tree,outgroup_list,tree_name,output_directory,min_dupl_tip_overlap,min_dupl_percentage_overlap,pruning,id2sp,keep_ss_dup)
    
    print("------------------------------------------------------------\n")
    print("Output trees write to " + output_directory + "\n")

    process_trees_count = len([f for f in os.listdir(output_directory) if f.endswith(".tre")])

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"Done with orthology inference. Total time elapsed: {elapsed}\n")
    print(f"Total homolog trees read: {all_tree_count}\n")
    print(f"Total homolog trees processed: {process_trees_count}\n")
    if pruning != True:
        unpruned_ortho_count = len([f for f in os.listdir(output_directory + "unpruned/") if "ortho" in f])
        print(f"Total unpruned ortholog trees produced: {unpruned_ortho_count}\n")
    if pruning != False:
        pruned_ortho_count = len([f for f in os.listdir(output_directory + "pruned/") if "ortho" in f])
        print(f"Total pruned ortholog trees produced: {pruned_ortho_count}")

    print("\n------------------------------------------------------------\n")
