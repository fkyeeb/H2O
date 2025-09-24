"""
This is the script for the subcommand infer_ortho
It roots the homolog trees
infers gene duplications
and outputs the ortholog trees
"""

from h2o import tree_reader as t
from h2o.node import Node
from h2o.utils import get_sister, check_folder
import copy
import os
import sys
import time

def list_minus_list(list1, list2):
    """
    This function takes two lists and returns a list of elements in list1 that are not in list2.
    If there are duplicates in list1, only the same number of copies in list2 will be removed in the output.

    :param list1: first list
    :param list2: second list
    :return: list of elements in list1 that are not in list2
    """
    out = list1[:]
    for i in list2:
        if i in list1:
            out.remove(i)
    return out

def check_outgroup_status(filename,tree,outgroup):
    """
    This function checks if the outgroups are in the tree
    and if the outgroup + ingroup both duplicated. 
    If so, skip the tree for now.

    :param filename: name of the tree file
    :param tree: root node of the tree
    :param outgroup: list of outgroups
    :return: True if the outgroup is not duplicated, False otherwise
    :rtype: bool
    """

    outgroup = set(outgroup)
    outgroup_count = {}
    all_tips = tree.lvsnms()

    # record outgroup number
    for tip in all_tips:
        if tip in outgroup:
            if tip not in outgroup_count:
                outgroup_count[tip] = 0
            outgroup_count[tip] += 1

    # check if outgroup exists in the tree
    if len(outgroup_count) == 0:
        print("None of the outgroups is in the tree, " + filename + " is skipped.\n")
        return False
    # elif len(outgroup_count) == 1:
    #     print("Only one outgroup is in the tree, " + filename + " is skipped.\n")
    #     return False
    
    # check if outgroup is duplicated
    limit = len(outgroup) * 2 // 3
    duplicated_outgroup = 0
    for count in outgroup_count.values():
        if count > 1:
            duplicated_outgroup += 1  
    if duplicated_outgroup !=0 and duplicated_outgroup >= limit:
        print("Outgroup is duplicated, " + filename + " is skipped.\n")
        return False
    return True

def locate_outgroup(file_name,tree,outgroup):
    """
    This function takes a tree and a list of outgroups and returns the node of the clade.

    :param file_name: name of the tree file
    :param tree: root node of the tree
    :param outgroup: list of outgroups
    :return: node of the clade and a boolean indicating if the tree needs to be rerooted
    :rtype: Node, bool
    """
    outgroup = set(outgroup)

    # make sure to only save the biggest outgroup clade
    all_tips = tree.lvsnms()
    to_root = True
    for node in tree.iternodes(order="postorder"):
        child_tips = node.lvsnms()
        if node.parent != None:
            other_tips = set(list_minus_list(all_tips,child_tips))
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

def corrent_parent(node):
    """
    After copying a node, this function makes sure that the children of the copied node has the correct parent.
    """
    for child in node.children:
        if child.parent != node:
            child.parent = node

def root_tree(file_name,tree,outgroup):
    """
    This function takes a tree and a list of outgroups and returns the tree rooted at the outgroup.

    :param tree: root node of the tree to be rooted
    :param outgroup: one outgroup for now
    """
    new_root = Node()
    dummy = new_root

    outgroup_node,to_root = locate_outgroup(file_name,tree,outgroup)

    if isinstance(outgroup_node,str):
        # when outgroup is polyphyletic
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
    for node in new_root.iternodes():
        if not node.istip:
            corrent_parent(node)
    return new_root

def prune_selected_taxa(node,selected_taxa):
    """
    This function takes a duplication node and prunes the taxa that are not duplicated.
    """
    node2prune = []
    for n in node.iternodes():
        if n.istip:
            if n.label in selected_taxa:
                node2prune.append(n)

    # make sure no tip is skipped
    for n in node2prune:
        n.prune()

def check_prune_dup(node,min_dup_overlap,bool):
    """
    This function checks if the node is a duplication node.
    And also saves any taxa that are not duplicated.
    """
    child1 = set(node.children[0].lvsnms())
    child2 = set(node.children[1].lvsnms())
    common = child1 & child2
    if len(child1) <= len(child2):
        child2prune = node.children[0]
    else:
        child2prune = node.children[1]

    if len(common) > 0:
        # prune overlap taxa if they are less than the minimum duplication overlap
        if len(common) < min_dup_overlap:
            # print(node)
            # print("Pruned taxa: ", common, "\n")
            prune_selected_taxa(child2prune,common)
            return
        # recognize duplication node if the overlap is big enough
        elif len(common) >= min_dup_overlap:
            node.duplication = True
            node.label = "D"
            if bool:
                node.missing_dup = child1 ^ child2
                # print(node)
                # print(node.missing_dup,"\n")
                prune_selected_taxa(node,node.missing_dup)

def label_duplication_node(tree,min_dup_overlap,bool):
    """
    This function takes a node and checks if it is a duplication node if it is not a tip.
    It will prune the duplications that are smaller than the specified minimum clade size.
    It also saves any taxa that are not duplicated in the duplication node.
    """
    for node in tree.iternodes():
        if not node.istip:
            # pruning happens at the same time as the duplication check
            check_prune_dup(node,min_dup_overlap,bool)
        
def get_orthologs(root):
    """
    This function splits all duplication nodes and into separate ortholog trees
    returns a list of ortholog trees
    """
    ortho_trees = []
    for node in root.iternodes(order="postorder"):
        if node.duplication:
            # split off the smaller clade
            clade1_size = len(set(node.children[0].lvsnms()))
            clade2_size = len(set(node.children[1].lvsnms()))
            if clade1_size >= clade2_size:
                ortho_trees.append(node.children[1])
                node.children[1].prune()
            else:
                ortho_trees.append(node.children[0])
                node.children[0].prune()
    if len(root.children) == 1:
        root = root.children[0]
    ortho_trees.append(root)

    return ortho_trees

def prune_or_not(tree,min_dup_overlap,output_directory,tree_name,bool):

    label_duplication_node(tree,min_dup_overlap,bool)
    
    # remove node when there is only one child left
    if len(tree.children) == 1:
        tree = tree.children[0]
    for n in tree.iternodes():
        if len(n.children) == 1:
            p = n.parent
            if p != None:
                n.children[0].length += n.length
                p.add_child(n.children[0])
                p.remove_child(n)

    # print("pruned tree: ", rooted_tree.get_newick_repr(showbl=True), "\n")
    with open(output_directory + tree_name + "_rooted_pruned.tre","w") as f:
        f.write(tree.get_newick_repr(showbl=True) + ";\n")
    # sys.exit(0)
    # get orthologs
    ortho_trees = get_orthologs(tree)

    # print(len(ortho_trees))
    l = len(ortho_trees)
    for i in range(l-1,-1,-1):
        with open(output_directory + tree_name + "_ortho" + str(l-i) + ".tre","w") as f:
            string = ortho_trees[i].get_newick_repr(showbl=True) + ";\n"
            # print("ortholog tree ",l-i," : ",string, "\n")
            f.write(string)

def process_trees(tree,outgroup_list,tree_name,output_directory,min_dup_overlap):
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

    if check_outgroup_status(tree_name,tree,outgroup_list):
        rooted_tree = root_tree(tree_name,tree,outgroup_list)
    else:
        # if less than 2 ourgroups or if outgroup is duplicated, skip the tree
        rooted_tree = None

    if rooted_tree:
        # print("rooted tree",rooted_tree.get_newick_repr(), "\n")
        with open(output_directory + tree_name + "_rooted.tre","w") as f:
            f.write(rooted_tree.get_newick_repr(showbl=True) + ";\n")

        # label duplication node and prune missing taxa
        rooted_tree_2prune = copy.deepcopy(rooted_tree)

        prune_or_not(rooted_tree_2prune,min_dup_overlap,output_directory + "pruned/",tree_name,True)
        prune_or_not(rooted_tree,min_dup_overlap,output_directory + "no_pruning/",tree_name,False)


def main(args):
    # read tree folder and tree file ending
    tree_folder = args.homolog_tree_dir
    tree_folder = check_folder(tree_folder)
    if not os.path.exists(tree_folder):
        print("Error: The folder " + tree_folder + " does not exist.")
        sys.exit(2)
    tree_file_ending = args.tree_file_ending

    # read outgroup list
    if hasattr(args, "outgroup_list"):
        outgroup_list = args.outgroup_list.split(",")
    elif hasattr(args, "outgroup_file"):
        with open(args.outgroup_file,"r") as f:
            outgroup_list = f.read().splitlines()

    # read min ingroup taxa
    if args.min_ingroup_taxa:
        min_dup_overlap = args.min_ingroup_taxa
    else:
        min_dup_overlap = 3

    # read output folder
    output_directory = args.output_directory
    output_directory = check_folder(output_directory)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if not os.path.exists(output_directory + "pruned/"):
        os.makedirs(output_directory + "pruned/")
    if not os.path.exists(output_directory + "no_pruning/"):
        os.makedirs(output_directory + "no_pruning/")

    # run the script
    print("------------------------------------------------------------\n\n")
    print(time.ctime() + "\n")
    print("Processing homolog trees\n")

    for tree_file in os.listdir(tree_folder):
        if tree_file.endswith(tree_file_ending):
            print("Processing tree: " + tree_file + "\n")
            with open(tree_folder + tree_file,"r") as f:
                tree = t.read_tree_string(f.readline().strip())
            tree_name = tree_file[:-len(tree_file_ending)]
            process_trees(tree,outgroup_list,tree_name,output_directory,min_dup_overlap)
    
    print("Output trees are saved in " + output_directory + "\n")
    print("Done with orthology inference at " + time.ctime())
    print("\n------------------------------------------------------------\n\n")
