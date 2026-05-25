"""
This is the script for the subcommand map_gene_loss
it maps gene losses to species tree
and output the results in a tsv  and a tree file
"""

from h2o.utils import (
    check_path,
    precompute_leaf_names_number_nodes,
    transform_elapsed_time,
    get_deepest_dup_parent,
    get_deepest_non_dup_parent,
    get_sister
)
from h2o import tree_reader as t
import sys
import time

def get_missing_tips(dup_tree,leaf_cache,gfe_tips_dict,gfe_nodes,ils,gfe_sis_dict):
    """
    if node is a dup node, if deepest_dup_parent is gfe, get missing tips of both child nodes compared to itself
    this is counting gene loss following parsimony, in case missing tips are counted twice

    :param [Node] dup_tree: root node of the dup tree
    :param [dict] leaf_cache: dictionary with node numbers as keys and set(tips) as values
    :param [dict] gfe_tips_dict: dictionary with gfe node numbers as keys and tips [set(ingroup_tips),set(outgroup_tips)] as values
    :param [list] gfe_nodes: list of gfe node numbers in this tree
    :param [bool] ils: flag indicating whether ILS corrected tree is used
    :param [dict] gfe_sis_dict: dictionary with gfe node numbers as keys and set(tips) as values for their sister nodes
    """
    all_dup_tips = leaf_cache["0"]
    missing_tips = []
    gfe_node_output = []

    for node in dup_tree.iternodes():
        if node.label == "D":
            deepest_dup_parent = get_deepest_dup_parent(node)
            deepest_non_dup_parent = get_deepest_non_dup_parent(deepest_dup_parent)
            dup_tips = set(leaf_cache[deepest_dup_parent.cache_label])
            deepest_non_dup_parent_tips = set(leaf_cache[deepest_non_dup_parent.cache_label])
            other_dup_tips = deepest_non_dup_parent_tips - dup_tips

            for gfe_node in gfe_nodes:
                gfe_tips = gfe_tips_dict[gfe_node][0]
                gfe_other_tips = gfe_tips_dict[gfe_node][1]
                if dup_tips <= gfe_tips and other_dup_tips <= gfe_other_tips:
                    for child in node.children:
                        child_tips = set(leaf_cache[child.cache_label])
                        if len(child_tips) >= 2:
                            node_missing_tips = gfe_tips - set(leaf_cache[child.cache_label])
                            if len(node_missing_tips) > 0:
                                missing_tips.append(node_missing_tips)
                                gfe_node_output.append(gfe_node)
                elif ils:
                    sis_tips = gfe_sis_dict[gfe_node]
                    dup_ils_tips = dup_tips - sis_tips
                    if dup_ils_tips <= gfe_tips and other_dup_tips <= gfe_other_tips:
                        for child in node.children:
                            child_ils_tips = set(leaf_cache[child.cache_label]) - sis_tips
                            if len(child_ils_tips) >= 2:
                                node_missing_tips = gfe_tips - child_ils_tips
                                if len(node_missing_tips) > 0:
                                    missing_tips.append(node_missing_tips)
                                    gfe_node_output.append(gfe_node)

    return missing_tips,gfe_node_output

def map_gene_loss(missing_tips,sp_tree_leaf_cache):
    """
    Mapping gene loss to species tree nodes

    :param [list] missing_tips: list of set(missing tips) for each gfe dup node
    :param [dict] sp_tree_leaf_cache: dictionary with node numbers as keys and set(tips) as values for the species tree
    :return [list]: list of list of gene loss node numbers for each gfe dup node
    :return [list]: list of set(missing tips) for each gfe dup node in the same order as gene loss node numbers
    """
    missing_nodes = []
    for one_clade_missing_tips in missing_tips:
        possible_nodes = []
        for node_num in sp_tree_leaf_cache:
            if sp_tree_leaf_cache[node_num] <= one_clade_missing_tips:
                possible_nodes.append([node_num,sp_tree_leaf_cache[node_num]])
        possible_nodes.sort(key=lambda x: len(x[1]),reverse=True)

        one_clade_missing_nodes = []
        for node_info in possible_nodes:
            node_num,node_tips = node_info
            if node_tips <= one_clade_missing_tips:
                one_clade_missing_tips -= node_tips
                one_clade_missing_nodes.append(node_num)
        missing_nodes.append(one_clade_missing_nodes)

    return missing_nodes,missing_tips

def main(args):
    """
    Main function to map gene losses
    """

    processed_tree_folder = check_path(args.processed_tree_dir) + "unpruned/"
    processed_tree_folder = check_path(processed_tree_folder,error_if_not_exists=True)

    default_output_folder = "other_output/"
    output_folder = check_path(args.output_directory,default_path=default_output_folder,create_if_not_exists=True)
    
    gfe_nodes = args.gfe_nodes.split(",")
    try:
        for node in gfe_nodes:
            int(node)
    except ValueError:
        print("Error: gfe node numbers must be integers.")
        sys.exit(2)

    ils = args.ils_correction
    default_dup_dir = "other_output/"
    duplication_counts_dir = check_path(args.duplication_counts_dir,default_path=default_dup_dir,error_if_not_exists=True)

    if ils:
        duplication_counts_file = check_path(duplication_counts_dir + "duplication_counts_ils_corrected.tsv",is_folder=False,error_if_not_exists=True)
    else:
        duplication_counts_file = check_path(duplication_counts_dir + "duplication_counts.tsv",is_folder=False,error_if_not_exists=True)
    numbered_tree_file = check_path(duplication_counts_dir + "consensus_tree_numbered.tre",is_folder=False,error_if_not_exists=True)

    if args.id2sp_file:
        id2sp = check_path(args.id2sp_file,error_if_not_exists=True,is_folder=False)
    else:
        id2sp = None

    print("------------------------------------------------------------\n")
    print(time.ctime() + "\n")
    start_time = time.time()
    
    with open(numbered_tree_file,"r") as f:
        numbered_tree = t.read_tree_string(f.readline().strip())
    
    gene_loss_node_counts = {}
    sp_tree_leaf_cache = precompute_leaf_names_number_nodes(numbered_tree,use_label=True,return_set=True)
    all_sp_tree_tips = sp_tree_leaf_cache["0"]
    gfe_tips_dict = {}
    gfe_sis_dict = {}
    for gfe_node in gfe_nodes:
        tips = sp_tree_leaf_cache[gfe_node]
        gfe_tips_dict[gfe_node] = [tips,all_sp_tree_tips - tips]
    if ils:
        for node in numbered_tree.iternodes():
            if node.label in gfe_nodes:
                sister = get_sister(node)
                if sister.istip:
                    gfe_sis_dict[node.label] = set([sister.label])
                else:
                    gfe_sis_dict[node.label] = set(sp_tree_leaf_cache[sister.label])


    for node in sp_tree_leaf_cache:
        gene_loss_node_counts[node] = 0
    
    gene_loss_per_tree = {}
    
    with open(duplication_counts_file,"r") as f:
        f.readline()
        if ils:
            f.readline() # use the third ils corrected tree
        for line in f:
            splt = line.strip().split("\t")
            tree_name = splt.pop(0)

            if any(int(node) >= len(splt) or int(node) < 0 for node in gfe_nodes):
                print("Error: gfe node numbers provided is out of range for the consensus tree.")
                sys.exit(2)
            
            gfe_node_in_this_tree = [node for node in gfe_nodes if int(splt[int(node)]) > 0]
            if len(gfe_node_in_this_tree) > 0:
                gene_loss_per_tree[tree_name] = {}
                gene_loss_per_tree[tree_name]["gfe_node_in_this_tree"] = gfe_node_in_this_tree
                gene_loss_per_tree[tree_name]["gfe_node_output"] = []
                gene_loss_per_tree[tree_name]["gene_loss_nodes"] = []
                gene_loss_per_tree[tree_name]["gene_loss_tips"] = []

    if id2sp:
        tree_file_ending = "_rooted_processed_id.tre"
    else:
        tree_file_ending = "_rooted_processed.tre"
    for tree_name in gene_loss_per_tree:
        with open(processed_tree_folder + tree_name + tree_file_ending,"r") as f:
            tree = t.read_tree_string(f.readline().strip())
        
        leaf_cache = precompute_leaf_names_number_nodes(tree,id2sp=id2sp)
        gfe_node_in_this_tree = gene_loss_per_tree[tree_name]["gfe_node_in_this_tree"]
        missing_tips,gfe_node_output = get_missing_tips(tree,leaf_cache,gfe_tips_dict,gfe_node_in_this_tree,ils,gfe_sis_dict)
        gene_loss_nodes,gene_loss_tips = map_gene_loss(missing_tips,sp_tree_leaf_cache)
        gene_loss_per_tree[tree_name]["gene_loss_nodes"] = gene_loss_nodes
        gene_loss_per_tree[tree_name]["gene_loss_tips"] = gene_loss_tips
        gene_loss_per_tree[tree_name]["gfe_node_output"] = gfe_node_output

        for l in gene_loss_nodes:
            for node in l:
                gene_loss_node_counts[node] += 1
    
    with open(output_folder + "gene_loss_counts.tsv","w") as f:
        f.write("tree\tdupl_clade_count\tgfe_node\tgene_loss_nodes\tgene_loss_tips\n")
        for tree_name in gene_loss_per_tree:
            data = gene_loss_per_tree[tree_name]
            for i in range(len(data["gfe_node_output"])):
                f.write(
                    tree_name + "\t" +
                    str(i+1) + "\t" +
                    str(data["gfe_node_output"][i]) + "\t" +
                    ",".join(data["gene_loss_nodes"][i]) + "\t" +
                    ",".join(data["gene_loss_tips"][i]) + "\n"
                )
    
    with open(output_folder + "gene_loss_counts.tre","w") as f:
        f.write(numbered_tree.get_newick_repr(showbl=True) + ";\n")
        
        for node in numbered_tree.iternodes():
            if node.label in gene_loss_node_counts:
                node.label = str(gene_loss_node_counts[node.label])
        f.write(numbered_tree.get_newick_repr(showbl=True) + ";\n")
    
    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"Done with gene loss mapping. Total time elapsed: {elapsed}")

    print("\n------------------------------------------------------------\n\n")
           

        