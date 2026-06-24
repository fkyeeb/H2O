'''
extract bp conflict result for ploting
'''

from h2o import tree_reader as t
from h2o.utils import check_path,precompute_leaf_names_number_nodes

def process_bp_tree(tree_line,output_file):
    tree = t.read_tree_string(tree_line)
    output_file.write(tree_line + "\n")
    for node in tree.iternodes():
        if not node.istip:
            if node.label == "":
                node.label = "0"
    return tree

def process_bp_file(bp_output_file,out,data,first=True):
    with open(bp_output_file,"r") as file:
        for i in file:
            if "TREES WITH CONFLICT" in i:
                conflict_line = file.readline().strip()
                conflict_tree = process_bp_tree(conflict_line,out)

                concord_line = file.readline().strip()
                concord_tree = process_bp_tree(concord_line,out)
                
                unsup_line = file.readline().strip()
                unsup_tree = process_bp_tree(unsup_line,out)
    
    bp_leaf_cache = precompute_leaf_names_number_nodes(conflict_tree,return_set=True)
    precompute_leaf_names_number_nodes(concord_tree,return_set=True)
    precompute_leaf_names_number_nodes(unsup_tree,return_set=True)

    for node in conflict_tree.iternodes():
        if not node.istip:
            if first:
                data[node.cache_label] = [int(node.label),0,0]
            else:
                data[node.cache_label][0] += int(node.label)
    for node in concord_tree.iternodes():
        if not node.istip:
            data[node.cache_label][1] += int(node.label)
    for node in unsup_tree.iternodes():
        if not node.istip:
            data[node.cache_label][2] += int(node.label)
    
    return data,bp_leaf_cache,conflict_tree

def main(args):
    bp_output_files = args.bp_output_file.split(",")
    for bp_output_file in bp_output_files:
        bp_output_file = check_path(bp_output_file,is_folder=False,error_if_not_exists=True)
    output_directory = check_path(args.output_directory,default_path="./",create_if_not_exists=True)

    if args.consensus_tree_file:
        consensus_tree_file = check_path(args.consensus_tree_file,is_folder=False,error_if_not_exists=True)
        with open(consensus_tree_file,"r") as file:
            consensus_tree = t.read_tree_string(file.readline().strip())
    if args.run_name:
        run_name = "_" + args.run_name
    else:
        run_name = ""

    with open(output_directory + "bp_output" + run_name + ".tre","w") as out:
        data = {}
        for index,bp_output_file in enumerate(bp_output_files):
            data,bp_leaf_cache,tree = process_bp_file(bp_output_file,out,data,first=(index==0))
    data = dict(sorted(data.items(), key=lambda x: int(x[0])))
    del data["0"]

    with open(output_directory + "bp_data" + run_name + ".tsv","w") as file:
        file.write("node_number\tconflict\tconcord\tuninformative\n")
        for node_number in data:
            nums = data[node_number]
            file.write(node_number + "\t" + str(nums[0]) + "\t" + str(nums[1]) + "\t" + str(nums[2]) + "\n")
    
    if args.consensus_tree_file:
        for node in consensus_tree.iternodes():
            if not node.istip:
                leaves = set(node.lvsnms())
                for node_number in bp_leaf_cache:
                    if bp_leaf_cache[node_number] == leaves:
                        if node_number != "0":
                            node.label = node_number
                        break
        output_tree = consensus_tree
    else:
        for node in tree.iternodes():
            if not node.istip:
                if node.cache_label != "0":
                    node.label = node.cache_label
        output_tree = tree
    with open(output_directory + "bp_consensus_tree_numbered" + run_name + ".tre","w") as out:
        out.write(output_tree.get_newick_repr(showbl=True) + ";\n")

    pies = {}
    for node_number in data:
        nums = data[node_number]
        if args.pie_option:
            pie = "[&pie=" + str(nums[0]) + "," + str(nums[1]) + "," + str(nums[2]) + "]"
        else:
            if nums[0] + nums[1] == 0:
                pie = ""
            else:
                pie = "[&pie=" + str(nums[0]) + "," + str(nums[1]) + "]"
        pies[node_number] = pie

    with open(output_directory + "gokstad_pie" + run_name + ".tre","w") as out:
        for node in output_tree.iternodes():
            if not node.istip:
                if node.label in pies:
                    node.label = pies[node.label]
        out.write(output_tree.get_newick_repr(showbl=True) + ";\n")