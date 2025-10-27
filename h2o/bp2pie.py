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

def main(args):
    bp_output_file = check_path(args.bp_output_file,is_folder=False,error_if_not_exists=True)
    output_directory = check_path(args.output_directory,default_path="./",create_if_not_exists=True)

    if args.summary_tree_file:
        summary_tree_file = check_path(args.summary_tree_file,is_folder=False,error_if_not_exists=True)
        with open(summary_tree_file,"r") as file:
            summary_tree = t.read_tree_string(file.readline().strip())

    with open(bp_output_file,"r") as file, open(output_directory + "bp_output.tre","w") as out:
        for i in file:
            if "TREES WITH CONFLICT" in i:
                conflict_line = file.readline().strip()
                conflict_tree = process_bp_tree(conflict_line,out)

                concord_line = file.readline().strip()
                concord_tree = process_bp_tree(concord_line,out)
                
                unsup_line = file.readline().strip()
                unsup_tree = process_bp_tree(unsup_line,out)
    
    bp_leaf_cache = precompute_leaf_names_number_nodes(conflict_tree)
    precompute_leaf_names_number_nodes(concord_tree)
    precompute_leaf_names_number_nodes(unsup_tree)

    data = {}
    for node in conflict_tree.iternodes():
        if not node.istip:
            data[node.cache_label] = [node.label,0,0]
    for node in concord_tree.iternodes():
        if not node.istip:
            data[node.cache_label][1] = node.label
    for node in unsup_tree.iternodes():
        if not node.istip:
            data[node.cache_label][2] = node.label
    data = dict(sorted(data.items(), key=lambda x: int(x[0])))
    del data["0"]

    with open(output_directory + "bp_data.tsv","w") as file:
        file.write("node_number\tconflict\tconcord\tunsupported\n")
        for node_number in data:
            nums = data[node_number]
            file.write(node_number + "\t" + nums[0] + "\t" + nums[1] + "\t" + nums[2] + "\n")
    
    if args.summary_tree_file:
        for node in summary_tree.iternodes():
            if not node.istip:
                leaves = set(node.lvsnms())
                for node_number in bp_leaf_cache:
                    if bp_leaf_cache[node_number] == leaves:
                        if node_number != "0":
                            node.label = node_number
                        break
        output_tree = summary_tree
    else:
        for node in conflict_tree.iternodes():
            if not node.istip:
                if node.cache_label != "0":
                    node.label = node.cache_label
        output_tree = conflict_tree
    with open(output_directory + "bp_summary_tree_numbered.tre","w") as out:
        out.write(output_tree.get_newick_repr(showbl=True) + ";\n")

    pies = {}
    for node_number in data:
        nums = data[node_number]
        if args.pie_option:
            pie = "[&pie=" + nums[0] +","+ nums[1] +","+ nums[2] + "]"
        else:
            if int(nums[0]) + int(nums[1]) == 0:
                pie = ""
            else:
                pie = "[&pie=" + nums[0] +","+ nums[1] + "]"
        pies[node_number] = pie

    with open(output_directory + "gokstad_pie.tre","w") as out:
        for node in output_tree.iternodes():
            if not node.istip:
                if node.label in pies:
                    node.label = pies[node.label]
        out.write(output_tree.get_newick_repr(showbl=True) + ";\n")