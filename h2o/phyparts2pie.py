'''
extract phyparts conflict result for ploting
'''

from h2o import tree_reader as t
from h2o.utils import check_path,precompute_leaf_names_number_nodes

def process_hist(hist_file,total_num_trees):
    conflict_data = {}
    with open(hist_file,"r") as file:
        for line in file:
            splt = line.strip().split(",")
            node_number = splt.pop(0).replace("Node","")
            total = int(splt.pop())
            splt = [float(i) for i in splt]
            conflict = max(splt)
            others = sum(splt) - conflict
            concord = total - conflict - others
            uninformative = total_num_trees - total
            conflict_data[node_number] = [str(concord),(str(conflict)),str(others),str(uninformative)]
    return conflict_data

def process_node_key(node_key_file):
    node_key = {}
    with open(node_key_file,"r") as file:
        for line in file:
            splt = line.strip().split(" ")
            node_number = splt[0]
            node_tree = t.read_tree_string(splt[1] + ";")
            node_key[node_number] = node_tree.lvsnms()
    return node_key

def match_node_with_data(tree,conflict_data,node_key):
    pies = {}
    for node_number in conflict_data:
        nums = conflict_data[node_number]
        pie = "[&pie=" + nums[0] + "," + nums[1] + "," + nums[2] + "," + nums[3] + "]"
        pies[node_number] = pie

    for node in tree.iternodes():
        if not node.istip:
            lvs = node.lvsnms()
            for node_number in node_key:
                if lvs == node_key[node_number]:
                    node.label = pies[node_number]

def main(args):
    node_key_file = check_path(args.phyparts_node_key_file, is_folder=False, error_if_not_exists=True)
    hist_file = check_path(args.phyparts_hist_file, is_folder=False, error_if_not_exists=True)
    total_num_trees = int(args.total_tree_number)
    consensus_tree_file = check_path(args.consensus_tree_file, is_folder=False, error_if_not_exists=True)
    output_directory = check_path(args.output_directory,default_path="./",create_if_not_exists=True)

    with open(consensus_tree_file,"r") as file:
        consensus_tree = t.read_tree_string(file.readline().strip())

    conflict_data = process_hist(hist_file,total_num_trees)
    node_key = process_node_key(node_key_file)

    with open(output_directory + "phyparts_summary.tsv","w") as file:
        file.write("node_number\tconcord\tmain_conflict\tother_conflict\tunimformative\n")
        for node_number in conflict_data:
            nums = conflict_data[node_number]
            file.write(node_number + "\t" + nums[0] + "\t" + nums[1] + "\t" + nums[2] + "\t" + nums[3] + "\n")

    match_node_with_data(consensus_tree,conflict_data,node_key)
    with open(output_directory + "gokstad_pie.tre","w") as out:
        out.write(consensus_tree.get_newick_repr(showbl=True) + ";\n")