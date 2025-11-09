from h2o import tree_reader as t
from h2o.utils import check_path
import random as r

def randomly_select_tip_in_clade(tree,nodes2keep):
    tips2keep = []
    for node in tree.iternodes():
        if node.label in nodes2keep:
            if node.istip:
                tips2keep.append(node.label)
            else:
                tip = r.choice(node.lvsnms())
                tips2keep.append(tip)
    return tips2keep

def extract_constraint_tree_by_tips(tree,tips2keep):
    for node in tree.iternodes():
        if node.istip:
            if node.label not in tips2keep:
                node.prune()
    
    for n in tree.iternodes():
        if len(n.children) == 1:
            p = n.parent
            if p != None:
                n.children[0].length += n.length
                p.add_child(n.children[0])
                p.remove_child(n)
            

def main(args):
    consensus_tree_file = check_path(args.consensus_tree_file,is_folder=False,error_if_not_exists=True)
    output_directory = check_path(args.output_directory,default_path="./",create_if_not_exists=True)
    
    with open(consensus_tree_file,"r") as f:
        consensus_tree = t.read_tree_string(f.readline().strip())
    
    # option 1: enter nodes
    if args.nodes:
        nodes2keep = args.nodes.split(",")
        tips2keep = randomly_select_tip_in_clade(consensus_tree,nodes2keep)
    elif args.tips_file:
        with open(args.tips_file,"r") as f:
            tips2keep = f.readlines()
            tips2keep = [line.strip() for line in tips2keep]
    
    extract_constraint_tree_by_tips(consensus_tree,tips2keep)

    with open(output_directory + "constraint_tree.tre","w") as f:
        f.write(consensus_tree.get_newick_repr(showbl=True) + ";\n")