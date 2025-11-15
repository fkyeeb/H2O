"""
Module for removing specified tips from trees in a directory.
Or removing all tips except specified ones.
"""
from h2o import tree_reader as t
from h2o.utils import check_path
import os

def main(args):
    tree_dir = check_path(args.tree_dir, error_if_not_exists=True)
    output_directory = check_path(args.output_directory, create_if_not_exists=True)
    tree_file_ending = args.tree_file_ending
    if not args.minimum_taxa:
        minimum_taxa = 5
    else:
        minimum_taxa = int(args.minimum_taxa)

    if args.tips2remove:
        tips_file = check_path(args.tips2remove, is_folder=False, error_if_not_exists=True)
        with open(tips_file, "r") as file:
            tips = set(line.strip() for line in file)
            remove = True
    else:
        tips_file = check_path(args.tips2save, is_folder=False, error_if_not_exists=True)
        with open(tips_file, "r") as file:
            tips = set(line.strip() for line in file)
            remove = False
    
    # run the script
    print("------------------------------------------------------------\n")

    # Filter tree files more efficiently using list comprehension
    tree_files = [f for f in os.listdir(tree_dir) if f.endswith(tree_file_ending)]
    
    all_tree_count = len(tree_files)
    for tree_file in tree_files:
        print("Processing tree: " + tree_file + "\n")
        with open(tree_dir + tree_file,"r") as f:
            tree = t.read_tree_string(f.readline().strip())
        
        for node in tree.iternodes():
            if node.istip:
                if remove:
                    if node.label in tips:
                        node.prune()
                else:
                    if node.label not in tips:
                        node.prune()
        
        if len(tree.lvsnms()) < minimum_taxa:
            print("Skipping tree " + tree_file + " because it has less than " + str(minimum_taxa) + " taxa after tip removal.\n")
            continue

        for n in tree.iternodes(order="postorder"):
            if len(n.children) == 1:
                p = n.parent
                if p != None:
                    n.children[0].length += n.length
                    p.add_child(n.children[0])
                    p.remove_child(n)
                else:
                    tree = n.children[0]
        
        with open(output_directory + tree_file,"w") as out:
            out.write(tree.get_newick_repr(showbl=True) + ";\n")

    print("------------------------------------------------------------\n")
    print("Output trees are saved in " + output_directory + "\n")