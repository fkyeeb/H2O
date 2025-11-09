"""
Utility functions for H2O
"""
import subprocess
import time
import sys
import os

def get_sister(node):
    """
    This function takes a node and returns its sister nodes.
    """
    for child in node.parent.children:
        if child != node:
            return child

def check_path(path,is_folder=True,default_path=None,error_if_not_exists=False,create_if_not_exists=False):
    """
    check given path
    
    :param path: path
    :param is_folder: whether the path is a folder
    :param error_if_not_exists: error if folder does not exist
    :param create_if_not_exists: create folder if it does not exist
    :param default_path: default path if path is not provided
    """
    
    if path == None:
        path = default_path
    
    if is_folder:
        if path[-1] != "/":
            path = path + "/"
    
    if error_if_not_exists:
        if not os.path.exists(path):
            print("Error: " + path + " does not exist.")
            sys.exit(2)
    
    if create_if_not_exists:
        if not os.path.exists(path):
            os.makedirs(path)
    
    return path

def precompute_leaf_names_number_nodes(tree,use_label=False,label=False,return_set=True):
    """
    Pre-compute leaf names for all nodes and number all nodes
    node number is stored in node.cache_label

    :param [Node] tree: root node of the tree
    :param bool label: whether to put numbers as node.label
    :return [dict]: dictionary with node numbers as keys and tips [list] or set() as values
    """
    leaf_cache = {}
    num = 0
    for node in tree.iternodes():
        if not node.istip:
            if not use_label:
                node.cache_label = str(num)
                if label:
                    node.label = node.cache_label
            else:
                node.cache_label = node.label
            if return_set:
                leaf_cache[node.cache_label] = set(node.lvsnms())
            else:
                leaf_cache[node.cache_label] = node.lvsnms()
            num += 1
    return leaf_cache

def transform_elapsed_time(start_time,end_time):
    """
    Record the time elapsed
    """
    elapsed = end_time - start_time
    if elapsed < 60:
        time_str = f"{elapsed:.2f} seconds"
    else:
        minutes = int(elapsed // 60)
        seconds = elapsed % 60
        time_str = f"{minutes} minutes {seconds:.2f} seconds"
    
    return time_str

def run_shell_command(cmd):
    """
    run a bash command
    
    :param [str] cmd: bash command
    """
    start_time = time.time()
    print("\n------------------------------------------------------------")
    print("\nRunning command: " + cmd + "\n")

    subprocess.run([os.environ["SHELL"], "-c",cmd], check=True)

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"\nCommand finished in {elapsed}\n")
    print("------------------------------------------------------------\n")

def get_deepest_dup_parent(node):
    """
    This function takes a node and returns its deepest duplication parent.
    If node has no parent, return node
    """
    current_node = node

    while current_node.parent != None:
        if current_node.parent.label != "D":
            return current_node
        else:
            current_node = current_node.parent

    return current_node