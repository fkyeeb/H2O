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

def check_folder(folder,default_folder=None,error_if_not_exists=False,create_if_not_exists=False):
    """
    check if folder has a trailing slash
    
    :param folder: folder path
    :param error_if_not_exists: error if folder does not exist
    :param create_if_not_exists: create folder if it does not exist
    :param default_folder: default folder if folder is not provided
    """
    
    if folder == None:
        folder = default_folder
    
    if folder[-1] != "/":
        folder = folder + "/"
    
    if error_if_not_exists:
        if not os.path.exists(folder):
            print("Error: The folder " + folder + " does not exist.")
            sys.exit(2)
    
    if create_if_not_exists:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    return folder

def precompute_leaf_names_number_nodes(tree,label=False):
    """
    Pre-compute leaf names for all nodes and number all nodes
    node number is stored in node.cache_label

    :param [Node] tree: root node of the tree
    :param bool label: whether to put numbers as node.label
    :return [dict]: dictionary with node numbers as keys and tips [list] as values
    """
    leaf_cache = {}
    num = 0
    for node in tree.iternodes():
        if not node.istip:
            node.cache_label = str(num)
            if label:
                node.label = node.cache_label
            leaf_cache[node.cache_label] = set(node.lvsnms())
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

    subprocess.run([os.environ["SHELL"], "-i", "-c",cmd], check=True)

    end_time = time.time()
    elapsed = transform_elapsed_time(start_time,end_time)
    print(f"\nCommand finished in {elapsed}\n")
    print("------------------------------------------------------------\n")