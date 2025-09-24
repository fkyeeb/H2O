"""
Utility functions for H2O
"""

def get_sister(node):
    """
    This function takes a node and returns its sister nodes.
    """
    for child in node.parent.children:
        if child != node:
            return child

def check_folder(arg):
    """
    check if folder has a trailing slash
    
    :param arg: folder path
    :return: folder path with trailing slash
    """
    if arg[-1] != "/":
        return arg + "/"
    else:
        return arg