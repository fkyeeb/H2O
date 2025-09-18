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