from h2o import utils
from h2o import tree_reader as t

def test_check_folder():
    """Test the check_folder function, see if it only adds a trailing slash when there is none"""

    assert utils.check_folder("tests/test_data/") == "tests/test_data/"

def test_check_folder2():
    """Test the check_folder function, see if it only adds a trailing slash when there is none"""

    assert utils.check_folder("tests/test_data") == "tests/test_data/"

def test_get_sister():
    """Test the get_sister function, see if it returns the sister node"""

    tree = t.read_tree_string("((a,b),c);")
    assert utils.get_sister(tree.children[0]) == tree.children[1]
