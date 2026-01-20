from h2o import utils
from h2o import tree_reader as t
import pytest
import os
import shutil

def test_check_path():
    """Test the check_path function, see if it only adds a trailing slash when there is none"""

    assert utils.check_path("tests/test_data/") == "tests/test_data/"
    assert utils.check_path("tests/test_data") == "tests/test_data/"
    assert utils.check_path(None,default_path="tests/test_data") == "tests/test_data/"

def test_check_path_error(capsys):
    """Test the check_path function, see if it only adds a trailing slash when there is none"""

    with pytest.raises(SystemExit) as exc_info:
        utils.check_path("tests/test_error",error_if_not_exists=True)
    
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert "Error: tests/test_error/ does not exist." in captured.out

def test_check_path_create():
    """Test the check_path function, see if it creates a folder if it does not exist"""
    
    utils.check_path("tests/test_create",create_if_not_exists=True)
    assert os.path.exists("tests/test_create")

    shutil.rmtree("tests/test_create")

def test_get_sister():
    """Test the get_sister function, see if it returns the sister node"""

    tree = t.read_tree_string("((a,b),c);")
    assert utils.get_sister(tree.children[0]) == tree.children[1]

def test_precompute_leaf_names_number_nodes():
    """Test the precompute_leaf_names_number_nodes function"""

    tree = t.read_tree_string("((a,b),c);")
    leaf_cache = utils.precompute_leaf_names_number_nodes(tree,return_set=True)

    assert leaf_cache == {"0": {"a", "b", "c"}, '1': {'a', 'b'}}

def test_run_shell_command(capsys):
    """Test the run_shell_command function"""

    utils.run_shell_command("echo 'Hello from H2O test'")

    captured = capsys.readouterr()
    assert "Hello from H2O test" in captured.out

def test_get_tips_in_ascending_order():
    """Test the get_tips_in_ascending_order function"""

    tree = t.read_tree_string("(o,(((a,b),(c,(d,e))),(f,g)));")
    tips = utils.get_tips_in_ascending_order(tree)

    assert tips == ['o', 'f', 'g', 'a', 'b', 'c', 'd', 'e']