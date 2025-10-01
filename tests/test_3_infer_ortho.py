from h2o import infer_orthology
import argparse
import shutil
import os

def test_infer_ortho():
    """Test the infer_ortho function, see if it runs correctly"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees",
        outgroup_list="o",
        tree_file_ending=".tre",
        min_ingroup_taxa=3,
        output_directory=None,
        no_pruning=False,
        just_pruning=False
    )
    infer_orthology.main(args)
    
    with open("tests/test_data/processed_trees/unpruned/dup_ortho1.tre", "r") as f:
        assert f.read() == "(o:0.0,((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees/dup_difficult_rooted.tre", "r") as f:
        assert f.read() == "(o:0.0,(a:0.0,((a:0.0,b:0.0):0.0,((c:0.0,d:0.0):0.0,(((c:0.0,d:0.0):0.0,e:0.0):0.0,e:0.0):0.0):0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees/unpruned/dup_difficult_ortho1.tre", "r") as f:
        assert f.read() == "(o:0.0,((a:0.0,b:0.0):0.0,((c:0.0,d:0.0):0.0,e:0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees/pruned/dup_bl_rooted_processed.tre", "r") as f:
        assert f.read() == "((((((a:1.0,b:1.0):1.0,c:2.0):1.0,((a:1.0,b:1.0):1.0,c:2.0):1.0)D:1.0,(((a:1.0,b:1.0):1.0,d:2.0):1.0,((a:1.0,b:1.0):1.0,d:2.0):1.0)D:1.0)D:1.0,(e:3.0,f:3.0):2.0):1.0,o:6.0):0.0;\n"   
    with open("tests/test_data/processed_trees/unpruned/dup_bl_ortho1.tre", "r") as f:
        assert f.read() == "(((e:3.0,f:3.0):2.0,((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0):3.0):1.0,o:6.0):0.0;\n"
    
    # remove all the output files created by unit test
    # shutil.rmtree('tests/test_data/processed_trees')

def test_infer_ortho_rooting(capsys):
    """Test the infer_ortho function, see if it roots correctly"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees_rooting",
        outgroup_file="tests/test_data/outgroups.txt",
        tree_file_ending=".tre",
        min_ingroup_taxa=3,
        output_directory="tests/test_data/processed_trees_rooting",
        no_pruning=False,
        just_pruning=False
    )
    infer_orthology.main(args)

    with open("tests/test_data/processed_trees_rooting/outgroup_rooted.tre", "r") as f:
        assert f.read() == "(((x:0.0,y:0.0):0.0,z:0.0):0.0,(d:0.0,(c:0.0,((a:0.0,b:0.0):0.0,e:0.0):0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees_rooting/pruned/outgroup2_ortho1.tre", "r") as f:
        assert f.read() == "(((x:0.0,y:0.0):0.0,z:0.0):0.0,((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees_rooting/unpruned/outgroup2_ortho1.tre", "r") as f:
        assert f.read() == "(((x:0.0,y:0.0):0.0,z:0.0):0.0,(((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0,(e:0.0,f:0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees_rooting/outgroup_para_rooted.tre", "r") as f:
        assert f.read() == "((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0;\n"
    assert not os.path.exists("tests/test_data/processed_trees_rooting/outgroup_poly_rooted.tre")

    captured = capsys.readouterr()
    assert "Outgroup is polyphyletic, outgroup_poly is skipped.\n" in captured.out
    assert "None of the outgroups is in the tree, no_outgroup is skipped.\n" in captured.out
    
    # remove all the output files created by unit test
    shutil.rmtree('tests/test_data/processed_trees_rooting')