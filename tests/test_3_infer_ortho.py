from h2o import infer_orthology, cli
import argparse
import shutil
import os
import sys

def test_parse_arguments(monkeypatch):
    """Test the parse_arguments function for the infer_ortho command. See if arguments are parsed correctly."""

    # Store the original infer_orthology.main function
    original_main = infer_orthology.main
    
    # Create a mock to capture the arguments passed to main()
    captured_args = None
    
    def mock_main(args):
        nonlocal captured_args
        captured_args = args
        # Don't actually run the main function, just capture args
    
    # Replace infer_orthology.main with our mock
    monkeypatch.setattr(infer_orthology, 'main', mock_main)

    test_args = ["h2o", "infer_ortho", "-t", "tests/test_data/homolog_trees", "-of", "tests/test_data/outgroups.txt", "-e", ".tre", "-m", "3", "-mp", "0.5", "-od", "tests/test_data/ortholog_trees"]
    monkeypatch.setattr(sys, "argv", test_args)
    
    args = cli.parse_arguments()
    args.func(args)
    assert captured_args.homolog_tree_dir == "tests/test_data/homolog_trees"
    assert captured_args.outgroup_file == "tests/test_data/outgroups.txt"
    assert captured_args.tree_file_ending == ".tre"
    assert captured_args.min_dupl_tip_overlap == 3
    assert captured_args.min_dupl_percentage_overlap == 0.5
    assert captured_args.output_directory == "tests/test_data/ortholog_trees"
    assert captured_args.no_pruning == False
    assert captured_args.just_pruning == False
    assert captured_args.id2sp_file == None
    assert captured_args.single_sample_duplications == False

    # Restore the original function
    monkeypatch.setattr(infer_orthology, 'main', original_main)

def test_infer_ortho():
    """Test the infer_ortho function, see if it runs correctly"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees",
        outgroup_list="o",
        tree_file_ending=".tre",
        min_dupl_tip_overlap=None,
        min_dupl_percentage_overlap=None,
        output_directory="tests/test_data/processed_trees",
        no_pruning=False,
        just_pruning=False,
        id2sp_file=None,
        single_sample_duplications=False
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
    with open("tests/test_data/processed_trees/pruned/dup_reroot_rooted_processed.tre", "r") as f:
        assert f.read() == "(o:0.0,(((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0,((a:0.0,b:0.0):0.0,(c:0.0,d:0.0):0.0):0.0)D:0.0):0.0;\n"
    
    # remove all the output files created by unit test
    # shutil.rmtree('tests/test_data/processed_trees')

def test_infer_ortho_rooting(capsys):
    """Test the infer_ortho function, see if it roots correctly"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees_rooting",
        outgroup_list=None,
        outgroup_file="tests/test_data/outgroups.txt",
        tree_file_ending=".tre",
        min_dupl_tip_overlap=None,
        min_dupl_percentage_overlap=None,
        output_directory="tests/test_data/processed_trees_rooting",
        no_pruning=False,
        just_pruning=False,
        id2sp_file=None,
        single_sample_duplications=False
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

def test_infer_ortho_single_sample_duplications():
    """Test the single_sample_duplications"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees_id",
        outgroup_list="o",
        tree_file_ending=".tre",
        min_dupl_tip_overlap=None,
        min_dupl_percentage_overlap=None,
        output_directory="tests/test_data/processed_trees_id",
        no_pruning=False,
        just_pruning=False,
        id2sp_file="tests/test_data/id2sp.txt",
        single_sample_duplications=True
    )
    infer_orthology.main(args)

    with open("tests/test_data/processed_trees_id/unpruned/dup_ssd_rooted_processed.tre", "r") as f:
        assert f.read() == "(o:0.0,(((a:0.0,a:0.0):0.0,b:0.0):0.0,c:0.0):0.0):0.0;\n"

    # remove all the output files created by unit test
    shutil.rmtree('tests/test_data/processed_trees_id')

def test_infer_ortho_id():
    """Test id2sp procedure"""

    args = argparse.Namespace(
        homolog_tree_dir="tests/test_data/homolog_trees_id",
        outgroup_list="o",
        tree_file_ending=".tre",
        min_dupl_tip_overlap=None,
        min_dupl_percentage_overlap=None,
        output_directory="tests/test_data/processed_trees_id",
        no_pruning=False,
        just_pruning=False,
        id2sp_file="tests/test_data/id2sp.txt",
        single_sample_duplications=False
    )
    infer_orthology.main(args)

    with open("tests/test_data/processed_trees_id/dup_rooted.tre", "r") as f:
        assert f.read() == "(o:0.0,(((a:0.0,b:0.0):0.0,d:0.0):0.0,((a:0.0,b:0.0):0.0,c:0.0):0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees_id/unpruned/dup_rooted_processed_id.tre", "r") as f:
        assert f.read() == "(o@12:0.0,(((a@123:0.0,b@123:0.0):0.0,d@123:0.0):0.0,((a@345:0.0,b@345:0.0):0.0,c@345:0.0):0.0)D:0.0):0.0;\n"
    with open("tests/test_data/processed_trees_id/unpruned/dup_ortho1_id.tre", "r") as f:
        assert f.read() == "(o@12:0.0,((a@345:0.0,b@345:0.0):0.0,c@345:0.0):0.0):0.0;\n"
    with open("tests/test_data/processed_trees_id/unpruned/dup_ssd_rooted_processed.tre", "r") as f:
        assert f.read() == "(o:0.0,((b:0.0,a:0.0):0.0,c:0.0):0.0):0.0;\n"

    # remove all the output files created by unit test
    # shutil.rmtree('tests/test_data/processed_trees_id')