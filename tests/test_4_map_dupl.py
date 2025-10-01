import argparse
import os
import shutil
from h2o import map_duplications


def test_map_dupl():
    """Test the map_dupl function with advance"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        output_directory="tests/test_data/map_dupl_output",
        species_tree_file="tests/test_data/species.tre"
    )
    map_duplications.main(args)


    with open('tests/test_data/map_dupl_output/duplication_counts.tsv', 'r') as f:
        lines = f.readlines()
        assert lines[3] == 'dup_bl\t0\t0\t0\t3\t0\t0\t0\t\n'

    with open('tests/test_data/map_dupl_output/summary_tree_numbered.tre', 'r') as f:
        assert f.read() == '(o:0.0,((e:0.0,f:0.0)2:0.0,((a:0.0,b:0.0)4:0.0,(c:0.0,d:0.0)5:0.0)3:0.0)1:0.0)0:0.0;\n(o:0.0,((e:0.0,f:0.0)0:0.0,((a:0.0,b:0.0)0:0.0,(c:0.0,d:0.0)0:0.0)4:0.0)0:0.0)0:0.0;\n'

    shutil.rmtree('tests/test_data/map_dupl_output')

def test_map_dupl_no_output_directory():
    """Test the map_dupl function with no output directory"""

    args = argparse.Namespace(
        processed_tree_dir="tests/test_data/processed_trees",
        species_tree_file="tests/test_data/species.tre",
        output_directory=None
    )
    map_duplications.main(args)

    assert os.path.exists('tests/test_data/other_output/duplication_counts.tsv')

    # shutil.rmtree('tests/test_data/processed_trees')
    # shutil.rmtree('tests/test_data/other_output')
