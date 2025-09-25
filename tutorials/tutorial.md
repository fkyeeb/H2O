Table of Contents
- [Overview](#overview)
- [1. Orthology and Gene duplication Inference - `h2o infer_ortho`](#1-orthology-and-gene-duplication-inference---h2o-infer_ortho)
  - [1.1 Input](#11-input)
  - [1.2 How to run](#12-how-to-run)
    - [1.2.1 Command Options](#121-command-options)
    - [Running the command with example dataset](#running-the-command-with-example-dataset)
- [under construction](#under-construction)

# Overview
`h2o` uses the command line interface to run.

`h2o` is designed for datasets that have putative WGDs, but is also compatible with datasets without putative ancient whole-genome duplications (WGD). `h2o` employs two approaches to reduce gene tree conflicts induced by WGDs:
1. Remove tips that have lost one of the gene copy from WGD - subsequently refer as *pruning*  - implemented in `infer_ortho`
2. Select ortholog trees that are from homolog trees that show the gene duplication from WGD, and only use them for species tree inference - implemented in ``

There is a very small subset of Ericales phylotranscriptomic dataset (Carruthers et al. 2024) in the `example_data` folder. The following tutorial ran on this example dataset.

Carruthers, T., D. J. P. Gonçalves, P. Li, A. S. Chanderbali, C. W. Dick, P. W. Fritsch, D. A. Larson, et al. 2024. Repeated shifts out of tropical climates preceded by whole genome duplication. *New Phytologist* 244: 2561–2575.

# 1. Orthology and Gene duplication Inference - `h2o infer_ortho`

This is the most important step for `h2o`. It takes homolog trees as input and output cleaned ortholog trees for species tree inference. Output includes both ortholog trees experienced *pruning* and the same set of tree without *pruning*.

! make them able to disable pruning

## 1.1 Input

Homolog trees in individual newick tree files in one folder. 

**We recommend users not to input homolog trees with branches longer than 0.5 substitution/site.** If homolog trees contain branches that is longer than 0.5 substitution/site, then this tree should be separated into two homolog trees at this branch. 0.5 substitution/site means on average half of the sites have substitutions. This branch is too long for a real homologous gene family. `h2o` does not give warnings about branches longer than 0.5, and is likely to take it as a gene duplication and break it into orthologs. In reality this may not interfere with species tree inference. It may inflate gene duplication numbers, but the inflation should not be significant if the dataset is reasonably good.

## 1.2 How to run

### 1.2.1 Command Options

You can also find this information by using `h2o infer_ortho -h`

| Options  | Longer Name | Description
| ------------- | ------------- | ------------- |
| `-d` | `--homolog_tree_dir` | Folder containing homolog trees (required)
| `-o` | `--outgroup_list` | List of outgroup taxa, separated by commas, no spaces (required, mutually exclusive with `-of`)
| `-of` | `--outgroup_file` | File containing the outgroup taxa, each line is a taxon (required, mutually exclusive with `-o`)
| `-t` | `--tree_file_ending` | File ending of the homolog trees (required)
| `-m` | `--min_ingroup_taxa` | Minimum number of ingroup taxa, default is 3
| `-od` | `--output_directory` | Output directory (required)
| `-p` | `--just_pruning` | Only produce pruned ortholog trees (optional)
| `-np` | `--no_pruning` | Only produce unpruned ortholog trees (optional)


### Running the command with example dataset

The folders may vary based on your current directory and where you downloaded the dataset. Before running `h2o`, I do not know whether *pruning* is going to be helpful to for my dataset, so I decided to run the default setting, which is producing both the *pruned* and *unpruned* orthologs for comparison. To do this, I do not flag `-p` or `-np`. If I want to run `infer_ortho` again to only produce pruned ortholog trees, I will flag `-p`. To only unpruned ortholog trees, I will flag `-np`.

```console
$ h2o infer_ortho -d example_data/ERIC_homolog_subset/ -of example_data/ERIC_outgroup.txt -t .subtree -od example_data/ERIC_ortholog
```
# under construction

There will be 2 folders inside the output folder:
- `unpruned/` - contains trees with less tip pruning - I will explain more later
- `pruned/` - contain trees with more tip pruning

Each of the two folders contain these following processed trees for each homolog tree file. If the input tree file name is `cluster1.subtree`, the script will write:
- `cluster1_rooted.tre`, `cluster1_rooted_pruned.tre`
- `cluster1_ortho1.tre`...