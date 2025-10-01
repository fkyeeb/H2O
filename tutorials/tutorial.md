# Table of Contents <!-- omit in toc -->
- [impletement no ortholog production? and some all\_in\_1 options](#impletement-no-ortholog-production-and-some-all_in_1-options)
- [Overview](#overview)
- [1. Orthology Inference - `h2o infer_ortho`](#1-orthology-inference---h2o-infer_ortho)
  - [1.1 Input](#11-input)
  - [1.2 How to run](#12-how-to-run)
  - [1.3 Output](#13-output)
- [2. Gene duplication summary statistics - `h2o map_dupl`](#2-gene-duplication-summary-statistics---h2o-map_dupl)
  - [2.1 Input](#21-input)
  - [2.2 How to run](#22-how-to-run)
  - [2.3 Output](#23-output)
- [3. Extracting trees that shows gene duplications at WGD events - `h2o extract_wgd_trees`](#3-extracting-trees-that-shows-gene-duplications-at-wgd-events---h2o-extract_wgd_trees)
  - [3.1 Input](#31-input)
  - [3.2 How to run](#32-how-to-run)
  - [3.3 Output](#33-output)

# impletement no ortholog production? and some all_in_1 options

# Overview
`h2o` uses the command line interface to run.

`h2o` is designed for datasets that have putative WGDs, but is also compatible with datasets without putative ancient whole-genome duplications (WGD). If your dataset is not susceptible to gene tree conflict potentially caused by WGDs, the orthology inference command [`h2o infer_ortho`](#1-orthology-inference---h2o-infer_ortho) can produce cleaned ortholog trees very fast and easy. If your dataset is susceptible to gene tree conflict potentially caused by WGDs, I recommend going through the whole tutorial and try all the commands.

`h2o` employs two approaches to reduce gene tree conflicts induced by WGDs:
1. Remove tips that have lost one of the gene copy from WGD - subsequently refer as *pruning*  - implemented in [`infer_ortho`](#1-orthology-inference---h2o-infer_ortho)
2. Select ortholog trees that are from homolog trees that show the gene duplication from WGD, and only use them for species tree inference - implemented in [`extract_wgd_trees`](#3-extracting-trees-that-shows-gene-duplication-at-wgd-events---h2o-extract_wgd_trees)

Most of the subcommands are not independent of each other. They usually require some results from the step before it. **To run through the whole pipeline, please go through the tutorial in order.**

There is a very small subset of Ericales phylotranscriptomic dataset (Carruthers et al. 2024) in the `example_data` folder. The following tutorial ran on this example dataset.

Carruthers, T., D. J. P. Gonçalves, P. Li, A. S. Chanderbali, C. W. Dick, P. W. Fritsch, D. A. Larson, et al. 2024. Repeated shifts out of tropical climates preceded by whole genome duplication. *New Phytologist* 244: 2561–2575.

# 1. Orthology Inference - `h2o infer_ortho`

This is the most important step for `h2o`. It takes homolog trees as input and output cleaned ortholog trees for species tree inference. Output includes both ortholog trees experienced *pruning* and the same set of tree without *pruning*.

## 1.1 Input

Homolog trees in individual newick tree files in one folder. 

**We recommend users not to input homolog trees with branches longer than 0.5 substitution/site.** If homolog trees contain branches that is longer than 0.5 substitution/site, then this tree should be separated into two homolog trees at this branch. 0.5 substitution/site means on average half of the sites have substitutions. This branch is too long for a real homologous gene family. `h2o` does not give warnings about branches longer than 0.5, and is likely to take it as a gene duplication and break it into orthologs. In reality this may not interfere with species tree inference. It may inflate gene duplication numbers, but the inflation should not be significant if the dataset is reasonably good.

## 1.2 How to run

### 1.2.1 Command Options <!-- omit in toc -->

You can also find this information by using `h2o infer_ortho -h`

| Option  | Long Option Name | Required | Description
| ------------- | ------------- | ------------- | ------------- |
| `-d` | `--homolog_tree_dir` | Yes | Folder containing homolog trees
| `-o` | `--outgroup_list` | Yes | List of outgroup taxa, separated by commas, no spaces (mutually exclusive with `-of`)
| `-of` | `--outgroup_file` | Yes | File containing the outgroup taxa, each line is a taxon (mutually exclusive with `-o`)
| `-t` | `--tree_file_ending` | Yes | File ending of the homolog trees
| `-od` | `--output_directory` | No | Output directory, default is creating a `processed_trees/` directory in the parent directory of `homolog_tree_dir` 
| `-m` | `--min_ingroup_taxa` | No | Minimum number of ingroup taxa, default is 3
| `-p` | `--just_pruning` | No | Only produce pruned ortholog trees
| `-np` | `--no_pruning` | No | Only produce unpruned ortholog trees

### 1.2.2 Running the command with example dataset <!-- omit in toc -->

The folders may vary based on your current directory and where you downloaded the dataset. Before running `h2o`, I do not know whether *pruning* is going to be helpful to for my dataset, so I decided to run the default setting, which is producing both the *pruned* and *unpruned* orthologs for comparison. To do this, I do not flag `-p` or `-np`. 

```console
$ h2o infer_ortho -d example_data/ERIC_homolog_subset/ -of example_data/ERIC_outgroup.txt -t .subtree
```
If the output directory does not exist, `h2o` will create the folder. Because all output trees are going to be created inside the output directory, this directory, if exists, is recommended to be empty before running the command.

If I want to run `infer_ortho` again to only produce pruned ortholog trees, I will add `-p` to the command. To only unpruned ortholog trees, I will add `-np`. For example, **if you are only looking for some fast ortholog trees without _pruning_**, do:

```console
$ h2o infer_ortho -d example_data/ERIC_homolog_subset/ -of example_data/ERIC_outgroup.txt -t .subtree -np
```

## 1.3 Output

For the default setting, there will be some output and 2 folders write to the output folder `processed_trees/`:
- Rooted homolog trees
- `unpruned/` - contains processed* homolog trees and *pruned* ortholog trees
- `pruned/` - contains processed* homolog trees and *unpruned* ortholog trees

*_processing is a necessary tree cleaning process, involving labeling duplication nodes and removing tip duplications that are less than `min_ingroup_taxa`._

If the input tree file name is `cluster1.subtree` and `.subtree` is entered as the `tree_file_ending`:
- Rooted homolog trees  - `cluster1_rooted.tre`
- Rooted processed trees - `cluster1_rooted_processed.tre`
- Rooted ortholog trees - `cluster1_ortho1.tre`, `cluster1_ortho2.tre` ...


# 2. Gene duplication summary statistics - `h2o map_dupl`
Now that we processed all the trees, we can count the gene duplications at each node to support the inference of putative WGDs. This command has a similar functionality as [phyparts](https://bitbucket.org/blackrim/phyparts/src/master/) duplication command. It provides slightly more detailed results for gene duplication counts as a tsv file.

## 2.1 Input
This command requires that you already ran `h2o infer_ortho` because all the gene duplication inference for individual homolog trees were done in the `infer_ortho` command. `map_dupl` just summarizes all the information. 

This command requires the *unpruned* `*_rooted_processed.tre` from the previous command and also a rooted summary "species" tree computed from all the orthologs. We are only mapping gene duplications with *unpruned* homolog trees because *pruned* trees are meant for summary tree inference and not for gene duplication mapping. If *pruned* trees are used, the gene duplications will be mapped to more nested nodes compared to the correct one.

I provided a rooted summary tree computed from the *unpruned* full Ericales dataset from Carruthers et al. (2024) in `example_data/`. A summary tree made from the subset dataset can be different from the one I provided. To compute the summary tree for the example dataset, you can do the following commands, which uses two other packages:
```console
$ cat example_data/ERIC_ortholog/unpruned/*_ortho[0-9].tre > example_data/ERIC_ASTRAL_in_unpruned.tre
$ cat example_data/ERIC_ortholog/pruned/*_ortho[0-9].tre > example_data/ERIC_ASTRAL_in_pruned.tre
```
This combines all the ortholog trees in one file, one for the *unpruned* trees, one for the *pruned* ones.

Then computes the summary tree with [astral4](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md) and rerooted the tree with [phyx](https://github.com/FePhyFoFum/phyx).  `astral4` is the newest version of ASTRAL as I write this tutorial. The old versions are also fine. I used `-t` to use more threads with `astral4` to make the analyses go faster, you can add `-t` as you see fit for your machine. `phyx` has a lot of useful commands for phylogenetics, highly recommend!

```console
$ astral4 -i example_data/ERIC_ASTRAL_in_unpruned.tre -o example_data/ERIC_ASTRAL_out_unpruned.tre
$ pxrr -t example_data/ERIC_ASTRAL_out_unpruned.tre -f example_data/ERIC_outgroup.txt -o example_data/ERIC_ASTRAL_rooted_unpruned.tre
```
Do the same for *pruned* trees, commands omitted here.

## 2.2 How to run

### 2.2.1 Command Options <!-- omit in toc -->

You can also find this information by using `h2o map_dupl -h`

| Option  | Long Option Name | Required | Description |
| ------------- | ------------- | ------------- | ------------- |
| `-d` | `--processed_tree_dir` | Yes | Folder containing processed trees, is the output folder from `infer_ortho` |
| `-t` | `--species_tree_file` | Yes | Species tree file |
| `-od` | `--output_directory` | No | Output directory, default is creating an `other_output/` directory in the parent directory of `processed_tree_dir` |

### 2.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o map_dupl -d example_data/processed_trees -t example_data/ERIC_ASTRAL_rooted_unpruned.tre
```

## 2.3 Output
Without specifying a different output directory, a folder named `other_output/` will be created inside `example_data/`. Inside the folder, there will be:
- `duplication_counts.tsv` - gene duplication counts for *unpruned* and *processed* homolog clusters. The tab-delimited file will look something like this:

  | tree | 0 | 1 | 2 | 3 | 4 | 5 | ... | n/a |
  |------|---|---|---|---|---|---|---|---|
  | cluster1 | 0 | 0 | 0 | 1 | 0 | 0 | ... | 0 |
  | cluster2 | 0 | 0 | 0 | 0 | 0 | 0 | ... | 0 |
  | cluster3 | 0 | 0 | 0 | 3 | 0 | 0 | ... | 0 |
  | ... |
  
  It summarizes all the gene duplications at each node of each homolog tree. The first column is the name of the homolog tree without `.tre`. The first row is the node numbers in the summary tree, the same node numbers in `summary_tree_numbered.tre`. The last column "n/a" lists the number of gene duplications that is not matched with any node in the summary tree, likely due to gene tree conflict.
- `summary_tree_numbered.tre` - there are two trees inside this newick tree file. The first tree has node numbers as node labels; the second tree has the sum of gene duplications of all homolog trees at each node as node labels.

# 3. Extracting trees that shows gene duplications at WGD events - `h2o extract_wgd_trees`
This subcommand will gather the names of homolog trees that shows gene duplications at the node of *known putative WGD events**. Then it will extract the ortholog trees inferred from these homolog trees and concatenate them into a file that is ready for ASTRAL input.

**Known putative WGD events should be supported by multiple pieces of evidence, not simply from elevated gene duplication counts at nodes. As examples in Yang et al. (2018) and Feng et al. (2024), putative WGDs are supported by Ks plots, chromosome counts, and elevated gene duplication counts.*

> [!TIP]
> If your dataset has multiple WGD events, the events are far from each other in the species tree, and they correlate with elevated gene tree conflicts in separate parts of the tree, we do not recommend running the pipeline in this standardized way. This limitation of `h2o` is specified in the publication. `h2o` will not run into error if this is the case, but the program output may not be as helpful. Feel free to contact KF for help if your dataset has this problem.

Feng, K., J. F. Walker, H. E. Marx, Y. Yang, S. F. Brockington, M. J. Moore, R. K. Rabeler, and S. A. Smith. 2024. The link between ancient whole‐genome duplications and cold adaptations in the Caryophyllaceae. *American Journal of Botany* e16350.

Yang, Y., Moore, M. J., Brockington, S. F., Mikenas, J., Olivieri, J., Walker, J. F., & Smith, S. A. (2018). Improved transcriptome sampling pinpoints 26 ancient and more recent polyploidy events in Caryophyllales, including two allopolyploidy events. *New Phytologist*, 217(2), 855-870.

## 3.1 Input
This subcommand requires the processed ortholog trees from `infer_ortho`, the `duplication_counts.tsv` from `map_dupl`, and node numbers corresponding to WGD events. 

If one of the `pruned/` or `unpruned/`directories doesn't exist in `processed_tree_dir`, `h2o` will only extract ortholog trees from the existing folder. If both of those directories do not exist. This command will not run. The corresponding node numbers have to match with those in the first tree of `summary_tree_numbered.tre`.

## 3.2 How to run

### 3.2.1 Command Options <!-- omit in toc -->

| Option | Long Option Name | Required | Description |
|--------|-------------|----------|-------------|
| `-d` | `--processed_tree_dir` | Yes | Folder containing processed trees, is the output folder from `infer_ortho` |
| `-n` | `--wgd_nodes` | Yes | List of WGD node numbers, separated by commas, no spaces |
| `-f` | `--duplication_counts_file` | No | Duplication counts file, is `other_output/duplication_counts.tsv` from `map_dupl` output* |
| `-od` | `--output_directory` | No | Output directory, default is creating an `other_output/` directory in the parent directory of `processed_tree_dir`  |

**If the directory `duplication_counts.tsv` is currently in is not `other_output/`, please specify the location of `duplication_counts.tsv` with `-f`*

### 3.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o extract_wgd_trees -d example_data/processed_trees -n idk
```

## 3.3 Output
Unless specified otherwise, `other_output/` is the default output folder. Inside the folder, there will be:
- `cat_unpruned_wgd_trees.sh` and `cat_pruned_wgd_trees.sh` - bash script to concatenate the *unpruned* and *pruned* ortholog trees into one file
  - `extract_wgd_trees` already ran this script for you. If needed, to run these script again, do `bash cat_unpruned_wgd_trees.sh` or `bash cat_pruned_wgd_trees.sh`
- `ASTRAL_in_unpruned_wgd_trees.tre` and `ASTRAL_in_pruned_wgd_trees.tre` - concatenated ortholog tree files
  - One only contains *unpruned* orthologs that are from homolog trees that shows gene duplications at the node of WGD events; the other one only contains those that are *pruned*.
  - These tree files can be used directly as ASTRAL inputs.