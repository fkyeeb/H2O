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
- [4. Inferring gene loss after WGD events - `h2o gene_loss`](#4-inferring-gene-loss-after-wgd-events---h2o-gene_loss)
  - [4.1 Input](#41-input)
  - [4.2 How to run](#42-how-to-run)
  - [4.3 Output](#43-output)
- [5. Extracting `bp` conflict result for ploting - `h2o bp2pie`](#5-extracting-bp-conflict-result-for-ploting---h2o-bp2pie)
  - [5.1 Input](#51-input)
  - [5.2 How to run](#52-how-to-run)
  - [5.3 Output](#53-output)
- [6. Extracting constraint tree - `h2o constraint`](#6-extracting-constraint-tree---h2o-constraint)
  - [6.1 input](#61-input)
  - [6.3 Output](#63-output)

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
| `-t` | `--homolog_tree_dir` | Yes | Folder containing homolog trees
| `-o` | `--outgroup_list` | Yes | List of outgroup taxa, separated by commas, no spaces (mutually exclusive with `-of`)
| `-of` | `--outgroup_file` | Yes | File containing the outgroup taxa, each line is a taxon (mutually exclusive with `-o`)
| `-e` | `--tree_file_ending` | Yes | File ending of the homolog trees
| `-od` | `--output_directory` | No | Output directory, default is creating a `processed_trees/` directory in the parent directory of `homolog_tree_dir` 
| `-m` | `--min_ingroup_taxa` | No | Minimum number of ingroup taxa, default is 3
| `-p` | `--just_pruning` | No | Only produce pruned ortholog trees
| `-np` | `--no_pruning` | No | Only produce unpruned ortholog trees

### 1.2.2 Running the command with example dataset <!-- omit in toc -->

The folders may vary based on your current directory and where you downloaded the dataset. Before running `h2o`, I do not know whether *pruning* is going to be helpful to for my dataset, so I decided to run the default setting, which is producing both the *pruned* and *unpruned* orthologs for comparison. To do this, I do not flag `-p` or `-np`. 

```console
$ h2o infer_ortho -t example_data/ERIC_homolog_subset/ -of example_data/ERIC_outgroup.txt -e .subtree
```
If the output directory does not exist, `h2o` will create the folder. Because all output trees are going to be created inside the output directory, this directory, if exists, is recommended to be empty before running the command.

If I want to run `infer_ortho` again to only produce pruned ortholog trees, I will add `-p` to the command. To only unpruned ortholog trees, I will add `-np`. For example, **if you are only looking for some fast ortholog trees without _pruning_**, do:

```console
$ h2o infer_ortho -t example_data/ERIC_homolog_subset/ -of example_data/ERIC_outgroup.txt -e .subtree -np
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
| `-t` | `--processed_tree_dir` | Yes | Folder containing processed trees, is the output folder from `infer_ortho` |
| `-s` | `--species_tree_file` | Yes | Species tree file |
| `-od` | `--output_directory` | No | Output directory, default is creating an `other_output/` directory in the parent directory of `processed_tree_dir` |

### 2.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o map_dupl -t example_data/processed_trees -s example_data/ERIC_ASTRAL_rooted_unpruned.tre
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
| `-t` | `--processed_tree_dir` | Yes | Folder containing processed trees, is the output folder from `infer_ortho` |
| `-n` | `--wgd_nodes` | Yes | List of WGD node numbers, separated by commas, no spaces |
| `-d` | `--duplication_counts_dir` | No | Duplication counts directory, default is `other_output` directory in the parent directory of `processed_tree_dir`* |
| `-od` | `--output_directory` | No | Output directory, default is creating an `other_output/` directory in the parent directory of `processed_tree_dir`  |

**If the directory `duplication_counts.tsv` is currently not in `other_output/`, please specify the location of the folder with `-f`.*

### 3.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o extract_wgd_trees -t example_data/processed_trees -n idk
```

## 3.3 Output
Unless specified otherwise, `other_output/` is the default output folder. Inside the folder, these new files will be created:
- `cat_unpruned_wgd_trees.sh` and `cat_pruned_wgd_trees.sh` - bash script to concatenate the *unpruned* and *pruned* ortholog trees into one file
  - `extract_wgd_trees` already ran this script for you. If needed, to run these script again, do `bash cat_unpruned_wgd_trees.sh` or `bash cat_pruned_wgd_trees.sh`
- `ASTRAL_in_unpruned_wgd_trees.tre` and `ASTRAL_in_pruned_wgd_trees.tre` - concatenated ortholog tree files
  - One only contains *unpruned* orthologs that are from homolog trees that shows gene duplications at the node of WGD events; the other one only contains those that are *pruned*.
  - These tree files can be used directly as ASTRAL inputs.

# 4. Inferring gene loss after WGD events - `h2o gene_loss`

This subcommand records and maps gene loss after known putative WGD events. 

## 4.1 Input
This command requires the *unpruned* `*_rooted_processed.tre` from `infer_ortho` command, `duplication_counts.tsv` and `summary_tree_numbered.tre` from `map_dupl` command. It also requires known putative WGD events.

## 4.2 How to run

### 4.2.1 Command Options <!-- omit in toc -->

| Option | Long Option Name | Required | Description |
|--------|-------------|----------|-------------|
| `-t` | `--processed_tree_dir` | Yes | Folder containing processed trees, is the output folder from `infer_ortho` |
| `-n` | `--wgd_nodes` | Yes | List of WGD node numbers, separated by commas, no spaces |
| `-d` | `--duplication_counts_dir` | No | Duplication counts directory, default is `other_output` directory in the parent directory of `processed_tree_dir`* |
| `-od` | `--output_directory` | No | Output directory, default is creating an `other_output/` directory in the parent directory of `processed_tree_dir`  |

**If the directory `duplication_counts.tsv` and `summary_tree_numbered.tre` are currently not in `other_output/`, please specify the location of the folder with `-f`.*

### 4.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o gene_loss -t example_data/processed_trees -n idk
```

## 4.3 Output

Unless specified otherwise, `other_output/` is the default output folder. Inside the folder, these new files will be created:
- `gene_loss_counts.tsv` - gene loss data after given WGD events. The tab-delimited file will look something like this:

  | tree | dupl_clade_count | wgd_node | gene_loss_nodes | gene_loss_tips |
  |------|---|---|---|---|
  | cluster1 | 1 | 3 | 4 | tip1 |
  | cluster1 | 2 | 3 | 5 | tip2,tip4 |
  | cluster2 | 1 | 3 | 4 |  |
  | ... |

  For example, the 3rd line of this example file means: it's data for tree <u>cluster1</u>, it's the <u>second duplicated clade</u> in this tree for <u>WGD event at node 3</u> of the specis tree. For this duplicated clade, there is a <u>gene loss at node 5</u> of the species tree, there is also single tip gene loss for <u>tip2</u> and <u>tip4</u>.
- `gene_loss_counts.tre` - gene loss data summarized on the summary species tree. There are two trees inside this newick tree file. The first tree has node numbers as node labels; the second tree has the number of gene loss events of all homolog trees with the given WGD gene duplication at each node as node labels.

# 5. Extracting `bp` conflict result for ploting - `h2o bp2pie`

`h2o` is built to reduce gene tree conflict induced by WGDs, then evaluating gene tree conflict throughout the process is important. We tend to use [bellerophon](https://git.sr.ht/~hms/bellerophon) (`bp`) to infer gene tree conflict, so we implemented a subcommand in `h2o` to take `bp` output directly and digest it into files that is easy to plot in R and with[`gokstad`](https://git.sr.ht/~hms/gokstad).

## 5.1 Input
`bp` output file is required. To produce the appropriate `bp` result file, `-tv` has to be flagged in the command:

```console
$ bp -c example_data/ERIC_ASTRAL_rooted_unpruned.tre -t example_data/ERIC_ASTRAL_in_unpruned.tre -tv > bp_output_unpruned.txt
```

We also recommend:
- flag `-v` for more verbose results if you would like to see if there is a dominant conflicting topology
- flag `-scut` for support cutoff if you do not want to consider relationships with lower support
- flag `-rng` if you have too many trees and would like to "parallelize" `bp` by hand. `bp` can be a bit slow with a great number of trees and does not have parallel options. `h2o bp2pie` do not support summarizing multiple `bp` output for now, but can be implemented in the future.

## 5.2 How to run

### 5.2.1 Command Options <!-- omit in toc -->

| Option | Long Option Name | Required | Description |
|--------|-------------|----------|-------------|
| `-f` | `--bp_output_file` | Yes | bp output file, `-tv` has to be flagged when running bp |
| `-s` | `--summary_tree_file` | No | Summary tree file*, provide if branch length different from bp tree |
| `-p` | `--pie_option` | No | Flag to include unsupported** counts in the gokstad pie tree |
| `-od` | `--output_directory` | No | Output directory, default is the current directory  |

**In case users want to plot `bp` results with a different branch length, you can provide a tree with the same topology but different branch length, compared to the tree you ran `bp` with.*

***In `bp` results, unsupported means that the tree provided does not contain information about this specific bipart/relationship in the summary tree. It is usually due to missing taxa and having a large number of unsupported trees for each node is normal in large genomic datasets. However, as most folks will interpret unsupported as low support at first glance, the default of `bp2pie` does not include unsupported counts for `gokstad` plotting.*

### 5.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o bp2pie -f bp_output_unpruned.txt
```

## 5.3 Output

Unless specified otherwise, the current directory is the default output folder. Inside the folder, these new files will be created:
- `bp_output.tre` - Counts for conflict, concordance, and unsupported as node labels on each tree. This is part of `bp` output. We provide it as a way for visualization and data exploration.
- `bp_data.tsv` - Counts for conflict, concordance, and unsupported listed with corresponding node number in the summary tree. This is for plotting pie charts in the tree with R. The tab-delimited file will look something like this:

  | node_number | conflict | concord | unsupported |
  |------|---|---|---|
  | 1 | 1 | 3 | 0 |
  | 2 | 0 | 2 | 2 |
  | 3 | 4 | 0 | 0 |
  | ... |
-  `bp_summary_tree_numbered.tre` - Node number as node label, to correpond with `bp_data.tsv`. The default tree is the summary tree used in `bp` analysis. If a summary tree with different branch length is provided with `-s`, then this corresponding tree with `bp_data.tsv` will be the provided tree.
-  `gokstad_pie.tre` - The input tree for `gokstad` plotting, cannot be opened with any tree visualizing application, such as figtree. Example usage with `gokstad`:
```console
$ gokstad -s -d -b -pie gokstad_pie.tre -o gokstad.svg
```

# 6. Extracting constraint tree - `h2o constraint`
`h2o` removes a lot of data from phylogenomic datasets. Although it offers a topology with more support for relationships right after WGDs, it can lead to less support for more embedded and well-defined clades. To resolve this dilemma, we offer an option in `h2o` to extract the "wanted" relationships from a summary tree to use as constraint for phylogenetic analysis.

## 6.1 input
This command requires a summary tree file and a list of nodes or tips to keep in the constraint tree. 

### 6.2.1 Command Options <!-- omit in toc -->

| Option | Long Option Name | Required | Description |
|--------|-------------|----------|-------------|
| `-s` | `--summary_tree_file` | Yes | Summary tree file*, provide if branch length different from bp tree |
| `-od` | `--output_directory` | No | Output directory, default is the current directory  |
| `-n` | `--nodes` | Yes | List** of nodes to keep, node labels separated by commas, no spaces (mutually exclusive with `-t`) |
| `-t` | `--tips_file` | Yes | File*** containing the tips to keep, each line is a tip (mutually exclusive with `-n`) |

**The summary tree file needs to have unique node labels for each node, e.g. in `summary_tree_numbered.tre`. Only the first tree will be read as summary tree. Other trees will be ingnored.*

***Only one tip is going to be kept for each node in the constraint tree. This one tip is going to be selected by random. If there are single tips that you would like to keep aside from the list of nodes, you can list the tip name with the node numbers, e.g. 2,3,tip_name or 6,tip_name,20*

****Only the tips listed in this file are going to show up in the constraint tree.*

### 6.2.2 Running the command with example dataset <!-- omit in toc -->

```console
$ h2o constraint -s summary_tree_numbered.tre -n idk
```
["141","132","115","103","77","101","87","2"]
["Cornales_Nyssaceae_Nyssa_sinensis","Ericales_Balsaminaceae_Impatiens_hawkeri","Ericales_Ericaceae_Gaultheria_nummularioides","Ericales_Lecythidaceae_Lecythis_congestiflora","Ericales_Polemoniaceae_Linanthus_californicus","Ericales_Primulaceae_Primula_veris","Ericales_Sapotaceae_Sarcosperma_laurinum","Ericales_Sapotaceae_Manilkara_sapota","Ericales_Ebenaceae_Diospyros_lotus"]

## 6.3 Output
Unless specified otherwise, the current directory is the default output folder. Inside the folder, this new file will be created:
- `constraint_tree.tre` - the constraint tree, example usage with ASTRAL:
```console
$ astral4 -o ERIC_ASTRAL_out_constraint.tre -c constraint_tree.tre ERIC_ASTRAL_in.tre
```