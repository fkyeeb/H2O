# H2O  <!-- omit in toc -->
![version](https://img.shields.io/badge/version-v.0.9-blue)
[![CI](https://github.com/fkyeeb/H2O/actions/workflows/ci.yml/badge.svg)](https://github.com/fkyeeb/H2O/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fkyeeb/H2O/graph/badge.svg?token=ZHPI95S70O)](https://codecov.io/gh/fkyeeb/H2O)

- [Overview](#overview)
- [Installation](#installation)
  - [Install the most recent release with PyPI](#install-the-most-recent-release-with-pypi)
  - [Build and install from the source](#build-and-install-from-the-source)
  - [Dependency](#dependency)
- [Tutorial](#tutorial)
  - [Example workflows](#example-workflows)
    - [1. *Main workflow* - reducing gene tree discordance created from erroneous orthology inference associated with WGDs](#1-main-workflow---reducing-gene-tree-discordance-created-from-erroneous-orthology-inference-associated-with-wgds)
    - [2. Exploring gene duplications and gene copy losses after WGD](#2-exploring-gene-duplications-and-gene-copy-losses-after-wgd)
    - [3. Simple orthology inference](#3-simple-orthology-inference)
- [References](#references)


# Overview
`h2o` is an ancient-whole-genome-duplication-aware homolog to ortholog trees command-line tool, for plant phylogenomics. It can infer orthology from homologous phylogenetic trees, calculate gene duplication counts, and most importantly, **reduce erroneous orthology inference in clades with putative ancient whole genome duplications (WGD)**. Ancient WGDs tend to associate with high levels of gene tree discordance, and `h2o` is designed to reduce the discordance created from erroneous orthology inference associated with WGDs.

The publication of `h2o` is still in prep.

# Installation
Installation should work the same way on macOS, Linux, and Windows. You can either install with PyPI or from the source.

## Install the most recent release with PyPI

```console
pip install h2o-phy
```

## Build and install from the source

`git clone` the repo to your local directory
```console
git clone https://github.com/fkyeeb/H2O.git
```

Navigate inside the `H2O/` folder and build the package
```console
cd H2O
python -m build
```
Then install the package
```console
pip install dist/h2o-*.tar.gz
```

## Dependency
This package only requires Python version >= 3.8 to run, and no extra Python libraries installation is required.

# Tutorial
**The detailed tutorial of each subcommand of `h2o` is [here](tutorials/tutorial.md).** The content below provides some example workflows for different user scenarios.

> [!TIP]
> `h2o` orthology inference requires homolog trees to have **monophyletic outgroup(s)** to be processed. It will skip the tree and print out error messages if a homolog tree does not have monophyletic outgroup(s). Consider cleaning your dataset if most of your outgroups are not monophyletic.

## Example workflows

### 1. *Main workflow* - reducing gene tree discordance created from erroneous orthology inference associated with WGDs

**User scenario**: I have a plant phylogenomic dataset with known WGD event(s), and these events are correlated with gene tree discordance. It is difficult to resolve the relationships right after WGD events, and I would like to explore alternative relationships.

`h2o` employs two approaches to reduce gene tree conflicts induced by WGDs:
1. Remove taxa that no longer retain both gene copies generated from WGD - modified trees are subsequently referred to as <u>*pruned*</u>
2. Select homolog trees that show the gene duplication from WGD, and only use ortholog trees from these homologs for species tree inference - selected ortholog trees are subsequently referred to as <u>*WGD*</u> trees

![drawio](tutorials/main_workflow.drawio.svg)

This workflow will produce 4 different sets of ortholog trees (purple box) as in the table below. **Both approaches to reduce gene tree conflict remove a lot of data from the trees.** We highly recommend users compare the summary "species" tree produced from each set of orthologs and the gene tree conflict results to select the best supported summary tree inferred from your dataset, as shown in the flowchart below. The steps enclosed by dashed lines are optional.

|  | unpruned |  pruned |
|------|---|---|
| No WGD selection | Ortholog Trees | Pruned Ortholog Trees |
| WGD | WGD Ortholog Trees | WGD Pruned Ortholog Trees |

The purpose is to explore a better hypothesis for the relationships after WGD events. While `h2o` removes data from the trees and reduces gene tree conflict, it can cause more nested relationships to be less supported. 

**In your "best supported summary tree", if you found more support in the relationships right after WGD events but less support in more nested relationships, please proceed with the optional steps. Subcommand `h2o constraint` can easily extract a constraint tree (for the relationships right after WGD events) from a summary tree, and users can use the constraint tree with the full ortholog dataset to produce a constrained summary tree. Then the "best summary tree hypothesis" will contain the better supported relationships for both the relationships right after WGD and also more nested relationships.*

![drawio](tutorials/main_workflow2.drawio.svg)

Details of the subcommands and example external software mentioned are in this [tutorial](tutorials/tutorial.md).

### 2. Exploring gene duplications and gene copy losses after WGD

**User scenario**: I would like to explore the patterns of gene duplications in my plant phylogenomic dataset and gene copy losses after putative ancient WGD events.

> [!TIP]
> The summary tree topology is going to affect the counts. If you are exploring different phylogenetic hypotheses with the main workflow, you may want to use the best-supported topology here.

![drawio](tutorials/dup_workflow.drawio.svg)

Details of the subcommands are in this [tutorial](tutorials/tutorial.md).

### 3. Simple orthology inference

**User scenario**: I have a plant phylogenomic dataset with good outgroups, and I would like to produce ortholog trees in one simple step.

This workflow is the first step of the main workflow. Although `infer_ortho` is similar to ortholog inference methods in Yang and Smith (2014), it is not the same as any of the 4 methods. `h2o infer_ortho` requires monophyletic outgroups in the tree and does not keep any tree without outgroups.

![drawio](tutorials/ortho_workflow.drawio.svg)

# References

Yang, Y., and Smith, S. A. (2014). Orthology inference in nonmodel organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. *Molecular biology and evolution*, 31(11), 3081-3092.