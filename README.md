# H2O  <!-- omit in toc -->
![version](https://img.shields.io/badge/version-v.0.9-blue)
[![CI](https://github.com/fkyeeb/H2O/actions/workflows/ci.yml/badge.svg)](https://github.com/fkyeeb/H2O/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fkyeeb/H2O/graph/badge.svg?token=ZHPI95S70O)](https://codecov.io/gh/fkyeeb/H2O)

- [Overview](#overview)
- [Installation](#installation)
- [Tutorial](#tutorial)
  - [Who should use `h2o`?](#who-should-use-h2o)
    - [I have phylogenomic datasets with known WGD events and these events are correlated with gene tree conflict](#i-have-phylogenomic-datasets-with-known-wgd-events-and-these-events-are-correlated-with-gene-tree-conflict)
    - [I have phylogenomic datasets (homolog trees) and I want to produced cleaned ortholog trees and/or gene duplication counts.](#i-have-phylogenomic-datasets-homolog-trees-and-i-want-to-produced-cleaned-ortholog-trees-andor-gene-duplication-counts)
    - [Other useful subcommands](#other-useful-subcommands)


# Overview

Ancient WGD-aware homolog to ortholog trees command-line tool

# Installation

upload to PyPI?

This package only requires python version >= 3.6 to run and no extra python libraries installation required.

# Tutorial
## Who should use `h2o`?

### I have phylogenomic datasets with known WGD events and these events are correlated with gene tree conflict
![drawio](tutorials/package_workflow.drawio.svg)

### I have phylogenomic datasets (homolog trees) and I want to produced cleaned ortholog trees and/or gene duplication counts.

### Other useful subcommands

gene loss, bp2pie, extract constraint tree

The detailed tutorial of each subcommand of `h2o` is [here](tutorials/tutorial.md).