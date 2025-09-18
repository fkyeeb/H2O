Table of Contents
- [Overview](#overview)
- [Orthology and Gene duplication Inference](#orthology-and-gene-duplication-inference)

# Overview
H2O uses the command line interface to run.

# Orthology and Gene duplication Inference

```console
$ h2o infer_ortho
```

There will be 2 folders inside the output folder:
- `no_pruning/` - contains trees with less tip pruning - I will explain more later
- `pruned/` - contain trees with more tip pruning

Each of the two folders contain these following processed trees for each homolog tree file. If the input tree file name is `cluster1.subtree`, the script will write:
- `cluster1_rooted.tre`, `cluster1_rooted_pruned.tre`
- `cluster1_ortho1.tre`...