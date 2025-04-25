[![Build Status](https://travis-ci.com/evogytis/baltic.svg?branch=master)](https://travis-ci.com/evogytis/baltic)
[![downloads](https://anaconda.org/bioconda/baltic/badges/downloads.svg)](https://anaconda.org/bioconda/baltic)

---

### The most recent development version of BALTIC can be found on branch `dev-main`

### For documentation, check out the GitHub [wiki](https://github.com/evogytis/baltic/wiki). For full examples of usage, you can search through the [notebooks](https://github.com/evogytis/baltic/tree/dev-main/docs/notebooks):

- [austechia](https://github.com/evogytis/baltic/blob/dev-main/docs/notebooks/austechia.ipynb): most common BALTIC use cases such as plotting trees, collapsing clades, plotting tanglegrams, etc,
- [curonia](https://github.com/evogytis/baltic/blob/dev-main/docs/notebooks/curonia.ipynb): animations and geospatial integrations,
- [galindia](https://github.com/evogytis/baltic/blob/dev-main/docs/notebooks/galindia.ipynb): integrations with trees generated using the nextstrain toolkit.

### The [BALTIC gallery](https://phylo-baltic.github.io/baltic-gallery/) has several examples of figures generated using BALTIC, and example code for each.

---

# baltic

baltic is a Python library for parsing phylogenetic trees. It takes newick, Nexus or nextstrain JSON trees and allows you to manipulate, explore and visualise them. baltic stands for Backronymed Adaptable Lightweight Tree Import Code if you like that sort of thing.

---

## Installation

Use package manager `pip` to install baltic:

```
pip install baltic
```

---

## Usage

```python
import baltic as bt

# When called with a tree string the `make_tree()` function return a baltic tree object:

treeString='((A:1.0,B:2.0):1.0,C:3.0);'
myTree = bt.make_tree(treeString)

# Otherwise you can import trees from newick, nexus or nextstrain JSON files

newickPath='/Users/myUsername/tree.newick'
myTree = bt.loadNewick(newickPath)

nexusPath='/Users/myUsername/tree.nex'
myTree = bt.loadNexus(nexusPath, absoluteTime = False)

nextstrainPath='https://nextstrain.org/charon/getDataset?prefix=/dengue/denv1'
myTree, myMeta = bt.loadJSON(nextstrainPath)

```

---

Copyright 2016 [Gytis Dudas](https://twitter.com/evogytis). Licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).
