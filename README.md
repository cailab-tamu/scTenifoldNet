scTenifoldNet
=============

A R/MATLAB package to construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on machine learning methods. scTenifoldNet uses principal component regression, tensor decomposition, and manifold alignment, to accurately identify even subtly shifted gene expression programs.

Install:
-------
**scTenifoldNet** requires *numpy* and *scipy* from ***Python***, to install them, we recommend to install *Miniconda Python*: https://docs.conda.io/en/latest/miniconda.html.

After install the ***Phyton*** dependencies, install **scTenifoldNet**, using the following command:

```{R}
library(devtools)
install_github('cailab-tamu/scTenifoldNet')
library(scTenifoldNet)
```

Available functions:
--------------------

|Code| Function |
|:-|:-|
|installPyDependencies||
|scQC||
|cpmNormalization||
|pcNet||
|makeNetworks||
|tensorDecomposition||
|manifoldAlignment||
|dCoexpression||
|scTenifoldNet||
