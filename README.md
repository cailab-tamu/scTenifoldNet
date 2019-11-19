scTenifoldNet
=============

A R/MATLAB package to construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on machine learning methods. scTenifoldNet uses principal component regression, tensor decomposition, and manifold alignment, to accurately identify even subtly shifted gene expression programs.

Install:
-------
**scTenifoldNet** requires *numpy* and *scipy* from the ***Python*** programming language, to install them, we recommend to install *miniconda Python*: https://docs.conda.io/en/latest/miniconda.html.

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
|installPyDependencies|Install the *Phyton* dependencies using the *reticulate* package|
|scQC|Performs single-cell data quality control|
|cpmNormalization|Performs counts per million (CPM) data normalization|
|pcNet|Computes a gene regulatory network based on principal component regression|
|makeNetworks|Computes gene regulatory networks for subsamples of cells based on principal component regression|
|tensorDecomposition|Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition|
|manifoldAlignment|Performs non-linear manifold alignment of two gene regulatory networks|
|dCoexpression|Evaluates gene differential coexpression based on manifold alignment distances|
|scTenifoldNet|Gonstruct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on principal component regression, tensor decomposition, and manifold alignment.|
