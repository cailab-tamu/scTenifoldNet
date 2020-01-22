scTenifoldNet
=============

A R/MATLAB package to construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data collected from different conditions based on machine learning methods. scTenifoldNet uses principal component regression, tensor decomposition, and manifold alignment, to accurately identify even subtly shifted gene expression programs.

Install:
-------
**scTenifoldNet** requires *numpy* and *scipy* from the ***Python*** programming language, to install them, we recommend to install *miniconda Python*: https://docs.conda.io/en/latest/miniconda.html.

**scTenifoldNet** is available through the CRAN repositories, after install the ***Phyton*** dependencies, you can install **scTenifoldNet**, using the following command:
```{R}
install.packages('scTenifoldNet')
library(scTenifoldNet)
```
Or if you are interested in the version in development, you can install **scTenifoldNet**, using the following command:
```{R}
library(devtools)
install_github('cailab-tamu/scTenifoldNet')
library(scTenifoldNet)
```

Available functions:
--------------------

|Code| Function |
|:-|:-|
|scQC|Performs single-cell data quality control|
|cpmNormalization|Performs counts per million (CPM) data normalization|
|pcNet|Computes a gene regulatory network based on principal component regression|
|makeNetworks|Computes gene regulatory networks for subsamples of cells based on principal component regression|
|tensorDecomposition|Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition|
|manifoldAlignment|Performs non-linear manifold alignment of two gene regulatory networks|
|dRegulation|Evaluates gene differential regulation based on manifold alignment distances|
|scTenifoldNet|Construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on principal component regression, tensor decomposition, and manifold alignment.|

Example:
--------
#### Loading scTenifoldNet
Once installed, **scTenifoldNet** can be loaded typing:
```{r}
library(scTenifoldNet)
```

#### Simulating of a dataset 
Here we simulate a dataset of 2000 cells (columns) and 100 genes (rows) following the negative binomial distribution with high sparsity (~67%). We label the last ten genes as mitochondrial genes ('mt-') to perform single-cell quality control.
```{r}
nCells = 2000
nGenes = 100
set.seed(1)
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- round(X)
X <- matrix(X, ncol = nCells)
rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
```

#### Generating a perturbed network 
We generate a perturbed network modifying the expression of genes 10, 2, and 3 and replacing them with the expression of genes 50, 11, and 5.
```{r}
Y <- X
Y[10,] <- Y[50,]
Y[2,] <- Y[11,]
Y[3,] <- Y[5,]
```
#### scTenifoldNet
Here we run **scTenifoldNet** under the H0 (there is no change in the regulation of the gene) using the same matrix as input and under the HA (there is a change in the regulation of the genes) using the control and the perturbed network.
```{r}
outputH0 <- scTenifoldNet(X = X, Y = X,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minFC = 0)

outputHA <- scTenifoldNet(X = X, Y = Y,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minFC = 0)
```
#### Differential regulation based on manifold alignment distances
As is shown below, under the H0, none of the genes shown a significative difference in regulatory profiles using an FDR cut-off of 0.1, but under the HA, the 6 genes involved in the perturbation (50, 11, 2, 10, 5, and 3) are identified as perturbed.
```
head(outputH0$diffRegulation, n = 10)
#    gene     distance        Z       FC    p.value     p.adj
# 23 ng23 3.070131e-15 2.322136 3.043773 0.08104736 0.5707127
# 2   ng2 3.028738e-15 2.262939 2.962251 0.08522920 0.5707127
# 19 ng19 2.912207e-15 2.094281 2.738690 0.09794501 0.5707127
# 34 ng34 2.769999e-15 1.884258 2.477752 0.11546719 0.5707127
# 9   ng9 2.617076e-15 1.652876 2.211725 0.13696542 0.5707127
# 79 ng79 2.602655e-15 1.630744 2.187418 0.13914237 0.5707127
# 61 ng61 2.582711e-15 1.600042 2.154021 0.14219705 0.5707127
# 71 ng71 2.489774e-15 1.455541 2.001789 0.15711364 0.5707127
# 20 ng20 2.479350e-15 1.439182 1.985062 0.15885815 0.5707127
# 16 ng16 2.461000e-15 1.410309 1.955788 0.16196447 0.5707127

head(outputHA$diffRegulation, n = 10)
#    gene   distance        Z        FC      p.value      p.adj
# 50 ng50 0.03122694 3.015237 13.196102 0.0002805319 0.01518175
# 11 ng11 0.03027085 2.954575 12.400413 0.0004292389 0.01518175
# 2   ng2 0.03013548 2.945872 12.289751 0.0004554525 0.01518175
# 10 ng10 0.02759308 2.776751 10.303556 0.0013277411 0.03319353
# 5   ng5 0.02657637 2.705876  9.558250 0.0019905282 0.03760933
# 3   ng3 0.02625458 2.683024  9.328183 0.0022565601 0.03760933
# 31 ng31 0.01309219 1.490233  2.319592 0.1277535062 0.91141393
# 96 mt-6 0.01099588 1.223243  1.636240 0.2008421562 0.91141393
# 6   ng6 0.01067100 1.178555  1.540980 0.2144719537 0.91141393
# 59 ng59 0.01063989 1.174223  1.532010 0.2158110783 0.91141393
```

#### Plotting the results
Results can be easily displayed using quantile-quantile plots. Here we labeled in red the identified perturbed genes with FDR < 0.05.
![Example](https://raw.githubusercontent.com/cailab-tamu/scTenifoldNet/master/inst/readmeExample.png)
```{r}
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
geneColor <- ifelse(outputH0$diffRegulation$p.adj < 0.05, 'red', 'black')
qqnorm(outputH0$diffRegulation$Z, pch = 16, main = 'H0', col = geneColor)
qqline(outputH0$diffRegulation$Z)
legend('bottomright', legend = c('FDR < 0.05'), pch = 16, col = 'red', bty='n', cex = 0.7)

geneColor <- ifelse(outputHA$diffRegulation$p.adj < 0.05, 'red', 'black')
qqnorm(outputHA$diffRegulation$Z, pch = 16, main = 'HA', col = geneColor)
qqline(outputHA$diffRegulation$Z)
legend('bottomright', legend = c('FDR < 0.05'), pch = 16, col = 'red', bty='n', cex = 0.7)
```

Citation
--------
To cite **scTenifoldNet** in publications use:

  Daniel Osorio, Yan Zhong, Guanxun Li, Jianhua Huang and James Cai (2019). scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data. R package version 1.1.0.
  https://CRAN.R-project.org/package=scTenifoldNet

A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data},
    author = {Daniel Osorio and Yan Zhong and Guanxun Li and Jianhua Huang and James Cai},
    year = {2019},
    note = {R package version 1.1.0},
    url = {https://CRAN.R-project.org/package=scTenifoldNet},
  }
  ```
