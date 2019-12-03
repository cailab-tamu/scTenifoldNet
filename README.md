scTenifoldNet
=============

A R/MATLAB package to construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data collected from different conditions based on machine learning methods. scTenifoldNet uses principal component regression, tensor decomposition, and manifold alignment, to accurately identify even subtly shifted gene expression programs.

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
|scQC|Performs single-cell data quality control|
|cpmNormalization|Performs counts per million (CPM) data normalization|
|pcNet|Computes a gene regulatory network based on principal component regression|
|makeNetworks|Computes gene regulatory networks for subsamples of cells based on principal component regression|
|tensorDecomposition|Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition|
|manifoldAlignment|Performs non-linear manifold alignment of two gene regulatory networks|
|dCoexpression|Evaluates gene differential coexpression based on manifold alignment distances|
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
Here we run **scTenifoldNet** under the H0 (there is no change in the coexpression profiles) using the same matrix as input and under the HA (there is a change in the coexpression profiles of some genes) using the control and the perturbed network.
```{r}
outputH0 <- scTenifoldNet(X = X, Y = X,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minDist = 0)

outputHA <- scTenifoldNet(X = X, Y = Y,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minDist = 0)
```
#### Differential coexpression based on manifold alignment distances
As is shown below, under the H0, none of the genes shown a significative difference in coexpression profiles using an FDR cut-off of 0.1, but under the HA, the 6 genes involved in the perturbation (50, 11, 2, 10, 5, and 3) are identified as perturbed.
```
head(outputH0$diffCoexpression, n = 10)
#   gene     distance        Z       FC   p.value     p.adj
#23 ng23 3.070131e-15 2.322136 1.828489 0.1763061 0.4731993
#2   ng2 3.028738e-15 2.262939 1.803837 0.1792493 0.4731993
#19 ng19 2.912207e-15 2.094281 1.734434 0.1878456 0.4731993
#34 ng34 2.769999e-15 1.884258 1.649739 0.1989945 0.4731993
#9   ng9 2.617076e-15 1.652876 1.558662 0.2118613 0.4731993
#79 ng79 2.602655e-15 1.630744 1.550073 0.2131246 0.4731993
#61 ng61 2.582711e-15 1.600042 1.538195 0.2148867 0.4731993
#71 ng71 2.489774e-15 1.455541 1.482844 0.2233301 0.4731993
#20 ng20 2.479350e-15 1.439182 1.476636 0.2243016 0.4731993
#16 ng16 2.461000e-15 1.410309 1.465707 0.2260242 0.4731993

head(outputHA$diffCoexpression, n = 10)
   gene   distance        Z       FC    p.value     p.adj
#50 ng50 0.03122694 3.015237 5.407703 0.02004808 0.5497305
#11 ng11 0.03027085 2.954575 5.242133 0.02204623 0.5497305
#2   ng2 0.03013548 2.945872 5.218690 0.02234538 0.5497305
#10 ng10 0.02759308 2.776751 4.778411 0.02881870 0.5497305
#5   ng5 0.02657637 2.705876 4.602345 0.03192826 0.5497305
#3   ng3 0.02625458 2.683024 4.546618 0.03298383 0.5497305
#31 ng31 0.01309219 1.490233 2.267231 0.13213580 0.6588093
#96 mt-6 0.01099588 1.223243 1.904204 0.16760853 0.6588093
#6   ng6 0.01067100 1.178555 1.847943 0.17402286 0.6588093
#59 ng59 0.01063989 1.174223 1.842556 0.17465166 0.6588093
```

#### Plotting the results
Results can be easily displayed using quantile-quantile plots. Here we labeled in red the identified perturbed genes with FDR < 0.1.
![Example](https://raw.githubusercontent.com/cailab-tamu/scTenifoldNet/master/inst/readmeExample.png)
```{r}
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
geneColor <- ifelse(outputH0$diffCoexpression$p.value < 0.05, 'red', 'black')
qqnorm(outputH0$diffCoexpression$Z, pch = 16, main = 'H0', col = geneColor)
qqline(outputH0$diffCoexpression$Z)
legend('bottomright', legend = c('P < 0.05'), pch = 16, col = 'red', bty='n', cex = 0.7)

geneColor <- ifelse(outputHA$diffCoexpression$p.value < 0.05, 'red', 'black')
qqnorm(outputHA$diffCoexpression$Z, pch = 16, main = 'HA', col = geneColor)
qqline(outputHA$diffCoexpression$Z)
legend('bottomright', legend = c('P < 0.05'), pch = 16, col = 'red', bty='n', cex = 0.7)
```
