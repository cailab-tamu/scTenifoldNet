scTenifoldNet
=============

A R/MATLAB package to construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data collected from different conditions based on machine learning methods. scTenifoldNet uses principal component regression, tensor decomposition, and manifold alignment, to accurately identify even subtly shifted gene expression programs.

Install:
-------
**scTenifoldNet** is available through the CRAN repositories, you can install **scTenifoldNet**, using the following command:
```{R}
install.packages('scTenifoldNet')
library(scTenifoldNet)
```
Or if you are interested in the version in development, you can install **scTenifoldNet**, using the following command:
```{R}
library(remotes)
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

Input:
--------
The required input for **scTenifoldNet** is an expression matrix with genes in the rows and cells (barcodes) in the columns. Data is expected to be _not normalized_ if the main **scTenifoldNet** function is used. Given the modular structure of the package, users are free to include modifications in each step to perform their analysis.

Output:
--------
The output of **scTenifoldNet** is a list with 3 slots as follows: 
  * **tensorNetworks**: The computed weight-averaged denoised gene regulatory networks after CANDECOMP/PARAFAC (CP) tensor decomposition. It includes two slots with:
    * **X**: The constructed network for the _X_ sample.
    * **Y**: The constructed network for the _Y_ sample.
  * **manifoldAlignment**: The generated low-dimensional features result of the non-linear manifold alignment. It is a data frame with _genes_ in the rows and _d_ (default= 30) dimensions in the columns
  * **diffRegulation**: The results of the differential regulation analysis.

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
                          td_K = 3, qc_minLibSize = 30)

outputHA <- scTenifoldNet(X = X, Y = Y,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30)
```
#### Differential regulation based on manifold alignment distances
As is shown below, under the H0, none of the genes shown a significative difference in regulatory profiles using an FDR cut-off of 0.1, but under the HA, the 6 genes involved in the perturbation (50, 11, 2, 10, 5, and 3) are identified as perturbed.
```
head(outputH0$diffRegulation, n = 10)
#    gene     distance        Z       FC    p.value     p.adj
# 59 ng59 1.138659e-15 2.497859 4.338493 0.03725989 0.6411916
# 32 ng32 9.970451e-16 2.081719 3.326448 0.06817394 0.6411916
# 72 ng72 9.452170e-16 1.917581 2.989608 0.08380044 0.6411916
# 89 ng89 8.924023e-16 1.742758 2.664849 0.10258758 0.6411916
# 12 ng12 8.851226e-16 1.718018 2.621550 0.10542145 0.6411916
# 17 ng17 8.784453e-16 1.695183 2.582145 0.10807510 0.6411916
# 31 ng31 8.621191e-16 1.638760 2.487057 0.11478618 0.6411916
# 95 mt-5 8.578947e-16 1.624022 2.462743 0.11657501 0.6411916
# 57 ng57 8.141541e-16 1.467912 2.218015 0.13640837 0.6411916
# 77 ng77 7.888569e-16 1.374548 2.082321 0.14901342 0.6411916

head(outputHA$diffRegulation, n = 10)
#    gene   distance        Z        FC      p.value      p.adj
# 2   ng2 0.023526702 2.762449 12.193413 0.0004795855 0.02414332
# 50 ng50 0.023514429 2.761550 12.180695 0.0004828665 0.02414332
# 11 ng11 0.022443941 2.681598 11.096894 0.0008647241 0.02882414
# 3   ng3 0.020263415 2.508478  9.045415 0.0026335445 0.06583861
# 10 ng10 0.019194561 2.417929  8.116328 0.0043868326 0.07711821
# 5   ng5 0.019079975 2.407977  8.019712 0.0046270923 0.07711821
# 31 ng31 0.013632541 1.865506  4.094085 0.0430335257 0.61476465
# 96 mt-6 0.011401177 1.589757  2.863536 0.0906081350 0.90977795
# 59 ng59 0.009835354 1.368238  2.130999 0.1443466682 0.90977795
# 62 ng62 0.007995812 1.067193  1.408408 0.2353209153 0.90977795
```

#### Plotting the results
Results can be easily displayed using quantile-quantile plots. Here we labeled in red the identified perturbed genes with FDR < 0.1.
![Example](https://raw.githubusercontent.com/cailab-tamu/scTenifoldNet/master/inst/readmeExample.png)
```{r}
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
set.seed(1)
qChisq <- rchisq(100,1)
geneColor <- rev(ifelse(outputH0$diffRegulation$p.adj < 0.1, 10,1))
qqplot(qChisq, outputH0$diffRegulation$FC, pch = 16, main = 'H0', col = geneColor, 
       xlab = expression(X^2~Quantiles), ylab = 'FC', xlim=c(0,8), ylim=c(0,13))
qqline(qChisq)
legend('bottomright', legend = c('FDR < 0.1'), pch = 16, col = 'red', bty='n', cex = 0.7)

geneColor <- rev(ifelse(outputHA$diffRegulation$p.adj < 0.1, 'red', 'black'))
qqplot(qChisq, outputHA$diffRegulation$FC, pch = 16, main = 'HA', col = geneColor, 
       xlab = expression(X^2~Quantiles), ylab = 'FC', xlim=c(0,8), ylim=c(0,13))
qqline(qChisq)
legend('bottomright', legend = c('FDR < 0.1'), pch = 16, col = 'red', bty='n', cex = 0.7)
```

Citation
--------
To cite **scTenifoldNet** in publications use:

  Daniel Osorio, Yan Zhong, Guanxun Li, Jianhua Huang and James Cai (2019). scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data. R package version 1.2.0.
  https://CRAN.R-project.org/package=scTenifoldNet

A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data},
    author = {Daniel Osorio and Yan Zhong and Guanxun Li and Jianhua Huang and James Cai},
    year = {2019},
    note = {R package version 1.2.0},
    url = {https://CRAN.R-project.org/package=scTenifoldNet},
  }
  ```
