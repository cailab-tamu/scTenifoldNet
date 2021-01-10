# ScTenifoldNet.jl
 ScTenifoldNet.jl: constructing and comparing scGRNs from single-cell RNAseq (scRNAseq) data

This package provides a Julia implementation of ScTenifoldNet.
See [bioRxiv - scTenifoldNet: a machine learning workflow for constructing and comparing transcriptome-wide gene regulatory networks from single-cell data](https://doi.org/10.1101/2020.02.12.931469)
for more information.

## Installation

```jl
] add https://github.com/jamesjcai/ScTenifoldNet.jl
```

## Usage

Here is a simple example using randomly generated data.

```jl
using ScTenifoldNet
X0=rand(100,1000);
X1=rand(100,1000);
Z0=tenrnet(X0)
Z1=tenrnet(X1)
d,aln0,aln1=manialn(Z0,Z1)
fc,p,adjp=drgenes(d)
```
## Exported functions
|Code| Function |
|:-|:-|
|pcnet|Computes a gene regulatory network based on principal component regression|
|tensordecomp|Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition|
|manialn|Performs non-linear manifold alignment of two gene regulatory networks|
|drgenes|Evaluates gene differential regulation based on manifold alignment distances|
|tenrnet|Subsamples cells, constructs single-cell gene regulatory networks (scGRNs) using principal component regression (pcnet), and denoises scGRNs using tensor decomposition (tensordecomp).|

## Example
#### Loading ScTenifoldNet
Once installed, **ScTenifoldNet.jl** can be loaded typing:
```julia
using ScTenifoldNet
```

#### Simulating of a dataset 
Here we simulate a dataset of 2000 cells (columns) and 100 genes (rows) following the negative binomial distribution with high sparsity (~67%).
```julia
d=NegativeBinomial(20,0.98)
X=rand(d,100,2000)
```

#### Generating a perturbed network 
We generate a perturbed network modifying the expression of genes 10, 2, and 3 and replacing them with the expression of genes 50, 11, and 5.
```julia
Y=copy(X)
Y[10,:]=Y[50,:]
Y[2,:]=Y[11,:]
Y[3,:]=Y[5,:]

X=X[:,vec(sum(X,dims=1).>30)]
Y=Y[:,vec(sum(Y,dims=1).>30)]
```
#### ScTenifoldNet
Here we run **ScTenifoldNet** under the H0 (there is no change in the regulation of the gene) using the same matrix as input and under the HA (there is a change in the regulation of the genes) using the control and the perturbed network.
```julia
Z0=ScTenifoldNet.tenrnet(X, donorm=true)
Z1=ScTenifoldNet.tenrnet(Y, donorm=true)
```
#### Differential regulation based on manifold alignment distances
As is shown below, under the H0, none of the genes shown a significative difference in regulatory profiles using an FDR cut-off of 0.1, but under the HA, the 6 genes involved in the perturbation (50, 11, 2, 10, 5, and 3) are identified as perturbed.
```julia
d,aln0,aln1=ScTenifoldNet.manialn(Z0,Z1)
fc,p,adjp=ScTenifoldNet.drgenes(d)
```

#### Plotting the results
Results can be easily displayed using quantile-quantile plots. <br /><br />
![qqplot](https://raw.githubusercontent.com/jamesjcai/ScTenifoldNet.jl/master/qq.png)
```julia
using StatsPlots, Distributions
x=rand(Chisq(1), length(fc))
qqplot(x, fc)
```

## Citation
To cite **ScTenifoldNet.jl** in publications use:

Daniel Osorio, Yan Zhong, Guanxun Li, Jianhua Z. Huang, James J. Cai. 
scTenifoldNet: a machine learning workflow for constructing and comparing transcriptome-wide gene regulatory networks from single-cell data. bioRxiv 2020.02.12.931469; doi: https://doi.org/10.1101/2020.02.12.931469

## Author
James Cai -  @jamesjcai - jcai@tamu.edu
