# MATLAB Tensor Toolbox, Version 3.1

Tensors (also known as multidimensional arrays or N-way arrays) are used in a 
variety of applications ranging from chemometrics to network analysis. 
The Tensor Toolbox provides the following classes for manipulating dense, 
sparse, and structured tensors using MATLAB's object-oriented features:

* `tensor` - Dense tensors, extending MATLAB's native mutlidimensional array capabilities.
* `sptensor` - Sparse tensors, only stores the nonzeros and their indices.
* `symtensor` - Symmetric tensor, only stores the unique entries.
* `ttensor` - Tucker decomposed tensor, stored as a core and factor matrices.
* `ktensor` - Kruskal decomposed tensor, stored as weight and factor matrices.
* `symktensor` - Kruskal decomposed _symmetric_ tensor, stored as weight and factor matrix.
* `sumtensor` - Sum of different types of tensors, never formed explicitly.
* `tenmat` - Tensor as a matrix, with extra information so that it can be converted back into a tensor.
* `sptenmat` - Store an sptensor as a sparse matrix in coordinate format, with extra information so that it can be converted back into an sptensor.

The Tensor Toolbox for MATLAB is open source, but we ask that you please cite 
the appropriate references ([listed below](#how-to-cite)) so that we can continue to show the 
relevance of this work. Your contributions are warmly welcomed as well; please 
see the [contribution guide](CONTRIBUTION_GUIDE.md). 
Previous contributors are listed in [contributors](CONTRIBUTORS.md). 
Full details of the license can be found in [license](LICENSE.txt).

### What's new in Version 3.0?

Version 3.0 adds

* New classes and functions for symmetric tensors: `symtensor`, `symktensor`, `cp_sym`
* New class for sums of different tensor types: `sumtensor`
* Function to compute HOSVD and ST-HOSVD: `hosvd`

We have also fixed many bugs. View the [RELEASE_NOTES.txt](RELEASE_NOTES.txt) file for details.

## How to cite

If you have used the Tensor Toolbox in your work in any way, 
please cite the software itself along with at least one publication or preprint. 
Thanks very much for your support.

__General software reference, should always be cited:__
Brett W. Bader, Tamara G. Kolda and others. 
MATLAB Tensor Toolbox Version 3.1, 
Available online, June 2019. 
URL: https://gitlab.com/tensors/tensor_toolbox. 
_Consider adding the short hash for the exact version that was used. 
If you clone the repository, use the command 
`git log --pretty=format:'%h' -n 1`. 
If you download, the long hash is baked into the filename, but you need only use 
the first 8 characters._ 

``` bibtex
@misc{TTB_Software,
  author = {Brett W. Bader and Tamara G. Kolda and others},
  title = {MATLAB Tensor Toolbox Version 3.0-dev},
  howpublished = {Available online},
  month = aug,
  year = {2017},
  url = {https://gitlab.com/tensors/tensor_toolbox}
}
```

__Default citation for dense computations:__
B. W. Bader and T. G. Kolda.
Algorithm 862: MATLAB tensor classes for fast algorithm prototyping,
ACM Transactions on Mathematical Software 32(4):635-653, December 2006.
DOI: 10.1145/1186785.1186794. 

``` bibtex
@article{TTB_Dense,
  author = {Brett W. Bader and Tamara G. Kolda},
  title = {Algorithm 862: {MATLAB} tensor classes for fast algorithm prototyping},
  journal = {ACM Transactions on Mathematical Software},
  month = dec,
  year = {2006},
  volume = {32},
  number = {4},
  pages = {635--653},
  doi = {10.1145/1186785.1186794}
}
```

__Default citation for sparse computations:__
B. W. Bader and T. G. Kolda.
Efficient MATLAB computations with sparse and factored tensors,
SIAM Journal on Scientific Computing 30(1):205-231, December 2007. DOI: 10.1137/060676489. 

``` bibtex
@article{TTB_Sparse,
  author = {Brett W. Bader and Tamara G. Kolda},
  title = {Efficient {MATLAB} computations with sparse and factored tensors},
  journal = {SIAM Journal on Scientific Computing},
  month = dec,
  year = {2007},
  volume = {30},
  number = {1},
  pages = {205--231},
  doi = {10.1137/060676489}
}
```

__Citation for all-at-once optimization for CP factorization (`cp_opt`):__
E. Acar, D. M. Dunlavy and T. G. Kolda.
A Scalable Optimization Approach for Fitting Canonical Tensor Decompositions,
Journal of Chemometrics 25(2):67-86, February 2011. DOI: 10.1002/cem.1335. 

``` bibtex
@article{TTB_CPOPT,
  author = {Evrim Acar and Daniel M. Dunlavy and Tamara G. Kolda},
  title = {A Scalable Optimization Approach for Fitting Canonical Tensor Decompositions},
  journal = {Journal of Chemometrics},
  month = feb,
  year = {2011},
  volume = {25},
  number = {2},
  pages = {67--86},
  doi = {10.1002/cem.1335}
}
```

__Citation for CP factorization with missing data (`cp_wopt`):__
E. Acar, D. M. Dunlavy, T. G. Kolda and M. Mørup.
Scalable Tensor Factorizations for Incomplete Data,
Chemometrics and Intelligent Laboratory Systems 106(1):41-56, March 2011.
DOI: 10.1016/j.chemolab.2010.08.004. 

``` bibtex
@article{TTB_CPWOPT,
  author = {Evrim Acar and Daniel M. Dunlavy and Tamara G. Kolda and Morten M{\o}rup},
  title = {Scalable Tensor Factorizations for Incomplete Data},
  journal = {Chemometrics and Intelligent Laboratory Systems},
  month = mar,
  year = {2011},
  volume = {106},
  number = {1},
  pages = {41--56},
  doi = {10.1016/j.chemolab.2010.08.004}
}
```

__Citation for Shifted Symmetric Higher-Order Power Method for Tensor Eigenvalues (`eig_sshopm`):__
T. G. Kolda and J. R. Mayo.
Shifted Power Method for Computing Tensor Eigenpairs,
SIAM Journal on Matrix Analysis and Applications 32(4):1095-1124, October 2011. DOI: 10.1137/100801482. 

``` bibtex
@article{TTB_SSHOPM,
  author = {Tamara G. Kolda and Jackson R. Mayo},
  title = {Shifted Power Method for Computing Tensor Eigenpairs},
  journal = {SIAM Journal on Matrix Analysis and Applications},
  month = oct,
  year = {2011},
  volume = {32},
  number = {4},
  pages = {1095-1124},
  doi = {10.1137/100801482}
}
```

__Citation for Method for Generalized Tensor Eigenvalues (`eig_geap`):__
T. G. Kolda and J. R. Mayo.
An Adaptive Shifted Power Method for Computing Generalized Tensor Eigenpairs,
SIAM Journal on Matrix Analysis and Applications 35(4):1563-1581, December 2014.
DOI: 10.1137/100801482. 

``` bibtex
@Article{TTB_EIGGEAP,
  title                    = {An Adaptive Shifted Power Method for Computing Generalized Tensor Eigenpairs},
  author                   = {Tamara G. Kolda and Jackson R. Mayo},
  doi                      = {10.1137/140951758},
  journal                  = {SIAM Journal on Matrix Analysis and Applications},
  number                   = {4},
  volume                   = {35},
  year                     = {2014},
  month                    = dec,
  pages                    = {1563--1581},
  url                      = {http://epubs.siam.org/toc/sjmael/35/4},
}
```

__Citations for on Poisson Tensor Factorization (`cp_apr`)__
1. E. C. Chi and T. G. Kolda.
On Tensors, Sparsity, and Nonnegative Factorizations,
SIAM Journal on Matrix Analysis and Applications 33(4):1272-1299, December 2012.
2. S. Hansen, T. Plantenga and T. G. Kolda.
Newton-Based Optimization for Kullback-Leibler Nonnegative Tensor Factorizations,
Optimization Methods and Software 30(5):1002-1029, April 2015. 

``` bibtex
@Article{TTB_CPAPR,
  title                    = {On Tensors, Sparsity, and Nonnegative Factorizations},
  author                   = {Eric C. Chi and Tamara G. Kolda},
  doi                      = {10.1137/110859063},
  journal                  = {SIAM Journal on Matrix Analysis and Applications},
  number                   = {4},
  volume                   = {33},
  year                     = {2012},
  month                    = dec,
  pages                    = {1272-1299},
}

@Article{TTB_CPAPRB,
  author = {Samantha Hansen and Todd Plantenga and Tamara G. Kolda}, 
  title = {Newton-Based Optimization for {Kullback-Leibler} Nonnegative Tensor Factorizations}, 
  journal = {Optimization Methods and Software}, 
  volume = {30}, 
  number = {5}, 
  pages = {1002-1029},
  month = {April}, 
  year = {2015},
  doi = {10.1080/10556788.2015.1009977},
} 
```
__Citation for Symmetric CP (`cp_sym`):__
T. G. Kolda, 
Numerical Optimization for Symmetric Tensor Decomposition, 
Mathematical Programming B 151(1):225-248, April 2015, doi:10.1007/s10107-015-0895-0
``` bibtex
@article{TTB_CPSYM,  
author = {Tamara G. Kolda}, 
title = {Numerical Optimization for Symmetric Tensor Decomposition}, 
journal = {Mathematical Programming B}, 
volume = {151}, 
number = {1}, 
pages = {225-248}, 
month = apr, 
year = {2015},
doi = {10.1007/s10107-015-0895-0},
}
```

__Citation for CP with randomized least squares (`cp_rals`):__
C. Battaglino, G. Ballard and T. G. Kolda, 
A Practical Randomized CP Tensor Decomposition, 
arXiv:1701.06600, January 2017

```bibtex
@misc{TTB_CPRALS,  
author = {Casey Battaglino and Grey Ballard and Tamara G. Kolda}, 
title = {A Practical Randomized {CP} Tensor Decomposition}, 
month = jan, 
year = {2017},
eprint = {1701.06600},
eprintclass = {cs.NA},
}
```
__Citation for Memory-Efficient Tucker (`tucker_me` and `ttm_me`):__
T. G. Kolda and J. Sun.
Scalable Tensor Decompositions for Multi-aspect Data Mining,
ICDM 2008: Proceedings of the 8th IEEE International Conference on Data Mining,
pp. 363-372, December 2008. DOI: 10.1109/ICDM.2008.89.
_This code is no longer included with the toolbox but can be found in 
Version 2.6._

``` bibtex
@inproceedings{TTB_MET,
  author = {Tamara G. Kolda and Jimeng Sun},
  title = {Scalable Tensor Decompositions for Multi-aspect Data Mining},
  booktitle = {ICDM 2008: Proceedings of the 8th IEEE International Conference on Data Mining},
  month = dec,
  year = {2008},
  pages = {363--372},
  doi = {10.1109/ICDM.2008.89}
}
```


## Getting started and using the software

### Download

The latest release can be obtained by clicking 
[Releases](https://gitlab.com/tensors/tensor_toolbox/releases) at left.
The latest development version can be obtained here by cloning or 
downloading using the buttons above. 
Version 2.6 and earlier can be obtained 
[here](http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html).

### Installation
1. Unpack the files, if necessary
2. Start MATLAB
3. Within MATLAB, navigate to the `tensor_toolbox` directory and execute the following commands:
    1. `addpath(pwd)`
    2. `savepath`

### Getting help
At any time, type `help tensor_toolbox` for help on classes or functions. 
You can also find a getting started guide via MATLAB's help system. Launch help
by pressing the question mark button and look for Tensor Toolbox under supplemental
software, as highlighted in the image below.

![Navigating MATLAB Help Screen](doc/images/helpscreen.PNG "Navigating MATLAB Help Screen")

Copyright 2019, Sandia Corporation.

