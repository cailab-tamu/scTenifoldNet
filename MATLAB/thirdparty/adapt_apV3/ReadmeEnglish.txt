Help file of Adaptive Affinity Propagation clustering (Version 3.0)

  Affinity propagation clustering (AP) is a clustering algorithm proposed in "Brendan J. Frey and Delbert Dueck. Clustering by Passing Messages Between Data Points. Science, 2007, 315(5814), 972-976". It has some advantages: speed, general applicability, and suitable for large number of clusters. AP has two limitations: it is hard to known what value of parameter ¡®preference¡¯ can yield optimal clustering solutions, and oscillations cannot be eliminated automatically if occur.

  Adaptive AP improves AP in these items: adaptive adjustment of the damping factor to eliminate oscillations (called adaptive damping), adaptive escaping oscillations, and adaptive searching the space of preference parameter to find out the optimal clustering solution suitable to a data set (called adaptive preference scanning). 

  Your tests and comments are helpful and welcome (E-mail: wangkjun@yahoo.com). The citation information is:
K. Wang, J. Zhang, D. Li, X. Zhang and T. Guo. Adaptive Affinity Propagation Clustering. Acta Automatica Sinica, 33(12):1242-1246, 2007. (in Chinese)
(The English version of this paper is available from: http://arxiv.org/abs/0805.1096)

  The codes of Adaptive AP are available from: 
http://www.mathworks.com/matlabcentral/fileexchange/18244


(1) Contents of "Main_adaptAP_demo.m"
  This program includes six function parts: 
Part 1: Selecting a data set & initialization. 
Part 2: Loading a data file or similarity matrix by "data_load.m"
Part 3: Running Adaptive AP by "adapt_apcluster.m", or original AP by "adapt_apcluster.m"
Part 4: Computing Silhouette indices for clustering results by "solution_evaluation.m"
Part 5: Finding an optimal clustering solution based on Silhouette indices by "solution_findK.m"
Part 6: Computing Fowlkes-Mallows validity index (compared with true answers) by "valid_external.m"
           Computing error rates of a clustering solution (compared with true answers) by "valid_errorate.m"

Note 1: The programs are tested under Matlab 7.2.
Note 2: Statistics Toolbox of Matlab needs to be installed, since it contains routines such as Silhouette index.

(2)  Pearson similarity/distance
  Pearson similarity R(i,j) is the linear correlation coefficient between two vectors and has its value range from -1 to 1, and it is commonly used to measure the similarity/distances between genes. 
  As negative value of R(i,j) in [-1,1] will lead to confusion of energy computation, in actual computation the R(i,j) is transformed to positive P(i,j) in [0,1], i.e., P(i,j)=0.5-0.5R(i,j), and then P(i,j) gives a positive value between two elements. It is important to remark that P(i,j) is adopted only for the convenience of energy computation (but not a new measure), and it is easy to go back: R(i,j)=1-2P(i,j).
   Euclidean distances (type = 1) or Pearson similarity (type = 2) as the similarity measure may be selected in initialization part from row 62.

(3) Input: a data file like "yourdata.txt"
  The input data file is the tab delimited text file with numeric tabular data or similar Matlab file format (e.g. rows denote data points/elements and columns denote dimensions). Please note that all the data should be numeric values and without missing values. Different data files should be loaded in different rows of the programs (row 19 - row 58).
 
(4) Output
  The class labels (clustering solutions) at every number of clusters (NC) are stored in variable "labels" corresponding to every NC in variable "NCs". The optimal NC is found and stored in variable "NCopt", and the optimal clustering solution is the class labels in variable "labels" corresponding to NCopt.

  The Silhouette indices corresponding to every NC in variable "NCs" are stored in variable "Sil", and in "Silmin" for the minimal value of Silhouette indices of any pair of clusters.
  
(5) Parameters
  nrun --- max iteration times of the algorithm, default 50000 is enough for most cases.
  nconv --- convergence condition of the algorithm, default 50. The condition is more strict if bigger.  
  pstep --- decreasing step of parameter "preferences": pstep*pmedian, default 0.01. The searching of "preferences" is finer and the algorithm runs slower if pstep is smaller.
  lam --- initial damping factor, default 0.5.
  cut --- after clustering, drop an cluster with number of samples < cut.
  splot --- observing clustering processes of AP when it is on: = 'plot'.

(6) Error rates compared with true labels by "valid_errorate.m"
   The true class labels in 1st column of a data file rank in nature order (e.g. [1 1 1 ... 1 2 2 2 ... 2 3 3 3 ...3] for 3 classes). The error rate will be accurate and useful if  the clustering solution under true NC is relatively accurate (e.g. error rate < 20%), otherwise it might be inaccurate (since "valid_errorate.m" designed here can not deal with complex cases). In inaccurate case, Fowlkes-Mallows validity index (FM) is recommended.

  FM with values in [0,1] measures agreement between clustering solutions and true class labels, and bigger the value better the agreement, e.g., FM value is 1 when the solution and true labels are the same. FM index is used to evaluate clustering quality when the NC of a clustering solution is different from true NC.

(7) Demo running of "Main_adaptAP_demo.m"
- starting the demo programs of "Main_adaptAP_demo.m" under default parameters, and set id =5;
- the data file of 'ionosphere.txt' is loaded by "data_load.m"
- Adaptive AP is running by "adapt_apcluster.m", and gives class labels in variable
   "labels" (n rows x d columns) correspond to NC in variable "NCs" (2, 3, ...).
  This means that n class labels in the first column are the clustering solution at 
  NC=2, and so on.
- Silhouette indices are computed by "solution_evaluation.m", which gives Silhouette
   indices and the minimal value of Silhouette indices for every clustering solution.
- the optimal NC=2 (optimal clustering solution) is found based on Silhouette indices 
  by "solution_findK.m", i.e., the Sil at NC=2 is the maximum among all the Sil values.
- Fowlkes-Mallows index and error rate are computed by "valid_external.m" &
   "valid_errorate.m".

(8) Demo data sets
  The following data sets are included:
data sets                  NC	Number of samples	Dimension
- Ionosphere	2	351                  	4
- Wine                        3	178                	3

Note 1. The following data sets are from "Brendan J. Frey and Delbert Dueck. Clustering
by Passing Messages Between Data Points. Science, 2007, 315(5814), 972-976.", and are
available from: http://www.psi.toronto.edu/affinitypropagation/
FaceClusteringSimilarities.txt
DocumentSummarization.mat
TravelRouting.mat
GeneFindingProblem.mat

Note 2. The pre-processing methods of some data sets used in this program are in the following:
# dataset 14k10close is from dataset 14k10 added with the random Gaussian noise of variance 1.8.
# dataset 22k10far consists of dataset 14k10 and the additional data of its 8 adjusted clusters:
    Let "data" be dataset 14k10, then dataset 22k10far is [data; data(1:40,:)+2; data(191:210,:)-4; data(211:260,:)-2;data(41:90,:)-3; data(91:140,:)+3; data(301:350,:)+3; data(351:400,:)-3].
# dataset Ionosphere here is obtained by feature selection, i.e. it comprises the second, 4th, 6th and 8th columns of original dataset (the class labels are stored in the first column).
# dataset Wine here is obtained by feature selection, i.e. it comprises the second, 8th and 13th columns of original dataset (the class labels are stored in the first column).
# dataset NCI60 here is obtained by removing two samples of class PROSTATE and then selecting 20 genes from 6830 genes according to their good distinguishable effects on 8 classes. Thus, the dataset has 58  20-dimensional samples from 8 classes.

---------------------------------------------------------------------------------------------------
This software is distributed under the BSD license. 
Copyright (C) 2007-2008.
Last modified: July 26, 2009