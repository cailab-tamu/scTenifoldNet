自适应仿射传播聚类（版本3.0） 帮助文件

  仿射传播聚类（Affinity propagation clustering, AP）是在Science上提出的一种新聚类算法，参见："Brendan J. Frey and Delbert Dueck. Clustering by Passing Messages Between Data Points. Science, 2007, 315(5814), 972-976."。其优势体现在处理类数很多的情况时运算速度快和适合大类数。仿射传播聚类有两个尚未解决的问题：一是很难确定偏向参数取何值能够使算法产生最优的聚类结果，另一是当震荡发生后算法不能自动消除震荡并收敛。
  自适应仿射传播聚类（adAP）方法为了解决这两个问题，设计了扫描偏向参数空间来搜索聚类个数空间以寻找最优聚类结果（称为自适应扫描）、调整阻尼因子来消除震荡（称为自适应阻尼）以及降低参数preference的值以逃离震荡（称为自适应逃离）。

  欢迎使用和评论自适应仿射传播聚类方法，您的意见是对我们工作的支持（E-mail: wangkjun@yahoo.com）。引用信息是：
  王开军,张军英,李丹,张新娜,郭涛.自适应仿射传播聚类. 自动化学报, 33(12):1242-1246, 2007.
  （此文章可由自动化学报的网站获得: http://www.aas.net.cn/）

  本程序的下载网址是：
http://www.mathworks.com/matlabcentral/fileexchange/18244


(1) 主文件"Main_adaptAP_demo.m"的内容
第1部分:  选择数据集、初始化 
第2部分:  加载数据或相似度矩阵文件（由程序文件"data_load.m"执行）
第3部分:  运行自适应AP（由程序文件"adapt_apcluster.m"执行），AP则是"adapt_apcluster.m"
第4部分:  计算聚类结果的Silhouette指标值（由程序文件"solution_evaluation.m"执行）
第5部分:  由Silhouette指标值找出最优聚类结果（由程序文件"solution_findK.m"执行）
第6部分:  计算聚类结果的Fowlkes-Mallows指标（由程序文件"valid_external.m"执行）
                 和错误率（由程序文件"valid_errorate.m"执行）

注1: 本程序在 Matlab 7.2下经过测试。
注2: 本程序的运行需要安装Matlab的统计工具箱（Statistics Toolbox）, 因需要Silhouette指标等的子程序 .

(2)  Pearson相似度和距离
    基因表达数据聚类分析中基因之间的相似性测度定义为Pearson相关系数，即两个样本/向量i和j之间的线性相关系数R(i,j)，其取值范围[-1,1]。
    为了便于各指标的计算，将R(i,j)进行这样的转换：R(i,j)=0.5-0.5R(i,j)，从而R(i,j)取值范围[0,1]。这样R(i,j)表示正的Pearson距离和不相似度，而相似度则为1-R(i,j)。 例如最远Pearson距离的两个基因 g1和 g2有R(g1,g2)=1 以及自身最近R(g1,g1)=0 。
    选择欧式距离（type = 1）或Pearson相似度（type = 2）作为计算相似度的测度是由初始设置（第62行开始的 initialization）确定的。

(3) 输入：数据文件 "yourdata.txt"
    文件的每一行表示一个数据样本或基因，列表示维数，数据应为数值并无缺失数据。若类标已知，类标存放于第一列；若类标未知，第一列则放数据。不同的数据文件应放置在程序中的不同位置，请按照程序中的放置说明调用数据文件。
 
(4) 输出
  对应于不同类数的聚类结果（类标）存放在变量"labels"中，而对应的类数存放在变量"NCs"中。若求得最优类数（存放于变量"NCopt"），则最优聚类结果是"labels"中对应于最优类数NCopt的那一列类标。
  对应于变量"NCs"中不同类数（下聚类结果）的Silhouette指标存放在变量"Sil"中，每一个聚类结果中任两个聚类的Silhouette指标中最小值存放在变量"Silmin"中。
  
(5) 参数
  nrun --- 算法的最大循环次数，默认值 50000对大多数情况已足够。
  nconv --- 达到收敛时的循环次数，默认值 50。值越大，收敛条件越严格。
  pstep --- 参数"preferences"的下降步长: pstep*pmedian，默认值 0.01。值越小，搜索越细致、耗时越多。
  lam --- 初始的阻尼系数，默认值0.5。
  cut --- 若一个聚类包含样本数目少于cut，则忽略此聚类。
  splot --- 启用状态（= 'plot'）时可观察算法的运行过程。

(6) 与正确类标相比，聚类结果的错误率（由"valid_errorate.m"计算）
   放置于数据文件第一列的类标是用自然数标识的，例如3个聚类的类标： [1 1 1 ... 1 2 2 2 ... 2 3 3 3 ...3] 。如果聚类结果的错误率>20%，错误率指标可能不准确。因为程序"valid_errorate.m"比较简单，不能处理复杂的情况。此时可以采用Fowlkes-Mallows指标。 
  FM表示测量聚类结果与正确类标一致性的 Fowlkes-Mallows指标（外部有效性指标）。FM的值处于0与1之间且越大表示一致性越好，例如当聚类结果与正确类标完全一致时FM=1。当聚类结果的类数与正确类数不同时，可采用FM指标对聚类结果的质量进行评价。

(7) "Main_adaptAP_demo.m"的示例
- 使用默认参数，选id = 5，启动 "Main_adaptAP_demo.m"
- "data_load.m"加载数据文件'ionosphere.txt'
- "adapt_apcluster.m"运行自适应AP算法，输出聚类结果（类标）并存放于变量
   "labels" (n行 x d列)中，对应变量"NCs" (2, 3, ...)中的类数。这意味着
  "labels"中第一列的n个类标是类数为2时的聚类结果，依次类推。
-  由"solution_evaluation.m"计算每个聚类结果的Silhouette指标和所有任两个聚类的
   Silhouette指标中的最小值。 
- "solution_findK.m"依据Silhouette指标值找出最优类数NC=2 (对应最优聚类结果) 。
  也就是类数NC=2时Sil是所有Silhouette值中的最大值。
- 由"valid_external.m" &  "valid_errorate.m"计算出Fowlkes-Mallows值FM和聚类
  结果的错误率。

(8) 示例数据集 (数据文件的第一列为正确的类标，其余列为数据)
  数据集                  类数	样本数	维数
 Ionosphere	2	351           4
 Wine                           3	178           3

注1. 如下4个数据集：
  FaceClusteringSimilarities.txt
  DocumentSummarization.mat
  TravelRouting.mat
  GeneFindingProblem.mat
来自于 "Brendan J. Frey and Delbert Dueck. Clustering by Passing Messages Between Data Points. Science, 2007, 315(5814), 972-976."。http://www.psi.toronto.edu/affinitypropagation/

注2. 本程序中使用的数据集的预处理方法如下：
＃数据集14k10close是将数据集14k10加上方差为1.8的随机高斯噪声后产生的。
＃数据集22k10far是由数据集14k10再与其中8个聚类的调整数据合并产生的：
    设数据集14k10为data，则数据集22k10far为[data; data(1:40,:)+2; data(191:210,:)-4; data(211:260,:)-2;data(41:90,:)-3; data(91:140,:)+3; data(301:350,:)+3; data(351:400,:)-3]。
＃这里的数据集Ionosphere是经过特征选择后的数据，即由原数据集的第2、4、6、8列（第1列为类标）组成。
＃这里的数据集Wine是经过特征选择后的数据，即由原数据集的第2、8、13列（第1列为类标）组成。
＃这里的数据集NCI60是去掉了2个样本的PROSTATE类后剩下8个类的58个样本，同时从
6830个基因中选出对区分这8个类效果好的20个基因，则数据集成为20维的数据。

---------------------------------------------------------------------------------------------------
本程序在BSD许可证下分发
Copyright (C) 2007-2008.
最后修改: 2009.7.26