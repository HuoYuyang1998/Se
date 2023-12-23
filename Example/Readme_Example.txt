*This folder contains codes of an example to implement SCOP (Selection conditional conformal prediction) procedure.


*Code Description:

SCOP+Example.R: The R code of an example to use SCOP and obtain FCP value and length of constructed prediction intervals.

Function_scop.R: The R code for functions used in the simulations.

AlgorithmClass_scop.R: The R code for algorithm (model) classes used in the simulations.


*Description of main function SCOP:

Inputs:
R_cal: n*1 vector of residuals of calibration data
R_test: m*1 vector of responses of test data
T_cal: n*1 vector of selection scores of calibration data
T_test: m*1 vector of selection scores of calibration data
Select_test: the index vector of selected samples in test data
L: threshold for selection
alpha: pre-specified FCR level
Selection_rule: a list containing the type of selection rule and other parameters
exchangeable_adapted: "TRUE" means we will use the Algorithm 1+ in the paper, "FALSE" means we will use the original one, Algorithm 1.

Outputs:
FCP : the false coverage proportion for the constructed prediction intervals
Length : the average length of the constructed prediction intervals
Method : the method name


*Description of other functions in the Function_scop.R:

DataSplit: splitting the generated data into training, calibrating and unlabeled datasets

POS_threshold: computing the threshold by selection rule T-pos(b_0,FDRlevel)

Clustering_threshold: computing the threshold by selection rule T-clu

Selection_thresholding: computing the threshold by assigning a selection rule

OCP: using ordinary conformal prediction to give prediction intervals after selection

ACP:  using adjusted conformal prediction to give prediction intervals after selection. If the parameter "simplified=TRUE", we will use |\hat{\mathcal{S}}_u| as approximation of \mathcal{M}^j_{\min}.
  