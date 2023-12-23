*This folder contains codes for all simulation results presented in the article and the supplementary material.


*Code Description:

Simulation_Figure 1.R: The R code for reproducing Figure 1. The replications are implemented by parallelization.

Simulation_Table 1.R: The R code for reproducing Table 1. The replications are implemented by parallelization.

Simulation_Figure S1.R: The R code for reproducing Figure S1 of the density plots for residuals. 

Simulation_Figure S2.R: The R code for reproducing Figure S2 of prediction interval plots.

Simulation_Figure S3.R: The R code for reproducing Figure S3 which considers LORD-CI. The replications are implemented by parallelization.

Simulation_Table S1.R: The R code for reproducing Table S1 which utilizes quantile regression nonconformity scores. The replications are implemented by parallelization.

Simulation_Table Supp.R: The R code for reproducing other tables in Supplementary Material. The replications are implemented by parallelization.

Function_scop.R: The R code for functions used in the simulations.

AlgorithmClass_scop.R: The R code for algorithm (model) classes used in the simulations.


*Guidelines:

To reproduce the figures and tables, just use the default parameter settings of the corresponding parts.

To reproduce other tables in the supplementary material, just use Simulation_Table Supp.R and change:

1. Table S2 and S3: the DataGen function for each algorithm
2. Table S4: the parameter "exchangeable_adapted" for function SCOP
3. Table S5: the value of the variable "n" and "m" for different tuples (n,m). Assign "T-clu" to "Thresholding_type_array"  and "RF" to "algo_array".
