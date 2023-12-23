*Data Description:

Davis3.txt: The original DAVIS dataset for predicting drug-target pair's binding affinity. It contains 25, 772 drug-target pairs. Each pair includes the binding affinity, the structural information of the drug compound, and the amino acid sequence of the target protein.


*Code Description:
Drug discovery main.py: The python code for processing the DAVIS data. It requires the DeepPurpose library to encode the structural information. The results for SCOP and other two methods are saved in file plotdata.csv.

Drug discovery plot.R: The R code for reproducing Figure 2. It uses the data in plotdata.csv for plotting. 

plotdata.csv: The processed results by Drug discovery main.py, containing the FCP and length for each method, selection rule and each repetition.





