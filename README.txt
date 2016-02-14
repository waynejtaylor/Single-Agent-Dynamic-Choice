This folder contains three groups of self-contained files:

First Group: Nested Fixed Point and Nested Pseudo Likelihood------------------
1) "Rust Data Generating Process.R” illustrates a simple DGP for the Rust bus dataset.
   It generates the data and estimates the parameters using a nested fixed point algorithm
   The output file “bus_df_in.csv” is used in the Nested Pseudo Likelihood algorithm

2) "Import Data and Estimate" implements Aguirregabiria and Mira’s (2002) NPL algorithm
   It imports “bus_df_in.csv” generated above and calls "npl_sing.R" and "clogit.R"

Second Group: Arcidiacono and Miller (2011)------------------------------------
1) “AM2011Table1cols2356.R” recreates columns 2,3,5 and 6 from Arcidiacono and Miller (2011)
   Data is simulated within the program
   It calls the following support functions:
   -xgrid.R,wlogitd.R,wlogit.R,likebusML4.R,genbus4.cpp,fvdataBOTH.cpp,intcond.R,intcondP.R

Third Group: Bayesian DDC (2009)------------------------------------------------
*) "BayesianDDCEstimateDataRustvEmax.R" uses the output from "Rust Data Generating Process.R” to
replicate the Imai, Jain, and Ching (2009) method, but does not allow for random effects. Use this
file to get a base understanding of the model before proceeding to the hierarchical version.

1) "RustDGPwithHierarchicalRE.R" generates Rust data with a hierarchical mixing

2) "EstimateBayesianHierarchicalDDC.R" follows a method similar to Imai, Jain, and Ching (2009)
for estimation. See the code for options on variations of the estimation procedure.

3) "EstimateBayesianHierarchicalDDC cpp.R" is similar to #2 but uses c++ for speed improvements.
The main c++ program is in "bddcMCMCloop.cpp"

IMPORTANT: this is a stylized example meant to highlight the mechanics of the process. Care must be
taken when selecting the priors and scaling parameters. Here convergence is achieved, but only because
of how I set up the problem.