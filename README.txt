This folder contains two groups of self-contained files:

First Group: Nested Fixed Point and Nested Pseudo Likelihood—————————————
1) "Rust Data Generating Process.R” illustrates a simple DGP for the Rust bus dataset.
   It generates the data and estimates the parameters using a nested fixed point algorithm
   The output file “bus_df_in.csv” is used in the Nested Pseudo Likelihood algorithm

2) "Import Data and Estimate" implements Aguirregabiria and Mira’s (2002) NPL algorithm
   It imports “bus_df_in.csv” generated above and calls "npl_sing.R" and "clogit.R"

Second Group: Arcidiacono and Miller (2011)—————————————————————————————-
1) “AM2011Table1cols2356.R” recreates columns 2,3,5 and 6 from Arcidiacono and Miller (2011)
   Data is simulated within the program
   It calls the following support functions:
   -xgrid.R,wlogitd.R,wlogit.R,likebusML4.R,genbus4.cpp,fvdataBOTH.cpp,intcond.R,intcondP.R