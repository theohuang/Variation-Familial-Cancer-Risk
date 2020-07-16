# Variation in Familial Cancer Risk

## Description

This respository contains code to run simulations (both "sample samples" and "large samples") to obtain O/E ratios to analyze variation in familial cancer risk. All code is written by Theodore Huang except where specified.

## Running the code
In order to reproduce the results for the main ("small samples") analysis, run the following four files that correspond to the four types of models:

* OE_Sim_Frailty_CNC_main.R : frailty on both mutation carriers and non-carriers
* OE_Sim_Frailty_C_main.R : frailty on only mutation carriers
* OE_Sim_Frailty_NC_main.R : frailty on only mutation non-carriers
* OE_Sim_Frailty_none_main.R : no frailty

In order to reprodue the results for the "one-family-per-frailty" analysis, run "OE_Sim_Frailty_1fam.R".

In order to reproduce the results for the "large samples" analysis, run the following three files that correspond to the three types of frailty models:

* OE_Sim_Frailty_CNC_LargeSamples.R : frailty on both mutation carriers and non-carriers
* OE_Sim_Frailty_C_LargeSamples.R : frailty on only mutation carriers
* OE_Sim_Frailty_NC_LargeSamples.R : frailty on only mutation non-carriers
    
## Other files are described below:

* OE Functions.R : functions to calculate the O/E ratios
* OE Simulation Frailty Analysis.R : obtaining the main ("small samples") simulation results
* OE Simulation Large Samples Frailty Analysis.R : obtaining the "large samples" simulation results
* oe_sim_cnc.job, oe_sim_c.job, oe_sim_nc.job : files to run the simulations on the Harvard FAS cluster
* pp.peelingParing.R : peeling-paring function to calculate carrier probabilities
* estLik.R : function to evaluate the likelihood
* penet.mmr.net.RData : colorectal and endometrial cancer penetrances from the BayesMendel R package
* death.othercauses.RData : hazard function for death from other causes from the BayesMendel R package
