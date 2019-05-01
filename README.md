# Variation in Familial Cancer Risk

## Description

This respository contains code to run simulations to obtain O/E ratios to analyze variation in familial cancer risk. All code is written by Theodore Huang except where specified.

## Running the code
In order to reproduce the results, run the following three files that correspond to the three types of frailty models:

* OE_Sim_Frailty_CNC.R : frailty on both mutation carriers and non-carriers
* OE_Sim_Frailty_C.R : frailty on only mutation carriers
* OE_Sim_Frailty_NC.R : frailty on only mutation non-carriers
    
## The files are all described below:

### General files

* OE Functions.R : functions to calculate the O/E ratios
* OE Simulation Frailty Analysis.R : obtaining the simulation results
* oe_sim_cnc.job, oe_sim_c.job, oe_sim_nc.job : files to run the simulations on the Harvard FAS cluster
* pp.peelingParing.R : peeling-paring function to calculate carrier probabilities
* estLik.R : function to evaluate the likelihood
* penet.mmr.net.RData : colorectal and endometrial cancer penetrances from the BayesMendel R package
* death.othercauses.RData : hazard function for death from other causes from the BayesMendel R package
