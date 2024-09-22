# Boreal_flowpaths

_Introduction_
This repository contains data and code for the analyses and figures for the manuscript titled "Permafrost and rain influence hydrologic flowpaths in boreal catchments" by Karen L. Jorgenson, Thomas A. Douglas, M. Torre Jorgenson, Neal J. Pastick, and Tamara K. Harms.

_Raw data files_
All_mixing_data.csv: All data for source-waters and stream water samples included in the mixing models
Precip_all.csv: Precipitation data collected by rain gages
Site_characteristics: Catchment area (km^2), slope (degrees), and the percent of the catchment area covered by south facing slopes, deciduous vegetation, permafrost, bedrock, loess and yedoma

_Output data files_
All_summary_data.csv: Estimated source proportions for source-waters at all sites and years
Cumulative_precip_calc.csv: Cumulative precipitation (cm), cumulative fall precipitation (cm), and the proportion of precipitation that fell in the fall for all sites and years
Model_sources_raw.csv: Chloride and magnesium concentrations for all source-water samples included in the mixing models
Slopes.csv: Linear slopes and intercepts fitted to curves for each source-water over days since snowmelt for each site and year. in the "Par" column, "alph" indicates the intercept for a single model draw, "beta" indicated the slope for a single model draw, "mu_i" indicates the estimate of the average intercept, "mu_s" indicates the estimate of the average slope, "s_in" indicates the error associated with the average intercept, "s_slope" indicated the error associated with the average slope, "s_res" indicates the residual error for the model, and "devi" indicates the model deviance. 


The "MixSIAR_model_scripts" folder contains a Rscript for each model, and the shell scripts used to run them on Chinook, University of Alaska's HPC cluster. 

The "MixSIAR_model_output" folder contains .RData files that contain the output from the mixing models.

_Site abbreviations_
MOOS = Moose Creek
FRCH = French Creek
POKE = Poker Creek
STRT = Stuart Creek
VAUL = Vault Creek


