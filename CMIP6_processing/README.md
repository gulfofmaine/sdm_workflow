
# CMIP6 Processing

**About:**

This section of the repository contains code detailing how the CMIP6 data was processed to 
prepare any bias-corrections. This readme exists to orient users through what files are connected 
and need to be run to maintain upstream-downstream connections.

## Support Function Sources:

For both python and R Workflows the data processing steps were broken down into discrete steps
and written as individual functions. 

R functions for processing bias correction steps can be found [here](https://github.com/gulfofmaine/sdm_workflow/blob/main/CMIP6_processing/R/sdm_workflow_funs.R)

Python functions used to download and process cmip data can be found [here](https://github.com/gulfofmaine/sdm_workflow/blob/main/CMIP6_processing/fcts.py)


## Processing Order:

**Downloading CMIP6 Data**

CMIP6 data was downloaded using climate data operator (CDO) command line tools using the following 
specifications:

**Downloading Reference Data for Bias Corrections**

 1. OISSTv2
 
 OISSTv2 Data was downloaded as annual NetCDF files from the [Physical Science Laboratory](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html).
 
 2. SODA
 
 SODA data was accessed via
 
 **Preparing Date Keys for CMIP6 Models**

For some reason the R packages that support NetCDF file handling do not like the calendar structure
used by climate models. To still use R as a processing tool the dates for individual. To prepare
lookup tables that include the source model, the variable of interest, and the proper dates that 
go with them we used a python3 script [CMIP6 Variable Date Keys](https://github.com/gulfofmaine/sdm_workflow/blob/main/CMIP6_processing/date_key.py). This file 
was stepped for each of the four variables to create lookup tables for dates and also to flag 
models that had inconsistent structures/dates/variables.

**Preparing Climatologies**

For all data sources the reference period for the climatologies was 1985-2014. This
period was chosen because it is the last 30-year period available in the CMIP6 historical
period of 1950-2014.

The Daily Climatology for OISSTv2 data was performed separately as part of the [OISST Mainstays Workflow](https://github.com/adamkemberling/oisst_mainstays) using the 1985-2014 reference period. 
This was then converted to a monthly climatology as part of the processing for bias-correction.

SODA monthly climatologies were done for bottom temperature,  bottom salinity, surface temperature, 
and surface salinity in the python3 script: [soda climatologies](https://github.com/gulfofmaine/sdm_workflow/blob/main/CMIP6_processing/soda_climatology.py).

All CMIP6 climatologies were done using historical model runs, using the same reference period of 
1985-2014. These were done for each variable as part of the two rmarkdown workflows:
[OISSTv2 Bias Corrections](https://gulfofmaine.github.io/sdm_workflow/CMIP6_processing/R/CMIP_OISST_bias_corrections.html) & [SODA Bias Corrections](https://gulfofmaine.github.io/sdm_workflow/CMIP6_processing/R/CMIP_SODA_bias_corrections.html)

**Performing Bias Corrections**

Bias-corrections were performed by stepping through one of two rmarkdowns. One for SODA bias-corrections
and the other for OISSTv2 bias-corrections. The SODA rmarkdown has a parameter option that lets you 
select which variable to process, outputting everything according to the selection.

Details on individual steps can be found in the Rmarkdown report links:

[OISSTv2 Bias Corrections](https://gulfofmaine.github.io/sdm_workflow/CMIP6_processing/R/CMIP_OISST_bias_corrections.html) & [SODA Bias Corrections](https://gulfofmaine.github.io/sdm_workflow/CMIP6_processing/R/CMIP_SODA_bias_corrections.html)
