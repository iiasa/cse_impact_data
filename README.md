# Overview

This repository contains the scripts used to create the data for the [Climate Solutions Explorer](www.climate-solutions-explorer.eu). 

![CSE](https://github.com/iiasa/cse_impact_data/assets/91878469/1e6ea586-3be3-45d5-9928-78024d4b86f5)

The data are publicly available and can be downloaded from [https://zenodo.org/doi/10.5281/zenodo.7971429](https://zenodo.org/doi/10.5281/zenodo.7971429) (v0.4) (Werning et al., 2023). 

# Data and data storage

All scripts in this repository are designed to use ISIMIP2b or ISIMIP3b data as input for the calculation of the indicators. ISIMIP data can be downloaded from https://data.isimip.org/. The scripts for the indicator calculation expect the ISIMIP data to be stored using the following folder structure:

*period/GCM/ or GHM/period/GCM*

Period can be either *historical*, *piControl* or one of the RCPs, e.g. *ssp126/GFDL-ESM4* for ISIMIP3b or *rcp26/GFDL-ESM2M* for ISIMIP2b. Paths to the folders where the ISIMIP data are stored can be set once in cse_functions.py (in the set_protocol function) for both ISIMIP2b and ISIMIP3b and are then used automatically in the indicator calculation scripts (with the exception of the hydrology indicators where the input directory is set manually in the scripts). 

# General notes

Path to folder that contains the repository needs to be set in scripts to make sure that all modules are imported correctly. Required files are stores in *required_files* and paths in the scripts need to be adjusted accordingly.

# Scripts

Scripts are numbered in the order that they should be executed in. The following scripts are available:

### 00_cse_calculate_gcm_gmt: 
- Calculate temperature anomalies for ISIMIP2b/3b protocols using the pre-industrial control data

### 01_cse_calculate_cdd.py 
- Calculate indicator 'Consecutive dry days'

### 01_cse_calculate_drought_intensity.py
- Calculate indicator 'Drought intensity'
- Can be run using either runoff (qtot) or discharge (dis)

### 01_cse_calculate_heatwave.py
- Calculate indicator 'Heatwave events' (all variants)

### 01_cse_calculate_precipitation.py
- Calculate indicators 'Heavy precipitation days', 'Very heavy precipitation days', 'Wet days', 'Very wet days', and 'SDII'

### 01_cse_calculate_sdd.py 
- Calculate indicator 'Cooling degree days'

### 01_cse_calculate_seas_iavar.py 
- Calculate indicators 'Seasonality' and 'Inter-annual variability'
- Can be run using either runoff (qtot) or discharge (dis)

### 01_cse_calculate_tropical_nights.py
- Calculate indicator 'Tropical nights'

### 01_cse_calculate_water_stress_I.py 
- Calculate indicator 'Water stress index'
- First part of calculation
- Written by Yusuke Satoh and adapted by Michaela Werning
- Water demand from ISIMIP3b is used for calculation

### 01_cse_calculate_water_stress_II.py
- Calculate indicator 'Water stress index'
- Second part of calculation

### 02_cse_calculate_multi_model_means.py
- Calculate multi-model ensemble statistics for all indicators

### 03_cse_postprocessing_bivariate_scores.py 
- Calculate relative differences and bivariate scores for all indicators except from heatwave
- Can be run for any of the multi-model ensemble statistics (set in variable *stat*)

### 04_cse_postprocessing_bivariate_scores_heatwave.py
- Calculate relative differences and bivariate scores for indicator 'Heatwave events'
- Can be run for any of the multi-model ensemble statistics (set in variable *stat*)

### 04_cse_postprocessing_split_indicators.py 
- Split multi-model ensemble, difference, and score datasets into separate files 
- Can be run for any of the multi-model ensemble statistics (set in variable *stat*)

### 05_cse_table_output_scaling.py 
- Create the table output for all indicators
- Can be run in two modes: COUNTRIES (almost 200 countries) or R10 (R10 regions & EU)
- Can be run for any of the multi-model ensemble statistics (set in variable *stat*)

# Required files

### cse_params.yaml
- Yaml file with indicator-specific parameters, such as unit, min/max values, ISIMIP protocol, etc. 

### ddm30_flowdir_cru_neva.nc4
- Flow direction used for 'Water stress index' indicator
- File provide by Yusuke Satoh

### grd_ara.hlf
- File with grid area
- File provided by Yusuke Satoh

### gridarea05.nc
- NetCDF file with grid area
- Based on file from Yusuke Satoh, adapted by Edward Byers

### ISIMIP2b_GCM_GMT_1661_2099.xlsx & ISIMIP3b_GCM_GMT_1601_2100.xlsx
- Excel files containing temperature anomalies for all GCM/RCP combinations for ISIMIP2b and ISIMIP3b data
- Created using 00_isimip_calculate_gcm_gmt.py

### kg_class.nc
- NetCDF file with definitions of Koeppen-Geiger regions
- Based on Beck, Hylke E., Niklaus E. Zimmermann, Tim R. McVicar, Noemi Vergopolan, Alexis Berg, and Eric F. Wood. 2018. “Present and Future Köppen-Geiger Climate Classification Maps at 1-Km Resolution.” Scientific Data 5 (1): 180214. https://doi.org/10.1038/sdata.2018.214. 
- Adapted by Daniel Hooke

### landareamaskmap0.nc
- NetCDF file with grid area
- Based on file from Yusuke Satoh, adapted by Edward Byers

### region_classification.xlsx
- Excel file containing region classifications and mapping from countries to R10 regions
- Based on IPCC. 2022a. “Annex II: Definitions, Units and Conventions.” In Climate Change 2022: Mitigation of Climate Change. Contribution of Working Group III to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change, edited by A Al Khourdaje, R. van Dieman, W.F. Lamb, M. Pathak, A. Reisinger, S. de la Rue du Can, J. Skea, R. Slade, S. Some, and L. Steg, 1st ed., 1821–40. Cambridge, UK and New York, NY, USA: Cambridge University Press. https://doi.org/10.1017/9781009157926.021.


