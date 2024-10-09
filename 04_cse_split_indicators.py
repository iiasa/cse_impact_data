# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:18:28 2022

@author: werning
"""

import sys
sys.path.append('H:\\git\\climate_impacts_processing')
import os
import xarray as xr
import itertools as it
import cse_functions as cf
import cse_functions_pp as cfp
import cse_functions_attributes as cfa

# %% Set up
# -----------------------------------------------------------------------------

# Specify input/output directories and path to landmask
input_dir = 'multi_model'
output_dir = 'split_files'
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'
land_cover_path = 'path_to_folder_with_repo\\required_files\\landcovermasks.nc'

# Set parameters
indicators = ['cdd', 'dri', 'dri_qtot', 'heatwave', 'iavar', 'iavar_qtot', 
              'precip', 'seas', 'seas_qtot', 'sdd', 
              'sdd_24p0', 'sdd_20p0', 'sdd_18p3', 'tr20', 'wsi']
stat = 'median'
ssps = ['ssp1', 'ssp2', 'ssp3']
types = {'MM': {'short': 'abs', 'desc': 'absolute values'}, 
          'Diff': {'short': 'diff', 'desc': 'difference to historical'}, 
          'Scores': {'short': 'score', 'desc': 'risk score'}}
percentiles = [95, 97, 99]
dt = [3, 5, 7, 10]

# Load required files
land_mask = cfp.load_landmask(land_mask_path)
landcover_mask = cfp.load_landcover_mask(land_cover_path)
params = cfp.load_parameters(yaml_path)
    
# %% Split files
# -----------------------------------------------------------------------------

# Loop through indicators, ssps and file types
for ind, ssp, t in it.product(indicators, ssps, types):    
    
    # Loop through variables
    for var in params['indicators'][ind]:  
        
        # Split files for heatwave indicator
        if ind == 'heatwave':
                
            for p, d in it.product(percentiles, dt):
                
                # Create variable name based on percentiles and days over threshold
                var_name = f'{params["indicators"][ind][var]["short_name"]}_{p}_{d}'
                
                # Load multi-model input data, apply land mask and select stat
                input_data = xr.open_dataset(os.path.join(input_dir, 
                                   f'ISIMIP{params["protocol"][ind]}_{t}_{var_name}.nc4'))            
                input_data = cfp.apply_land_mask(input_data, land_mask)
                input_data = input_data.sel({'stats': stat})
                
                # For every GWL, save data as separate file
                for GWL in input_data.threshold.values:
                
                    data = input_data[var].sel({'threshold': GWL}).to_dataset(
                            name=f'{var_name}')        
                    data[f'{var_name}'].attrs = cfa.split_ind_attrs(
                        ind, var, var_name, ssp, params, types, t, GWL)                    
                    output_file = os.path.join(output_dir, ind, stat,   
                              f'cse_{var_name}_{ssp}_{(str(GWL)).replace(".", "p")}_{types[t]["short"]}.nc4')
                    cf.write_output(data, output_file, var)
        
        # Split files for other indicators
        else:
        
            # Load multi-model input data, apply land mask and select stat
            input_data = xr.open_dataset(os.path.join(input_dir, 
                               f'ISIMIP{params["protocol"][ind]}_{t}_{params["indicators"][ind][var]["short_name"]}.nc4'))            
            input_data = cfp.apply_land_mask(input_data, land_mask)
            input_data = input_data.sel({'stats': stat})
            
            # For every GWL, save data as separate file
            for GWL in input_data.threshold.values:
                
                # Mask out Greenland ice sheet and desert areas                    
                if ind in ['dri', 'dri_qtot', 'iavar', 'iavar_qtot', 'seas', 'seas_qtot', 'wsi', 'wsi_qtot']:                
                    input_data = input_data.where(landcover_mask.grid_index < 11)
                
                var_name = f'{params["indicators"][ind][var]["short_name"]}'                 
                data = input_data[var].sel({'threshold': GWL}).to_dataset(
                        name=f'{var_name}')
                
                if ind in ['dri', 'dri_qtot']:
                    data = data.drop('quantile')
    
                data[f'{var_name}'].attrs = cfa.split_ind_attrs(
                    ind, var, var_name, ssp, params, types, t, GWL)
                output_file = os.path.join(output_dir, ind, stat,   
                          f'cse_{var_name}_{ssp}_{(str(GWL)).replace(".", "p")}_{types[t]["short"]}.nc4')
                cf.write_output(data, output_file, var)

