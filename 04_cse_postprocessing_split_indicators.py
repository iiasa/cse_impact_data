# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:18:28 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import cse_functions_pp as cfp
import cse_functions as cf
import cse_functions_attributes as cfa
import os
import xarray as xr
import itertools as it

# %% Set up -------------------------------------------------------------------

# Specify input/output directories and path to landmask
input_dir = 'multi-model'
output_dir = 'split_files'
landmask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'

# Customise this section according to needs
indicators = ['cdd']
ssps = ['ssp1', 'ssp2', 'ssp3']
stat = 'mean'

if stat != 'mean':
    addon = f'_{stat}'
else:
    addon = ''

# Abbreviations and descriptions for different file types
types = {'MM': {'short': 'abs', 'desc': 'absolute values'}, 
          'Diff': {'short': 'diff', 'desc': 'difference to reference period'}, 
          'Scores': {'short': 'score', 'desc': 'risk score'}}

# %% Preparation --------------------------------------------------------------

land_mask = cfp.load_landmask(landmask_path)
params = cfp.load_parameters(yaml_path)
    
# %% Computation --------------------------------------------------------------

for ind, ssp, t in it.product(indicators, ssps, types):    
    
    # Make input specifiable
    if t == 'MM':
        input_data = xr.open_dataset(os.path.join(input_dir, 
                           f'ISIMIP{params["protocol"][ind]}_{t}_{ind}.nc4'))
    else: 
        input_data = xr.open_dataset(os.path.join(input_dir, 
                           f'ISIMIP{params["protocol"][ind]}_{t}_{ind}{addon}.nc4'))
    
    if t != 'Scores':
        input_data = input_data.sel({'stats': stat})
    
    # Apply landmask again in case it hasn't been applied yet
    input_data = cfp.apply_land_mask(input_data, land_mask)

    for var, threshold in it.product(params['indicators'][ind], input_data.threshold.values):
                                          
        if ind == 'heatwave':
            
            for d, q in it.product(input_data.dt.values, input_data.percentile.values):
                
                var_name = f'{params["indicators"][ind][var]["short_name"]}_{int(q*100)}_{d}'                
                data = input_data[var].sel({'threshold': threshold, 'dt': d, 'percentile': q}).to_dataset(
                        name=f'{var_name}')
                
                data[var_name].attrs = cfa.split_ind_attrs(
                    ind, var, var_name, ssp, params, types, t, threshold, q, d)
                
                # Write output
                output_file = os.path.join(output_dir, addon[1:], ind, 
                          f'ISIMIP{params["protocol"][ind]}_{var_name}_{types[t]["short"]}.nc4')
                
                if not os.path.exists(os.path.join(output_dir, addon[1:], ind)):
                    os.makedirs(os.path.join(output_dir, addon[1:], ind))
                
                cf.write_output(data, output_file, var)
              
        else: 
            
            var_name = f'{params["indicators"][ind][var]["short_name"]}'                 
            data = input_data[var].sel({'threshold': threshold}).to_dataset(
                    name=f'{var_name}')
            
            if ind in ['dri', 'dri_qtot']:
                data = data.drop('quantile')

            # Create new dataset attributes
            data[var_name].attrs = cfa.split_ind_attrs(
                ind, var, var_name, ssp, params, types, t, threshold)
            
            # Write output
            output_file = os.path.join(output_dir, addon[1:], ind,  
                      f'ISIMIP{params["protocol"][ind]}_{var_name}_{types[t]["short"]}.nc4')
            
            if not os.path.exists(os.path.join(output_dir, addon[1:], ind)):
                os.makedirs(os.path.join(output_dir, addon[1:], ind))
                
            cf.write_output(data, output_file, var)

