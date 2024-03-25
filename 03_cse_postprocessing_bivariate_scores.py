# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:52:11 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import cse_functions_pp as cfp
import cse_functions as cf
import cse_functions_attributes as cfa
import xarray as xr
import os

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths for yaml and land area mask
input_dir = 'multi-model'
ref_dir = 'output_dir_with_indicator_data'
output_dir = ''
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'
landmask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'
kg_class_path = 'path_to_folder_with_repo\\required_files\\kg_class.nc'

# Set parameters
GHMs = ['H08', 'LPJmL', 'MATSIRO']
indicators = ['cdd', 'dri', 'dri_qtot', 'iavar', 'iavar_qtot', \
              'precip', 'seas', 'seas_qtot', 'sdd', \
              'sdd_24p0', 'sdd_20p0', 'sdd_18p3', 'tr20', 'wsi']
indicators = ['cdd']
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]

# Set score parameters and multi-model ensemble statistic to use
score_bins = [3, 2, 1, 0]
stat = 'mean'

if stat != 'mean':
    addon = f'_{stat}'
else:
    addon = ''

# Set interpolation & load landarea mask and parameters
land_mask = cfp.load_landmask(landmask_path)
params = cfp.load_parameters(yaml_path)
with xr.open_dataset(kg_class_path, engine="netcdf4") as kg_class:
    kg_class.load()  
kg_class = kg_class.sel(band=0)

#%% Calculate difference and scores
# -----------------------------------------------------------------------------

for ind in indicators:
            
    isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(params['protocol'][ind])

    # Load future data and only select relevant variable
    future = cfp.load_multi_model_files(input_dir, params['protocol'][ind], 
                                        'MM', ind, 'future')
    future = cfp.apply_land_mask(future, land_mask)
    future = future[list(params['indicators'][ind])]
    
    # Load data for reference period and only select relevant variable
    ref = cfp.load_multi_model_files(input_dir, params['protocol'][ind], 'MM', 
                                     ind, 'historical')
    ref = cfp.apply_land_mask(ref, land_mask)
    ref = ref[list(params['indicators'][ind])]
    ref = ref.assign_coords({'threshold': GWLs})
    
    # Calculate relative difference between GWLs and reference period
    diff = cfp.calculate_differences(future, ref, stat, ind, 
                                     params['indicators'][ind].keys(), 
                                     params['indicators'][ind])    
    
    # Save difference data
    attributes = cfa.cust_attrs_diff(list(diff.threshold.values), 
                                     params['protocol'][ind], ind, stat)
    output_file = os.path.join(output_dir, 
                               f'ISIMIP{params["protocol"][ind]}_Diff_{ind}{addon}.nc4')
    cf.write_output(diff, output_file, ind, attributes)  
    
    # Calculate standard deviation
    std = cfp.calculate_standard_deviations(ind, params['protocol'][ind], 
                                            land_mask, ref_dir, GCMs, GHMs)
    
    # Calculate scores
    scores = cfp.calculate_bivariate_scores(future.sel({'stats': stat}), 
                                            ref.sel({'stats': stat}), std, 
                                            ind, params['protocol'][ind], kg_class)
    
    # Save scores
    attributes = cfa.cust_attrs_scores(list(scores.threshold.values), 
                                       params['protocol'][ind], ind, stat)
    output_file = os.path.join(output_dir, 
                               f'ISIMIP{params["protocol"][ind]}_Scores_{ind}{addon}.nc4')
    cf.write_output(scores, output_file, ind, attributes)
