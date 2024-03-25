# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 08:43:31 2023

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import cse_functions as cf
import cse_functions_pp as cfp
import cse_functions_attributes as cfa
import xarray as xr
import os
import numpy as np

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths for yaml and land area mask
input_dir = 'multi-model'
ref_dir = 'output_dir_with_indicator_data' 
output_dir = 'multi-model'
yaml_path = 'path_to_folder_with_repo\\required_data\\cse_params.yml'
landmask_path = 'path_to_folder_with_repo\\required_data\\landareamaskmap0.nc'

# Set parameters
ind = 'heatwave'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]

# Set score parameters
score_bins = [3, 2, 1, 0]
stat = 'mean'

if stat != 'mean':
    addon = f'_{stat}'
else:
    addon = ''

# Set interpolation & load landarea mask and parameters
land_mask = cfp.load_landmask(landmask_path)
params = cfp.load_parameters(yaml_path)

#%% Calculate difference and scores
# -----------------------------------------------------------------------------

isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(params['protocol'][ind])
 
# Load future data and only select relevant variable     
future = cfp.load_multi_model_files(input_dir, params['protocol'][ind], 'MM', 
                                    ind, 'future')
future = cfp.apply_land_mask(future, land_mask)
future = future[list(params['indicators'][ind])]

# Load data for reference period and only select relevant variable
ref = cfp.load_multi_model_files(input_dir, params['protocol'][ind], 'MM', 
                                 ind, 'historical')
ref = cfp.apply_land_mask(ref, land_mask)
ref = ref[list(params['indicators'][ind])]
ref = ref.assign_coords({'threshold': GWLs})

# Load wet-bulb temperature data
twb = cfp.load_multi_model_files(input_dir, '3b', 'MM', 'twb', 'future')
twb = cfp.apply_land_mask(twb, land_mask)

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

# Bin wet-blub temperature
future_binned = cfp.bin_future_data(twb.sel({'stats': stat}), {'twb': [0, 18, 25, 30]})

# Calculate standard deviation
std = cfp.calculate_standard_deviations(ind, params['protocol'][ind], land_mask, 
                                        ref_dir, GCMs, '')

# Create datasets to save data in 
future_bivariate = xr.Dataset(data_vars=dict(no_of_events=(["lat", "lon", "threshold", "percentile", 'dt'], 
                                             np.full([len(future.lat),len(future.lon), 
                                             len(GWLs), len(future.percentile),
                                             len(future.dt)], np.nan))), 
                              coords={'lon': future.lon, 'lat': future.lat, 
                                      'threshold': GWLs, 
                                      'percentile': future.percentile, 'dt': future.dt}) 
std_binned = xr.Dataset(data_vars=dict(no_of_events=(["lat", "lon", "threshold", "percentile", 'dt'], 
                                       np.full([len(future.lat),len(future.lon), 
                                       len(GWLs), len(future.percentile),
                                       len(future.dt)], np.nan))), 
                                coords={'lon': future.lon, 'lat': future.lat, 
                                        'threshold': GWLs, 
                                        'percentile': future.percentile, 'dt': future.dt}) 

for p in range(0, len(future.percentile)):
        
    future_data = future.sel({'percentile': future.percentile[p]}).sel({'stats': stat})
    ref_data = ref.sel({'percentile': future.percentile[p]}).sel({'stats': stat})
    std_data = std.sel({'percentile': future.percentile[p]})
    
    for d in range(0, len(future.dt)):
        
        future_sel = future_data.sel({'dt': future.dt[d]})
        ref_sel = ref_data.sel({'dt': future.dt[d]})
        std_sel = std_data.sel({'dt': future.dt[d]})

        std_change = (future_sel - ref_sel) / std_sel
        std_change = xr.where((ref_sel == 0) & (future_sel == 0), 0, std_change)
        std_change_binned = cfp.bin_std_data(std_change, [0, 1, 2, 3])       
        
        std_binned.no_of_events[:,:,:,p,d] = std_change_binned.no_of_events
        future_bivariate.no_of_events[:,:,:,p,d] = std_change_binned.no_of_events + future_binned.twb
        
# Save scores
attributes = cfa.cust_attrs_scores(list(future_bivariate.threshold.values), 
                                   params['protocol'][ind], ind, stat)
output_file = os.path.join(output_dir, 
                           f'ISIMIP{params["protocol"][ind]}_Scores_{ind}{addon}.nc4')
cf.write_output(future_bivariate, output_file, ind, attributes)