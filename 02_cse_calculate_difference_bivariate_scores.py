# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:22:22 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import numpy as np
import xarray as xr
import re
import cse_functions as cf
import cse_functions_pp as cfp
import cse_functions_attributes as cfa

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths for other files
input_dir = 'output_dir_of_indicator_scripts' 
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'
kg_class_path = 'path_to_folder_with_repo\\required_files\\kg_class.nc'
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'

# Set parameters
indicators =  ['cdd', 'dri', 'dri_qtot', 'heatwave', 'iavar', 'iavar_qtot', 
              'precip', 'seas', 'seas_qtot', 'sdd', 
              'sdd_24p0', 'sdd_20p0', 'sdd_18p3', 'tr20', 'wsi']
hist_years = [1974, 2004]
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
SOCs = ''
GHMs = ''

# Load files
land_mask = cfp.load_landmask(land_mask_path)
params = cfp.load_parameters(yaml_path)
kg_class = cfp.load_kg_class(kg_class_path)

#%% Calculate differences and scores
# -----------------------------------------------------------------------------

# Loop through all indicators
for ind in indicators:
    
    # Set protocol, GCMs, RCPs, GHMs, SOCs, etc.
    isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(params['protocol'][ind])        
    if params['protocol'][ind] == '2b':
        GHMs = ['H08', 'LPJmL', 'MATSIRO']
        SOCs = ['rcp26soc', 'rcp60soc']
    
    # Loop through GWLs
    for GWL in GWLs:
        
        # Load all required files
        files = cfp.create_esm_file_list(input_dir, GCMs, RCPs, ind, GWL, GHMs, SOCs)

        # Loop through all files    
        for f in files:
                            
            # Load files for future periods, apply land mask and select variables
            future = cf.load_netcdf(f)
            future = future.where(land_mask['land area'] > 0)
            future = future[list(params['indicators'][ind])]
            attrs = future.attrs                
            
            # Clip values for seasonality and interannual variability
            if ind in ['seas', 'seas_qtot', 'iavar', 'iavar_qtot']:
                future = xr.where(future > 5, np.nan, future)
                future = xr.where(future < 0, np.nan, future)                
            
            # Set path for reference period
            if params['protocol'][ind] == '2b':
                ref_file = re.sub('rcp[0-9]+soc', 'histsoc', f.replace(str(GWL).replace('.', 'p'), 'historical')[:-13] + f'{hist_years[0]}_{hist_years[1]}.nc4')
                ref_all_file = re.sub('rcp[0-9]+soc', 'histsoc', f.replace(str(GWL).replace('.', 'p'), 'historical')[:-20] + f'all_global_{hist_years[0]}_{hist_years[1]}.nc4').replace('rcp60', 'rcp26')      
            else:
                ref_file = f.replace(str(GWL).replace('.', 'p'), 'historical')[:-13] + f'{hist_years[0]}_{hist_years[1]}.nc4'
                ref_all_file =  f.replace(str(GWL).replace('.', 'p'), 'historical')[:-20] + f'all_global_{hist_years[0]}_{hist_years[1]}.nc4'               
            
            # Load files for reference period, apply land mask and select variables
            ref = cf.load_netcdf(ref_file)
            ref = ref.where(land_mask['land area'] > 0)
            ref = ref[list(params['indicators'][ind])]                
            ref_all = cf.load_netcdf(ref_all_file)            
            ref_all = ref_all.where(land_mask['land area'] > 0)
            ref_all = ref_all[list(params['indicators'][ind])]
            
            # Calculate difference and deal with infinite values
            diff = (future - ref) / ref * 100 
            diff = diff.where(diff.apply(np.isfinite)).fillna(np.nan)

            # Calculate standard deviation for all indicators apart from seasonality
            # and interannual variability
            if ind not in ['seas', 'seas_qtot', 'iavar', 'iavar_qtot']:
                std = ref_all.std(dim='year', skipna=True)
            else:
                std = ref_all 
                
            # Calculate scores for heatwave indicator
            if ind == 'heatwave':
                
                # Load twb values to bin future data
                twb = cf.load_netcdf(f.replace('heatwave', 'twb'))
                twb = twb.mean(dim='time')
                twb = twb.where(land_mask['land area'] > 0)                    
                future_binned = cfp.bin_future_data(twb, {'twb': [0, 18, 25, 30]})

                # Calculate relative change and bin it 
                std_change = (future - ref) / std
                std_change = xr.where((ref == 0) & (future == 0), 0, std_change)
                std_change_binned = cfp.bin_std_data(std_change, [0, 1, 2, 3])   
                
                # Calculate scores
                scores = std_change_binned + future_binned.twb
                
                # Set attributes and write out files for each percentile/days combination
                for var in params["indicators"][ind]:
                
                    for p, d in it.product(range(0, len(future.percentile)), range(0, len(future.dt))):
                                                 
                          attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'abs')
                          attributes['quantiles']=str(future.percentile[p].values)
                          attributes['no of days']=str(future.dt[d].values)
                          output_file = f'{f[:-4]}_abs.nc4'.replace(ind, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}')
                          cf.write_output(future[var][:,:,p,d].to_dataset(), output_file, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p]*100))}_{str(future.dt[d])}', '')
                         
                          attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'diff')
                          attributes['quantiles']=str(future.percentile[p].values)
                          attributes['no of days']=str(future.dt[d].values)
                          output_file = f'{f[:-4]}_diff.nc4'.replace(ind, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}')
                          cf.write_output(diff[var][:,:,p,d].to_dataset(), output_file, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p]*100))}_{str(future.dt[d])}', '') 
                         
                          attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'score')
                          attributes['quantiles']=str(future.percentile[p].values)
                          attributes['no of days']=str(future.dt[d].values)
                          output_file = f'{f[:-4]}_score.nc4'.replace(ind, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}')
                          cf.write_output(scores[var][:,:,p,d].to_dataset(), output_file, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}', '')
                         
                          attributes = cfa.amend_esm_attributes(ref.attrs, params["indicators"][ind][var]["long_name"], 'abs')
                          attributes['quantiles:']=str(future.percentile[p].values)
                          attributes['no of days']=str(future.dt[d].values)
                          output_file = f'{ref_file[:-4]}_abs.nc4'.replace(ind, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}')
                          cf.write_output(ref[var][:,:,p,d].to_dataset(), output_file, f'{params["indicators"][ind][var]["short_name"]}_{str(int(future.percentile[p].values*100))}_{str(future.dt[d].values)}', attributes)
                
            else:
                
                # Calculate scores
                scores = cfp.calculate_bivariate_scores(future, ref, std, ind, 
                                                        params['protocol'][ind], 
                                                        kg_class)               
            
                # Set attributes and write out files
                for var in list(future):
                    
                    attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'abs')
                    output_file = f'{f[:-4]}_abs.nc4'.replace(ind, params["indicators"][ind][var]["short_name"])
                    cf.write_output(future[var].to_dataset(), output_file, params["indicators"][ind][var]["short_name"], attributes)
                    
                    attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'diff')
                    output_file = f'{f[:-4]}_diff.nc4'.replace(ind, params["indicators"][ind][var]["short_name"])
                    cf.write_output(diff[var].to_dataset(), output_file, params["indicators"][ind][var]["short_name"], attributes) 
                    
                    attributes = cfa.amend_esm_attributes(attrs, params["indicators"][ind][var]["long_name"], 'score')
                    output_file = f'{f[:-4]}_score.nc4'.replace(ind, params["indicators"][ind][var]["short_name"])
                    cf.write_output(scores[var].to_dataset(), output_file, params["indicators"][ind][var]["short_name"], attributes)
                    
                    attributes = cfa.amend_esm_attributes(ref.attrs, params["indicators"][ind][var]["long_name"], 'abs')
                    output_file = f'{ref_file[:-4]}_abs.nc4'.replace(ind, params["indicators"][ind][var]["short_name"])
                    cf.write_output(ref[var].to_dataset(), output_file, params["indicators"][ind][var]["short_name"], attributes)
                