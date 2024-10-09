# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 14:29:15 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import os
import xarray as xr
import cse_functions_pp as cfp
import cse_functions as cf
import cse_functions_attributes as cfa

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths to required files
input_dir = 'output_dir_with_indicator_data'
output_dir = 'multi_model'
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'

# Set parameters
indicators = ['cdd', 'dri', 'dri_qtot', 'heatwave', 'iavar', 'iavar_qtot', 
              'precip', 'seas', 'seas_qtot', 'sdd', 
              'sdd_24p0', 'sdd_20p0', 'sdd_18p3', 'tr20', 'wsi']
thresholds = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5] 
GHMs = ''
SOCs = ''
percentiles = [95, 97, 99]
dt = [3, 5, 7, 10]
types = {'MM': {'short': 'abs', 'desc': 'absolute values'}, 
          'Diff': {'short': 'diff', 'desc': 'difference to historical'}, 
          'Scores': {'short': 'score', 'desc': 'risk score'}}

# Load required files
params = cfp.load_parameters(yaml_path)

#%% Calculate multi-model mean
#------------------------------------------------------------------------------    

data = xr.Dataset()

# Loop through indicators
for ind in indicators:    
    
    # Set protocol, GCMs, RCPs, GHMs, SOCs, etc.
    isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(params['protocol'][ind])    
    if params['protocol'][ind] == '2b':
        GHMs = ['H08', 'LPJmL', 'MATSIRO']
        SOCs = ['rcp26soc', 'rcp60soc']

    # Calculate multi-model ensemble statistics for heatwave indicators    
    if ind == 'heatwave':
        
        # Loop through variables
        for var in params['indicators'][ind]: 
        
            # Loop through percentiles and day over threshold combinations
            for p, d in it.product(percentiles, dt):
                
                # Create name for output variables and description based on
                # percentiles and days over threshold
                output_var = f'{params["indicators"][ind][var]["short_name"]}_{p}_{d}'
                hw_desc = f'{params["indicators"][ind][var]["long_name"]} - {p}th percentile / {d} days over threshold'
        
                # Loop through file types
                for ftype in types:                
                             
                    # Calculate multi model statistics and write to file
                    multi_model_stats = cfp.create_esm_multi_model_stats(
                                            input_dir, GCMs, RCPs, thresholds, 
                                            types[ftype]['short'], output_var, 
                                            GHMs, SOCs)
                    attributes = cfa.cust_attrs_mm_esm(
                                            list(multi_model_stats.threshold.values), 
                                            hw_desc, params['protocol'][ind], 
                                            types[ftype]['desc'])
                    output_file = os.path.join(output_dir, f'ISIMIP{params["protocol"][ind]}_{ftype}_{output_var}.nc4')
                    cf.write_output(multi_model_stats, output_file, var, attributes)
                
                # Calculate multi model statistics for reference period and write to file
                multi_model_stats_ref = cfp.create_esm_multi_model_stats(
                                            input_dir, GCMs, RCPs, ['historical'], 
                                            'abs', output_var, GHMs, SOCs)
                attributes = cfa.cust_attrs_mm_esm(
                                            'historical', hw_desc, params['protocol'][ind], 
                                            'absolute values')
                output_file = os.path.join(output_dir, f'ISIMIP{params["protocol"][ind]}_MM_historical_{output_var}.nc4')      
                cf.write_output(multi_model_stats_ref, output_file, var, attributes)
    
    # Calculate multi-model ensemble statistics for other indicators
    else:
    
        # Loop through variables
        for var in params['indicators'][ind]:

            # Loop through file types
            for ftype in types:
            
                # Calculate multi model statistics and write to file
                multi_model_stats = cfp.create_esm_multi_model_stats(
                                        input_dir, GCMs, RCPs, thresholds, 
                                        types[ftype]['short'], 
                                        params["indicators"][ind][var]["short_name"], 
                                        GHMs, SOCs)
                attributes = cfa.cust_attrs_mm_esm(
                                        list(multi_model_stats.threshold.values), 
                                        params['indicators'][ind][var]['long_name'], 
                                        params['protocol'][ind], types[ftype]['desc'])
                output_file = os.path.join(output_dir, f'ISIMIP{params["protocol"][ind]}_{ftype}_{params["indicators"][ind][var]["short_name"]}.nc4')      
                cf.write_output(multi_model_stats, output_file, params["indicators"][ind][var]["short_name"], 
                                        attributes)
            
            # Calculate multi model statistics for reference period and write to file
            if params['protocol'][ind] == '2b':     
                multi_model_stats_ref = cfp.create_esm_multi_model_stats(
                                            input_dir, GCMs, RCPs, ['historical'], 
                                            'abs', params["indicators"][ind][var]["short_name"], 
                                            GHMs, 'histsoc')
                attributes = cfa.cust_attrs_mm_esm(
                                            'historical', 
                                            params['indicators'][ind][var]['long_name'],
                                            params['protocol'][ind], 'absolute values')
                output_file = os.path.join(output_dir, f'ISIMIP{params["protocol"][ind]}_MM_historical_{params["indicators"][ind][var]["short_name"]}.nc4')      
                cf.write_output(multi_model_stats_ref, output_file, 
                                            params["indicators"][ind][var]["short_name"], 
                                            attributes)
            else:
                multi_model_stats_ref = cfp.create_esm_multi_model_stats(
                                            input_dir, GCMs, RCPs, ['historical'], 
                                            'abs', params["indicators"][ind][var]["short_name"], 
                                            GHMs, SOCs)
                attributes = cfa.cust_attrs_mm_esm(
                                            'historical', 
                                            params['indicators'][ind][var]['long_name'], 
                                            params['protocol'][ind], 'absolute values')
                output_file = os.path.join(output_dir, f'ISIMIP{params["protocol"][ind]}_MM_historical_{params["indicators"][ind][var]["short_name"]}.nc4')      
                cf.write_output(multi_model_stats_ref, output_file, params["indicators"][ind][var]["short_name"], 
                                            attributes)
