# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 14:29:15 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import cse_functions as cf
import cse_functions_attributes as cfa
import os
import xarray as xr

#%% Settings
#------------------------------------------------------------------------------

# Set input and output directories
input_dir = 'output_dir_of_indicator_scripts' 
output_dir = 'multi-model'

# Set netCDF format
netcdf4_format = 'NETCDF4_CLASSIC' 

# Set protocol
protocol = '3b' 
isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)

# Manually overwrite GCMs/RCP/timestep here if required
# GCMs = ['GFDL-ESM2M']
# RCPs = ['rcp26', 'rcp60', 'rcp85']
# timestep = ''

# Set variable
indicators = ['twb']
# indicators = ['cdd', 'heatwave', 'dri', 'dri_qtot', etc.]
 
# Set GHMs and socs if doing multi-model statistics for hydrology indicators
# GHMs = ['H08', 'LPJmL', 'MATSIRO']
GHMs = ''
# soc = ['rcp26soc', 'rcp60soc']
soc = ''
# histsoc = 'histsoc'
histsoc = ''

# Specify thresholds and year range
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]

#%% Calculate multi-model mean

data = xr.Dataset()

for ind in indicators:
    
    print(f'{ind}')
    
    # Calculate multi-model ensemble statistics
    multi_model_stats = cf.create_multi_model_stats(input_dir, GCMs, RCPs, GWLs, 
                                                   ind, GHMs, soc)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Save data
    attributes = cfa.cust_attrs_mm(list(multi_model_stats.threshold.values), protocol, ind)
    output_file = os.path.join(output_dir, f'ISIMIP{protocol}_MM_{ind}.nc4')        
    cf.write_output(multi_model_stats, output_file, ind, attributes)
    
    # Calculate multi-model ensemble statistics for reference period
    multi_model_stats_ref = cf.create_multi_model_stats(input_dir, GCMs, RCPs, 
                                                        ['historical'], ind, 
                                                        GHMs, histsoc)
    
    # Save data for reference period 
    attributes = cfa.cust_attrs_mm('historical', protocol, ind)
    output_file = os.path.join(output_dir, f'ISIMIP{protocol}_MM_historical_{ind}.nc4')        
    cf.write_output(multi_model_stats_ref, output_file, ind, attributes)
    