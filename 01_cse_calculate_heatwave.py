# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:30:23 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_indicators as cfi
import cse_functions_attributes as cfa
import xarray as xr

#%% Settings
#------------------------------------------------------------------------------

# Set protocol
protocol = '3b' 
input_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)

# Manually overwrite GCMs/RCP/timestep here if required
GCMs = ['MRI-ESM2-0']
RCPs = ['ssp126']
# timestep = ''

# Set output directory
output_dir = ''

# Choose required tas and set variable
input_var = 'tas'
output_var = 'heatwave'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_repo\\required_files\\ISIMIP3b_GCM_GMT_1601_2100.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
reference_period = [1974, 2004]
year_range = 30
quantiles = [0.95, 0.97, 0.99]
dt = [3, 5, 7, 10] # no of consecutive days

# Load land mask
landmask_path = '..\\required_files\\landareamaskmap0.nc'
with xr.open_dataset(landmask_path, engine="netcdf4") as land_mask:
    land_mask.load()

#%% Calculate heatwaves
#------------------------------------------------------------------------------

 # Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)

for GCM, RCP in it.product(GCMs, RCPs):
    
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
        
        print(f'{GCM}_{RCP}')      
            
        # Calculate wet bulb temperature, save as file and reload (to speed up computation)
        print('twb')
        twb_ref = cfi.calculate_wetbulb(input_dir, GCM, RCP, input_var, timestep, 
                                        reference_period, land_mask)
        attributes = cfa.cust_attrs_twb(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, 'twb')
        twb_ref.name = 'twb'
        twb_ref.attrs = attributes            
        twb_ref.to_netcdf(output_file)
        twb_ref = xr.open_dataset(output_file)
        
        # Calculate quantiles, save as file and reload (to speed up computation)
        print('twb_qnts')
        twb_ref_qnts = cf.calculate_quantiles(twb_ref, quantiles, 50)
        attributes = cfa.cust_attrs_twb_qnts(GCM, RCP, 'historical', reference_period, 
                                             protocol, quantiles)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, 'twb_qnts')
        twb_ref_qnts.attrs = attributes
        twb_ref_qnts.to_netcdf(output_file)
        twb_ref_qnts = xr.open_dataarray(output_file)
        
        # Calculate indicator for all 31 years and average over all years
        hw_ref = cfi.calculate_annual_heatwave(twb_ref.twb, twb_ref_qnts, quantiles, 
                                               dt, reference_period)
        hw_ref_mean = hw_ref.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_heatwave(GCM, RCP, 'historical', reference_period, 
                                             protocol, dt, quantiles)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}')
        cf.write_output(hw_ref_mean, output_file, output_var, attributes)     
        
        # Save averaged data
        attributes = cfa.cust_attrs_heatwave(GCM, RCP, 'historical', reference_period, 
                                             protocol, dt, quantiles)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}_all')
        cf.write_output(hw_ref, output_file, output_var, attributes)        
        
    else:
        continue
       
        
    for gwl in GWLs:
        
        combination = f'{GCM}_{RCP}_{gwl}'
        
        # Check if the current threshold exists for the GCM/RCP combination
        if not years[combination]:
            continue
           
        else:            
            
            # Calculate wet bulb temperature, save as file and reload (to speed up computation)
            print('twb')
            twb = cfi.calculate_wetbulb(input_dir, GCM, RCP, input_var, timestep, 
                                        years[combination], land_mask)
            attributes = cfa.cust_attrs_twb(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, years[combination], 'twb')
            twb.name = 'twb'
            twb.attrs = attributes
            twb.to_netcdf(output_file)
            twb = xr.open_dataset(output_file)            
            
            # Calculate indicator for all 31 years and average over all years
            hw = cfi.calculate_annual_heatwave(twb.twb, twb_ref_qnts, quantiles, 
                                                dt, years[combination])
            hw_mean = hw.mean(dim='year')
           
            # Save annual data
            attributes = cfa.cust_attrs_heatwave(GCM, RCP, gwl, years[combination], 
                                                  protocol, dt, quantiles)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, years[combination], 
                                                f'{output_var}')
            cf.write_output(hw_mean, output_file, output_var, attributes)     
    
            # Save averaged data
            attributes = cfa.cust_attrs_heatwave(GCM, RCP, gwl, years[combination], 
                                                  protocol, dt, quantiles)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, years[combination],
                                                f'{output_var}_all')
            cf.write_output(hw, output_file, output_var, attributes)
