# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 12:02:49 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import os
import pandas as pd
import sys
import cse_functions_tables as cft


#%% Settings
#------------------------------------------------------------------------------

# Specify paths to directories and files
output_dir = ''
pop_dir = 'population_data'
land_area_path = 'path_to_folder_with_repo\\required_files\\gridarea05.nc'
raster_path = 'isimip_fractional_country_mask'
regions_path = 'path_to_folder_with_repo\\required_files\\region_classification.xlsx'
R5_path = 'un_r5.yaml' 
ipcc_raster_path = 'path_to_folder_with_repo\\required_files\\IPCC-WGI-reference-regions_fractional_mask.nc4'

# Set parameters
modes = ['R10', 'R5', 'IPCC', 'ISIMIP']  

# Load files
land_area = cft.load_netcdf(land_area_path)
land_area = land_area / 1000000
pop_data = cft.load_population_data(pop_dir)    

#%% Create land area and population weighted rasters
#------------------------------------------------------------------------------

# Loop through all modes
for mode in modes:

    # Create raster files for countries/regions
    if mode == 'R10':  
        regions = pd.read_excel(regions_path,)
        r10_regions = regions[['ISO', 'region_r10_db']].rename(columns={'region_r10_db': 'region'})
        raster = cft.create_raster(raster_path, r10_regions, add_EU=True, add_world=True)  
        print('R10')        
    elif mode == 'R5':
        regions = pd.read_excel(regions_path,)
        r5_regions = regions[['ISO', 'region_unep_r5']].rename(columns={'region_unep_r5': 'region'})
        raster = cft.create_raster(raster_path, r5_regions, add_world=True)  
    elif mode == 'IPCC':
        raster = cft.create_raster(ipcc_raster_path)
        raster = raster[[i for i in list(raster) if 'Ocean' not in i and 'Antarctica' not in i and 'Sea' not in i and 'Bay' not in i]]
        print('IPCC')        
    else:
        raster = cft.create_raster(raster_path, add_world=True)
        print('ISIMIP')
        
    # Weight raster by land area and population and save to file
    raster_land = (raster * land_area)
    raster_land.to_netcdf(os.path.join(output_dir, f'{mode}_land_raster.nc4'))
    
    land_per_cntry = raster_land.sum(dim=['lat', 'lon'])
    land_per_cntry.to_netcdf(os.path.join(output_dir, f'{mode}_land_per_country.nc4')) 
    
    weighted_land = raster_land / (raster_land.sum(dim=['lat', 'lon']))
    weighted_land.to_netcdf(os.path.join(output_dir, f'{mode}_weighted_land.nc4'))
    
    raster_pop = raster * pop_data
    raster_pop.to_netcdf(os.path.join(output_dir, f'{mode}_population_raster.nc4'))
    
    pop_per_cntry = raster_pop.sum(dim=['lat', 'lon'])
    pop_per_cntry.to_netcdf(os.path.join(output_dir, f'{mode}_population_per_country.nc4')) 
    
    weighted_population = raster_pop / (raster_pop.sum(dim=['lat', 'lon']))
    weighted_population.to_netcdf(os.path.join(output_dir, f'{mode}_weighted_population.nc4'))
       