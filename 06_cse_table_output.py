# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:20:52 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import numpy as np
import os
import glob
import pandas as pd
import re
import cse_functions_pp as cfp
import cse_functions as cf
import cse_functions_tables as cft

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths to files
input_dir = 'output_dir_of_indicator_scripts'
output_dir = ''
pop_dir = 'population_data'
raster_dir = 'output_dir_of_weighted_rasters'
raster_path = 'ISIMIP_fractional_mask'
regions_path = 'path_to_folder_with_repo\\required_files\\region_classification.xlsx'
ipcc_raster_path = 'path_to_folder_with_repo\\required_files\\IPCC-WGI-reference-regions_fractional_mask.nc4'
R5_path = 'UN_R5_yaml' 
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'
land_area_path = 'path_to_folder_with_repo\\required_files\\gridarea05.nc'
land_cover_path = 'path_to_folder_with_repo\\required_files\\landcovermasks.nc4'

# Set parameters
indicators = ['cdd', 'dri', 'dri_qtot', 'heatwave', 'iavar', 'iavar_qtot', 
              'precip', 'sdd', 'sdd_18p3', 'sdd_20p0', 'sdd_24p0', 'seas', 
              'seas_qtot', 'tr20', 'wsi']
years = np.arange(2020,2101,10)
types = {'abs': 'Absolute', 'diff': 'Difference', 'score': 'Hazard score'}
high = 5 
medium = 3
low = 1
SSPs = ['ssp1', 'ssp2', 'ssp3', 'ssp4', 'ssp5']
GHMs = ''
SOCs = ''
quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]

# Select output mode - COUNTRIES/R10/R5/IPCC
mode = 'IPCC'

# Load files
params = cfp.load_parameters(yaml_path)
land_area = cft.load_netcdf(land_area_path)
land_area = land_area / 1000000       
landcover_mask = cfp.load_landcover_mask(land_cover_path)
pop_data = cft.load_population_data(pop_dir)
    
#%% Create table output - prepare rasters and load weighted rasters
#------------------------------------------------------------------------------

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
regions = list(raster)

raster_land, land_per_cntry, weighted_land, raster_pop, pop_per_cntry, weighted_pop = cft.load_weighted_rasters(raster_dir, mode)

#%% Create table output - prepare rasters
#------------------------------------------------------------------------------

# Loop through indicators
for ind in indicators:
    
    print(f'ind: {ind}')
    
    # Set protocol, GCMs, RCPs, GHMs, SOCs, etc.
    isimip_dir, GCMs, RCPs, timestep = cf.set_protocol(params['protocol'][ind])    
    if params['protocol'][ind] == '2b':
        GHMs = ['H08', 'LPJmL', 'MATSIRO']
        SOCs = ['rcp26soc', 'rcp60soc']       
     
    # Loop through variables                 
    for var in params['indicators'][ind]:
        
        print(var)
        
        # Create dataframe and short variable name
        df = pd.DataFrame(columns=['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist())                     
        short = params['indicators'][ind][var]['short_name']       
        
        # Create file list with all required files        
        if SOCs and GHMs:        
            files = sum([glob.glob(os.path.join(input_dir, GHM, RCP, GCM, 
                                                f'{GHM}_{GCM.lower()}_{RCP}_*_{SOC}_{short}_global_[1-9][0-9][0-9][0-9]_[1-9][0-9][0-9][0-9]_*')) \
                              for GHM, GCM, RCP, SOC in it.product(GHMs, GCMs, RCPs, SOCs)], [])  
        else:            
            files = sum([glob.glob(os.path.join(input_dir, RCP, GCM, 
                                                f'{GCM.lower()}_{RCP}_*{short}_global_[1-9][0-9][0-9][0-9]_[1-9][0-9][0-9][0-9]_*')) 
                         for GCM, RCP in it.product(GCMs, RCPs)], [])            

        # Loop through all files
        for f in files:           
           
            # Load data
            data = cft.load_netcdf(f)
            
            # Mask out Greenland ice sheet and deserts for hydrology indicators
            if ind in ['dri', 'dri_qtot', 'iavar', 'iavar_qtot', 'seas', 
                       'seas_qtot', 'wsi', 'wsi_qtot']:                
                data = data.where(landcover_mask.grid_index < 11)
            
            # Get parameters from file name
            esm = f.split('\\')[-1].split('_')[0] + '_' + f.split('\\')[-1].split('_')[1] if params['protocol'][ind] == '3b' else f.split('\\')[-1].split('_')[0] + '_' + f.split('\\')[-1].split('_')[1] + '_' + f.split('\\')[-1].split('_')[2]
            gwl  = f.split('\\')[-1].split('_')[2] if params['protocol'][ind] == '3b' else f.split('\\')[-1].split('_')[3]
            ftype = f.split('\\')[-1].split('_')[-1][:-4]                         
            if SOCs and GHMs:
                output_var = re.search('soc_(.*)_global', f).group(1)
            else:
                output_var = re.search(f'{gwl}_(.*)_global', f).group(1)
            unit = cft.set_unit(ftype, params['indicators'][ind][var])           
            
            # Calculate average value per country/region
            raster_data = raster * data            
            dim=['lat', 'lon']
            mean_value = (raster_data.mean(dim=dim))
            df_mean = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', mean_value[i].name, f'{output_var}|Hazard|{types[ftype]}', unit] + [mean_value[i].item()] * len(years) for i, ssp in it.product( list(mean_value), SSPs)])
            df_mean.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()                     
            df = pd.concat([df, df_mean])
            
            # Calculate average value per country/region weighted by land area
            mean_land = (data * weighted_land).sum(dim=dim)            
            df_mean_land = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', mean_land[i].name, f'{output_var}|Hazard|{types[ftype]}|Land area weighted', unit] + [mean_land[i].item()] * len(years) for i, ssp in it.product( list(mean_land), SSPs)])
            df_mean_land.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()            
            df = pd.concat([df, df_mean_land])            
            
            if ftype == 'score':
               
                # Calculate exposed land area to at least medium hazard
                exp_land = (raster_land.where(data >= medium).sum(dim=dim))
                df_exp_land = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land[i].name, f'{output_var}|Exposure|Land area', 'km2'] + [exp_land[i].item()] * len(years) for i, ssp in it.product(list(exp_land), SSPs)])
                df_exp_land.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land])
               
                # Calculate exposed land area to at least medium hazard (%)
                exp_land_perc = (exp_land / land_per_cntry * 100)
                df_exp_land_perc = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land_perc[i].name, f'{output_var}|Exposure|Land area|%', '%'] + [exp_land_perc[i].item()] * len(years) for i,ssp in it.product(list(exp_land_perc), SSPs)])
                df_exp_land_perc.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land_perc])
           
                # Calculate exposed land area to at least low hazard
                exp_land_low = (raster_land.where(data >= low).sum(dim=dim))
                df_exp_land_low = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land_low[i].name, f'{output_var}|Exposure|Land area|Low', 'km2'] + [exp_land_low[i].item()] * len(years) for i,ssp in it.product(list(exp_land_low), SSPs)])
                df_exp_land_low.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land_low])
       
                # Calculate exposed land area to at least low hazard (%)
                exp_land_low_perc = (exp_land_low / land_per_cntry * 100 )
                df_exp_land_low_perc = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land_low_perc[i].name, f'{output_var}|Exposure|Land area|Low|%', '%'] + [exp_land_low_perc[i].item()] * len(years) for i,ssp in it.product(list(exp_land_low_perc), SSPs)])
                df_exp_land_low_perc.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land_low_perc])            
               
                # Calculate exposed land area to at least high hazard
                exp_land_high = (raster_land.where(data >= high).sum(dim=dim))
                df_exp_land_high = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land_high[i].name, f'{output_var}|Exposure|Land area|High', 'km2'] + [exp_land_high[i].item()] * len(years) for i,ssp in it.product(list(exp_land_high), SSPs)])
                df_exp_land_high.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land_high])
               
                # Calculate exposed land area to at least high hazard (%)
                exp_land_high_perc = (exp_land_high / land_per_cntry * 100 )
                df_exp_land_high_perc = pd.DataFrame([[f'{esm}', f'{ssp}_{gwl}', exp_land_high_perc[i].name, f'{output_var}|Exposure|Land area|High|%', '%'] + [exp_land_high_perc[i].item()] * len(years) for i, ssp in it.product(list(exp_land_high_perc), SSPs)])
                df_exp_land_high_perc.columns = ['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()
                df = pd.concat([df, df_exp_land_high_perc])  
            
            # Calculate average value per country/region weighted by population
            mean_pop = (data * weighted_pop).sum(dim=dim)   
            df_mean_pop = mean_pop.to_array('Region').to_dataframe('value')
            df_mean_pop = df_mean_pop.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
            df_mean_pop['Scenario'] = df_mean_pop['Scenario'] + f'_{gwl}'
            df_mean_pop['Variable'] = f'{output_var}|Hazard|{types[ftype]}|Population weighted'
            df_mean_pop['Unit'] = unit
            
            if ftype == 'score': 
               
                # Calculate exposed population to at least medium hazard
                exp_pop = (raster_pop.where(data >= medium).sum(dim=dim) )
                df_exp_pop = exp_pop.to_array('Region').to_dataframe('value')
                df_exp_pop = df_exp_pop.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop['Scenario'] = df_exp_pop['Scenario'] + f'_{gwl}'
                df_exp_pop['Variable'] = f'{output_var}|Exposure|Population'
                df_exp_pop['Unit'] = 'people'

                # Calculate exposed population to at least medium hazard (%)
                exp_pop_perc = (exp_pop / pop_per_cntry * 100)
                df_exp_pop_perc = exp_pop_perc.to_array('Region').to_dataframe('value')
                df_exp_pop_perc = df_exp_pop_perc.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop_perc['Scenario'] = df_exp_pop_perc['Scenario'] + f'_{gwl}'
                df_exp_pop_perc['Variable'] = f'{output_var}|Exposure|Population|%'
                df_exp_pop_perc['Unit'] = '%'
                
                # Calculate exposed population to at least low hazard
                exp_pop_low = (raster_pop.where(data >= low).sum(dim=dim))
                df_exp_pop_low = exp_pop_low.to_array('Region').to_dataframe('value')
                df_exp_pop_low = df_exp_pop_low.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop_low['Scenario'] = df_exp_pop_low['Scenario'] + f'_{gwl}'
                df_exp_pop_low['Variable'] = f'{output_var}|Exposure|Population|Low'
                df_exp_pop_low['Unit'] = 'people'
               
                # Calculate exposed population to at least low hazard (%)
                exp_pop_low_perc = (exp_pop_low / pop_per_cntry * 100)
                df_exp_pop_low_perc = exp_pop_low_perc.to_array('Region').to_dataframe('value')
                df_exp_pop_low_perc = df_exp_pop_low_perc.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop_low_perc['Scenario'] = df_exp_pop_low_perc['Scenario'] + f'_{gwl}'
                df_exp_pop_low_perc['Variable'] = f'{output_var}|Exposure|Population|Low|%'
                df_exp_pop_low_perc['Unit'] = '%'

                # Calculate exposed population to at least high hazard
                exp_pop_high = (raster_pop.where(data >= high).sum(dim=dim))
                df_exp_pop_high = exp_pop_high.to_array('Region').to_dataframe('value')
                df_exp_pop_high = df_exp_pop_high.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop_high['Scenario'] = df_exp_pop_high['Scenario'] + f'_{gwl}'
                df_exp_pop_high['Variable'] = f'{output_var}|Exposure|Population|High'
                df_exp_pop_high['Unit'] = 'people'

                # Calculate exposed population to at least high hazard (%)
                exp_pop_high_perc = (exp_pop_high / pop_per_cntry * 100)
                df_exp_pop_high_perc = exp_pop_high_perc.to_array('Region').to_dataframe('value')
                df_exp_pop_high_perc = df_exp_pop_high_perc.pivot_table(values='value', index=['Region', 'Scenario'], columns='time').reset_index()
                df_exp_pop_high_perc['Scenario'] = df_exp_pop_high_perc['Scenario'] + f'_{gwl}'
                df_exp_pop_high_perc['Variable'] = f'{output_var}|Exposure|Population|High|%'
                df_exp_pop_high_perc['Unit'] = '%'
               
                df_pop = pd.concat([df_mean_pop, df_exp_pop, df_exp_pop_low, df_exp_pop_high, df_exp_pop_perc, df_exp_pop_low_perc, df_exp_pop_high_perc])
                 
            else:
                df_pop = df_mean_pop
               
            df_pop['Model'] = f'{esm}'
            df = pd.concat([df, df_pop])    
            
        # Calculate quantiles for all ensemble members
        quantile_data = df.groupby(by=['Variable', 'Region', 'Scenario', 'Unit']).quantile(quantiles, numeric_only=True).reset_index()  
        quantile_data['level_4'] = (quantile_data['level_4'] * 100)
        quantile_data['Variable'] = quantile_data['Variable'].astype(str) + '|' +  quantile_data['level_4'].astype(str) + 'th Percentile'
        quantile_data = quantile_data.drop('level_4', axis=1)
        quantile_data['Model'] = 'Climate Solutions'
          
        # Calculate mean of all ensemble members
        mean_data = df.groupby(by=['Variable', 'Region', 'Scenario', 'Unit']).mean(numeric_only=True).reset_index()
        mean_data['Model'] = 'Climate Solutions'
        mean_data['Variable'] = mean_data['Variable'].astype(str) + '|Mean'    
        df = pd.concat([df, quantile_data, mean_data],ignore_index=True)
        
        # Save table output
        df.to_csv(os.path.join(output_dir, f'table_output_{output_var}_{mode}.csv'), mode='w', index=False)    


