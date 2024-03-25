 # -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 10:17:08 2022

@author: werning
"""

import numpy as np
import pandas as pd
import sys
sys.path.append('path_to_folder_with_repo')
import glob
import os
import cse_functions_pp as cfp
import cse_functions_tables as cft
import re
import xarray as xr

#%% Settings
#------------------------------------------------------------------------------

# Specify paths to directories and files
input_dir = 'split_files'
output_dir = 'table_output'
raster_path = 'countrymasks_fractional.nc' # https://doi.org/10.48364/ISIMIP.635131.1
regions_path = 'path_to_final_repo\\required_files\\region_classification.xlsx'
yaml_path = 'path_to_final_repo\\required_files\\cse_params.yml'
pop_dir = 'population'
land_area_path = 'path_to_final_repo\\required_files\\gridarea05.nc'
# country_scalars_path = 'path_to_final_repo\\required_files\\Scale_population.nc4' # for countries
country_scalars_path = 'path_to_final_repo\\required_files\\Scale_population_R10_EU27.nc4' # for R10 and EU

# Select output mode - R10 or COUNTRIES
mode = 'R10'

# Specify parameters
indicators = ['cdd']
years = np.arange(2010,2101,10)
types = {'abs': 'Absolute', 'diff': 'Difference', 'score': 'Risk score'}
high = 5 
medium = 3
low = 1
stat = 'mean'

# %% Preparation
# -----------------------------------------------------------------------------

if stat != 'mean':
    addon = f'_{stat}'
else:
    addon = ''

# Load files
params = cfp.load_parameters(yaml_path)
land_area = cft.load_netcdf(land_area_path)
country_scalars = xr.open_dataset(country_scalars_path)
raster = cft.load_netcdf(raster_path)
raster = raster.rename({i: i[2:] for i in list(raster)})

# Create raster for countries or R10 regions
if mode == 'R10':
    countries = cft.create_raster(regions_path, raster)
    raster = countries
    
    # EU27 REGION
    EU27 = ['AUT', 'BEL', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', 'FIN',
            'FRA', 'DEU', 'GRC', 'HUN', 'IRL', 'ITA', 'LVA', 'LTU', 'LUX',
            'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'ESP', 'SWE']

    # Create country raster for EU27
    region_mask = xr.Dataset(coords={'lon': raster.lon, 'lat': raster.lat})
    region_mask['EU27'] = 0  

    for c in countries:    
        region_mask['EU27'] = region_mask['EU27'] + raster[c]   

    raster['EU'] = region_mask['EU27']    

else:
    raster = raster[[iso for iso in list(country_scalars) if len(iso) == 3]]        

# Set all grid cells not belonging to countries to NaN
raster_nan = xr.where(raster == 0, np.nan, raster)
regions = list(raster)

# Unit conversion
land_area = land_area / 1000000  

#%% Calculation
# -----------------------------------------------------------------------------
 
for ind in indicators:
        
    df = pd.DataFrame(columns=[['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()])
    
    print(f'ind: {ind}')
    
    for var in params['indicators'][ind]:
        
        print(var)
               
        # Create file list with all files for indicator
        short = params['indicators'][ind][var]['short_name']        
        files = glob.glob(os.path.join(input_dir, addon[1:], ind, f'*{short}*.nc4'))
            
        for f in files:
            
            # Load data
            data = cft.load_netcdf(f)   
            data = data.where(data.pipe(np.isfinite))   
            
            # Get ssp, gwl, and type from file name
            ssp, gwl, ftype = f.split('_')[-3:]
            ftype = ftype.split('.')[0]
            print(f'{ssp}_{gwl}_{ftype}')                
            pop_data = cft.load_population_data(pop_dir, ssp)
            
            # Interpolate for land cover data and get output variable from file name
            if ind in ['lc']:
                raster_nan = raster_nan.interp(lon=data.lon.values, lat=data.lat.values, method='nearest')
                land_area = land_area.interp(lon=data.lon.values, lat=data.lat.values, method='nearest')
                pop_data = pop_data.interp(lon=data.lon.values, lat=data.lat.values, method='nearest')            
            if 'pm2p5' in ind:
                output_var = re.search(r'GAINS_(.*?)_ssp', f).group(1)
            else:
                output_var = re.search(r'ISIMIP.._(.*?)_ssp', f).group(1)
               
            # Set unit for indicator
            unit = cft.set_unit(ftype, params['indicators'][ind][var])
            
            # Average data for country
            raster_data = raster_nan * data
            dim=['lat', 'lon']
            mean_value = raster_data.mean(dim=dim)
            df_mean = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', mean_value[i].name, f'{output_var}|Hazard|{types[ftype]}', unit] + 
                                    [mean_value[i].item()] * len(years) for i in list(mean_value)])
            df_mean.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
            df = pd.concat([df, df_mean])
            
            # Average data for country and weight by land area
            raster_land = raster_nan * land_area
            data_land = (data * raster_land).sum(dim=dim)
            land_per_cntry = raster_land.sum(dim=dim)
            mean_land = data_land / land_per_cntry
            df_mean_land = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', mean_land[i].name, f'{output_var}|Hazard|{types[ftype]}|Land area weighted', unit] + [mean_land[i].item()] * len(years) for i in list(mean_land)])
            df_mean_land.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
            df = pd.concat([df, df_mean_land])            
            
            if ftype == 'score':
                
                # Calculate exposed land area for score >= medium
                exp_land = raster_land.where(data >= medium).sum(dim=dim)
                df_exp_land = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land[i].name, f'{output_var}|Exposure|Land area', 'km2'] + [exp_land[i].item()] * len(years) for i in list(exp_land)])
                df_exp_land.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land])
                
                # Calculate exposed land area as percentage score >= medium
                exp_land_perc = exp_land / land_per_cntry * 100 
                df_exp_land_perc = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land_perc[i].name, f'{output_var}|Exposure|Land area|%', '%'] + [exp_land_perc[i].item()] * len(years) for i in list(exp_land_perc)])
                df_exp_land_perc.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land_perc])
            
                # Calculate exposed land area for score >= low
                exp_land_low = raster_land.where(data >= low).sum(dim=dim)
                df_exp_land_low = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land_low[i].name, f'{output_var}|Exposure|Land area|Low', 'km2'] + [exp_land_low[i].item()] * len(years) for i in list(exp_land_low)])
                df_exp_land_low.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land_low])
        
                # Calculate exposed land area as percentage for score >= low
                exp_land_low_perc = exp_land_low / land_per_cntry * 100 
                df_exp_land_low_perc = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land_low_perc[i].name, f'{output_var}|Exposure|Land area|Low|%', '%'] + [exp_land_low_perc[i].item()] * len(years) for i in list(exp_land_low_perc)])
                df_exp_land_low_perc.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land_low_perc])            
                
                # Calculate exposed land area for score >= high
                exp_land_high = raster_land.where(data >= high).sum(dim=dim)
                df_exp_land_high = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land_high[i].name, f'{output_var}|Exposure|Land area|High', 'km2'] + [exp_land_high[i].item()] * len(years) for i in list(exp_land_high)])
                df_exp_land_high.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land_high])
                
                # Calculate exposed land area as percentage for score >= high
                exp_land_high_perc = exp_land_high / land_per_cntry * 100 
                df_exp_land_high_perc = pd.DataFrame([['Climate Solutions', f'{ssp}_{gwl}', exp_land_high_perc[i].name, f'{output_var}|Exposure|Land area|High|%', '%'] + [exp_land_high_perc[i].item()] * len(years) for i in list(exp_land_high_perc)])
                df_exp_land_high_perc.columns = [['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
                df = pd.concat([df, df_exp_land_high_perc])           
    
            df_mean_pop, df_exp_pop, df_exp_pop_perc, df_exp_pop_low, \
            df_exp_pop_low_perc, df_exp_pop_high, df_exp_pop_high_perc = pd.Series(dtype='float64'), \
            pd.Series(dtype='float64'), pd.Series(dtype='float64'), pd.Series(dtype='float64'), \
            pd.Series(dtype='float64'), pd.Series(dtype='float64'), pd.Series(dtype='float64')
            
            for year in years:
                yr = f'{year}'
                
                if year < 2010:
                    popyear = 2010                 
                else:
                    popyear = year
                    
                # Average data for country and weight by population
                raster_pop = raster_nan * pop_data.sel(time=popyear) * country_scalars.sel({'ssp': ssp.upper(), 'year': str(popyear)})
                data_pop = (data * raster_pop).sum(dim=dim)
                pop_per_cntry = raster_pop.sum(dim=dim)
                mean_pop = data_pop / pop_per_cntry
                df_mean_pop = pd.concat([df_mean_pop, mean_pop.to_pandas()], axis=1)
                
                if ftype == 'score':                           
                    
                    # Calculate exposed population for score >= medium
                    exp_pop = raster_pop.where(data >= medium).sum(dim=dim) 
                    df_exp_pop = pd.concat([df_exp_pop, exp_pop.to_pandas()], axis=1)

                    # Calculate exposed population as percentage for score >= medium
                    exp_pop_perc = exp_pop / pop_per_cntry * 100
                    df_exp_pop_perc = pd.concat([df_exp_pop_perc, exp_pop_perc.to_pandas()], axis=1)

                    # Calculate exposed population for score >= low
                    exp_pop_low = raster_pop.where(data >= low).sum(dim=dim)
                    df_exp_pop_low = pd.concat([df_exp_pop_low, exp_pop_low.to_pandas()], axis=1)

                    # Calculate exposed population as percentage for score >= low
                    exp_pop_low_perc = exp_pop_low / pop_per_cntry * 100
                    df_exp_pop_low_perc = pd.concat([df_exp_pop_low_perc, exp_pop_low_perc.to_pandas()], axis=1)

                    # Calculate exposed population for score >= high
                    exp_pop_high = raster_pop.where(data >= high).sum(dim=dim)
                    df_exp_pop_high = pd.concat([df_exp_pop_high, exp_pop_high.to_pandas()], axis=1)

                    # Calculate exposed population as percentage for score >= high
                    exp_pop_high_perc = exp_pop_high / pop_per_cntry * 100
                    df_exp_pop_high_perc = pd.concat([df_exp_pop_high_perc, exp_pop_high_perc.to_pandas()], axis=1)
                          
            df_mean_pop = cft.set_df_columns(df_mean_pop, f'{output_var}|Hazard|{types[ftype]}|Population weighted', unit, regions, years)
            
            if ftype == 'score':
                
                df_exp_pop = cft.set_df_columns(df_exp_pop, f'{output_var}|Exposure|Population', 'people', regions, years)
                df_exp_pop_low = cft.set_df_columns(df_exp_pop_low, f'{output_var}|Exposure|Population|Low', 'people', regions, years)
                df_exp_pop_high = cft.set_df_columns(df_exp_pop_high, f'{output_var}|Exposure|Population|High', 'people', regions, years)
                df_exp_pop_perc = cft.set_df_columns(df_exp_pop_perc, f'{output_var}|Exposure|Population|%', '%', regions, years)
                df_exp_pop_low_perc = cft.set_df_columns(df_exp_pop_low_perc, f'{output_var}|Exposure|Population|Low|%', '%', regions, years)
                df_exp_pop_high_perc = cft.set_df_columns(df_exp_pop_high_perc, f'{output_var}|Exposure|Population|High|%', '%', regions, years)        
                df_pop = pd.concat([df_mean_pop, df_exp_pop, df_exp_pop_low, df_exp_pop_high, df_exp_pop_perc, df_exp_pop_low_perc, df_exp_pop_high_perc])
                
            else:
                df_pop = df_mean_pop
                
            df_pop['Model'] = 'Climate Solutions'
            df_pop['Scenario'] = f'{ssp}_{gwl}'
            df = pd.concat([df, df_pop])
    
    # Save table output
    df.to_csv(os.path.join(output_dir, f'table_output_{ind}_{mode}_pop_scaled{addon}.csv'), mode='w', index=False)    
