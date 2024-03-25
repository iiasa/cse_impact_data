# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 13:44:53 2023

@author: werning
"""

import xarray as xr
import pandas as pd
import glob
import numpy as np
from pandas.tseries.offsets import YearEnd

# -----------------------------------------------------------------------------

def load_netcdf(file_path):
    
    ''' Load netcdf files
    
        Arguments:
            - file_path: path to netcdf file'''
            
    try:
        with xr.open_dataarray(file_path, engine="netcdf4") as ds:
            return ds.load()
    except:
        try:
            with xr.open_dataset(file_path, engine="netcdf4") as ds:
                return ds.load()
        except:
            print('Not a netcdf file')
            
#------------------------------------------------------------------------------

def create_raster(file_path, raster):
    
    ''' Create raster for R10 regions
    
        Arguments:
            - file_path: path to file with region definitions
            - raster: xarray dataset with rasters for all countries'''            
    
    # Read region definitions and get names of R10 regions
    regions = pd.read_excel(file_path,)
    r10_regions = regions.region_r10_db.unique()
    region_mask = xr.Dataset(coords={'lon': raster.lon, 'lat': raster.lat})        

    for r in r10_regions:
        
        region_mask[r] = 0            
        countries = regions.ISO[regions.region_r10_db==r]
        # countries = [f'm_{c}' for c in countries]
        
        # UDDATE
        # region_mask = xr.Dataset({name: raster[name] for name in raster.keys() if name in countries})
        
        # Add up all countries belonging to R10 regions
        for c in countries:
            
            if c in list(raster):
                region_mask[r] = region_mask[r] + raster[c]
                
    return region_mask
                
# -----------------------------------------------------------------------------

def set_df_columns(df, variable, unit, region, years):
    
    ''' Set column names for dataframe
    
        Arguments:
            - df: Pandas dataframe with data
            - variable: string with output variable
            - unit: string with unit for output variable
            - region: list with region names
            - years: list with years ''' 
    
    df.columns = [['Variable'] + years.tolist()]
    df['Variable'] = variable
    df['Unit'] = unit
    df['Region'] = region
    
    return df               

# -----------------------------------------------------------------------------

def load_population_data(data_dir, ssp, option='total'):
    
    
    ''' Load population data for ssp (either total, urban, or rural)
    
        Arguments:
            - data_dir: string with path to directory with population data
            - ssp: string with ssp
            - option: string with either 'total', 'rural', or 'urban' ''' 
    
    if option != 'urban' and option != 'rural' and option != 'total':
        print('No valid option detected. Please choose either urban, rural or total.')
    
    else:
        
        pop_files = f'{data_dir}\\*{ssp.lower()}*{option}*.nc*'
        pop_file_list = glob.glob(pop_files)

        with xr.open_mfdataset(pop_file_list, parallel=True, use_cftime=True) as pop_data:
            pop_data.load()
            pop_data = pop_data.to_array(dim='time')
            pop_data['time'] = pd.to_datetime(np.arange(2010, 2101, 10), format='%Y')+YearEnd(1)
            pop_data['time'] = np.arange(2010, 2101, 10)
            return pop_data
        
# -----------------------------------------------------------------------------

def set_unit(ftype, params):
    
    ''' Set unit (either indicator unit or risk score)
    
        Arguments:
            - ftype: string with file type
            - params: content of yaml file for indicator ''' 
    
    if ftype == 'score':
        return 'risk score'
    
    else:
        return params['unit']

