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
import os

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

def create_raster(raster_path, regions=None, country_list=None, add_EU=False,
                  add_world=False):
    
    ''' Create raster files for countries or regions
    
        Arguments:
            - raster_path: path to netcdf file with raster
            - regions: specify regions
            - country_list: list to filter countries
            - add_EU: add EU to raster
            - add_world: add world to raster
    '''
    
    EU = ['AUT', 'BEL', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', 'FIN',
            'FRA', 'DEU', 'GRC', 'HUN', 'IRL', 'ITA', 'LVA', 'LTU', 'LUX',
            'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'ESP', 'SWE']
    
    raster = load_netcdf(raster_path)
    
    if raster.attrs:
        if raster.attrs['repository'] and raster.attrs['repository'] == 'https://github.com/ISI-MIP/isipedia-countries':        
            raster = raster.rename({i: i[2:] for i in list(raster)})  
        raster = raster / raster['world']
    
    if regions is not None:
           
        region_list = regions['region'].unique()
        region_list = [x for x in region_list if type(x) == str]
        
        region_mask = xr.Dataset(coords={'lon': raster.lon, 'lat': raster.lat}) 
        
        for r in region_list:
            
            print(r)
            region_mask[r] = 0
            countries = regions['ISO'][regions['region']==r]
            
            for c in countries:
                if c in list(raster):
                    region_mask[r] = region_mask[r] + raster[c]          
       
        if add_EU == True:
            region_mask['EU'] = 0
            countries = regions['ISO'][regions['region'].isin(EU)]
         
            for c in EU:    
                region_mask['EU'] = region_mask['EU'] + raster[c]  
    
        if add_world == True:
            region_mask['World'] = raster.world
    
        raster = region_mask
        
    raster = xr.where(raster == 0, np.nan, raster)
         
    if country_list:        
        raster = raster[[iso for iso in country_list]]
            
    return raster 

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
    
    elif ftype == 'diff':
        return '%'
    
    else:
        return params['unit']

# -----------------------------------------------------------------------------

def load_weighted_raster(raster_dir, mode):
    
    ''' Loads pre-calculated rasters weighted by land area or population
    
        Arguments:
            - raster_dir: string with path to rasters
            - mode: string with mode (COUNTRIES/R10/R5/IPCC)
    '''
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_land_raster.nc4'), engine="netcdf4") as raster_land:
        raster_land.load()
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_land_per_country.nc4'), engine="netcdf4") as land_per_cntry:
        land_per_cntry.load()
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_weighted_land.nc4'), engine="netcdf4") as weighted_land:
        weighted_land.load()
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_population_raster.nc4'), engine="netcdf4") as raster_pop:
        raster_pop.load()
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_population_per_country.nc4'), engine="netcdf4") as pop_per_cntry:
        pop_per_cntry.load()
        
    with xr.open_dataset(os.path.join(raster_dir, f'{mode}_weighted_population.nc4'), engine="netcdf4") as weighted_pop:
        weighted_pop.load()
        
    return raster_land, land_per_cntry, weighted_land, raster_pop, pop_per_cntry, weighted_pop
    