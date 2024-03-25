# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:27:44 2021

@author: werning
"""

import pandas as pd
import numpy as np
import itertools as it
import xarray as xr
import os
import glob
from dask.diagnostics import ProgressBar
import dask
import re

# %% Settings
# -----------------------------------------------------------------------------
netcdf4_format = 'NETCDF4_CLASSIC' 

# %% Functions
# -----------------------------------------------------------------------------

def find_year_range(GCMs, RCPs, GMT_anomaly_file, thresholds, year_range):
    
    ''' Find min and max years for temperature thresholds and GCM/RCP combination 
    
        Arguments:
            - GCMs: list of GCMs used to look up values in the anomaly spreadsheet
            - RCPs: list of RCPs used to look up values in the anomaly spreadsheet
            - GMT_anomaly_file: Excel spreadsheet containing temperature anomalies
                                for all GCM/RCP combination
            - thresholds: list of desired thresholds, e.g. 1.0, 1.5, etc.
            - year_range: interval around threshold year, e.g. 30 years 
                          (15 years on either side)'''    
    
    threshold_years = {}
    GMT_anomalies = pd.read_excel(GMT_anomaly_file)
    
    for GCM, RCP in it.product(GCMs, RCPs):
    
        column_name = f'{GCM}_{RCP}'    
    
        for i in thresholds:    
 
            # Check if the temperature threshold doesn't get crossed and return empty threshold years in this case
            if GMT_anomalies[(GMT_anomalies[column_name] >= i) & 
                             (GMT_anomalies[column_name].shift() < i)].empty:        
                threshold_years[f'{column_name}_{i}'] = []
            
            # If there are years where the temperature threshold gets crossed, select the first one and calculate corresponding years
            else:        
                threshold_index = GMT_anomalies[(GMT_anomalies[column_name] >= i) & 
                                                (GMT_anomalies[column_name].shift() < i)].index[0]
                threshold_years[f'{column_name}_{i}'] = [GMT_anomalies.year[threshold_index] - int(year_range/2), 
                                                         GMT_anomalies.year[threshold_index] + int(year_range/2)]
            
    return threshold_years 

#------------------------------------------------------------------------------

def set_protocol(protocol):
    
    ''' Sets the input directory, list of GCMs, list of RCPs, and timestep
        automatically for the specified protocol
        
        Arguments:
            - protocol: string with the chosen ISMIP protocol (2b or 3b) '''
    
    if protocol == '2b':
        input_dir = "P:\\watxene\\ISIMIP\\ISIMIP2b\\input"
        GCMs = ['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5']
        # RCPs = ['rcp26', 'rcp60', 'rcp85']
        RCPs = ['rcp26', 'rcp60']
        timestep = 'day'
        return input_dir, GCMs, RCPs, timestep
    
    elif protocol == '3b':
        input_dir = "P:\\watxene\\ISIMIP\\ISIMIP3b\\InputData\\climate_updated\\bias-adjusted"
        GCMs = ['GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL']
        RCPs = ['ssp126', 'ssp370', 'ssp585']
        timestep = 'daily'
        return input_dir, GCMs, RCPs, timestep
    
    else:
        print('Protocol not recognised')

#------------------------------------------------------------------------------

def load_data(files, years):
    
    ''' Open datatsets based on the created file list and limited to the required
        years only
        
        Arguments:
            - files: list containing all files for the selected GHM/GCM/RCP 
                     combination
            - years: list consisting of start and end year'''
    
    file_list = select_files_by_year_range(files, years)
    
    with ProgressBar():
        
        data_all = xr.open_mfdataset(file_list, parallel=True, use_cftime=True, 
                                     chunks={'lon':400, 'lat':200, 'time': 400})

        # take care of other date formats
        if type(data_all.indexes['time']) is not pd.DatetimeIndex:
            data_all['time'] = data_all.indexes['time'].to_datetimeindex()
    
    data_all = data_all.sel(time=slice(str(pd.to_datetime(years[0], format='%Y'))[:10], 
                                       str(pd.to_datetime(years[1], format='%Y') 
                                           + pd.offsets.YearEnd(0))[:10]))
    
    # Check that all required years are included
    if ((pd.to_datetime(years[0], format='%Y').date() != pd.to_datetime(data_all.time[0].values).date()) or \
        (pd.to_datetime(years[1], format='%Y') + pd.offsets.YearEnd(0)).date() != pd.to_datetime(data_all.time[-1].values).date()):        
        print('Missing data')
    
    return data_all

#------------------------------------------------------------------------------                                  

def create_file_list(input_dir, GCM, RCP, var, timestep, GHM=''):
    
    ''' Create file list with all matching input files 
        
        Arguments:
            - input_dir: path to the directory where the input data is stored 
            - GCM: string with currently selected GCM
            - RCP: string with currently selected RCP
            - var: variable indentifier in the output file
            - timestep: timestep of the input data
            - GHM: name of currently selected GHM - OPTIONAL '''
    
    file_path_future = f'{os.path.join(input_dir, GHM, RCP, GCM)}\\*.nc*'
    file_list_future = glob.glob(file_path_future)
       
    # Create file path and list for historical years
    file_path_hist = f'{os.path.join(input_dir, GHM, "historical", GCM)}\\*.nc*'
    file_list_hist = glob.glob(file_path_hist)
      
    # Combine both and exclude all datasets that are land only
    all_files = file_list_future + file_list_hist
    file_list = [item for item in all_files if f'{var}_' in item and timestep in item 
                 and 'landonly' not in item and re.findall("[0-9]{4,}", item)[0] < '2260']
    # restriction to years before 2260 needed due to cftime limitations - dates after 2262-4-12 are not supported
        
    return file_list

#------------------------------------------------------------------------------

def create_soc_file_list(input_dir, GCM, RCP, var, timestep, histsoc, futuresoc, GHM=''):
    
    ''' Create file list with all matching input files in case of SOC
        
        Arguments:
            - input_dir: path to the directory where the input data is stored 
            - GCM: string with currently selected GCM
            - RCP: string with currently selected RCP
            - var: variable indentifier in the output file
            - timestep: timestep of the input data
            - histsoc: string with SOC to use for historical SOC
            - futuresoc: string with SOC to use for future SOC
            - GHM: name of currently selected GHM - OPTIONAL '''
    
    file_path_future = f'{os.path.join(input_dir, GHM, RCP, GCM)}\\*.nc*'
    all_files_future = glob.glob(file_path_future)
    file_list_future = [item for item in all_files_future if var in item and timestep in item and futuresoc in item and 'landonly' not in item]
       
    # Create file path and list for historical years
    file_path_hist = f'{os.path.join(input_dir, GHM, "historical", GCM)}\\*.nc*'
    all_files_hist = glob.glob(file_path_hist)
    file_list_hist = [item for item in all_files_hist if var in item and timestep in item and histsoc in item and 'landonly' not in item]
      
    # Combine both and exclude all datasets that are land only
    file_list = file_list_future + file_list_hist
    
    return file_list

#------------------------------------------------------------------------------

def select_files_by_year_range(files, years):
    
    ''' Shorten file list containing all files for the selected GHM/GCM/RCP combination
        to only those that are required for selected year range
        
        Arguments:
            - files: list containing all files for the selected GHM/GCM/RCP
                     combination
            - years: list consisting of start and end year'''
    
    file_list = []
    
    for i in files:
        
        years_in_file = re.findall("[0-9]{4,}", i)    
        if (years_in_file[-1][0:4] >= str(years[0]) and 
            years_in_file[-2][0:4] <= str(years[1])):
            file_list.append(i)
            
    return file_list

#------------------------------------------------------------------------------

def convert_to_Celsius(input):  
    
    ''' Convert temperatures from Kelvin to Celsius '''
    
    return input - 273.15

#------------------------------------------------------------------------------      

def create_output_file(output_dir, GCM, RCP, threshold, years, var, GHM=''):
    
    ''' Create output file name and path 
        
        Arguments:
            - output_dir: path to the directory where the output will be saved
            - GCM: string with currently selected GCM
            - RCP: string with currently selected RCP
            - threshold: string with currently selected temperature threshold
            - years: list consisting of start and end year
            - var: variable indentifier in the output file
            - GHM: name of currently selected GHM - OPTIONAL '''
      
    full_output_dir = os.path.join(output_dir, GHM, RCP, GCM)
    
    # Create output directories if they don't exist yet
    if not os.path.exists(full_output_dir):
        os.makedirs(full_output_dir)    
               
    if not GHM:
        output_formatted = f'{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{var}_global_{years[0]}_{years[1]}.nc4'
        data_file = os.path.join(output_dir, RCP, GCM, output_formatted)
        
    else:
        output_formatted = f'{GHM}_{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{var}_global_{years[0]}_{years[1]}.nc4'
        data_file = os.path.join(output_dir, GHM, RCP, GCM, output_formatted)
        
    return data_file

#------------------------------------------------------------------------------ 
    
def create_soc_output_file(output_dir, GCM, RCP, SOC, threshold, years, var, GHM=''):
    
    ''' Create output file name and path in case of soc
        
        Arguments:
            - output_dir: path to the directory where the output will be saved
            - GCM: string with currently selected GCM
            - RCP: string with currently selected RCP
            - SOC: string with currently selected SOC
            - threshold: string with currently selected temperature threshold
            - years: list consisting of start and end year
            - var: variable indentifier in the output file
            - GHM: name of currently selected GHM - OPTIONAL '''
    
    full_output_dir = os.path.join(output_dir, GHM, RCP, GCM)
    
    # Create output directories if they don't exist yet
    if not os.path.exists(full_output_dir):
        os.makedirs(full_output_dir)    
               
    if not GHM:
        output_formatted = f'{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{SOC}_{var}_global_{years[0]}_{years[1]}.nc4'
        data_file = os.path.join(output_dir, RCP, GCM, output_formatted)
        
    else:
        output_formatted = f'{GHM}_{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{SOC}_{var}_global_{years[0]}_{years[1]}.nc4'
        data_file = os.path.join(output_dir, GHM, RCP, GCM, output_formatted)
        
    return data_file

# -----------------------------------------------------------------------------    

def write_output(data, output_file, var, attributes=''):  
    
    ''' Write ouput files for the currently selected GHM/GCM/RCP combination
    
        Arguments:
            - data: data to be written, either xarray dataset or data array
            - output_file: string with the format of the output file names
            - var: string with the name of the output variable
            - attr: attributes to be written''' 
    
    with ProgressBar():
        
        comp = dict(zlib=True, complevel=9)
        
        if isinstance(data, xr.Dataset):
            
            data.attrs = attributes            
            encoding = {var: comp for var in data.data_vars}
            data.to_netcdf(output_file, format=netcdf4_format, encoding = encoding)
            data.close()
            
        else:  
            
            data.name = var
            data = data.to_dataset()
            data.attrs = attributes
            data.to_netcdf(output_file, format=netcdf4_format, encoding = {var: comp})
            data.close()
            
#------------------------------------------------------------------------------

def calculate_quantiles(data, quantile_value, chunksize):
    
    ''' Calculate quantiles 
    
        Arguments:
            - data: data to be used for calculation, xarray dataset or data array
            - quantile_value: float specifying the required quantile
            - chunksize: chunksize used for rechunking data prior to quantile 
                         calculation '''    
    
    quantiles = data.chunk({'time': -1, 'lat':200, 'lon':400}).quantile(quantile_value, skipna=True, dim="time")
    return quantiles

#------------------------------------------------------------------------------

def create_multi_model_stats(input_dir, GCMs, RCPs, thresholds, var, GHMs='', SOC=''):
    
    ''' Calculate quantiles 
    
        Arguments:
            - input_dir: path to directory containing all the raw data
            - GCMs: list with GCMs
            - RCPs: list with RCPs
            - thresholds: list with GWLs
            - var: string with variable name
            - GHMs: list with GHMs if hydrology data
            - SOC: string with soc if hydrology data '''    
    
    data_all = xr.Dataset()
    thrshld = []
    
    for threshold in thresholds:
        
        if SOC and GHMs: 
        
            # Load all files
            if type(SOC) == str:            
                files = sum([glob.glob(os.path.join(input_dir, GHM, RCP, GCM, f'{GHM}_{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{SOC}_{var}_global_*')) \
                          for GHM, GCM, RCP in it.product(GHMs, GCMs, RCPs)], [])
            
            elif type(SOC) == list:            
                files = sum([glob.glob(os.path.join(input_dir, GHM, RCP, GCM, f'{GHM}_{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{SOCS}_{var}_global_*')) \
                          for GHM, GCM, RCP, SOCS in it.product(GHMs, GCMs, RCPs, SOC)], [])                
                    
            else:
                print('Wrong input format for SOC')
                
        else:            
            files = sum([glob.glob(os.path.join(input_dir, RCP, GCM, f'{GCM.lower()}_{RCP}_{str(threshold).replace(".", "p")}_{var}_global_*')) \
                      for GCM, RCP in it.product(GCMs, RCPs)], [])
        
        if len(files) <= 2:    
            continue
        
        else:
            
            thrshld.append(threshold)            
            
            # Load data
            data = xr.open_mfdataset(files, combine='nested', concat_dim='GCM_RCP', chunks={'lon':400, 'lat':200, 'time': 400})
                
            if var == 'twb':
                data = data.mean(dim='time')

            # Exclude outlier values for hydrology indicators
            if var in ['seas', 'seas_qtot', 'iavar', 'iavar_qtot']:
                data = xr.where(data > 5, np.nan, data)
                data = xr.where(data < 0, np.nan, data)
            
            if var == 'twb_qnts':
                data = data.rename({'quantile': 'percentile'})
             
            # Calculate quantiles
            data_qnt = data.chunk(dict(GCM_RCP=-1)).quantile([0.05, 0.25, 0.5, 0.75, 0.95], dim='GCM_RCP') \
                        .assign_coords({'quantile': ['q5', 'q25', 'q50', 'q75', 'q95']}).rename({'quantile': 'stats'})
            
            # Calculate multi-model ensemble statistics
            data = xr.concat([data.min(dim='GCM_RCP'), data.mean(dim='GCM_RCP'), 
                              data.median(dim='GCM_RCP'), data.max(dim='GCM_RCP'), \
                              data.std(dim='GCM_RCP'), \
                              data.std(dim='GCM_RCP')/data.mean(dim='GCM_RCP')], dim='stats') \
                              .assign_coords({"stats": ['min', 'mean', 'median', 'max', 'stdev', 'rsd']})               
            data = xr.concat([data, data_qnt], dim='stats')
         
        data_all = xr.concat([data_all, data], dim='threshold')
      
    return data_all.assign_coords({'threshold': thrshld})

#------------------------------------------------------------------------------