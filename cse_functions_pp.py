# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:08:08 2022

@author: werning
"""

import os
import numpy as np
import itertools as it
import xarray as xr
from scipy.interpolate import interp1d
import yaml
import math 

# -----------------------------------------------------------------------------

def load_parameters(yaml_path):
    
    ''' Load parameters for score conversion and plots 
    
        Arguments:
            - yaml_path: path to yaml file containing parameters ''' 
    
    with open(yaml_path, 'r', encoding='utf8') as f:
        return yaml.full_load(f)
    
# -----------------------------------------------------------------------------
    
def load_landmask(landmask_path):
    
    ''' Load land mask 
    
        Arguments:
            - landmask_path: path to land mask file''' 
    
    with xr.open_dataset(landmask_path, engine="netcdf4") as ds:
        return ds.load()
    
#------------------------------------------------------------------------------

def apply_land_mask(data, land_mask):
    
    ''' Apply land mask 
    
        Arguments:
            - data: dataset
            - land_mask: land mask to only keep land area ''' 
    
    return data.where(land_mask['land area'] > 0)
   
#------------------------------------------------------------------------------

def set_interpolation(interpolation):
    
    ''' Set string to specify interpolation option if used
    
        Arguments:
            - interpolation: string specifying used interpolation'''

    if interpolation:
        return f'_{interpolation}'
    else:
        return ''
    
#------------------------------------------------------------------------------

def load_multi_model_means(input_dir, protocol, ind, version, SOC=''):
    
    ''' Load files containing multi-model means 
    
        Arguments:
            - input_dir: string with path to multi-model mean files
            - protocol: string specifying protocol (ISIMIP 2b/3b)
            - ind: string with indicator to load
            - version: future or historical
            - SOC: optional, only required for loading GHM results ''' 
    
    soc, hist = '', ''
    
    if version != 'future' and version != 'historical':
        print('No valid version detected. Please choose either future or historical.')
        
    else:
    
        if SOC:
            soc = f'_{SOC}'
        
        if version == 'historical':
            hist = '_historical'
            
        data_set = f'ISIMIP{protocol}_MM{hist}_{ind}{soc}.nc4'
           
        with xr.open_dataset(os.path.join(input_dir, data_set), engine="netcdf4") as ds:
            return ds.load()     

#------------------------------------------------------------------------------

def load_multi_model_files(input_dir, protocol, file_type, ind, version, interpolation='', SOC=''):
    
    ''' Load files containing multi-model means 
    
        Arguments:
            - input_dir: string with path to multi-model mean files
            - protocol: string specifying protocol (ISIMIP 2b/3b)
            - ind: string with indicator to load
            - version: future or historical
            - SOC: optional, only required for loading GHM results ''' 
    
    soc  = ''
    
    if version != 'future' and version != 'historical':
        print('No valid version detected. Please choose either future or historical.')
        
    else:
    
        if SOC:
            soc = f'_{SOC}'
        
        if version == 'historical':
            data_set = f'ISIMIP{protocol}_{file_type}_historical_{ind}{soc}.nc4'
        
        else: 
           data_set = f'ISIMIP{protocol}_{file_type}_{ind}{soc}{interpolation}.nc4' 
        
        
        
           
        with xr.open_dataset(os.path.join(input_dir, data_set), engine="netcdf4") as ds:
            return ds.load()

#------------------------------------------------------------------------------

def calculate_differences(dataset_1, dataset_2, stats, ind, variables, params):
    
    ''' Calculate absolute or relative difference between two datasets 
    
        Arguments:
            - dataset_1: first dataset
            - dataset_2: second dataset
            - stats: list with statistics to be plotted
            - option: string specifying if absolute or relative difference '''
                          
    difference = xr.Dataset()
    dataset_1 = dataset_1.sel({'stats': stats})
    dataset_2 = dataset_2.sel({'stats': stats})

    for var in variables:                  
     
        option = params[var]['diff']
    
        if option != 'absolute' and option != 'relative':
            print('No valid version detected. Please choose either absolute or relative.')
        
        elif option == 'absolute':
            difference[var] = dataset_1[var] - dataset_2[var]
                
        else:
            difference[var] = (dataset_1[var] - dataset_2[var]) / (dataset_2[var]) * 100
    
    return difference
    
#------------------------------------------------------------------------------    

def std_av(gcms, axis=None): 
    
    ''''
    Calculate the mean of a list of standard deviations (NaNs are ignored)
    
    Arguments:
        - gcms: standard deviations for all GCMs/GHMs
    
    '''
   
    len_of_non_nans = len([x for x in gcms if not math.isnan(x)])
    if len_of_non_nans == 0:
        av = np.nan
    else:
        mean_of_stds = np.nansum([x**2 for x in gcms])
        av = math.sqrt(mean_of_stds/len_of_non_nans)
    return av

#------------------------------------------------------------------------------

def calculate_standard_deviations(ind, protocol, land_mask, hist_dir, GCMs, GHMs):
    
    ''' Calculate standard deviations for reference period and each GCM/GHM 
    
    Arguments:
        - ind: string with indicator name
        - protocol: string with protocol (either 2b or 3b)
        - land_mask: xarray dataset with land mask
        - hist_dir: path to reference data
        - GCMs: list with GCMs
        - GHMs: list with GHMs ''' 
    
    hist_all = xr.Dataset()
    
    if protocol == '2b':
        
        rcp = 'rcp26'            
        for GHM, GCM in it.product(GHMs, GCMs):                
            with xr.open_dataset(os.path.join(hist_dir, GHM, rcp, GCM, \
            f'{GHM}_{GCM.lower()}_{rcp}_historical_histsoc_{ind}_all_global_1974_2004.nc4')) \
            as hist:  
                hist = apply_land_mask(hist, land_mask)
                if ind not in ['seas', 'seas_qtot', 'iavar', 'iavar_qtot']:
                    hist_std = hist.std(dim='year', skipna=True)
                    hist_all = xr.concat([hist_all, hist_std], dim='GCM_RCP')
                else:                    
                    hist_all = xr.concat([hist_all, hist], dim='GCM_RCP')            
        
    else:  
          
        rcp = 'ssp126'                        
        for GCM in GCMs:                
            with xr.open_dataset(os.path.join(hist_dir, rcp, GCM, \
            f'{GCM.lower()}_{rcp}_historical_{ind}_all_global_1974_2004.nc4')) \
            as hist:                    
                hist = apply_land_mask(hist, land_mask) 
                hist_std = hist.std(dim='year', skipna=True)
                hist_all = xr.concat([hist_all, hist_std], dim='GCM_RCP')
                
    hist_all = hist_all.assign_coords({'GCM_RCP': list(range(0, len(hist_all.GCM_RCP)))})       
    hist_all_mean = xr.apply_ufunc(std_av, hist_all, input_core_dims=[["GCM_RCP"]], kwargs={"axis": -1}, vectorize=True)
    hist_all_mean = apply_land_mask(hist_all_mean, land_mask)
    
    return hist_all_mean

#------------------------------------------------------------------------------

def bin_std_data(data, bins):
    
    ''' Bin and score relative component  
    
    Arguments:
        - data: data to bin and score
        - bins: thresholds to use for scoring''' 
       
    for var in list(data):
        
        data[var] = xr.where(data[var] > bins[3], bins[3], data[var])
        data[var] = xr.where(data[var] < bins[0], bins[0], data[var])
        
        m = interp1d(bins, [0,1,2,3])
        pt = m(data[var].values)
        data[var].values = pt    
    
    return data

#------------------------------------------------------------------------------

def bin_future_data(data, bins):
    
    ''' Bin and score absolute component  
    
    Arguments:
        - data: data to bin and score
        - bins: thresholds to use for scoring''' 
    
    if type(bins) != dict:    
        indic_bins = {var: list(bins[var].values) for var in list(data)}
    else:
        indic_bins = bins
    
    for var in list(data):
        
        data[var] = xr.where(data[var] > indic_bins[var][3], indic_bins[var][3], data[var])
        data[var] = xr.where(data[var] < indic_bins[var][0], indic_bins[var][0], data[var])
        
        m = interp1d(indic_bins[var], [0,1,2,3])
        pt = m(data[var].values)
        data[var].values = pt    
    
    return data       

#------------------------------------------------------------------------------

def calculate_kg_scores(future, hist, std_change_binned, ind, protocol, kg_class):  
    
    ''' Calculate scores using Koeppen-Geiger climate zones
    
    Arguments:
        - future: future data 
        - hist: reference period data
        - std_change_binned: binned relative component
        - ind: string with indicator name
        - protocol: string with protocol (either 2b or 3b)
        - kg_class: array with Koeppen-Geiger zones ''' 
                
    future_binned = xr.full_like(future, np.nan)
    future_binned_output = xr.full_like(future, np.nan)
    no_hist_binned = xr.full_like(future, np.nan)

    for k in range(1,6):
               
        kg_hist = hist.where(kg_class.kg_class == k)
        kg_future = future.where(kg_class.kg_class == k)
        kg_std_change_binned = std_change_binned.where(kg_class.kg_class == k)                
  
        kg_bins = kg_hist.quantile([0, 0.5, 0.75, 0.95],  dim=(['lat', 'lon']))
        kg_future_binned = bin_future_data(kg_future, kg_bins)
        future_binned_output = xr.where(kg_class.kg_class == k, kg_future_binned, future_binned_output)

        kg_future_bivariate = kg_std_change_binned + kg_future_binned    
        future_binned =  xr.where(kg_class.kg_class == k, kg_future_bivariate, future_binned)
        
        no_hist = kg_future.where((kg_hist == 0) & (kg_future !=0))
        no_hist_binned = bin_future_data(no_hist, kg_bins)   

    return future_binned, no_hist_binned

#------------------------------------------------------------------------------

def calculate_bivariate_scores(future, hist, std, ind, protocol, kg_class):
    
    ''' Calculate scores using Koeppen-Geiger climate zones
    
    Arguments:
        - future: future data 
        - hist: reference period data
        - std: standard deviation data
        - ind: string with indicator name
        - protocol: string with protocol (either 2b or 3b)
        - kg_class: array with Koeppen-Geiger zones ''' 
    
    # Bin standard deviation change to 0-3
    std_change = (future - hist) / std
    std_change = xr.where(hist == 0, 0, std_change)
    std_change_binned = bin_std_data(std_change, [0, 1, 2, 3])
    
    # For tropical nights, set relative score component to 0 if fewer than 
    # consecutive tropical nights for GWLs
    if ind == 'tr20':
        std_change_binned = xr.where(future <= 5, 0, std_change_binned)
        
    # Bin future data based on reference period percentiles for each Koeppen-Geiger
    # climate zone
    if ind in ['seas', 'seas_qtot', 'iavar', 'iavar_qtot']:
        future_bivariate, no_hist_binned = calculate_kg_scores(future, hist, std_change_binned, 
                                                               ind, protocol, kg_class)
    
    else:
        
        # Use established values to score absolute component for water stress index
        if ind == 'wsi':
            bins = {'wsi': [0.1, 0.2, 0.3, 0.4]}
            
        # Otherwise use percentiles of reference period
        else:
            bins = hist.quantile([0, 0.50, 0.75, 0.95], dim=(['lat', 'lon']))
        future_binned = bin_future_data(future, bins)            
            
        # Deal with historical data that is zero (but future is not)
        no_hist = future.where((hist == 0) & (future !=0))
        no_hist_binned = bin_future_data(no_hist, bins)
    
        # Create bivariate score
        future_bivariate = std_change_binned + future_binned
        
    future_bivariate = xr.where(no_hist_binned >= 0, no_hist_binned, future_bivariate) 

    return future_bivariate