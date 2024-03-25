# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 11:34:42 2022

@author: werning
"""

import xarray as xr
import numpy as np
import cse_functions as cf
import itertools as it
import statistics
import time
import calendar
import dask

# -----------------------------------------------------------------------------

def calculate_wetbulb(input_dir, GCM, RCP, input_var, timestep, years, land_mask):
     
    ''' Calculate wet bulb temperature using Stull equation 
    (http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1)
    
    Arguments:
        - input_dir: path to ISIMIP files
        - GCM: string with currently selected GCM
        - RCP: string with currently selected RCP
        - input_var: temperature data to be used - tas/tasmin/tasmax
        - timestep: specifying timestep of data in the input files
        - years: list consisting of the start and end year '''
    
    start = time.time()

        
    tas_file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
    tas = cf.load_data(tas_file_list, years)
    
    tas = tas.where(land_mask['land area'] > 0)   
    tas = cf.convert_to_Celsius(tas)
    
    hurs_file_list = cf.create_file_list(input_dir, GCM, RCP, 'hurs', timestep)
    hurs = cf.load_data(hurs_file_list, years)
    hurs = hurs.where(land_mask['land area'] > 0) 
    
    twb = tas[input_var] * np.arctan(0.151977 * np.power((hurs.hurs + 8.313659), 0.5)) + np.arctan(tas[input_var] + hurs.hurs) \
          - np.arctan(hurs.hurs - 1.676331) + 0.00391838 * np.power(hurs.hurs, 1.5) * np.arctan(0.023101 * hurs.hurs) - 4.686035
          
    end = time.time()
    print(end-start)         
    
    return twb  

# -----------------------------------------------------------------------------

def calc_dd_c(t_out_ave, t_bal_fix):
    
    ''' Calculate cooling degree days (fixed balance temperature)
    
    Arguments:
        - t_out_ave: xarray dataset with global daily temperature
        - t_bal_fix: balance temperature in Celsius '''
    
    ddc = t_out_ave.where(t_out_ave > t_bal_fix) - t_bal_fix
    return ddc

# -----------------------------------------------------------------------------

def calc_dd_h(t_out_ave, t_bal_fix):
    
    ''' Calculate heating degree days (fixed balance temperature)
    
        Arguments:
            - t_out_ave: xarray dataset with global daily temperature
            - t_bal_fix:  balance temperature in Celsius '''
    
    ddh = t_bal_fix - t_out_ave.where(t_out_ave < t_bal_fix)
    return ddh

# -----------------------------------------------------------------------------

def calculate_sdd(TAS, balance_temperature_cooling_low, balance_temperature_heating_low):
    
    ''' Calculate simple degree days 
    
        Arguments:
            - TAS: xarray dataset with global daily temperature
            - balance_temperature_cooling_low: balance temperature for cooling degree days in Celsius
            - balance_temperature_heating_low: balance temperature for heating degreee days in Celsius '''
    
    simple_degree_days_heating = calc_dd_h(TAS.tas, balance_temperature_heating_low)
    simple_degree_days_cooling = calc_dd_c(TAS.tas, balance_temperature_cooling_low)        
    simple_degree_days = xr.Dataset(data_vars={'sdd_h':simple_degree_days_heating, 'sdd_c':simple_degree_days_cooling})
    
    simple_degree_days_yearly = simple_degree_days.groupby('time.year').sum(dim='time', skipna=True)
    simple_degree_days_yearly['total'] = simple_degree_days_yearly.sdd_c + simple_degree_days_yearly.sdd_h
    simple_degree_days_yearly['iavar'] = simple_degree_days_yearly.total.std(dim='year', skipna=True) / simple_degree_days_yearly.total.mean(dim='year', skipna=True) 
    
    return simple_degree_days_yearly  

# -----------------------------------------------------------------------------

def calculate_annual_heatwave(twb, twb_quantiles, quantiles, dt, years):
    
    ''' Calculate heatwave  
    
        Arguments:
            - twb: xarray dataset with the wetbulb temperature
            - twb_quantiles: TWB quantiles to be used for the calculation
            - quantiles: list with quantiles to be used for the calculation
            - dt: list with days over threshold
            - years: list with start and end year'''

    lats=len(twb.lat)
    lons=len(twb.lon) 
    years = list(range(years[0], years[1]+1))
    
    data = xr.Dataset(data_vars=dict(no_of_events=(["lat", "lon", "year", "percentile", 'dt'], np.zeros([lats,lons,len(years), len(quantiles), len(dt)])),
                                     mean_dur=(["lat", "lon", "year", "percentile", 'dt'], np.zeros([lats,lons,len(years), len(quantiles), len(dt)])),
                                     median_dur=(["lat", "lon", "year", "percentile", 'dt'], np.zeros([lats,lons,len(years), len(quantiles), len(dt)])),
                                     dot=(["lat", "lon", "year", "percentile", 'dt'], np.zeros([lats,lons,len(years), len(quantiles), len(dt)]))), 
                      coords={'lon': twb.lon, 'lat': twb.lat, 'year': years, 'percentile': quantiles, 'dt': dt})    
   
    
    def heatwave(annual_data, dt):
        
        dot = 0
        no_of_events = 0
        mean_dur = 0
        median_dur = 0
        
        if annual_data.any()==True: 
                
            a = [ sum( 1 for _ in group ) for key, group in it.groupby( annual_data ) if key ]
                    
            if len([item for item in a if item >= dt]) > 0:
                dot = sum([item for item in a if item >= dt])
                mean_dur = statistics.mean([item for item in a if item >= dt])
                median_dur = statistics.median([item for item in a if item >= dt])
                no_of_events = len([item for item in a if item >= dt])

        return dot, no_of_events, mean_dur, median_dur
      
    for q in range(len(quantiles)):
        
        hot_days = (twb > twb_quantiles[q,:,:])
        
        for year, d in it.product(range(len(years)), range(len(dt))):
            
            selected_data = hot_days.sel(({'time': f'{years[year]}'}))
            data.dot[:,:,year,q,d], 
            data.no_of_events[:,:,year,q,d], 
            data.mean_dur[:,:,year,q,d], 
            data.median_dur[:,:,year,q,d] = xr.apply_ufunc(heatwave, selected_data, 
                                            dt[d], input_core_dims=[['time'],[]],  
                                            output_core_dims= [[] for _ in range(4)], 
                                            output_dtypes = [float, float, float, float], 
                                            vectorize = True, dask='allowed')
    
    return data

# -----------------------------------------------------------------------------

def calculate_ctr20(tas, temp_threshold, years):
    
    ''' Calculate consecutive tropical nights  
    
        Arguments:
            - tas: xarray dataset with temperature
            - temp_threshold: int with minimum nighttime temperature
            - years: list with start and end year'''
    
    lats=len(tas.lat)
    lons=len(tas.lon)
    years = list(range(years[0], years[1]+1))

    data = xr.Dataset(data_vars=dict(max_dur=(["lat", "lon", "year"], np.zeros([lats,lons,len(years)]))),
                                     coords={'lon': tas.lon, 'lat': tas.lat, 'year': years})
    
    def tr20(annual_data):
        
        if annual_data.any()==True:

            a = [ sum( 1 for _ in group ) for key, group in it.groupby( annual_data ) if key ]
            maximum = max(a)

        else:
            maximum= 0

        return maximum
    
    
    for y in range(0, len(years)):
        
        selected_data = ((tas >= temp_threshold)).sel(time=slice(f'{years[y]-1}-07-01', f'{years[y]+1}-06-30')).compute()
        data.max_dur[:,:,y] = xr.apply_ufunc(tr20, selected_data, input_core_dims=[['time']],  
                                             output_core_dims= [[]], output_dtypes = [float], 
                                             vectorize = True, dask='parallelized')
        
        # If more than 365/366 consecutive nights (due to two year period considered
        # during calculation), set to 365/366
        max_value = 366 if calendar.isleap(years[y]) == True else 365
        data.max_dur[:,:,y] = xr.where(data.max_dur[:,:,y] >= 366, max_value, data.max_dur[:,:,y])
                  
    return data

# -----------------------------------------------------------------------------
                
def calculate_annual_drought_intensity(data, quantiles, input_var):
    
    ''' Calculate drought intensity 
    
        Arguments:
            - data: xarray data array with global daily discharge data
            - quantiles: q90 quantile of the above data
            - input_var: string with input_var (dis/qtot)'''
    
    duration = data.where(data < quantiles)
    duration_annual_sum = duration.groupby('time.year').count(dim='time')
    
    deficit = (quantiles - data).where(quantiles > data)
    deficit_annual_sum = deficit.groupby('time.year').sum(dim='time')
            
    drought_intensity_annual = deficit_annual_sum / duration_annual_sum
    
    drought_intensity_annual = xr.where(duration_annual_sum < 182, np.nan, drought_intensity_annual)
    deficit_annual_sum = xr.where(duration_annual_sum < 182, np.nan, deficit_annual_sum)
    duration_annual_sum = xr.where(duration_annual_sum < 182, np.nan, duration_annual_sum)
    
    return xr.merge([drought_intensity_annual.rename({input_var: 'drought_intensity'}), \
                     deficit_annual_sum.rename({input_var: 'deficit'}), \
                     duration_annual_sum.rename({input_var: 'duration'})])   

# -----------------------------------------------------------------------------
    
def calculate_seasonality(data):
    
    ''' Calculate seasonality and inter-annual variability
    
        Arguments:
            - data: xarray data array with global daily discharge data '''
    
    monthly_mean = data.groupby('time.month').mean('time', skipna=True)
    annual_mean = data.groupby('time.year').mean('time', skipna=True)
    
    mean = annual_mean.mean('year', skipna=True)
    
    variance_seasonality = np.power((monthly_mean - mean),2).sum(dim='month') / len(monthly_mean.month)
    standard_deviation_seasonality = np.sqrt(variance_seasonality)
    coeff_of_variation_seasonality = standard_deviation_seasonality / mean
    
    variance_ia_variability = np.power((annual_mean - mean),2).sum(dim='year') / len(annual_mean.year)
    standard_deviation_ia_variability = np.sqrt(variance_ia_variability)
    coeff_of_variation_ia_variability = standard_deviation_ia_variability / mean
    
    return coeff_of_variation_seasonality, coeff_of_variation_ia_variability, monthly_mean

# -----------------------------------------------------------------------------

def calculate_annual_seasonality(data):
    
    ''' Calculate seasonality and inter-annual variability
    
        Arguments:
            - data: xarray data array with global daily discharge data '''
    
    monthly_mean = data.groupby('time.month').mean('time', skipna=True)
    annual_mean = data.groupby('time.year').mean('time', skipna=True)
    
    mean = annual_mean.mean('year', skipna=True)
    
    variance_seasonality = np.power((monthly_mean - mean),2).sum(dim='month') / len(monthly_mean.month)
    standard_deviation_seasonality = np.sqrt(variance_seasonality)
    coeff_of_variation_seasonality = standard_deviation_seasonality / mean
    
    variance_ia_variability = np.power((annual_mean - mean),2).sum(dim='year') / len(annual_mean.year)
    standard_deviation_ia_variability = np.sqrt(variance_ia_variability)
    coeff_of_variation_ia_variability = standard_deviation_ia_variability / mean
    
    return standard_deviation_seasonality, standard_deviation_ia_variability

# -----------------------------------------------------------------------------

def calculate_precipitation(pr, pr_qnts_hist, years):
    
    ''' Calculate precipitation indicators based on the ETCCDI Climate Change Indices
    
        Arguments:
            - pr: xarray data array with precipiation data in mm
            - pr_qnts_hist: xarray data array with the 0.95/0.99 quantiles of the historic data
            - years: list containing the start and end year '''
    
    pr = pr.pr * 86400
    
    lats=len(pr.lat)
    lons=len(pr.lon)
    year=list(range(years[0], years[1]+1))
    
    data = xr.Dataset(data_vars=dict(r10=(["lat", "lon", "year"], np.zeros([lats,lons,len(year)])), 
                                     r20=(["lat", "lon", "year"], np.zeros([lats,lons,len(year)])),
                                     r95p=(["lat", "lon", "year"], np.zeros([lats,lons,len(year)])),
                                     r99p=(["lat", "lon", "year"], np.zeros([lats,lons,len(year)])),
                                     sdii=(["lat", "lon", "year"], np.zeros([lats,lons,len(year)]))),
                                     coords={'lon': pr.lon, 'lat': pr.lat, 'year': year})
    
    # R10mm - Annual count of days when RCP >= 10mm
    data['r10'] = pr.where(pr >= 10).groupby('time.year').count()

    # R20mm - Annual count of days when RCP >= 20mm
    data['r20'] = pr.where(pr >= 20).groupby('time.year').count()
    
    # R95pTOT - Annual total PRCP when RR > 95p
    data['r95p'] = pr.where((pr >= 1) & (pr > pr_qnts_hist.pr[0,:,:])).groupby('time.year').sum()
    
    # R99pTOT - Annual total PRCP when RR > 99p
    data['r99p'] = pr.where((pr >= 1) & (pr > pr_qnts_hist.pr[1,:,:])).groupby('time.year').sum()
    
    # SDII - Simple precipitation intensity index
    data['sdii'] = (pr.where(pr >= 1).groupby('time.year').sum()) / (pr.where(pr >= 1).groupby('time.year').count())
    
    return data

# -----------------------------------------------------------------------------

def calculate_cdd(data, years):
    
    ''' Calculate consecutive dry days  
    
        Arguments:
            - data: xarray dataset with precipitation data
            - years: list with start and end year'''
    
    lats = len(data.lat)
    lons = len(data.lon)
    years = list(range(years[0], years[1]+1))    
    
    data_mm = data * 86400
    dry_days = (data_mm.pr < 1)
    
    cdd = xr.Dataset(data_vars=dict(cdd=(["lat", "lon", "year"], np.full([lats, lons, len(years)], np.nan))),
                                    coords={'year': years, 'lon': data_mm.lon, 'lat': data_mm.lat})
    
    def calc_cdd(annual_data):
        
        if annual_data.any()==True:
            a = [ sum( 1 for _ in group ) for key, group in it.groupby( annual_data ) if key ]
            maximum = max(a)
        else:
            maximum= np.nan
        return maximum
    
    start = time.time()
    
    for year in range(len(years)):  
        
        selected_data = dry_days.sel(time=slice(f'{years[year]-1}-07-01', f'{years[year]+1}-06-30')).compute()
        cdd.cdd[:,:,year] = xr.apply_ufunc(calc_cdd, selected_data, input_core_dims=[['time']], output_core_dims= [[]], output_dtypes = [float], vectorize = True, dask='parallelized')                
        
        # If more than 365/366 consecutive nights (due to two year period considered
        # during calculation), set to 365/366
        max_value = 366 if calendar.isleap(years[year]) == True else 365
        cdd.cdd[:,:,year] = xr.where(cdd.cdd[:,:,year] >= 366, max_value, cdd.cdd[:,:,year])
        
    return cdd    