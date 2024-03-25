# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 09:45:54 2022

@author: werning
"""

from datetime import datetime

netcdf4_format = 'NETCDF4_CLASSIC' 

def cust_attrs_sdd(GCM, RCP, threshold, years, t_bal_cooling, t_bal_heating, protocol):

        return {'title': 'Simple Degree Days for heating and cooling: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning, E. Byers & A. Mastrucci',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'degree days',
                'set-point temperature (cooling, heating)': f'{t_bal_cooling}, {t_bal_heating}',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }        
                       
#------------------------------------------------------------------------------

def cust_attrs_pr(GCM, RCP, threshold, years, protocol):

        return {'title': 'Precipitation: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'r10 & r20: days, r95p & r99p: days, sdii: mm day-1',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_cdd(GCM, RCP, threshold, years, protocol):

        return {'title': 'Consecutive dry days: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'days',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_tr20(GCM, RCP, threshold, years, protocol):

        return {'title': 'Tropical nights: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'days',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }

#------------------------------------------------------------------------------

def cust_attrs_di(GHM, GCM, RCP, SOC, threshold, years, input_var, protocol):

        return {'title': 'Drought intensity, drought duration & cumulative drought deficit volume : Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org', 
                'input data': input_var,
                'GHM': GHM, 
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'unit': 'm3 s-1 day-1, day & m3 s-1', 
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 } 
                          
#------------------------------------------------------------------------------

def cust_attrs_twb(GCM, RCP, threshold, years, protocol):

        return {'title': 'Wet bulb temperature: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'degrees Celsius',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }

#------------------------------------------------------------------------------
    
def cust_attrs_twb_qnts(GCM, RCP, threshold, years, protocol, quantiles):

        return {'title': 'Wet bulb temperature quantiles: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'quantiles:': str(quantiles),
                'unit': 'degrees Celsius',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
#------------------------------------------------------------------------------

def cust_attrs_heatwave(GCM, RCP, threshold, years, protocol, dt, quantiles):

        return {'title': 'Quantiles and frequency of wet-bulb temperature heatwave events: Post-processed ISIMIP GCM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'no of days': str(dt),
                'quantiles:': str(quantiles),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_seasonality(GHM, GCM, RCP, SOC, input_var, threshold, years, protocol):

        return {'title': 'Seasonality: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GHM': GHM,
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_iavar(GHM, GCM, RCP, SOC, input_var, threshold, years, protocol):

        return {'title': 'Inter-annual variability: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GHM': GHM,
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_wsi(GHM, GCM, RCP, SOC, input_var, threshold, years, protocol):

        return {'title': 'Water stress index: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning, Y. Satoh & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': f'{input_var} created with script by Y. Satoh',
                'GHM': GHM,
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }

#------------------------------------------------------------------------------

def cust_attrs_ari(GHM, GCM, RCP, SOC, input_var, threshold, years, protocol):

        return {'title': 'Aridity index: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning, Y. Satoh & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': f'{input_var} created with adapted script by Y. Satoh',
                'GHM': GHM,
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }

#------------------------------------------------------------------------------

def cust_attrs_lc(GCM, RCP, threshold, years, protocol):

        return {'title': 'Land cover: Post-processed GLOBIOM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': 'GLOBIOM data provided by Stefan Frank (frank@iiasa.ac.at)',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'percentage',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_yields(GCM, RCP, threshold, years, protocol):

        return {'title': 'Crop yield data: Post-processed GLOBIOM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': 'GLOBIOM data provided by Stefan Frank (frank@iiasa.ac.at)',          
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'unit': 'ton/ha',
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                }
                 
#------------------------------------------------------------------------------

def cust_attrs_monthly_mean(GHM, GCM, RCP, SOC, input_var, threshold, years, protocol):

        return {'title': 'Monthly mean: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GHM': GHM,
                'GCM': GCM,        
                'scenario': RCP,
                'social forcing': SOC, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }  
    
#------------------------------------------------------------------------------

def cust_attrs_seasonality_pr(GCM, RCP, input_var, threshold, years, protocol):

        return {'title': 'Seasonality: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }
    
#------------------------------------------------------------------------------

def cust_attrs_iavar_pr(GCM, RCP, input_var, threshold, years, protocol):

        return {'title': 'Inter-annual variability: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GCM': GCM,        
                'scenario': RCP,
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 }

#------------------------------------------------------------------------------

def cust_attrs_monthly_mean_pr(GCM, RCP, input_var, threshold, years, protocol):

        return {'title': 'Monthly mean: Post-processed ISIMIP GHM data', 
                'file created by': 'M. Werning & E. Byers',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'input data': input_var,
                'GCM': GCM,        
                'scenario': RCP, 
                'warming': str(threshold),
                'years': f'{years[0]}-{years[1]}',
                'netCDF format': netcdf4_format,
                'software': 'File written using Python xarray'
                 } 

#------------------------------------------------------------------------------
    
def cust_attrs_mm(thresholds, protocol, variable, SOC=''):
    
    if SOC:
    
        return {'title': 'Multi-model statistics (min, mean, median, max, stdev, rel stdev): Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'social forcing': str(SOC),
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }
    
    else:
        
        return {'title': 'Multi-model statistics (min, mean, median, max, stdev, rel stdev): Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }

# -----------------------------------------------------------------------------
    
def cust_attrs_mm_lc(thresholds, protocol, variable, scen):
    
    return {'title': 'Multi-model statistics (min, mean, median, max, stdev, rel stdev): Post-processed ISI-MIP GHM data', 
            'file created by': 'M. Werning, E.Byers & colleagues',
            'contact': 'werning@iiasa.ac.at',
            'institution': 'IIASA',
            'date': str(datetime.now()), 
            'original source': 'GLOBIOM data provided by Stefan Frank (frank@iiasa.ac.at)',
            'indicator': variable,
            'scenario': str(scen),
            'warming': str(thresholds),
            'netCDF format': str(netcdf4_format),
            'software': 'File written using Python xarray'
            }
    
#------------------------------------------------------------------------------

def cust_attrs_diff(thresholds, protocol, variable, stats, SOC=''):
    
    if SOC:
    
        return {'title': 'Difference between historical and future scenarios: Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'statistics': str(stats),
                'social forcing': SOC,
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }
    
    else:
        
        return {'title': 'Difference between historical and future scenarios: Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'statistics': str(stats),
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }

#------------------------------------------------------------------------------
    
def cust_attrs_diff_lc(thresholds, protocol, variable, stats, scen):
    
    
    return {'title': 'Difference between historical and future scenarios: Post-processed GLOBIOM data', 
            'file created by': 'M. Werning, E.Byers & colleagues',
            'contact': 'werning@iiasa.ac.at',
            'institution': 'IIASA',
            'date': str(datetime.now()), 
            'original source': 'GLOBIOM data provided by Stefan Frank (frank@iiasa.ac.at)',
            'indicator': variable,
            'statistics': str(stats),
            'scenario': str(scen),
            'warming': str(thresholds),
            'netCDF format': str(netcdf4_format),
            'software': 'File written using Python xarray'
            }

#------------------------------------------------------------------------------
    
def cust_attrs_diff_yields(thresholds, protocol, variable, stats, scen):
    
    
    return {'title': 'Difference between historical and future scenarios: Post-processed GLOBIOM data', 
            'file created by': 'M. Werning, E.Byers & colleagues',
            'contact': 'werning@iiasa.ac.at',
            'institution': 'IIASA',
            'date': str(datetime.now()), 
            'original source': 'GLOBIOM data provided by Stefan Frank (frank@iiasa.ac.at)',
            'indicator': variable,
            'statistics': str(stats),
            'scenario': str(scen),
            'warming': str(thresholds),
            'netCDF format': str(netcdf4_format),
            'software': 'File written using Python xarray'
            }
    
#------------------------------------------------------------------------------

def cust_attrs_interp(thresholds, interp_thresholds, interp_int, protocol, variable, stats, SOC=''):
       
    return {'title': 'Interpolated indicator data: Post-processed ISI-MIP GHM data', 
            'file created by': 'M. Werning, E.Byers & colleagues',
            'contact': 'werning@iiasa.ac.at',
            'institution': 'IIASA',
            'date': str(datetime.now()), 
            'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
            'indicator': variable,
            'statistics': str(stats),
            'warming': str(thresholds),
            'interpolation interval': str(interp_int),
            'interpolated warming': str(interp_thresholds),
            'netCDF format': str(netcdf4_format),
            'software': 'File written using Python xarray'
            }
       
#------------------------------------------------------------------------------

def cust_attrs_scores(thresholds, protocol, variable, stats, SOC=''):
    
    if SOC:
    
        return {'title': 'Risk scores for future scenarios: Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'statistics': str(stats),
                'social forcing': SOC,
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }
    
    else:
        
        return {'title': 'Risk scores for future scenarios: Post-processed ISI-MIP GHM data', 
                'file created by': 'M. Werning, E.Byers & colleagues',
                'contact': 'werning@iiasa.ac.at',
                'institution': 'IIASA',
                'date': str(datetime.now()), 
                'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
                'indicator': variable,
                'statistics': str(stats),
                'warming': str(thresholds),
                'netCDF format': str(netcdf4_format),
                'software': 'File written using Python xarray'
                }  

#------------------------------------------------------------------------------

def cust_attrs_scores_lc(thresholds, protocol, variable, stats, scen):
      
    return {'title': 'Risk scores for future scenarios: Post-processed ISI-MIP GHM data', 
            'file created by': 'M. Werning, E.Byers & colleagues',
            'contact': 'werning@iiasa.ac.at',
            'institution': 'IIASA',
            'date': str(datetime.now()), 
            'original source': f'ISIMIP{protocol} GCM data see www.isimip.org',
            'indicator': variable,
            'statistics': str(stats),
            'scenario': str(scen),
            'warming': str(thresholds),
            'netCDF format': str(netcdf4_format),
            'software': 'File written using Python xarray'
            }
    
#------------------------------------------------------------------------------

def split_ind_attrs(ind, var, var_name, ssp, params, types, t, threshold, q='', d=''):
    
    if q and d:
        
        return {
            'authors': 'M. Werning, E. Byers & colleagues',
            'institution': 'IIASA',
            'contact': 'werning@iiasa.ac.at',
            'long_name':  f'{params["indicators"][ind][var]["long_name"]}',
            'short_name': f'{params["indicators"][ind][var]["short_name"]}_{int(q*100)}_{d}',
            'description': f'{params["indicators"][ind][var]["long_name"]} {types[t]["desc"]}',
            'percentile': q,
            'threshold (days)': d,
            'scen': f'{str(threshold)}°C',
            'scenario': f'{ssp}_{(str(threshold).replace(".", "p"))}',
            'model': 'climate-solutions',
            'layerid': f'{var_name}_{types[t]["short"]}',
            'variable': f'{var_name}|{types[t]["short"]}',
            }
    else:
    
        return {
            'authors': 'M. Werning, E. Byers & colleagues',
            'institution': 'IIASA',
            'contact': 'werning@iiasa.ac.at',
            'long_name':  f'{params["indicators"][ind][var]["long_name"]}',
            'short_name': f'{params["indicators"][ind][var]["short_name"]}',
            'description': f'{params["indicators"][ind][var]["long_name"]} {types[t]["desc"]}',
            'scen': f'{str(threshold)}°C',
            'scenario': f'{ssp}_{(str(threshold).replace(".", "p"))}',
            'model': 'climate-solutions',
            'layerid': f'{var_name}_{types[t]["short"]}',
            'variable': f'{var_name}|{types[t]["short"]}',
            }
    