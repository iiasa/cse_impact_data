# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:04:36 2024

@author: werning
"""

import xarray as xr
import pandas as pd
import glob 
import os
import seaborn as sns
import matplotlib.pyplot as plt
from natsort import natsorted


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''

# Set parameters
indicators = ['precip', 'cdd', 'dri', 'dri_qtot', 'seas', 'seas_qtot', 'iavar',
              'iavar_qtot', 'wsi', 'heatwave', 'tr20', 'sdd', 'sdd_18p3', 
              'sdd_20p0', 'sdd_24p0']
ind_list = ['Heavy precipitation days', 'Very heavy precipitation days*', 
            'Wet days', 'Very wet days*', 'Precipitation intensity index', 
            'Consecutive dry days*', 'Drought intensity (Discharge)', 
            'Drought intensity (Runoff)*', 'Seasonality (Discharge)', 
            'Seasonality (Runoff)*', 'Inter-annual variability (Discharge)', 
            'Inter-annual variability (Runoff)*', 'Water stress index*', 
            'Heatwave events - 95th perc. - 3 days', 'Heatwave events - 95th perc. - 5 days', 
            'Heatwave events - 95th perc. - 7 days*', 'Heatwave events - 95th perc. - 10 days',
            'Heatwave events - 97th perc. - 3 days', 'Heatwave events - 97th perc. - 5 days',
            'Heatwave events - 97th perc. - 7 days', 'Heatwave events - 97th perc. - 10 days',
            'Heatwave events - 99th perc. - 3 days*', 'Heatwave events - 99th perc. - 5 days',
            'Heatwave events - 99th perc. - 7 days', 'Heatwave events - 95th perc. - 10 days',
            'Heatwave days - 95th perc. - 3 days', 'Heatwave days - 95th perc. - 5 days', 
            'Heatwave days - 95th perc. - 7 days', 'Heatwave days - 95th perc. - 10 days',
            'Heatwave days - 97th perc. - 3 days', 'Heatwave days - 97th perc. - 5 days',
            'Heatwave days - 97th perc. - 7 days', 'Heatwave days - 97th perc. - 10 days',
            'Heatwave days - 99th perc. - 3 days', 'Heatwave days - 99th perc. - 5 days',
            'Heatwave days - 99th perc. - 7 days', 'Heatwave days - 95th perc. - 10 days',
            'Tropical nights*', 'Cooling degree days*', 'Cooling degree days (18.3 °C)',
            'Cooling degree days (20 °C)', 'Cooling degree days (24 °C)']
gwl = 2.0
stat = 'median'

#%% Correlation plot
#------------------------------------------------------------------------------

count = 0
all_data = pd.DataFrame()

# Load data
for ind in indicators:
       
    files = natsorted(glob.glob(os.path.join(input_dir, ind, stat, 
                      f'*_ssp2_{str(gwl).replace(".", "p")}_score.nc4')))
            
    for f in files: 
        data = xr.open_dataset(f)
        data = data.rename({list(data)[0]: list(data)[0].replace(f'_ssp2_{str(gwl).replace(".", "p")}_score', '')})
        data_df = data.to_dataframe()[[list(data)[0]]]        
        
        if count == 0:
            all_data = data_df
        else: 
            all_data = all_data.join(data_df)            
        count += 1
     
# Calculate correlation
all_data.columns = ind_list
corrs = all_data.corr().reset_index()
corrs_melt = pd.melt(corrs, id_vars='index')
corrs_melt.columns = ['x', 'y', 'value']

# Prepare plot
n_colors = 256 
palette = sns.diverging_palette(20, 220, n=n_colors) 
color_min, color_max = [-1, 1] 
sns.set(rc={'figure.figsize':(12,12), 'axes.facecolor':'white'})
sns.set(font_scale=0.7)
sns.set_style("white")

# Plot data
p = sns.relplot(x=corrs_melt['x'], y=corrs_melt['y'],
                hue=corrs_melt['value'], palette='RdBu',
                hue_norm=(-1, 1), size_norm=(-1, 1),
                marker="s", legend=False, linewidth=0,
                aspect=1.75)
ghost = p.ax.scatter([], [], c=[], vmin=-1, vmax=1, cmap='RdBu')
p.fig.colorbar(ghost).outline.set_visible(False)
p.ax.tick_params(axis='x', rotation=90)
p.ax.set(xlabel=None)
p.ax.set(ylabel=None)
p.despine(left=True, bottom=True)
p.savefig(os.path.join(output_dir, f'Correlation.eps'), format='eps', dpi=300)