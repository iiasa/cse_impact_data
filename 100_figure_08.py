# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:16:47 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import xarray as xr
import matplotlib.pyplot as plt
import os
import glob
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import cse_functions_pp as cfp

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''
pop_dir = 'population_data'
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'

# Set parameters
indicators = {'cdd': ['cdd'], 'dri': ['dri_qtot'], 'heatwave': ['hw_95_7', 'hw_99_3'], 
              'iavar': ['iavar_qtot'], 'precip': ['pr_r20', 'pr_r99p'], 
              'sdd': ['sdd_c'], 'seas': ['seas_qtot'], 'tr20': ['tr20'], 
              'wsi': ['wsi']}
gwls = [1.2, 1.5, 2.0, 3.0]
impact = 3
no_of_ind = sum(len(v) for v in indicators.values())
no_of_impacts = list(range(0, no_of_ind+1))
ssps = ['ssp1', 'ssp2', 'ssp3', 'ssp4', 'ssp5']
years = ['2030', '2050', '2070']
stat = 'median'

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')
land_mask = cfp.load_landmask('H:\\git\\climate_impacts_processing\\landareamaskmap0.nc')
land_mask = land_mask / 1000000

# Prepare variables
df = pd.DataFrame(columns=['year','gwl','ssp','value'])

#%% Exposure plot
#------------------------------------------------------------------------------

# Calculate required data
for ssp in ssps:
    
    for gwl in enumerate(gwls):
                
        files = sum(sum([[glob.glob(os.path.join(input_dir, ind, stat, f'*{var}_ssp2_{str(gwl[1]).replace(".", "p")}_score.nc4')) 
                      for var in indicators[ind]] for ind in indicators], []), [])
        scores = xr.open_mfdataset(files, compat='override')
        impacts = xr.where(scores >= impact, 1, 0)
        sum_impacts = impacts.to_array().sum('variable')
        sum_impacts = cfp.apply_land_mask(sum_impacts, land_mask)
        
        for y in years:
            pop_data = xr.open_dataset(os.path.join(pop_dir, f'{ssp}_{y}_total_hd.nc4'))
            df.loc[len(df.index)] = [y, gwl[1], ssp, pop_data[f'{ssp}_{y}'].where(sum_impacts >= 4).sum().values.item()]

# Convert to billion and get min/max values
df['value'] = df['value'] / 1000000000
df_min = df[['year', 'gwl', 'value']].groupby(by=['year', 'gwl']).min().reset_index()
df_max = df[['year', 'gwl', 'value']].groupby(by=['year', 'gwl']).max().reset_index()          

# Create figure
plt.figure(figsize=(10, 5))
sns.set_theme()
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

# Plot bars for first subplot
bar1 = sns.barplot(ax=axes[0], x="gwl",  y="value", data=df[(df.ssp == 'ssp2') & (df.year == '2070')], color='#7CCBA2')
bar2 = sns.barplot(ax=axes[0], x="gwl",  y="value", data=df[(df.ssp == 'ssp2') & (df.year == '2050')], color='#089099')
bar3 = sns.barplot(ax=axes[0], x="gwl",  y="value", data=df[(df.ssp == 'ssp2') & (df.year == '2030')], color='#045275')

# Calculate error bars for first subplot
errors = [list(df[(df.ssp == 'ssp2')].sort_values(by='year', ascending=False)['value'].reset_index(drop=True) - 
                df_min.sort_values(by='year', ascending=False)['value'].reset_index(drop=True)),
          list(df_max.sort_values(by='year', ascending=False)['value'].reset_index(drop=True) - 
                          df[(df.ssp == 'ssp2')].sort_values(by='year', ascending=False)['value'].reset_index(drop=True))]
x_coords = [p.get_x() + 0.5 * p.get_width() for p in axes[0].patches]
x_coords = [0.25, 1.25, 2.25, 3.25, 0.0, 1.0, 2.0, 3.0, -0.25, 0.75, 1.75, 2.75]
y_coords = [p.get_height() for p in axes[0].patches]  
bar1.errorbar(x=x_coords, y=y_coords, yerr = errors, fmt='none', capsize=6, c='black')

# Set labels for first subplot
bar3.set_xlabel('Global warming level', labelpad=10)
bar3.set_ylabel('Exposed population [billion]', labelpad=10)
bar3.set_xticks(range(len(gwls)), labels=[f'{g} °C' for g in gwls])

# Add legend for first subplot
top_bar = mpatches.Patch(color='#045275', label='2030')
middle_bar = mpatches.Patch(color='#089099', label='2050')
bottom_bar = mpatches.Patch(color='#7CCBA2', label='2070')
axes[0].legend(handles=[top_bar, middle_bar, bottom_bar], loc='lower center',ncol=3, bbox_to_anchor=(0.5,-0.35))

# Plot bars for the second subplot
bar1 = sns.barplot(ax=axes[1], x="year",  y="value", data=df[(df.ssp == 'ssp2') & (df.gwl == 3.0)], color='#7C1D6F')
bar2 = sns.barplot(ax=axes[1], x="year",  y="value", data=df[(df.ssp == 'ssp2') & (df.gwl == 2.0)], color='#DC3977')
bar3 = sns.barplot(ax=axes[1], x="year",  y="value", data=df[(df.ssp == 'ssp2') & (df.gwl == 1.5)], color='#F0746E')
bar4 = sns.barplot(ax=axes[1], x="year",  y="value", data=df[(df.ssp == 'ssp2') & (df.gwl == 1.2)], color='#FCDE9C')

# Calculate error bars for second subplot
errors = [list(df[(df.ssp == 'ssp2')].sort_values(by='gwl', ascending=False)['value'].reset_index(drop=True) - 
                df_min.sort_values(by='gwl', ascending=False)['value'].reset_index(drop=True)),
          list(df_max.sort_values(by='gwl', ascending=False)['value'].reset_index(drop=True) - 
                          df[(df.ssp == 'ssp2')].sort_values(by='gwl', ascending=False)['value'].reset_index(drop=True))]
x_coords = [p.get_x() + 0.5 * p.get_width() for p in axes[0].patches]
x_coords = x_coords = [0.3, 1.3, 2.3, 0.1, 1.1, 2.1, -0.1, 0.9, 1.9, -0.2, 0.8, 1.8]
y_coords = [p.get_height() for p in axes[1].patches]  
bar1.errorbar(x=x_coords, y=y_coords, yerr = errors, fmt='none', capsize=6, c='black')

# Set labels for second subplot
bar4.set_xlabel('Year', labelpad=10)
bar4.set_ylabel('', labelpad=10)

# Add legend for second subplot
bar1_legend = mpatches.Patch(color='#FCDE9C', label='1.2 °C')
bar2_legend = mpatches.Patch(color='#F0746E', label='1.5 °C')
bar3_legend = mpatches.Patch(color='#DC3977', label='2.0 °C')
bar4_legend = mpatches.Patch(color='#7C1D6F', label='3.0 °C')
axes[1].legend(handles=[bar1_legend, bar2_legend, bar3_legend, bar4_legend], 
            loc='lower center',ncol=4, bbox_to_anchor=(0.5,-0.35))
plt.savefig(os.path.join(output_dir, f'Exposure_ssp_year_gwl_{stat}.png'), 
            format='png', bbox_inches='tight')
