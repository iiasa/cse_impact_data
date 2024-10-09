# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:26:42 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import xarray as xr
import matplotlib.pyplot as plt
import os
import glob
import seaborn as sns
import pandas as pd
import cse_functions as cfp


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''
pop_dir = 'population_data'
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'

# Set parameters
colours = sns.color_palette(['#FCDE9C', '#F0746E', '#DC3977', '#7C1D6F'])
indicators = {'cdd': ['cdd'], 'dri': ['dri_qtot'], 'heatwave': ['hw_95_7', 'hw_99_3'], 
              'iavar': ['iavar_qtot'], 'precip': ['pr_r20', 'pr_r99p'], 
              'sdd': ['sdd_c'], 'seas': ['seas_qtot'], 'tr20': ['tr20'], 
              'wsi': ['wsi']}
gwls = [1.5, 2.0, 3.0]
impact = 1
ssp = 'ssp2'
year = '2050'
stat = 'median'

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')            
pop_data = xr.open_dataset(os.path.join(pop_dir, f'{ssp}_{year}_total_hd.nc4'))
land_mask = cfp.load_landmask(land_mask_path)
land_mask = land_mask / 1000000

# Prepare variables
no_of_ind = sum(len(v) for v in indicators.values())
no_of_impacts = list(range(0, no_of_ind+1))
exposed_land_all = {}
exposed_pop_all = {}
df = pd.DataFrame(columns=['no_impacts','gwl', 'stat', 'exp_land', 'exp_pop'])

#%% Exposure to indicators plot
#------------------------------------------------------------------------------

# Load data and calculate number of indicators per grid cell
for gwl in enumerate(gwls):
    
    files = sum(sum([[glob.glob(os.path.join(input_dir, ind, stat, f'*{var}_ssp2_{str(gwl[1]).replace(".", "p")}_score.nc4')) 
                  for var in indicators[ind]] for ind in indicators], []), [])
    scores = xr.open_mfdataset(files, compat='override')
    impacts = xr.where(scores >= impact, 1, 0)
    sum_impacts = impacts.to_array().sum('variable')
    sum_impacts = cfp.apply_land_mask(sum_impacts, land_mask)

    for i in no_of_impacts:
        df.loc[len(df.index)] = [i, gwl[1], 'mean', land_mask['land area'].where(sum_impacts == i).sum().values.item(),
                              pop_data[f'{ssp}_{year}'].where(sum_impacts == i).sum().values.item()]

# Add percentage for exposed land area and exposed population       
df['exp_land_perc'] = df['exp_land'] / land_mask['land area'].sum().item() * 100
df['exp_pop_perc'] = df['exp_pop'] / pop_data[f'{ssp}_{year}'].sum().item() * 100

# Prepare plot
sns.set_theme()
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
fig.suptitle('Exposure to indicators with at least low hazard', y=1.05)

# Exposed land area subplot
sns.barplot(ax=axes[0], data=df[(df.no_impacts < 9)],  x='no_impacts', y='exp_land_perc', hue='gwl', palette=colours)
axes[0].set_xlabel('Number of indicators', labelpad=10)
axes[0].set_ylabel('%', labelpad=10)
axes[0].set_title('Exposed land area')
axes[0].get_legend().remove()

# Exposed popoulation subplot
sns.barplot(ax=axes[1], data=df[(df.no_impacts < 9)],  x='no_impacts', y='exp_pop_perc', hue='gwl', palette=colours)
axes[1].set_xlabel('Number of indicators', labelpad=10)
axes[1].set_ylabel('', labelpad=10)
axes[1].set_title('Exposed population')
axes[1].get_legend().remove()

handles, labels = axes[1].get_legend_handles_labels()
lgd = fig.legend(handles, [f'{gwl} °C' for gwl in gwls], loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.15))
fig.savefig(os.path.join(output_dir, 
            f'Exposure_{impact}_{stat}_fix_wsi.png'), 
            bbox_inches='tight', dpi=300)
