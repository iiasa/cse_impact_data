# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:59:42 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import xarray as xr
import glob 
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cmcrameri import cm
import cse_plotting_funcs as cpf
import cse_functions as cfp


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
no_of_ind = sum(len(v) for v in indicators.values())
gwls = [1.2, 1.5, 2.0, 3.0]
impact = 3
ssp = 'ssp2'
year = 2050
stat = 'mean'

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')
land_mask = cfp.load_landmask(land_mask_path)
land_mask = land_mask / 1000000
pop = xr.open_dataset(os.path.join(pop_dir, f'{ssp}_{year}_total_hd.nc4'))

# Prepare variables
pop_dens = pop[f'{ssp}_{year}'] / land_mask['land area']
pop_dens_sel = pop_dens.where(pop_dens >= 10)

#%% Population density plot
#------------------------------------------------------------------------------

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(13, 10), 
                       subplot_kw={"projection": ccrs.Robinson()},
                       layout = 'constrained')

for gwl in enumerate(gwls):

    files = sum(sum([[glob.glob(os.path.join(input_dir, ind, stat, f'*{var}_ssp2_{str(gwl[1]).replace(".", "p")}_score.nc4')) 
                  for var in indicators[ind]] for ind in indicators], []), [])
    scores = xr.open_mfdataset(files, compat='override')
    impacts = xr.where(scores >= impact, 1, 0)
    sum_impacts = impacts.to_array().sum('variable')
    sum_impacts = cfp.apply_land_mask(sum_impacts, land_mask)
    sum_impacts_sel = sum_impacts.where((sum_impacts >= 4) & (pop_dens >= 10))
    
    p = cpf.make_plot_dens_map(
            pop_dens_sel,
            ax=ax.flatten()[gwl[0]],
            title='',
            cmap='Greys',
            cbar=False,
            plotargs={'norm': colors.LogNorm()})
    
    if gwl[0] == (len(gwls)-1):
        fig.colorbar(p, ax=ax.ravel().tolist(), location='bottom',
                      label='People km$^{-2}$', shrink=0.35)
    
    p = cpf.make_plot_dens_map(
            sum_impacts_sel,
            ax=ax.flatten()[gwl[0]],
            title=f'{gwl[1]} °C',
            cmap=cm.lipari_r,
            cbar=False,
            add_borders=True,
            plotargs={'vmin': 0, 'vmax': 12, 'levels': 13, 'alpha': 0.75})
    
    if gwl[0] == (len(gwls)-1):
        fig.colorbar(p, ax=ax.ravel().tolist(), location='bottom',
                      label='Number of indicators', shrink=0.35)
        
fig.suptitle(f'Number of indicators with at least medium hazard in densely populated areas', 
             fontsize=16)
fig.savefig(os.path.join(output_dir, f'Pop_density_{impact}_{stat}.png'), 
            dpi=300, bbox_inches='tight')   
