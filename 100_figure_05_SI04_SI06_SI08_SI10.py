# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:58:25 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import xarray as xr
import os
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cmcrameri import cm
import cse_plotting_funcs as cpf
import cse_functions as cfp


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'

# Set parameters
indicators = {'cdd': ['cdd'], 'dri': ['dri_qtot'], 'heatwave': ['hw_95_7', 'hw_99_3'], 
              'iavar': ['iavar_qtot'], 'precip': ['pr_r20', 'pr_r99p'], 
              'sdd': ['sdd_c'], 'seas': ['seas_qtot'], 'tr20': ['tr20'], 
              'wsi': ['wsi']}
no_of_ind = sum(len(v) for v in indicators.values())
gwls = [1.2, 1.5, 2.0, 3.0]
stat = 'median'
impact = 5

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')            
land_mask = cfp.load_landmask(land_mask_path)

#%% Number of indicators plot
#------------------------------------------------------------------------------

fig, ax = plt.subplots(
      ncols=2, nrows=2, figsize=(12, 9), 
      subplot_kw={"projection": ccrs.Robinson()},
      layout = 'constrained')

for gwl in enumerate(gwls):

    files = sum(sum([[glob.glob(os.path.join(input_dir, ind, stat, 
                      f'*{var}_ssp2_{str(gwl[1]).replace(".", "p")}_score.nc4')) 
                      for var in indicators[ind]] for ind in indicators], []), [])
    scores = xr.open_mfdataset(files, compat='override')
    impacts = xr.where(scores >= impact, 1, 0)
    sum_impacts = impacts.to_array().sum('variable')
    sum_impacts = cfp.apply_land_mask(sum_impacts, land_mask)

    p = cpf.make_discrete_map_from_da(
                    sum_impacts,
                    levels = no_of_ind+2,
                    ax=ax.flatten()[gwl[0]],
                    title=f'{gwl[1]} °C',
                    cmap=cm.lipari_r,
                    cbar=True if gwl[0] == 4 else False,
                    color_lims=[0, no_of_ind+1],                    
                    )            
            
fig.colorbar(p, ax=ax.ravel().tolist(), location='bottom',
              label='Number of indicators', shrink=0.35)
fig.suptitle(f'Number of indicators with a high or very high hazard', fontsize=16)
fig.savefig(os.path.join(output_dir, f'No_of_indicators_impact_{impact}_{stat}.png'), 
            dpi=300, bbox_inches='tight')   
