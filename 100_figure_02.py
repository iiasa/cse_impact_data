# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:30:04 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import xarray as xr
import os
import matplotlib.pyplot as plt
import cse_plotting_functions as cpf
import cse_functions_pp as cfp


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
hist_path = 'path_to_historical_file_for_chosen_ensemble_member'
absolute_path = 'path_to_absolute_file_for_chosen_ensemble_member'
diff_path = 'path_to_difference_file_for_chosen_ensemble_member'
score_path = 'path_to_score_file_for_chosen_ensemble_member'
std_path = 'path_to_std_file_for_chosen_ensemble_member' # Can be written out when calculating the bivariate score
std_change_binned_path = 'path_to_std_change_binned_file_for_chosen_ensemble_member' # Can be written out when calculating the bivariate score
future_binned_path = 'path_to_future_binned_file_for_chosen_ensemble_member' # Can be written out when caluclating the bivariate score
output_dir = ''
land_mask_path = 'path_to_folder_with_repo\\required_files\\landareamaskmap0.nc'
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'      

# Set parameters
gwl = 2.0
fig_size = (16,12)

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')            
land_mask = cfp.load_landmask(land_mask_path)
params = cfp.load_parameters(yaml_path)

#%% Score methodology plot
#------------------------------------------------------------------------------ 
        
# Load datasets and apply land area mask   
hist = xr.open_dataset(hist_path)
hist = cfp.apply_land_mask(hist, land_mask = land_mask)
      
absolute = xr.open_dataset(absolute_path)
absolute = cfp.apply_land_mask(absolute, land_mask = land_mask)

diff = xr.open_dataset(diff_path)
diff = cfp.apply_land_mask(diff, land_mask = land_mask)

scores = xr.open_dataset(score_path)
std = xr.open_dataset(std_path)      
std_change = xr.open_dataset(std_change_binned_path)
future = xr.open_dataset(future_binned_path)      
           
# Create figure
fig, axs = cpf.create_fig_for_score_plot(fig_size)        
cpf.create_score_methodology_plot(hist['r20'], absolute['r20'], future['r20'], 
                              absolute['r20']-hist['r20'] , std['r20'], std_change['r20'],
                              scores['r20'], params, 'precip', 'r20', '', axs)
fig.suptitle(f'{params["indicators"]["precip"]["r20"]["long_name"]}', y=0.92)

# Save output
fig.savefig(os.path.join(output_dir, 
            f'Score_methodology_{params["indicators"]["precip"]["r20"]["short_name"]}_new.png'), 
            dpi=300, bbox_inches='tight')