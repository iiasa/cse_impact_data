# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:44:49 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import matplotlib.pyplot as plt
import itertools as it
import os
import glob
import numpy as np
import cse_functions_pp as cfp
import cse_plotting_funcs as cpf


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'
 
# Set parameters 
indicators = ['precip']
percentiles = [0.95, 0.97, 0.99]
dt = [3, 5, 7, 10]
fig_size = (16,12)
stat = 'median'

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')            
params = cfp.load_parameters(yaml_path)

#%% Indicator plot
#------------------------------------------------------------------------------

for ind in indicators:
    for var in params['indicators'][ind]:  
                
        if ind == 'heatwave':
            
            for p, d in it.product(percentiles, dt):
                
                files = cpf.get_data(ind, f"{params['indicators'][ind][var]['short_name']}_{int(p*100)}_{d}", input_dir)
                rows = len(files)       
                unit = '' if params['indicators'][ind][var]['unit'] == '-' else params['indicators'][ind][var]['unit']
                
                fig = cpf.make_multi_maps(files, rows,
                                    title= f'{params["indicators"][ind][var]["long_name"]} - {int(p*100)}th percentile - {d} days', 
                                    cmaps=[params['indicators'][ind][var]['ind_cmap'], 
                                            params['indicators'][ind][var]['diff_cmap'],
                                            'magma_r'], 
                                    clims=[[params['indicators'][ind][var]['ind_min'], 
                                            params['indicators'][ind][var]['ind_max']], 
                                            [params['indicators'][ind][var]['diff_min'], 
                                            params['indicators'][ind][var]['diff_max']], 
                                            [0, 6]],
                                    cticks=[[params['indicators'][ind][var]['ind_min'],  
                                              params['indicators'][ind][var]['ind_max']], 
                                            [params['indicators'][ind][var]['diff_min'], 
                                              params['indicators'][ind][var]['diff_max']], 
                                            [0, 6]],
                                    clabels=[unit, 'Difference (%)', 'Score']
                                    )
                
                fig.savefig(os.path.join(output_dir, 
                            f"Indicator_{params['indicators'][ind][var]['short_name']}_{int(p*100)}_{d}.png"), 
                            dpi=300)
                
        else:
        
            all_files = glob.glob(os.path.join(input_dir, ind, stat, f'*{params["indicators"][ind][var]["short_name"]}*ssp2*.nc4'))
            nthreshs = len(all_files) // 3
            files = np.array(sorted(all_files)).reshape(nthreshs, 3).tolist()
            rows = len(files)       
            unit = '' if params['indicators'][ind][var]['unit'] == '-' else params['indicators'][ind][var]['unit_SI']
            
            fig = cpf.make_multi_maps(files, rows,
                                title= params['indicators'][ind][var]['long_name'], 
                                cmaps=[params['indicators'][ind][var]['ind_cmap'], 
                                        params['indicators'][ind][var]['diff_cmap'],
                                        'magma_r'], 
                                clims=[[params['indicators'][ind][var]['ind_min'], 
                                        params['indicators'][ind][var]['ind_max']], 
                                        [params['indicators'][ind][var]['diff_min'], 
                                        params['indicators'][ind][var]['diff_max']], 
                                        [0, 6]],
                                cticks=[[params['indicators'][ind][var]['ind_min'],  
                                          params['indicators'][ind][var]['ind_max']], 
                                        [params['indicators'][ind][var]['diff_min'], 
                                          params['indicators'][ind][var]['diff_max']], 
                                        [0, 6]],
                                clabels=[unit, '%', '']
                                )
        
            fig.savefig(os.path.join(output_dir, 
                        f"Indicator_{params['indicators'][ind][var]['short_name']}.eps"), 
                        format='eps', dpi=300, bbox_inches='tight')