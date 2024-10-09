# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:36:59 2024

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
sys.path.append('path_to_folder_with_repo\\paper_figures')
import matplotlib.pyplot as plt
import os
import glob
import cse_functions as cfp
import cse_plotting_funcs as cpf


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'split_files'
output_dir = ''
yaml_path = 'path_to_folder_with_repo\\required_files\\cse_params.yml'

# Set parameters 
indicators = {'cdd': ['cdd'], 'dri_qtot': ['dri_qtot'], 'heatwave': ['hw_95_7',
              'hw_99_3'], 'iavar_qtot': ['iavar_qtot'], 'precip': 
              ['pr_r20', 'pr_r99p'], 'seas_qtot': ['seas_qtot'], 
              'sdd': ['sdd_c'], 'tr20': ['tr20'], 'wsi': ['wsi']}
gwl = 2.0
fig_size = (16,12)
stat = 'median'

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')            
params = cfp.load_parameters(yaml_path)

#%% All scores plot
# -----------------------------------------------------------------------------

all_scores_fig = cpf.make_all_scores(indicators, input_dir, stat, 
                                     sum(len(v) for v in indicators.values()), 
                                     gwl)
all_scores_fig.savefig(os.path.join(output_dir, 
                                    f'All_scores_{str(gwl).replace(".", "p")}_{stat}.eps'), 
                       format='eps', dpi=300, bbox_inches='tight')