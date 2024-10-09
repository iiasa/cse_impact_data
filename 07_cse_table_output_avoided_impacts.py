# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 10:17:08 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import numpy as np
import pandas as pd
import os
import itertools as it

#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths to files
input_dir = 'table_outputs'
output_dir = ''

# Set parameters
indicators = ['cdd', 'dri', 'dri_qtot', 'heatwave', 'iavar', 'iavar_qtot',
               'hw_95_3', 'hw_95_5', 'hw_95_7', 'hw_95_10',
               'hw_97_3', 'hw_97_5', 'hw_97_7', 'hw_97_10'
               'hw_99_3', 'hw_99_5', 'hw_99_7', 'hw_99_10',
               'hwd_95_3', 'hwd_95_5', 'hwd_95_7', 'hwd_95_10'
               'hwd_97_3', 'hwd_97_5', 'hwd_97_7', 'hwd_97_10',
               'hwd_99_3', 'hwd_99_5', 'hwd_99_7', 'hwd_99_10',
               'pr_r10', 'pr_r20', 'pr_r95p', 'pr_r99p', 'sdii',
               'seas', 'seas_qtot', 'sdd', 'sdd_24p0', 'sdd_20p0', 
               'sdd_18p3', 'tr20', 'wsi']
baseline = 1.5
years = np.arange(2020,2101,10)
ssps = ['ssp1', 'ssp2', 'ssp3', 'ssp4', 'ssp5']

# Select output mode - COUNTRIES/IPCC/R10/R5
mode = 'COUNTRIES'

#%% Calculation
# -----------------------------------------------------------------------------

previous = pd.DataFrame()

# Loop through indicators  
for ind in indicators:
    
    data_list = []
    
    # Read in table output for indicator and select the 50.0th percentile
    data = pd.read_csv(os.path.join(input_dir, f'table_output_{ind}_{mode}.csv'))
    data_sel = data[(data.Model == 'Climate Solutions') & data.Variable.str.contains('50.0th Percentile') ]
    
    # Select exposure data only
    exposure = data_sel.loc[data_sel.Variable.str.contains('Exposure', regex=False),:]
    
    # Loop through all variables, regions and ssps
    for var, r, ssp in it.product(exposure.Variable.unique(), exposure.Region.unique(), ssps):

        # Set GWLs 
        if ind == 'lc' and ssp == 'ssp1':
            GWLs = [1.5, 2.0, 2.5]
        elif ind in ['dri', 'dri_qtot', 'iavar', 'iavar_qtot', 'seas', 'seas_qtot', 'wsi']:
            GWLs = [1.5, 2.0, 2.5, 3.0]
        else:
            GWLs = [1.5, 2.0, 2.5, 3.0, 3.5]
        
        # Select base data for avoided impacts
        base_data = exposure[(exposure.Scenario == f'{ssp}_{str(baseline).replace(".", "p")}') & (exposure.Variable == var) & (exposure.Region == r)]
       
        # Loop through all GWLs
        for GWL in GWLs:            
            
            # Select data for GWL
            future_data = exposure[(exposure.Scenario == f'{ssp}_{str(GWL).replace(".", "p")}') & (exposure.Variable == var) & (exposure.Region == r)]            
            
            if future_data.empty:
                print(var, r, 'missing')
                continue
            
            else:
            
                # Save unit for variable
                unit = future_data.Unit.values[0]
                
                # If GWL > 1.5, then calculate avoided impacts compared to 1.5
                if GWL > 1.5:
        
                    avoided = (future_data.iloc[:, 5:]).sub(base_data.iloc[:, 5:].squeeze())
                    data_list.append(['Climate Solutions', f'{ssp}_{str(baseline).replace(".", "p")}', r, f'{var}|Avoided impacts|{str(GWL).replace(".", "p")}', unit] + avoided.values.tolist()[0])
                               
                    previous = exposure[(exposure.Scenario == f'{ssp}_{str(GWLs[GWLs.index(GWL)-1]).replace(".", "p")}') & (exposure.Variable == var) & (exposure.Region == r)]
                    diff = (future_data.iloc[:, 5:]).sub(previous.iloc[:, 5:].squeeze())
                    data_list.append(['Climate Solutions', f'{ssp}_{str(GWL).replace(".", "p")}', r, f'{var}|Difference|{str(GWLs[GWLs.index(GWL)-1]).replace(".", "p")}', unit] + diff.values.tolist()[0])
                
                # Else calculate difference between 1.5 and 1.2
                else: 
                    
                    data_1p2 = exposure[(exposure.Scenario == f'{ssp}_1p2') & (exposure.Variable == var) & (exposure.Region == r)]
                    diff = (future_data.iloc[:, 5:]).sub(data_1p2.iloc[:, 5:].squeeze())
                    data_list.append(['Climate Solutions', f'{ssp}_{str(GWL).replace(".", "p")}', r, f'{var}|Difference|1p2', unit] + diff.values.tolist()[0])
 
    # Set column names and save
    df = pd.DataFrame(data_list)
    df.columns=[['Model','Scenario', 'Region', 'Variable', 'Unit'] + years.tolist()]
    df.to_csv(os.path.join(output_dir, f'table_output_avoided_impacts_{ind}_{mode}.csv'), mode='w', index=False)

    # Remove negative values and save as separate file
    num = df._get_numeric_data()
    num[num < 0] = 0
    df.to_csv(os.path.join(output_dir, f'table_output_avoided_impacts_{ind}_{mode}_NoNeg.csv'), mode='w', index=False)    
    