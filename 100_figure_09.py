# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:55:42 2024

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


#%% Settings
#------------------------------------------------------------------------------

# Set input/output directories and paths
input_dir = 'table_output'
output_dir = ''

# Set parameters
indicators = ['pr_r99p', 'wsi', 'hw_95_7']
gmt = ['1p5', '2p0', '3p0']
year = '2050'
var = 'Exposure\|Population\|%\|50.0th Percentile'
var_p5 = 'Exposure\|Population\|%\|5.0th Percentile'
var_p95 = 'Exposure\|Population\|%\|95.0th Percentile'
ssp = 'ssp2'
r10_labels = ['Latin\nAmerica', 'Southern\nAsia', 'Africa', 'Eastern\nAsia', 'Middle\nEast', 'EU', 'Europe', 'North\nAmerica', 'South-East\nAsia', 'Asia-\nPacific\nDeveloped', 'Eurasia']
labels = ['Very wet days', 'Water stress index', 'Heatwave events - 95th percentile - 7 days']
colours = sns.color_palette(['#FCDE9C', '#F0746E', '#DC3977', '#7C1D6F'])

# Load files
plt.style.use('path_to_folder_with_repo\\required_files\\style.mplstyle')

# Prepare variables  
columns = 1
rows = int(len(indicators)/columns)

#%% Region exposure plot
#------------------------------------------------------------------------------

# Load and prepare data
file_list = glob.glob(os.path.join(input_dir, '*_R10.csv'))
file_list = [j for j in [i for i in file_list if 'lc' not in i] if 'pm2p5' not in j]
data = pd.concat((pd.read_csv(f) for f in file_list), ignore_index=True)
data[['Scenario', 'GMT']] = data['Scenario'].str.split('_', expand=True)
df = data.melt(id_vars=['Model', 'Scenario', 'Region', 'Variable', 'Unit', 'GMT'], var_name="Year", value_name="Value", value_vars=["2020", "2030", "2040", "2050", '2060', '2070', '2080', '2090', '2100'])
df = df[df.Region != 'World']

# Create plot
sns.set_theme()
fig, axes = plt.subplots(rows, columns, figsize=(15, 10), sharey=True)

ind_count = 0

for i in range(0, rows):
    
    ind = indicators[ind_count]  
    str_match = f'{ind}\|{var}'
    str_match_p5 = f'{ind}\|{var_p5}'
    str_match_p95 = f'{ind}\|{var_p95}'
    
    # Select datat for plot
    plot_data = df[(df.Scenario == ssp) &
                      (df.Year == year) &
                        (df.Variable.str.contains('^' + str_match + '$', regex=True)) & 
                      (df.GMT.isin(gmt))]  
    
    # Calculate error bars
    errors = [list(df[(df.Scenario == ssp) &
                (df.Year == year) &
                (df.Variable.str.contains('^' + str_match + '$', regex=True)) & 
                (df.GMT.isin(gmt))].sort_values(['GMT', 'Region'])['Value'].reset_index(drop=True) - 
                
                df[(df.Scenario == ssp) &
                (df.Year == year) &
                (df.Variable.str.contains('^' + str_match_p5 + '$', regex=True)) & 
                (df.GMT.isin(gmt))].sort_values(['GMT', 'Region'])['Value'].reset_index(drop=True)),
                
                list(df[(df.Scenario == ssp) &
                (df.Year == year) &
                (df.Variable.str.contains('^' + str_match_p95 + '$', regex=True)) & 
                (df.GMT.isin(gmt))].sort_values(['GMT', 'Region'])['Value'].reset_index(drop=True) - 
                
                df[(df.Scenario == ssp) &
                (df.Year == year) &
                (df.Variable.str.contains('^' + str_match + '$', regex=True)) & 
                (df.GMT.isin(gmt))].sort_values(['GMT', 'Region'])['Value'].reset_index(drop=True))]    
    errors = [[j if j >= 0 else 0 for j in i] for i in errors]
    
    # Plot data    
    sns.barplot(ax=axes[i], data = df[(df.Scenario == ssp) &
                      (df.Year == year) &
                      (df.Variable.str.contains('^' + str_match + '$', regex=True)) & 
                      (df.GMT.isin(gmt))].sort_values(['GMT', 'Region']), 
                  x='Region', y='Value', hue='GMT', palette=colours,
                  capsize=.1, errcolor='black', errwidth=1)
    axes[i].get_legend().remove()    
    
    # Plot error bars
    x_coords = [p.get_x() + 0.5 * p.get_width() for p in axes[i].patches]
    y_coords = [p.get_height() for p in axes[i].patches]  
    axes[i].errorbar(x=x_coords, y=y_coords, yerr = errors, fmt='none', capsize=4, c='black')
        
    # Set labels
    axes[i].set(title=labels[i], xlabel=None)
    axes[i].set_ylabel('Exposed population [%]', labelpad=10)
        
    if i == rows-1:
        axes[i].set_xticks(range(len(df.Region.unique())), labels=r10_labels)
    else:
        axes[i].set_xticks(range(len(df.Region.unique())), labels=[])
    
    ind_count += 1  

handles, plot_labels = axes[i].get_legend_handles_labels()
fig.legend(handles, [f'{g.replace("p", ".")}  °C' for g in gmt], loc='lower center',ncol=3)
plt.savefig(os.path.join(output_dir, f'R10_regions_median.eps'), format='eps', 
            bbox_inches='tight', dpi=300)
