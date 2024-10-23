import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import glob
import os
import numpy as np
import matplotlib.colors
import math
from cmcrameri import cm
import re

#%% Define colours etc.
#------------------------------------------------------------------------------

cmap_blues = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#e8e8e8","#ace4e4","#5ac8c8"])
cmap_pinks = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#e8e8e8","#dfb0d6","#be64ac"])

score_args = dict(extend="neither", cbar_kwargs=dict(orientation="horizontal",
                  shrink=0.5, pad=0.05, aspect=10, ticks=[0, 1, 2, 3, 4, 5, 6],
                  label="Score"))

abs_cbar_kwargs = dict(orientation="horizontal", shrink=0.5, pad=0.05,
                       aspect=10, label="Score")

diff_args = dict(extend="neither", cbar_kwargs=dict(orientation="horizontal",
                 shrink=0.5, pad=0.05, aspect=10, label="Difference (%)"))

diff_cbar_kwargs = dict(orientation="horizontal", shrink=0.5, pad=0.05,
                        spect=10, ticks=[-20, -10, 0, 10, 20], label="Difference (%)")

abs_args = dict(extend="neither", cbar_kwargs=dict(orientation="horizontal",
                shrink=0.5, pad=0.05, aspect=10, label="days/year"))

abs_cbar_kwargs = dict(orientation="horizontal", shrink=0.5, pad=0.05,
                       aspect=10, label="days year⁻¹")

JK1_cmap = ["#e8e8e8", "#ace4e4", "#5ac8c8", "#dfb0d6", "#a5add3", "#5698b9",
            "#be64ac", "#8c62aa", "#3b4994",]
JK1 = matplotlib.colors.ListedColormap(JK1_cmap)

#%% Functions
#------------------------------------------------------------------------------

def make_map(fp, ax=None, title="Add a better title", ylabel="", color_lims=[0, 6],
             savefolder='',
             add_coastlines=True, add_borders=False, add_gridlines=False,
             cmap="magma_r", cbar=True, crs=ccrs.Robinson(), plotargs={},
             gwl = ''):
    
    """
    Makes a map from a given file path for a single datarray.
    """
    
    data = xr.open_dataarray(fp)
    
    if gwl:
        data = data.sel({'threshold': gwl})    
    if ax is None:
        fig, ax = plt.subplots(subplot_kw=dict(projection=crs))
    elif not isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot):
        msg = "Must provide a cartopy axes object, not: {}"
        raise ValueError(msg.format(type(ax)))

    if add_coastlines:
        ax.add_feature(cartopy.feature.OCEAN, facecolor="#FFFFFF", zorder=0)
        coasts = cartopy.feature.NaturalEarthFeature(
            category="physical",
            scale="110m",
            name="coastline",
            facecolor="none",
            zorder=2)
        ax.add_feature(coasts, edgecolor="dimgrey", linewidth=0.2, zorder=2)
    if add_borders:
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="#778899")
    if add_gridlines:
        ax.gridlines(zorder=3, color="#f5f5f5")

    if data is not None:
        p = data.plot(
            transform=ccrs.PlateCarree(),
            ax=ax,
            vmin=color_lims[0],
            vmax=color_lims[1],
            cmap=cmap,
            zorder=1,
            add_colorbar=cbar,
            rasterized=True,
            **plotargs)
    ax.set_title(title, pad=15)
    return p

#------------------------------------------------------------------------------

def make_map_from_da(da, ax=None, title="Add a better title", ylabel="",
                     color_lims=[0, 6], savefolder='',
                     add_coastlines=True, add_borders=False,
                     add_gridlines=False, cmap="magma_r", cbar=True,
                     crs=ccrs.Robinson(), plotargs={}, gwl=''):
    
    """
    Makes a map from an Xarray dataarray without multiple dimensions (exclusing lon and lat).
    Used for making maps of historical data where the stats have already been selected.
    """
    
    if gwl:
        da = da.sel({'threshold': gwl})

    if add_coastlines:
        ax.add_feature(cartopy.feature.OCEAN, facecolor="#FFFFFF", zorder=0)
        coasts = cartopy.feature.NaturalEarthFeature(
            category="physical",
            scale="110m",
            name="coastline",
            facecolor="none",
            zorder=2)
        ax.add_feature(coasts, edgecolor="dimgrey", linewidth=0.2, zorder=2)
    if add_borders:
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="#778899")
    if add_gridlines:
        ax.gridlines(zorder=3, color="#f5f5f5")

    if da is not None:
        p = da.plot(
            transform=ccrs.PlateCarree(),
            ax=ax,
            vmin=color_lims[0],
            vmax=color_lims[1],
            cmap=cmap,
            zorder=1,
            add_colorbar=cbar,
            **plotargs)
    ax.set_title(title)
    return p

#------------------------------------------------------------------------------

def make_grid(nrows=6):
    
    """
    Make a grid of subplots with a Robinson projection for them multimap figure.
    """
    
    fig, ax = plt.subplots(
        ncols=3,
        nrows=nrows,
        figsize=(9,10),
        layout="constrained",
        gridspec_kw={"left": 0.15, "right": 0.95, "top": 1,
                     "bottom": 0.15, "wspace": 0.05, "hspace": -0.3},
        subplot_kw={"projection": ccrs.Robinson()},
    )
    return fig, ax

#------------------------------------------------------------------------------

def get_data(ind, var, stat='', folder=''):
    
    """
    Make an array with the abs, diff and score data for one indicator. Shapes the filenames into an array in the correct shape for the figure axes.
    """
    
    all_files = glob.glob(os.path.join(folder, ind, stat, f'*{var}*ssp2*.nc4'))
    nthreshs = len(all_files) // 3
    files = np.array(sorted(all_files)).reshape(nthreshs, 3).tolist()
    return files

#------------------------------------------------------------------------------

def make_multi_maps(files, nrows=6, title="", cmaps=["YlOrRd", "Reds", "magma_r"],
                    clims=[[0, 350], [0, 5170], [0, 6]], args=[abs_args, diff_args, score_args],
                    cticks=[[0, 30, 60, 90, 120, 150, 180],
                            [-20, -10, 0, 10, 20],
                            [0, 1, 2, 3, 4, 5, 6]],
                    clabels=["days per year", "Difference (%)", "Score"]):
    
    """
    Make a map showing the abs diff and score of a given indicator at different thresholds.
    """
    
    degrees = ["1.2 °C", "1.5 °C", "2.0 °C", "2.5 °C", "3.0 °C", "3.5 °C"]
    titles = [['Absolute values', 'Relative difference', 'Hazard score'], 
              ['', '', ''], ['', '', ''], ['', '', ''], 
              ['', '', ''], ['', '', ''], ['', '', '']]
    
    fig, axs = make_grid(nrows=nrows)
    for i, degree in zip(range(nrows), degrees):
        for j, cmap, clim, ticks, label in zip(range(3), cmaps, clims, cticks, clabels):
            arg = dict(
                extend="neither",
                cbar_kwargs=dict(
                    orientation="horizontal",
                    shrink=0.5,
                    pad=0.05,
                    aspect=10,
                    ticks=ticks,
                    label=label))

            make_map(
                files[i][j],
                ax=axs[i, j],
                color_lims=clim,
                title=titles[i][j],
                ylabel="Hello",
                cbar=True if i == (nrows - 1) else False,
                cmap=cmap,
                plotargs=arg if i == (nrows - 1) else dict(),
            )
            
            if j == 0:
                axs[i, j].text(-0.05, 0.5, f"{degree}",
                               transform=axs[i, j].transAxes,
                               ha="right", va="center")
                
        fig.suptitle(title, x=0.55, y=1.05, fontsize=16)
    return fig

#------------------------------------------------------------------------------

def get_score_data(ind, gwl, folder):
    
    """
    Get the score data for a given indicator, ssp, threshold and metric.
    """
    
    filename = f"{ind}_ssp2_{str(gwl).replace('.', 'p')}_score.nc4"
    fp = folder + filename
    files = glob.glob(fp, recursive=True)
    clean_files = []
    for file in files:
        if ("q5" not in file) and ("q95" not in file):
            clean_files.append(file)
    return clean_files

#------------------------------------------------------------------------------

def all_scores_grid(no_of_ind):
    
    """
    Make a grid of subplots with a Robinson projection for the all scores figure.
    """
    
    nrows = math.ceil(no_of_ind / 3)
    
    fig, ax = plt.subplots(
        ncols=3, nrows=nrows, figsize=(13, 10), 
        subplot_kw={"projection": ccrs.Robinson()},
        layout = 'constrained'
    )
    
    rm_ax = 3 - (no_of_ind % 3)
    
    for i in range(0, rm_ax):
        fig.delaxes(ax[nrows-1, 3-(i+1)])
   
    return fig, ax

#------------------------------------------------------------------------------

def make_all_scores(inds, folder, stat, no_of_ind, gwl):
    
    """
    Make a map of all scores for a given set of indicators at a given threshold.
    """
    
    fig, ax = all_scores_grid(no_of_ind)
    count = 0
    
    for ind, variables in inds.items():
        for var in variables:

            file = glob.glob(os.path.join(folder, ind, stat, f"*_{var}_ssp2_{str(gwl).replace('.', 'p')}_score.nc4"))
            data = xr.open_dataset(file[0])       
            
            long_name = data[list(data)[0]].attrs['long_name']
            
            if re.findall("\d{1}[C]", long_name):
                degree = re.findall("\d{1}[C]", long_name)
                plot_title = long_name.replace(degree[0], degree[0][0] + ' °' + degree[0][1])
            elif ind == 'heatwave':
                plot_title = f'Heatwave events - {var.split("_")[1]}th percentile - {var.split("_")[2]} days'
            else:
                plot_title = data[list(data)[0]].attrs['long_name']
                
            p = make_map_from_da(
                    data[list(data)[0]],
                    ax=ax.flatten()[count],
                    title=plot_title,
                    cmap="magma_r",
                    cbar=True if count == (no_of_ind) else False,
                    color_lims=[0, 6],
                    plotargs=score_args if count == (no_of_ind) else dict())
            
            count +=1
            
    fig.colorbar(p, ax=ax.ravel().tolist(), location='bottom',
                 label='', shrink=0.35)
    fig.suptitle(f'Hazard scores at {gwl} °C', y=1.05, fontsize=16)    
    return fig

#------------------------------------------------------------------------------

def bivar_bin(ds_future, ds_std):
    
    """
    Make an Xarray dataset with values 0-9 for bivariate scores.
    """
    
    ds_binned = xr.where((ds_future <= 1) & (ds_std <= 1), 1, np.nan)  # hist 0-1, temp 0-1 = 1
    ds_binned = xr.where((ds_future <= 1) & ((ds_std <= 2) & (ds_std > 1)),
                         4, ds_binned)  # hist 0-1, temp 1-2 = 2
    ds_binned = xr.where((ds_future <= 1) & ((ds_std <= 3) & (ds_std > 2)), 
                         7, ds_binned)  # hist 0-1, temp 2-3 = 3
    ds_binned = xr.where(((ds_future > 1) & (ds_future <= 2)) & (ds_std <= 1), 
                         2, ds_binned)  # hist 1-2, temp 0-1 = 4
    ds_binned = xr.where(((ds_future > 1) & (ds_future <= 2)) & ((ds_std > 1) & (ds_std <= 2)), 
                         5, ds_binned)  # hist 1-2, temp 1-2 = 5
    ds_binned = xr.where(((ds_future > 1) & (ds_future <= 2)) & ((ds_std > 2) & (ds_std <= 3)), 
                         8, ds_binned)  # hist 1-2, temp 2-3 = 6
    ds_binned = xr.where(((ds_future > 2) & (ds_future <= 3)) & (ds_std <= 1), 
                         3, ds_binned)  # hist 2-3, temp 0-1 = 7
    ds_binned = xr.where(((ds_future > 2) & (ds_future <= 3)) & (ds_std > 1) & (ds_std <= 2), 
                         6, ds_binned)  # hist 2-3, temp 1-2 = 8
    ds_binned = xr.where(((ds_future > 2) & (ds_future <= 3)) & ((ds_std > 2) & (ds_std <= 3)), 
                         9, ds_binned)  # hist 2-3, temp 2-3 = 9
    return ds_binned

#------------------------------------------------------------------------------

def make_bivar_map(ds_future, ds_std, ax=None, cbar=None, gwl=2, 
                   add_coastlines=True, title="new title", 
                   cmap=JK1, crs=ccrs.Robinson()):
    
    """
    Make a bivariate map of historic and future scores from two datasets.
    """

    if ax is None:
        fig, ax = plt.subplots(subplot_kw=dict(projection=crs))
    elif not isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot):
        msg = "Must provide a cartopy axes object, not: {}"
        raise ValueError(msg.format(type(ax)))
        
    if gwl:
        ds_future = ds_future.sel({'threshold': gwl})
        ds_std = ds_std.sel({'threshold': gwl})

    ds_binned = bivar_bin(ds_future, ds_std)
    p = ds_binned.plot(
        x="lon",
        y="lat",
        cmap=cmap,
        clim=(1, 9),
        vmin=1, vmax=9,
        transform=ccrs.PlateCarree(),
        add_colorbar=False,
        ax=ax)
    if add_coastlines:
        ax.add_feature(cartopy.feature.OCEAN, facecolor="#FFFFFF", zorder=0)
        coasts = cartopy.feature.NaturalEarthFeature(
            category="physical",
            scale="110m",
            name="coastline",
            facecolor="none",
            zorder=2)
        ax.add_feature(coasts, edgecolor="dimgrey", linewidth=0.2, zorder=2)
    ax.set_title(title)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    legend_numbers = np.array([[6.5, 7.5, 8.5], [3.5, 4.5, 5.5], [0.5, 1.5, 2.5]])
    cbar.imshow(legend_numbers, cmap=cmap, clim=(0, 9))
    cbar.set_ylabel("Relative score")
    cbar.set_xlabel("Absolute score")
    cbar.set_xticks([])
    cbar.set_yticks([-0.5, 0.5, 1.5, 2.5])
    cbar.set_yticklabels([3, 2, 1, 0])
    cbar.set_xticks([-0.5, 0.5, 1.5, 2.5])
    cbar.set_xticklabels([0, 1, 2, 3])
    return p, cbar

#------------------------------------------------------------------------------

def create_fig_for_score_plot(fig_size):   
    
    """
    Create figure for score methodology plot.
    """

    fig, axs = plt.subplot_mosaic([
        ['abs_hist', 'abs_hist', 'abs_gwl', 'abs_gwl', 'fut_bin', 'fut_bin'],
        ['diff', 'diff', 'std', 'std', 'std_bin','std_bin'],
        ['score', 'score', 'score', 'bivar_cmap', 'bivar','bivar']],
        figsize=(fig_size), 
        per_subplot_kw={('abs_hist', 'abs_gwl', 'fut_bin', 'diff',
                          'std', 'std_bin', 'bivar', 'bivar_cmap', 'score'): 
                        {'projection': ccrs.Robinson()}})            
    
    return fig, axs

#------------------------------------------------------------------------------
    
def create_score_methodology_plot(hist, absolute, future, diff, std, std_change,
                                  scores, params, ind, var, gwl, axs):
    
    """
    Make a plot explaining the score methodology.
    """
    
    unit = '' if params['indicators'][ind][var]['unit'] == '-' else params['indicators'][ind][var]['unit_SI']
    
    make_map_from_da(hist, ax=axs['abs_hist'], 
                         title='Reference period',
                         cmap=params['indicators'][ind][var]['ind_cmap'],
                         cbar=True,
                         color_lims=[params['indicators'][ind][var]['ind_min'], 
                                     params['indicators'][ind][var]['ind_max']],
                         plotargs= dict(
                           extend="neither",
                           cbar_kwargs=dict(
                             orientation="horizontal",
                             shrink=0.5,
                             pad=0.05,
                             aspect=10,
                             ticks=[params['indicators'][ind][var]['ind_min'], 
                                    params['indicators'][ind][var]['ind_max']],
                             label = unit)))

    make_map_from_da(absolute, ax=axs['abs_gwl'],
                       title=f'2.0 °C',
                       cmap=params['indicators'][ind][var]['ind_cmap'],
                       color_lims=[params['indicators'][ind][var]['ind_min'], 
                                   params['indicators'][ind][var]['ind_max']],
                       plotargs= dict(
                         extend="neither",
                         cbar_kwargs=dict(
                            orientation="horizontal",
                            shrink=0.5,
                            pad=0.05,
                            aspect=10,
                            ticks=[params['indicators'][ind][var]['ind_min'], 
                                   params['indicators'][ind][var]['ind_max']],
                            label = unit)),
                           gwl=gwl)

    make_map_from_da(future, ax=axs['fut_bin'],
                       title='Absolute score component',
                       cmap=cmap_blues,
                       color_lims=[0, 3],
                       plotargs= dict(
                         extend="neither",
                         cbar_kwargs=dict(
                            orientation="horizontal",
                            shrink=0.5,
                            pad=0.05,
                            aspect=10,
                            ticks=[0, 3],
                            label="")),
                        gwl=gwl)

    make_map_from_da(diff, ax=axs['diff'],
                       title='Difference',
                       cmap=params['indicators'][ind][var]['ind_cmap'],
                       color_lims=[0, 
                                   10],
                       plotargs= dict(
                         extend="neither",
                         cbar_kwargs=dict(
                            orientation="horizontal",
                            shrink=0.5,
                            pad=0.05,
                            aspect=10,
                            ticks=[0,
                                   10],
                            label=unit)),
                         gwl=gwl)


    make_map_from_da(std, ax=axs['std'], 
                       title='Standard deviation',
                       cmap=cm.batlow_r,
                       color_lims=[math.floor(std.min()), 10],
                       cbar=True,
                       plotargs= dict(
                         extend="neither",
                         cbar_kwargs=dict(
                            orientation="horizontal",
                            shrink=0.5,
                            pad=0.05,
                            aspect=10,
                            ticks=[math.floor(std.min()), 10],
                            label="")))

    make_map_from_da(std_change, ax=axs['std_bin'],
                 title='Relative score component',
                 cmap=cmap_pinks,
                 color_lims=[0, 3],
                 plotargs= dict(
                 extend="neither",
                 cbar_kwargs=dict(
                     orientation="horizontal",
                     shrink=0.5,
                     pad=0.05,
                     aspect=10,
                     ticks=[0, 3],
                     label="")),
                  gwl=gwl)

    make_bivar_map(future, std_change, 
                          ax=axs['bivar'], cbar=axs['bivar_cmap'],
                          gwl=gwl, 
                          title='Bivariate hazard score')

    make_map_from_da(scores, ax=axs['score'],
               title='Hazard score',
               cmap='magma_r',
               color_lims=[0, 6],
               plotargs= dict(
                   extend="neither",
                   cbar_kwargs=dict(
                       orientation="horizontal",
                       shrink=0.5,
                       pad=0.05,
                       aspect=10,
                       ticks=[0, 1, 2, 3, 4, 5, 6],
                       label="")),
               gwl=gwl)

#------------------------------------------------------------------------------
    
def make_discrete_map_from_da(da, levels, ax=None, title="Add a better title",
                              ylabel="", color_lims=[0, 6], add_coastlines=True,
                              add_borders=False, add_gridlines=False, cmap="magma_r",
                              cbar=True, crs=ccrs.Robinson(), plotargs={},
                              gwl=''):
    
    """
    Makes a map from an Xarray dataarray without multiple dimensions (exclusing lon and lat).
    Used for making maps of historical data where the stats have already been selected.
    """
    
    if gwl:
        da = da.sel({'threshold': gwl})

    if add_coastlines:
        ax.add_feature(cartopy.feature.OCEAN, facecolor="#FFFFFF", zorder=0)
        coasts = cartopy.feature.NaturalEarthFeature(
            category="physical",
            scale="110m",
            name="coastline",
            facecolor="none",
            zorder=2)
        ax.add_feature(coasts, edgecolor="dimgrey", linewidth=0.2, zorder=2)
    if add_borders:
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="#778899")
    if add_gridlines:
        ax.gridlines(zorder=3, color="#f5f5f5")

    if da is not None:
        p = da.plot(
            transform=ccrs.PlateCarree(),
            ax=ax,
            vmin=color_lims[0],
            vmax=color_lims[1],
            cmap=cmap,
            levels=levels,
            zorder=1,
            add_colorbar=cbar,
            **plotargs)
    ax.set_title(title)
    return p

#------------------------------------------------------------------------------

def make_plot_dens_map(da, ax=None, title="Add a better title", ylabel="",
                       add_coastlines=True, add_borders=False,
                       add_gridlines=False, cmap=cm.lipari_r, cbar=True,
                       crs=ccrs.Robinson(), plotargs={}):
    
    """
    Makes population density plot.
    """
   
    if add_coastlines:
        ax.add_feature(cartopy.feature.OCEAN, facecolor="#FFFFFF", zorder=0)
        coasts = cartopy.feature.NaturalEarthFeature(
            category="physical",
            scale="110m",
            name="coastline",
            facecolor="none",
            zorder=2)
        ax.add_feature(coasts, edgecolor="dimgrey", linewidth=0.2, zorder=2)
    if add_borders:
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="#778899")
    if add_gridlines:
        ax.gridlines(zorder=3, color="#f5f5f5")

    if da is not None:
        p = da.plot(
            transform=ccrs.PlateCarree(),
            ax=ax,
            cmap=cmap,
            zorder=1,
            add_colorbar=cbar,
            **plotargs)
    ax.set_title(title)
    return p