#! /usr/local/bin/python
# To calculate water scarcity index from ISIMIP2b data
# By Yusuke Satoh (satoh@iiasa.ac.at) (adapted by Michaela Werning)
# On 23 Nov 2021
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
from netCDF4 import Dataset

# ---(Edit these for your process)--------------------------------------------------------------------------------------
CHECKFIG = True

scenarios = [
    ('rcp26', '2005soc'),
    ('rcp60', '2005soc')]

temporal_target = 'Yearly'
gcms = ['GFDL-ESM2M', 'HADGem2-ES', 'IPSL-CM5A-LR', 'MIROC5']
ghms = ['H08', 'LPJmL', 'MATSIRO']  
syear, eyear = 1971, 2099

"""
Caution!! 
syear has to be the first year of the first file. i.e., 1??1, 20?1, or 2006
eyear has to be the last year of the final file. i.e, 1??0, 20?0, 2005, or 2099
"""

if np.mod(syear, 10) != 1 and syear != 2006:
    raise ValueError('syear has to be 1??1, 20?1, or 2006.')
if np.mod(eyear, 10) != 0 and eyear != 2005 and eyear != 2099:
    raise ValueError('eyear has to be 1??0, 20?0, 2005, or 2099.')

# your directories...  (Edit here for your process)
isimip_directory_main = '\\ISIMIP\\ISIMIP2b\\output'
waterdemand_directory_main = '\\ISIMIP3b\\data\\water_abstraction'
output_directory_main = ''

# paths for some input data...
gridarea_path = 'path_to_folder_with_repo\\required_files\\grd_ara.hlf'
flow_direction_path = 'path_to_folder_with_repo\\required_files\\ddm30_flowdir_cru_neva.nc4'
landsea_mask_path = 'ISIMIP2b_landseamask_generic.nc4' # (https://data.isimip.org/datasets/871048b3-e771-42c0-ba14-28463452efff/)


# ---(Global paras basically you don't have to touch)-------------------------------------------------------------------
supply_type = 'Inflow'  # this is default
#supply_type = 'nonInflow'  # consider only local runoff

# parameters related to year
years = range(syear, eyear+1)
nyear = len(years)
year_chunks_of_ncfiles = []
if eyear <= 2005:  # only historical period
    syears = range(syear, eyear, 10)
elif syear < 2006 and 2010 <= eyear:  # both historical and future
    syears = list(range(syear, eyear, 10)) + [2006]
    syears.sort()
elif 2006 == syear and 2010 < eyear:  # only future starting 2006
    syears = [2006] + list(range(2011, eyear, 10))
elif 2006 < syear:  # only future
    syears = range(syear, eyear, 10)
for _syear in syears:
    if _syear == 2001: year_chunks_of_ncfiles.append((2001, 2005))
    elif _syear == 2006: year_chunks_of_ncfiles.append((2006, 2010))
    elif _syear == 2091: year_chunks_of_ncfiles.append((2091, 2099))
    else: year_chunks_of_ncfiles.append((_syear, _syear+9))
print(f'year_chunks_of_ncfiles: {year_chunks_of_ncfiles}')

base_year = 1901

# other...
nt,ny,nx = nyear, 360, 720

# ----------------------------------------------------------------------------------------------------------------------
def mk_inflow(src, flowdirection):  # src (nyear, 12, 280, 720)

    mask = src.mask
    YY, XX = np.where(mask[0,0]==False)  # pick up land grid cells
    src = src.filled(0)
    inflow = np.zeros(src.shape, 'float32')
    
    for iy, ix in zip(YY, XX):
        if   flowdirection[iy,ix] == 1 and ix != 719: inflow[:, :, iy,   ix+1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 1 and ix == 719: inflow[:, :, iy,      0] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 2 and ix != 719: inflow[:, :, iy+1, ix+1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 2 and ix == 719: inflow[:, :, iy+1,    0] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 3              : inflow[:, :, iy+1, ix  ] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 4 and ix !=   0: inflow[:, :, iy+1, ix-1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 4 and ix ==   0: inflow[:, :, iy+1,  719] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 5 and ix !=   0: inflow[:, :, iy  , ix-1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 5 and ix ==   0: inflow[:, :, iy  ,  719] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 6 and ix !=   0: inflow[:, :, iy-1, ix-1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 6 and ix ==   0: inflow[:, :, iy-1,  719] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 7              : inflow[:, :, iy-1, ix  ] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 8 and ix != 719: inflow[:, :, iy-1, ix+1] += src[:, :, iy, ix]
        elif flowdirection[iy,ix] == 8 and ix == 719: inflow[:, :, iy-1,    0] += src[:, :, iy, ix]

    return np.ma.masked_array(inflow, mask=mask)  # (nyear, 12, 280, 720) <numpy.ma.core.MaskedArray>

# ----------------------------------------------------------------------------------------------------------------------
def read_netcdfs(variable, paths, dim, decode_times=True):

    def load_dataset(_path):
        if not os.path.isfile(_path):
            raise FileNotFoundError(_path)
        else:
            print(f'loading... {_path}')
            dataset = xr.open_dataset(_path, decode_times=decode_times)
            if '2099' in _path: print(f'{dataset[variable]}\n')
            return dataset

    if variable in ['dis', 'qtot']:  # daily input --> (nyear*12, 360, 720) <numpy.ndarray>
        datasets = [load_dataset(path).resample(time='M').mean() for path in paths]  # day2month
    elif variable in ['pirrww',           # monthly input  --> (_nyear*12, 360, 720)
                      'domww', 'indww']:  # yearly input   --> (_nyear, 280, 720)
        datasets = [load_dataset(path) for path in paths]  # decode_times=False
    src = xr.concat(datasets, dim)[variable].fillna(1e+20).data 
    return np.ma.masked_equal(src, 1e+20)


def load_modeloutput(gcm, ghm, scenario, soc, variable):

    def gen_path_daily(_syear, _eyear):
        _ghm = ghm.lower().replace('_', '')
        if _eyear < 2006: _period, _scenario, _soc = 'historical', 'historical', 'histsoc'
        # else:             _period, _scenario, _soc = 'future', scenario, soc
        else:             _period, _scenario, _soc = scenario, scenario, soc
        if _ghm == 'clm45' and _period == 'historical': _soc = '2005soc'
        file_name = f'{_ghm}_{gcm}_ewembi_{_scenario}_{_soc}_co2_{variable}_global_daily_{_syear}_{_eyear}.nc4'
        return os.path.join(isimip_directory_main, ghm, _period, gcm, file_name)

    def gen_path_monthly(_scenario, _soc, _syear, _eyear):
        _ghm = ghm.lower().replace('_', '')
        _period = 'historical' if _scenario == 'historical' else scenario
        file_name = f'{_ghm}_{gcm}_ewembi_{_scenario}_{_soc}_co2_{variable}_global_monthly_{_syear}_{_eyear}.nc4'
        return os.path.join(isimip_directory_main, ghm, _period, gcm, file_name)

    if not variable in ['dis', 'qtot', 'pirrww']:
        raise ValueError(f'this function is for dis, qtot, and pirrww, but varialbe = {variable}.')
    # load src
    
    if variable in ['dis', 'qtot']:  # daily2monthly
        srcpaths = [gen_path_daily(_syear, _eyear) for _syear, _eyear in year_chunks_of_ncfiles]
        src = read_netcdfs(variable, srcpaths, 'time')             # (nyear*12,  360, 720)
        src = src.reshape(-1, 12, 360, 720)[:, :, 12:12 + 280, :]  # (nyear, 12, 280, 720)
        src = np.ma.masked_less(src, 0)  # This is an exceptional process for minus runoff in WaterGAP data...
    elif variable == 'pirrww':  # monthly
        if eyear <= 2005:  # only historical
            _years = range(1861, 2005 + 1)
            srcpaths = [gen_path_monthly('historical', 'histsoc', 1861, 2005)]
        elif syear < 2006 and 2005 < eyear:  # historical+future
            _years = range(1861, 2099 + 1)
            srcpaths = [gen_path_monthly('historical', 'histsoc', 1861, 2005), gen_path_monthly(scenario, soc, 2006, 2099)]
        elif 2006 <= syear:  # only future
            _years = range(2006, 2099 + 1)
            srcpaths = [gen_path_monthly(scenario, soc, 2006, 2099)]
        src = read_netcdfs(variable, srcpaths, 'time', False)                                            # (_nyear*12, 360, 720)
        src = src.reshape(-1, 12, 360, 720)[_years.index(syear):_years.index(eyear)+1, :, 12:12+280, :]  # (nyear, 12, 280, 720)
    # post process. Unit conversion.
    if variable == 'dis':
        src = src * 3600 * 24  # (nyear, 12, 280, 720)  UnitConv: [m3/s] -> [m3/dy]
        flowdirection = xr.open_dataset(flow_direction_path)['flowdirection'].data[::-1]  # (280,720) <numpy.ndarray>
        src = mk_inflow(src, flowdirection)  # (nyear, 12, 280, 720) <numpy.ma.core.MaskedArray>
    elif variable in ['qtot', 'pirrww']:
        area = np.fromfile(gridarea_path, 'float32').byteswap().reshape(360, 720)[12:12+280, :]  # (280, 720) [m2]
        src = src * area * 86400 * 1e-3  # (nyear, 12, 280, 720)  Unit Conv: [kg/m2/s] -> [kg/s] -> [kg/dy] -> [m3/dy]  <numpy.ma.core.MaskedArray>
    return src.filled(1e+20)  # (nyear, 12, 280, 720) <numpy.ndarray>


def load_demandinput(soc, variable):

    def gen_path(_soc, _syear, _eyear):
        # file_name = f'{variable}_{_soc}_annual_{_syear}-{_eyear}.nc'
        file_name = f'{variable}_{_soc}_annual_{_syear}_{_eyear}.nc'
        # return os.path.join(waterdemand_directory_main, _soc, file_name)
        return os.path.join(waterdemand_directory_main, file_name)
        
    if eyear <= 2005:  # only historical
        _years = range(1901, 2014+1)
        srcpaths = [gen_path('histsoc', 1901, 2014)]
    elif syear < 2015 and 2014 < eyear:  # historical+future
        _years = range(1901, 2100+1)
        srcpaths = [gen_path('histsoc', 1901, 2014), gen_path('2015soc', 2015, 2100)]
    elif 2015 <= syear:  # only future
        _years = range(2014, 2100+1)
        # srcpaths = [gen_path(soc, 2015, 2099)]
        srcpaths = [gen_path('2015soc', 2014, 2100)]        
        
    src = read_netcdfs(variable, srcpaths, 'time', False)  # (_nyear, 280, 720)   decode_times=False
    src = src[_years.index(syear):_years.index(eyear)+1]   # (nyear, 280, 720)
    # post process. year2month & unit conv  [m3/yr] -> [m3/dy] 
    src = np.array([[src[iyear]/((datetime(year,12,31)-datetime(year,1,1)).days+1) for imon in range(12)] 
                                                                                   for iyear, year in enumerate(years)])  # (nyear, 12, 280, 720)
    return np.ma.masked_equal(src, 1e+20)

# ----------------------------------------------------------------------------------------------------------------------
def draw_a_map(variables, srcs, year, gcm, rcp, soc, output_directory):

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 4, figure=fig)
    gs.update(left=0.05, right=0.95, bottom=0.04, top=0.9, wspace=0.06, hspace=0.04)

    positions = [(0,0), (0,1), (0,2), (0,3), (1,0), (1,1), (1,2), (1,3)]
    for (irow, icol), variable, src in zip(positions, variables, srcs):

        ax = fig.add_subplot(gs[irow,icol])
        ax.axis('off')
        ax.set_title(variable)
        im = ax.imshow(src)
        pos = ax.get_position()
        cax = fig.add_axes([pos.x0, pos.y0, pos.x1-pos.x0, 0.01])
        plt.colorbar(im, cax, orientation='horizontal', pad=0.01)

    figure_path = os.path.join(output_directory, f'wsi_{gcm}_{rcp}_{soc}_{year}.png')
    plt.savefig(figure_path)
    print(f'savefig: {figure_path}')
    plt.close()
    del fig

# ----------------------------------------------------------------------------------------------------------------------
def write_nc(variable_name, src, ghm, gcm, rcp, soc, output_directory):

    if variable_name == 'Stress':
        long_name = 'Water Stress Index'
        acronym = 'wsi'
        unit = '-'
    else:
        raise ValueError(f'This script is currently for WSI, but variable_name is {variable_name}')

    mask = Dataset(landsea_mask_path)['LSM'][0].mask  # (360, 720)

    outSrc = np.ma.masked_equal(np.zeros((nt, ny, nx), 'float32'), 0).filled(1e+20)
    if   src.shape[1] == 280: outSrc[:,12:12+280,:] = src.filled(1e+20)
    elif src.shape[1] == 360: outSrc = src.filled(1e+20)
    else: raise ValueError(f'something is wrong with shape of src... {src.shape}')
    outSrc = np.ma.masked_equal(outSrc, 1e+20)  # (nyear, 360, 720)
    outSrc = np.ma.masked_array(outSrc, mask=np.resize(mask, outSrc.shape))

    # --- create a file
    filename = f'{acronym}_{ghm}_{gcm}_{rcp}_{soc}_{syear}_{eyear}.nc4'
    output_path = os.path.join(output_directory, filename)
    rootgrp = Dataset(output_path, 'w', format='NETCDF4')
    lon  = rootgrp.createDimension('lon', nx)
    lat  = rootgrp.createDimension('lat', ny)
    time = rootgrp.createDimension('time', nt)
    longitudes = rootgrp.createVariable('lon', 'f8', ('lon'))
    latitudes  = rootgrp.createVariable('lat', 'f8', ('lat'))
    times      = rootgrp.createVariable('time','f8', ('time'))
    srcs       = rootgrp.createVariable(acronym, 'f4', ('time','lat','lon'), fill_value=1e+20)
    # --- attribute
    import time
    rootgrp.title            = f'{variable_name} {ghm} {rcp} {soc} {gcm}'
    rootgrp.description      = f'{variable_name} {syear}-{eyear} under {rcp}_{soc} projected by {gcm} with {", ".join(ghms)}'
    rootgrp.history          = 'Created ' + time.ctime(time.time())
    rootgrp.source           = 'ISI-MIP2b'
    rootgrp.contact          = 'xxxxxxx@iiasa.ac.at'
    rootgrp.institution      = 'IIASA'
    longitudes.long_name     = 'longitude'
    longitudes.units         = 'degrees east'
    longitudes.standard_name = 'longitude'
    longitudes.axis          = 'X'
    latitudes.long_name      = 'latitude'
    latitudes.units          = 'degrees north'
    latitudes.standard_name  = 'latitude'
    latitudes.axis           = 'Y'
    times.units              = f'years since {base_year}-01-01 00:00:00.0'
    times.calendar           = 'gregorian'
    srcs.long_name           = long_name
    srcs.unit                = unit
    # --- allocate data
    times[:] = range(syear-base_year, eyear-base_year+1)
    latitudes[:] = np.arange(89.75, -90, -0.5)
    longitudes[:] = np.arange(-179.75, 180, 0.5)
    srcs[:] = outSrc
    # --- close
    rootgrp.close()
    print(f'write out: {output_path}\n\n')


# ----------------------------------------------------------------------------------------------------------------------
def main():

    seamask = Dataset(landsea_mask_path)['LSM'][0,12:12+280,:].mask  # (280, 720)

    for rcp, soc in scenarios:
        for ghm in ghms:
            for gcm in gcms:

                output_directory = os.path.join(output_directory_main, supply_type, f'{rcp}_{soc}_{temporal_target}')
                if not os.path.isdir(output_directory): os.makedirs(output_directory)
    
                # load supply side
                qtot = np.array(load_modeloutput(gcm, ghm, rcp, soc, 'qtot'))          # (nghm, nyear, 12, 280, 720)  [m3/dy]
                if supply_type == 'Inflow':
                    inflow = np.array(load_modeloutput(gcm, ghm, rcp, soc, 'dis'))     # (nghm, nyear, 12, 280, 720)  [m3/dy]
                # load water demand (potential water withdrawal)
                irr_demand = np.array(load_modeloutput(gcm, ghm, rcp, soc, 'pirrww'))  # (nghm, nyear, 12, 280, 720)  [m3/dy]
                ind_demand = load_demandinput(soc, 'indww')                                              # (      nyear, 12, 280, 720)
                dom_demand = load_demandinput(soc, 'domww')                                              # (      nyear, 12, 280, 720
                print(f'qtot.shape: {qtot.shape}')
                print(f'inflow.shape: {inflow.shape}')
                print(f'irr_demand.shape: {irr_demand.shape}')
                print(f'ind_demand.shape: {ind_demand.shape}')
                print(f'dom_demand.shape: {dom_demand.shape}')
    
                # mask missing value
                qtot       = np.ma.masked_array(np.ma.masked_greater(qtot,       1e+19), mask=np.resize(seamask, qtot.shape))
                if supply_type == 'Inflow':
                    inflow = np.ma.masked_array(np.ma.masked_greater(inflow,     1e+19), mask=np.resize(seamask, inflow.shape))
                irr_demand = np.ma.masked_array(np.ma.masked_greater(irr_demand, 1e+19), mask=np.resize(seamask, irr_demand.shape))
                ind_demand = np.ma.masked_array(np.ma.masked_greater(ind_demand, 1e+19), mask=np.resize(seamask, ind_demand.shape))
                dom_demand = np.ma.masked_array(np.ma.masked_greater(dom_demand, 1e+19), mask=np.resize(seamask, dom_demand.shape))
    
                if supply_type == 'Inflow':
                    water_supply = qtot + inflow.filled(0)                                              # (nyear,12,280,720)  [m3/dy]
                else:
                    water_supply = qtot                                                                 # (nyear,12,280,720)  [m3/dy]
                water_demand = irr_demand.filled(0) + ind_demand.filled(0) + dom_demand.filled(0)       # (nyear,12,280,720)  [m3/dy]
                water_supply = np.ma.masked_array(water_supply, mask=np.resize(seamask, water_supply.shape))
                water_demand = np.ma.masked_array(water_demand, mask=np.resize(seamask, water_demand.shape))
    
                # calculate WSI at the temporal_target
                if temporal_target == 'Yearly':
                    water_stress = water_demand.mean(axis=1) / water_supply.mean(axis=1)                # (nyear,280,720)  [-]
                    water_stress[np.where(water_supply.mean(axis=1)==0)] = 1.
                elif temporal_target == 'DryMonth':
                    water_stress = (water_demand / water_supply).max(axis=1)                            # (nyear,280,720)  [-]
                    water_stress[np.where(water_supply.min(1)==0)] = 1.
    
                water_stress[np.where(water_stress>1)] = 1.
                water_stress = np.ma.masked_array(water_stress, mask=np.resize(seamask, water_stress.shape))
                write_nc('Stress', water_stress, ghm, gcm, rcp, soc, output_directory)
    

if __name__=='__main__':
    main()
