# This script contains code to regrid gridded climate data from 
# PRISM and gridMET, as well as county-level WNV case data, into a
# hexagonal grid of specified width.
#
# Initialized by RDH on 7/7/2023

#%% Rewriting code to assign data to a hexagonal grid because I forgot to save before exiting Vim........

## start by converting all .bil PRISM files to .nc files
# these two steps convert .bil PRISM files to .nc including
from osgeo import gdal
import glob
all_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/PRISM_vpdmin_stable_4kmM3_198101_202303_bil/*.bil')
for ifile in all_files:
    ofile = ifile[:-3] + 'nc'
    ds = gdal.Translate(ofile, ifile, format='NetCDF')

combined_ds = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/PRISM_vpdmin_stable_4kmM3_198101_202303_bil/*.nc', combine='nested', concat_dim='time', join='override')
time_array = pd.date_range(start='1981-01', end='2023-03', freq='MS')
combined_ds['time'] = time_array
combined_ds = combined_ds.drop_vars('crs')
combined_ds = combined_ds.rename({'Band1': 'mon_min_vpd'})
combined_ds['mon_min_vpd'].attrs['long_name'] = 'Monthly Min Vapor Deficit (hPa)'
combined_ds.to_netcdf('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_vpdmin_198101-202303.nc')

#
ds = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/soil_moisture/*.nc4')


# pasting code to create a hexagonal grid from github
# pulled from https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html
from shapely.geometry import Polygon
from shapely.geometry import Point
import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import time

ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/soil_moisture/PRISM_vpdmax_198101-202303.nc')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

unit = 2
xmin = -125 - 2 * unit
xmax = -66 + 2 * unit
ymin = 24 - 2 * unit
ymax = 50 + 2 * unit

a = np.sin(np.pi / 3)
cols = np.arange(np.floor(xmin), np.ceil(xmax), 3 * unit)
rows = np.arange(np.floor(ymin) / a, np.ceil(ymax) / a, unit)

hexagons = []
for x in cols:
  for i, y in enumerate(rows):
    if (i % 2 == 0):
      x0 = x
    else:
      x0 = x + 1.5 * unit
    hexagons.append(Polygon([
      (x0, y * a),
      (x0 + unit, y * a),
      (x0 + (1.5 * unit), (y + unit) * a),
      (x0 + unit, (y + (2 * unit)) * a),
      (x0, (y + (2 * unit)) * a),
      (x0 - (0.5 * unit), (y + unit) * a),
    ]))

grid = gpd.GeoDataFrame({'geometry': hexagons})
# grid["grid_area"] = grid.area
grid = grid.reset_index().rename(columns={"index": "grid_id"})

grid['min_lon'] = grid.bounds['minx']
grid['max_lon'] = grid.bounds['maxx']
grid['min_lat'] = grid.bounds['miny']
grid['max_lat'] = grid.bounds['maxy']

grid.to_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')


## use previously developed code to create a mask which assigns each point to a given hex grid cell
# identifying which hex a given county falls within
def find_which_shapefile(lon, lat, shapefile_list):
    pt = Point([lon, lat])
    for i in np.arange(len(shapefile_list)):
        if pt.within(shapefile_list[i]):
            break
    return i

# looping through all grid points to create mask
grid = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan
c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_hexes = (lat > grid['min_lat']) & (lat < grid['max_lat']) & (lon > grid['min_lon']) & (lon < grid['max_lon'])
    if np.sum(possible_hexes) == 1:
        grid_mask[ind] = grid['grid_id'][possible_hexes]
    elif np.sum(possible_hexes) == 0:
        print('Failed: ' + str(ind))
    else:
        hex_list = grid['grid_id'][possible_hexes]
        hex_shapefiles = list()
        for hex_ind in hex_list:
            hex_shapefiles.append(grid[grid['grid_id'] == hex_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, hex_shapefiles)
        grid_mask[ind] = hex_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')
grid['num_hexes'] = np.nan
for g in grid['grid_id']:
    grid['num_hexes'].loc[g] = np.sum(grid_mask==g)

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_grid_2_deg.npy', grid_mask)


# getting mean of all points in a given grid cell and transforming to a dataframe for saving
grid_array = np.zeros([len(ds.time), np.shape(grid)[0]]); grid_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
# grid_array = np.zeros([len(ds.day), np.shape(grid)[0]]); grid_array[:,:] = np.nan
start_time = time.time()
for g in grid['grid_id']:
    if g % 10 == 0:
        print('Grid ID: ' + str(g) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    if grid['num_hexes'][grid['grid_id'] == g].values[0] == 0:
        continue
    grid_array[:, g] = ds['SOILM'].where(np.transpose(grid_mask) == g).mean(dim=['lat', 'lon'], skipna=True).compute()  # transposed for PRISM, not transposed for gridMET
ds_processed = pd.DataFrame(data=grid_array, index=ds.time.values, columns = grid['grid_id'].values)
ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/NLDAS_soil_temp_' + str(depth) + '_199901-202212.csv')

grid_array = np.zeros([len(ds.time), np.shape(grid)[0]]); grid_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
# grid_array = np.zeros([len(ds.day), np.shape(grid)[0]]); grid_array[:,:] = np.nan
start_time = time.time()
for d in ds.depth.values:
    ds_level = ds.sel(depth=d)
    for g in grid['grid_id']:
        if g % 10 == 0:
            print('Grid ID: ' + str(g) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
        if grid['num_hexes'][grid['grid_id'] == g].values[0] == 0:
            continue
        grid_array[:, g] = ds_level['SOILM'].where(np.transpose(grid_mask) == g).mean(dim=['lat', 'lon'], skipna=True).compute()  # transposed for PRISM, not transposed for gridMET
    ds_processed = pd.DataFrame(data=grid_array, index=ds_level.time.values, columns = grid['grid_id'].values)
    ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/NLDAS_soil_moisture_' + str(int(d)) + '_199901-202212.csv')


#%% Prepping for WNV data assignment
# loading in hex grid
side_len = 2
gdf = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_' + str(side_len) + '_deg.shp')
bounds = gdf.bounds
gdf = gdf.join(bounds)


# loading in WNV data and supporting county info
df_wnv = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/NeuroWNV_by_county_2000-2022_monthly.csv')
df_pop_cens = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_centers_of_pop.txt', dtype={'STATEFP': str, 'COUNTYFP': str})
df_pop_cens['FIPS'] = df_pop_cens['STATEFP'] + df_pop_cens['COUNTYFP']
excluded_states = ['02', '15', '72']

