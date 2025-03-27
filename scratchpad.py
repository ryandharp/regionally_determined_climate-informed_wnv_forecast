# This script is a scratchpad for code development related to the WNV spatiotemporal pattern project.
# Initial plans are to run an EOF analysis on the interannual variability of WNV. 
# 
# Initialized by RDH on 1/3/2023

## Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd

## Loading data
df = pd.read_csv('../raw_data/NeuroWNV_by_county_2000-2020_monthly.csv')
df = df.drop(columns='PERCENT')
min_year = 2003
df = df[df['Year'] > min_year]

## Setting up time series array
max_year = df['Year'].max(axis=0)
num_counties = df['County'].nunique()
WNV_array = np.empty([(max_year-min_year)*12, num_counties])
# WNV_array[:] = np.nan

## Transforming dataframe to county-month time series array
df['time_step'] = (df['Year']-(min_year+1))*12+df['month']
df = df[~np.isnan(df['time_step'])]
df['time_step'] = df['time_step']-1
df['time_step'] = df['time_step'].astype('int')

WNV_County = df['County'].unique()
for index, row in df.iterrows():
    WNV_array[int(row['time_step']), (WNV_County == row['County']).nonzero()] = row['COUNT']

## Running the EOF analysis
from eofs.standard import Eof
from eofs.examples import example_data_path

solver = Eof(WNV_array)  # this actually runs the EOF
eof_results = solver.eofsAsCorrelation(neofs=2)  # grabbing the first EOF
pc_results = solver.pcs(npcs=2, pcscaling=1)  # grabbing the first PC
varfrac = solver.varianceFraction()

## Merging county locations for plotting
county_locs = pd.read_csv('../supplemental_files/2021_Gaz_counties_national.txt', sep='\t', header=0)
county_locs.columns.values[-1] = 'INTPTLON'
WNV_lat = np.empty(np.shape(WNV_County))
WNV_lat[:] = np.nan
WNV_lon = np.empty(np.shape(WNV_County))
WNV_lon[:] = np.nan
for index, fips in enumerate(WNV_County):
    if fips == 19999:
       print(index)
       continue
    WNV_lat[index] = county_locs['INTPTLAT'][county_locs['GEOID']==fips].values
    WNV_lon[index] = county_locs['INTPTLON'][county_locs['GEOID']==fips].values
WNV_lat = np.delete(WNV_lat, 522)
WNV_lon = np.delete(WNV_lon, 522)

## Trying to plot eof1
eof1 = np.squeeze(eof_results[0])
eof1 = np.delete(eof1, 522)
plt.scatter(WNV_lon, WNV_lat, c=eof1, s=5, cmap='Blues')
plt.colorbar()
plt.show()

eof2 = np.squeeze(eof_results[1])
eof2 = np.delete(eof2, 522)
plt.scatter(WNV_lon, WNV_lat, c=eof2, s=5, cmap='Blues')
plt.colorbar()
plt.show()

plt.plot(pc_results[:,0])
plt.xticks(np.arange(0, (max_year - min_year + 1)*12, 12), np.arange(min_year, max_year + 1) )
plt.show()

plt.plot(pc_results[:,1])
plt.xticks(np.arange(0, (max_year - min_year + 1)*12, 12), np.arange(min_year, max_year + 1) )
plt.show()


## Attempting to Create a Hexagonal Grid
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import Point
import time

N = 121
ratio = np.sqrt(3)/2
N_X = int(np.sqrt(N)/ratio)
N_Y = N // N_X
xv, yv = np.meshgrid(np.arange(N_X), np.arange(N_Y), sparse = False, indexing = 'xy')

xv = xv * ratio
xv[::2, :] = xv[::2, :] + ratio/2

lat_point_list = [yv[i-1, j-1], yv[i-1, j], yv[i, j+1], yv[i+1, j], yv[i+1, j-1], yv[i, j-1], yv[i-1, j-1]]
lon_point_list = [xv[i-1, j-1], xv[i-1, j], xv[i, j], xv[i+1, j], xv[i+1, j-1], xv[i, j-2], xv[i-1, j-1]]
fig, ax = plt.subplots()
ax.scatter(xv, yv, s=1)
ax.scatter(lon_point_list, lat_point_list, s = 5)
plt.show()


polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
polygon = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon_geom])

lon_min = -125 - 1
lon_max = -67 + 2
lat_min = 25 - 2
lat_max = 53 + 1
side_len = 1
ratio = np.sqrt(3)/2 * side_len

lat_array = np.arange(lat_max, lat_min, -ratio)
lon_array = np.arange(lon_min, lon_max, side_len) 

xv, yv = np.meshgrid(lon_array, lat_array, sparse = False, indexing = 'xy')
xv = xv.astype('float')
xv[1::2] = xv[1::2, :] + side_len/2


# recombobulating
polygon_geom = list()
lon_min = -125 - 1
lon_max = -67 + 2
lat_min = 25 - 2
lat_max = 53 + 1
side_len = 2
ratio = np.sqrt(3)/2 * side_len

lat_array = np.arange(lat_max, lat_min, -ratio)
lon_array = np.arange(lon_min, lon_max, side_len) 

xv, yv = np.meshgrid(lon_array, lat_array, sparse = False, indexing = 'xy')
xv = xv.astype('float')
xv[1::2] = xv[1::2, :] + side_len/2

for j in np.arange(1, len(lat_array)-1, 2):
    for i in np.arange(1, len(lon_array)-1, 3):
        lat_point_list = [yv[j-1, i], yv[j-1, i+1], yv[j, i+1], yv[j+1, i+1], yv[j+1, i], yv[j, i-1], yv[j-1, i]]
        lon_point_list = [xv[j-1, i], xv[j-1, i+1], xv[j, i+1], xv[j+1, i+1], xv[j+1, i], xv[j, i-1], xv[j-1, i]]
        polygon_geom.append(Polygon(zip(lon_point_list, lat_point_list)))

for j in np.arange(2, len(lat_array)-1, 2):
    for i in np.arange(3, len(lon_array)-1, 3):
        lat_point_list = [yv[j-1, i-1], yv[j-1, i], yv[j, i+1], yv[j+1, i], yv[j+1, i-1], yv[j, i-1], yv[j-1, i-1]]
        lon_point_list = [xv[j-1, i-1], xv[j-1, i], xv[j, i+1], xv[j+1, i], xv[j+1, i-1], xv[j, i-1], xv[j-1, i-1]]
        polygon_geom.append(Polygon(zip(lon_point_list, lat_point_list)))

polygon = gpd.GeoDataFrame(crs='epsg:4326', geometry=polygon_geom)
polygon.to_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/hex_grid_' + str(side_len) + '_deg.shp')



# loading modules
import pandas as pd
import numpy as np
import geopandas as gpd




## Loading data
df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_time_series_2_deg.csv', index_col=0)
df.index = pd.to_datetime(df.index)
df['year'] = pd.DatetimeIndex(df.index).year
min_year = 2005
df = df[df['year'] >= min_year]
df = df.groupby('year').sum()
df = df.transpose()
wnv_array = np.array(df)


## Running the EOF analysis
from eofs.standard import Eof
from eofs.examples import example_data_path

solver = Eof(np.transpose(wnv_array))  # this actually runs the EOF
eof_results = solver.eofsAsCorrelation(neofs=3)  # grabbing the first EOF
pc_results = solver.pcs(npcs=3, pcscaling=1)  # grabbing the first PC
varfrac = solver.varianceFraction()


# grabbing hex grid info for plotting
gdf = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/hex_grid_2_deg.shp')
gdf['centroid'] = gdf['geometry'].centroid
gdf['hex_lats'] = gdf['centroid'].y
gdf['hex_lons'] = gdf['centroid'].x
df.index = df.index.astype(int)
df = df.merge(gdf, how='left', left_on=df.index, right_on='FID')
df = df.drop(columns=['geometry', 'centroid'])


## Merging county locations for plotting
county_locs = pd.read_csv('../supplemental_files/2021_Gaz_counties_national.txt', sep='\t', header=0)
county_locs.columns.values[-1] = 'INTPTLON'
WNV_lat = np.empty(np.shape(WNV_County))
WNV_lat[:] = np.nan
WNV_lon = np.empty(np.shape(WNV_County))
WNV_lon[:] = np.nan
for index, fips in enumerate(WNV_County):
    if fips == 19999:
       print(index)
       continue
    WNV_lat[index] = county_locs['INTPTLAT'][county_locs['GEOID']==fips].values
    WNV_lon[index] = county_locs['INTPTLON'][county_locs['GEOID']==fips].values
WNV_lat = np.delete(WNV_lat, 522)
WNV_lon = np.delete(WNV_lon, 522)

## Trying to plot eof1
eof1 = np.squeeze(eof_results[0])
eof1 = np.delete(eof1, 522)
plt.scatter(df['hex_lons'], df['hex_lats'], c=eof1, s=5, cmap='Blues')
plt.colorbar()
plt.show()

eof2 = np.squeeze(eof_results[1])
eof2 = np.delete(eof2, 522)
plt.scatter(WNV_lon, WNV_lat, c=eof2, s=5, cmap='Blues')
plt.colorbar()
plt.show()

plt.plot(pc_results[:,0])
plt.xticks(np.arange(0, (max_year - min_year + 1)*12, 12), np.arange(min_year, max_year + 1) )
plt.show()

plt.plot(pc_results[:,1])
plt.xticks(np.arange(0, (max_year - min_year + 1)*12, 12), np.arange(min_year, max_year + 1) )
plt.show()


# plotting eof1 results
import matplotlib
fig, ax = plt.subplots(1, 1)
gdf.plot(column='eof1', ax=ax, legend=True, vmin=-1, vmax=1, cmap='bwr')
plt.show()

plt.plot(pc_results[:,2])
plt.xticks(np.arange(18), ['2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'])

plt.plot(np.arange(1, 11, 1), varfrac[:10])
plt.show()


# testing with preprocessing on the WNV data
df = df[df.sum(axis=1)>=18]
df2 = df.subtract(df.mean(axis=1), axis=0).divide(df.std(axis=1), axis=0)
wnv_array = np.array(df2)
df2 = df2.merge(gdf, how='left', left_on='hex', right_on='FID')
df2['eof1'] = eof1
df2['eof2'] = eof2
df3['eof3'] = eof3
df2.index = df2.index.astype(int)
gdf = gdf.merge(df2, how='left', left_on='FID', right_on='hex')

fig, ax = plt.subplots(1,1)
gdf.plot(column='eof1', ax=ax, legend=True, vmin=-1, vmax=1, cmap='bwr_r')


import cartopy.crs as ccrs
import cartopy.feature as cfeature

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

gdf['corr_map'] = df3.corr(method='kendall')[corr_num]
gdf['corr_map_sq'] = gdf['corr_map'].pow(2) * np.sign(gdf['corr_map'])
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column='corr_map', ax=ax, legend=True, vmin=-1, vmax=1, cmap='bwr_r')
plt.scatter(gdf['hex_lons_x'][gdf['FID'] == corr_num], gdf['hex_lats_y'][gdf['FID'] == corr_num], c='k', s=100)
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=18)
plt.yticks(np.arange(30, 60, 10), fontsize=18)
plt.show()

# corr_nums:
    # TX: 129
    # OR/WA: 21
    # ND/SD: 92
    # Maricopa: 52
    # Denver: 43


## plotting WNV data
# import everything
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt


# loading WNV data
df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_time_series_2_deg.csv', index_col=0)
df.index = pd.to_datetime(df.index)
df['year'] = pd.DatetimeIndex(df.index).year
min_year = 2001
df = df[df['year'] >= min_year]
df = df.groupby('year').sum()
df = df.transpose()
df.index = df.index.rename('FID')

# loading hex grid data
gdf = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/hex_grid_2_deg.shp')
gdf['centroid'] = gdf['geometry'].centroid
gdf['hex_lats'] = gdf['centroid'].y
gdf['hex_lons'] = gdf['centroid'].x
df.index = df.index.astype(int)

gdf = gdf.merge(df, how='left', left_on='FID', right_on='FID')

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

for yr in np.arange(2001, 2023, 1):
# yr = 2021
    gdf[yr][gdf[yr]>0] = np.log10(gdf[yr])
    gdf[yr][np.isnan(gdf[yr])] = 0
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
    ax.set_aspect('equal')
    gdf.plot(column=yr, ax=ax, legend=True, vmin=0, vmax=2.5, cmap='Blues')
    # gdf.plot(column=yr, ax=ax, legend=True, cmap='Blues')
    plt.xlim([-129, -62])
    plt.ylim([23, 52])
    plt.xticks(np.arange(-125, -55, 10), fontsize=18)
    plt.yticks(np.arange(30, 60, 10), fontsize=18)
    plt.title(str(yr))
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/' + str(yr) + '_log.png')
    plt.close()


# loading WNV data
df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_time_series_2_deg.csv', index_col=0)
df.index = pd.to_datetime(df.index)
df['year'] = pd.DatetimeIndex(df.index).year
min_year = 2001
df = df[df['year'] >= min_year]
df = df.groupby('year').sum()
df = df[df>0].transpose().rank(pct=True, axis=1)
df.index = df.index.rename('FID')

# loading hex grid data
gdf = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/hex_grid_2_deg.shp')
gdf['centroid'] = gdf['geometry'].centroid
gdf['hex_lats'] = gdf['centroid'].y
gdf['hex_lons'] = gdf['centroid'].x
df.index = df.index.astype(int)
gdf = gdf.merge(df, how='left', left_on='FID', right_on='FID')

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

for yr in np.arange(2001, 2023, 1):
# yr = 2021
    gdf[yr][np.isnan(gdf[yr])] = 0
    # gdf[yr][gdf[yr]<.8] = 0
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
    ax.set_aspect('equal')
    gdf.plot(column=yr, ax=ax, legend=True, vmin=0, vmax=1, cmap='Blues')
    # gdf.plot(column=yr, ax=ax, legend=True, cmap='bwr_r')
    plt.xlim([-129, -62])
    plt.ylim([23, 52])
    plt.xticks(np.arange(-125, -55, 10), fontsize=18)
    plt.yticks(np.arange(30, 60, 10), fontsize=18)
    plt.title(str(yr))
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/' + str(yr) + '_percentile.png')
    plt.close()

# using imagemagick at the terminal: convert -delay 100 -loop 0 *.png wnv_percentile_map.gif




#%% Starting preliminary analysis on climate/environmental factors
# starting with some composites of outbreak years
# loading modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd

# loading data
hex_grid = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/updated_hex_time_series_2_deg.csv')  # for hex-grid analysis
temp_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_tmean_198101-202211_gridded-2.csv')
precip_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_precip_198101-202211_gridded-2.csv')
dew_point_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_td_198101-202211_gridded-2.csv')
pdsi_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridMET_PDSI_199901-202212_gridded-2.csv')

# wnv_data_test = pd.DataFrame(index = time_index, columns=np.sort(wnv_data['County'].unique().astype(int)))
# wnv_data_test = wnv_data_test.fillna(value=0)
# for i, row in wnv_data.iterrows():
#     if (np.isnan(row['month']) | np.isnan(row['County'])):
#         print(row)
#         continue
#     wnv_data_test.loc[row['date']][int(row['County'])] = wnv_data_test.loc[row['date']][int(row['County'])] + row['COUNT']
# wnv_data = wnv_data_test.copy()

# identifying highest case load hexes
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_sums = wnv_data.sum(axis=0)
wnv_sums = wnv_sums.sort_values(ascending=False)

wnv_sums[:31]  # this gives you 80% of all caseload for hexes
wnv_ann = wnv_data.resample('YS').sum()
wnv_data = wnv_data.resample('Y').sum()

# testing some time series hypotheses
test_rank = wnv_data['27'].rank().astype(int)
test_value = wnv_data['27']

comp_array_zero = np.array([])
comp_array_one = np.array([])
comp_array_two = np.array([])
comp_array_zero_two = np.array([])
for hex in wnv_sums[:31].index:
    test_rank = wnv_data[hex].rank().astype(int)
    test_value = wnv_data[hex]
    max_lag = 2
    lag_array_rank = np.zeros([max_lag+1, len(test_rank)+2])
    lag_array_value = np.zeros([max_lag+1, len(test_value)+2])
    for lag in np.arange(0, max_lag + 1):
        for i in np.arange(len(test_rank)) + lag:
            lag_array_rank[lag, i] = test_rank[i-lag]
            lag_array_value[lag, i] = test_value[i-lag]
    lag_array_value = lag_array_value[:, 2:-2]
    # lag_array_value[lag_array_value == 0] = 1
    lag_array_value = lag_array_value + 1
    lag_array_log_value = np.log10(abs(lag_array_value))
    norm_value = lag_array_log_value - lag_array_log_value[0,:].mean()
    plt.scatter(norm_value[2,:], norm_value[0,:], s=2)
    comp_array_zero = np.append(comp_array_zero, norm_value[0, :])
    comp_array_one = np.append(comp_array_one, norm_value[1, :])
    comp_array_two = np.append(comp_array_two, norm_value[2, :])
    plt.title(str(hex) + ' two year lag')




## redone for county-level calculations
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/NeuroWNV_by_county_2000-2022_monthly.csv')
wnv_data['date'] = pd.to_datetime(dict(year=df_wnv['Year'], month=df_wnv['month'], day=1))
time_index = pd.date_range(wnv_data['date'].min(), wnv_data['date'].max(), freq='MS')

wnv_data_test = pd.DataFrame(index = time_index, columns=np.sort(wnv_data['County'].unique().astype(int)))
wnv_data_test = wnv_data_test.fillna(value=0)
for i, row in wnv_data.iterrows():
    if (np.isnan(row['month']) | np.isnan(row['County'])):
        print(row)
        continue
    wnv_data_test.loc[row['date']][int(row['County'])] = wnv_data_test.loc[row['date']][int(row['County'])] + row['COUNT']

del wnv_data
wnv_data = wnv_data_test.copy()

# identifying highest case load hexes
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_sums = wnv_data.sum(axis=0)
wnv_sums = wnv_sums.sort_values(ascending=False)

wnv_sums[:49]  # this gives you 50% of all caseload for counties
wnv_ann = wnv_data.resample('YS').sum()
wnv_data = wnv_data.resample('Y').sum()

# testing some time series hypotheses
comp_array_zero = np.array([])
comp_array_one = np.array([])
comp_array_two = np.array([])
comp_array_zero_two = np.array([])
county_correlations_one = np.array([])
county_correlations_two = np.array([])
for hex in wnv_sums[:49].index:
    test_rank = wnv_data[hex].rank().astype(int)
    test_value = wnv_data[hex]
    max_lag = 2
    lag_array_rank = np.zeros([max_lag+1, len(test_rank)+2])
    lag_array_value = np.zeros([max_lag+1, len(test_value)+2])
    for lag in np.arange(0, max_lag + 1):
        for i in np.arange(len(test_rank)) + lag:
            lag_array_rank[lag, i] = test_rank[i-lag]
            lag_array_value[lag, i] = test_value[i-lag]
    lag_array_value = lag_array_value[:, 2:-2]
    lag_array_value = lag_array_value + 1
    lag_array_log_value = np.log10(abs(lag_array_value))
    norm_value = lag_array_log_value - lag_array_log_value[0,:].mean()
    lag_array_value = norm_value
    plt.scatter(lag_array_value[2,:], lag_array_value[0,:], s=2)
    comp_array_zero = np.append(comp_array_zero, lag_array_value[0, :])
    comp_array_one = np.append(comp_array_one, lag_array_value[1, :])
    comp_array_two = np.append(comp_array_two, lag_array_value[2, :])
    plt.title(str(hex) + ' two year lag')
    county_correlations_one = np.append(county_correlations_one, np.corrcoef(lag_array_value[1,:], lag_array_value[0,:])[0,1])
    county_correlations_two = np.append(county_correlations_two, np.corrcoef(lag_array_value[2,:], lag_array_value[0,:])[0,1])







df_wnv = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/NeuroWNV_by_county_2000-2022_monthly.csv')

temp_data = temp_data.dropna(axis=1, how='all')
temp_data.index = temp_data['Unnamed: 0'].values
temp_data.index = pd.to_datetime(temp_data.index)
temp_data = temp_data.drop('Unnamed: 0', axis=1)
temp_data = temp_data[temp_data.index >= '2003-01-01']
temp_data['days_in_month'] = temp_data.index.days_in_month
temp_data['year'] = temp_data.index.year
# temp_data[temp_data.sum(axis=0) > 0]
temp_data['date'] = temp_data.index
temp_data.index = np.arange(np.shape(temp_data)[0])


cols = ['47']
col_label = '47'

# transforming data to seasonal mean
def seasonal_mean(data, cols):
    return pd.Series(np.average(data[cols], weights=data['days_in_month'], axis=0), cols)

seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
temp_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
temp_seasonal['seasonal_date'] = seasonal_date

for index, row in temp_seasonal.iterrows():
    if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
        i = temp_data[temp_data['date'] == row['seasonal_date']].index[0]
        seasonal_data = temp_data.iloc[i:i+3]
        test = seasonal_mean(seasonal_data, cols)
        for c in np.arange(len(cols)):  # this is awful but whatever it works
            temp_seasonal.loc[temp_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]


hex_df = pd.DataFrame(wnv_data[col_label])
hex_df = hex_df.rename(columns={col_label: 'wnv'})
hex_df['spring_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
hex_df['spring_temp_anom'] = (hex_df['spring_temp'] - hex_df['spring_temp'].mean())/hex_df['spring_temp'].std().astype(float)
hex_df['summer_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
hex_df['summer_temp_anom'] = (hex_df['summer_temp'] - hex_df['summer_temp'].mean())/hex_df['summer_temp'].std().astype(float)
hex_df['fall_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
hex_df['fall_temp_anom'] = (hex_df['fall_temp'] - hex_df['fall_temp'].mean())/hex_df['fall_temp'].std().astype(float)
hex_df['winter_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
hex_df['winter_temp_anom'] = (hex_df['winter_temp'] - hex_df['winter_temp'].mean())/hex_df['winter_temp'].std().astype(float)


hex_df['wnv'].nlargest(3)



pdsi_data = pdsi_data.dropna(axis=1, how='all')
pdsi_data.index = pdsi_data['Unnamed: 0'].values
pdsi_data.index = pd.to_datetime(pdsi_data.index)
pdsi_data = pdsi_data.drop('Unnamed: 0', axis=1)
pdsi_data = pdsi_data[pdsi_data.index >= '2003-01-01']
pdsi_data['days_in_month'] = pdsi_data.index.days_in_month
pdsi_data['year'] = pdsi_data.index.year
# temp_data[temp_data.sum(axis=0) > 0]
pdsi_data['date'] = pdsi_data.index
pdsi_data.index = np.arange(np.shape(pdsi_data)[0])

seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
pdsi_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
pdsi_seasonal['seasonal_date'] = seasonal_date

for index, row in pdsi_seasonal.iterrows():
    if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
        i = pdsi_data[pdsi_data['date'] == row['seasonal_date']].index[0]
        seasonal_data = pdsi_data.iloc[i:i+3]
        test = seasonal_mean(seasonal_data, cols)
        for c in np.arange(len(cols)):  # this is awful but whatever it works
            pdsi_seasonal.loc[pdsi_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]

hex_df['spring_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
hex_df['spring_pdsi_anom'] = (hex_df['spring_pdsi'] - hex_df['spring_pdsi'].mean())/hex_df['spring_pdsi'].std().astype(float)
hex_df['summer_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
hex_df['summer_pdsi_anom'] = (hex_df['summer_pdsi'] - hex_df['summer_pdsi'].mean())/hex_df['summer_pdsi'].std().astype(float)
hex_df['fall_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
hex_df['fall_pdsi_anom'] = (hex_df['fall_pdsi'] - hex_df['fall_pdsi'].mean())/hex_df['fall_pdsi'].std().astype(float)
hex_df['winter_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
hex_df['winter_pdsi_anom'] = (hex_df['winter_pdsi'] - hex_df['winter_pdsi'].mean())/hex_df['winter_pdsi'].std().astype(float)



precip_data = precip_data.dropna(axis=1, how='all')
precip_data.index = precip_data['Unnamed: 0'].values
precip_data.index = pd.to_datetime(precip_data.index)
precip_data = precip_data.drop('Unnamed: 0', axis=1)
precip_data = precip_data[precip_data.index >= '2003-01-01']
precip_data['days_in_month'] = precip_data.index.days_in_month
precip_data['year'] = precip_data.index.year
# temp_data[temp_data.sum(axis=0) > 0]
precip_data['date'] = precip_data.index
precip_data.index = np.arange(np.shape(precip_data)[0])

seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
precip_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
precip_seasonal['seasonal_date'] = seasonal_date

for index, row in precip_seasonal.iterrows():
    if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
        i = precip_data[precip_data['date'] == row['seasonal_date']].index[0]
        seasonal_data = precip_data.iloc[i:i+3]
        test = seasonal_mean(seasonal_data, cols)
        for c in np.arange(len(cols)):  # this is awful but whatever it works
            precip_seasonal.loc[precip_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]

hex_df['spring_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
hex_df['spring_precip_anom'] = (hex_df['spring_precip'] - hex_df['spring_precip'].mean())/hex_df['spring_precip'].std().astype(float)
hex_df['summer_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
hex_df['summer_precip_anom'] = (hex_df['summer_precip'] - hex_df['summer_precip'].mean())/hex_df['summer_precip'].std().astype(float)
hex_df['fall_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
hex_df['fall_precip_anom'] = (hex_df['fall_precip'] - hex_df['fall_precip'].mean())/hex_df['fall_precip'].std().astype(float)
hex_df['winter_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
hex_df['winter_precip_anom'] = (hex_df['winter_precip'] - hex_df['winter_precip'].mean())/hex_df['winter_precip'].std().astype(float)


# creating some composites
hex_df['wnv'].nlargest(3)  # pulling three largest outbreaks
wnv_threshold = hex_df['wnv'].nlargest(3).values[-1]

hex_df[hex_df['wnv'] >= wnv_threshold].mean(axis=0)
hex_df.corr()['wnv']



# looping the above
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# identifying highest case load hexes
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/updated_hex_time_series_2_deg.csv')
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_sums = wnv_data.sum(axis=0)
wnv_sums = wnv_sums.sort_values(ascending=False)


def seasonal_mean(data, cols):
    return pd.Series(np.average(data[cols], weights=data['days_in_month'], axis=0), cols)


num_hexes = len(wnv_sums[wnv_sums>=19])
cols = wnv_sums[:num_hexes].index 
# pearson_corrs = np.zeros([num_hexes, 12]); pearson_corrs[:] = np.nan
# spearman_corrs = np.zeros([num_hexes, 12]); spearman_corrs[:] = np.nan
# outbreak_anoms = np.zeros([num_hexes, 24]); outbreak_anoms[:] = np.nan

counter = 0
for col in cols:
    print(counter)
    col_label = col
    col = [col]
    hex_grid = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
    wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/updated_hex_time_series_2_deg.csv')
    wnv_data.index = wnv_data['Unnamed: 0'].values
    wnv_data.index = pd.to_datetime(wnv_data.index)
    wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
    wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
    wnv_ann = wnv_data.resample('YS').sum()
    temp_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_tmean_198101-202211_gridded-2.csv')
    precip_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_precip_198101-202211_gridded-2.csv')
    dew_point_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/PRISM_td_198101-202211_gridded-2.csv')
    pdsi_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridMET_PDSI_199901-202212_gridded-2.csv')
    temp_data = temp_data.dropna(axis=1, how='all')
    temp_data.index = temp_data['Unnamed: 0'].values
    temp_data.index = pd.to_datetime(temp_data.index)
    temp_data = temp_data.drop('Unnamed: 0', axis=1)
    temp_data = temp_data[temp_data.index >= '2003-01-01']
    temp_data['days_in_month'] = temp_data.index.days_in_month
    temp_data['year'] = temp_data.index.year
    # temp_data[temp_data.sum(axis=0) > 0]
    temp_data['date'] = temp_data.index
    temp_data.index = np.arange(np.shape(temp_data)[0])
    seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
    temp_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
    temp_seasonal['seasonal_date'] = seasonal_date
    for index, row in temp_seasonal.iterrows():
        if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
            i = temp_data[temp_data['date'] == row['seasonal_date']].index[0]
            seasonal_data = temp_data.iloc[i:i+3]
            test = seasonal_mean(seasonal_data, cols)
            for c in np.arange(len(cols)):  # this is awful but whatever it works
                temp_seasonal.loc[temp_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]
    hex_df = pd.DataFrame(wnv_ann[col_label])
    hex_df = hex_df.rename(columns={col_label: 'wnv'})
    hex_df['spring_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
    hex_df['spring_temp_anom'] = (hex_df['spring_temp'] - hex_df['spring_temp'].mean())/hex_df['spring_temp'].std().astype(float)
    hex_df['summer_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
    hex_df['summer_temp_anom'] = (hex_df['summer_temp'] - hex_df['summer_temp'].mean())/hex_df['summer_temp'].std().astype(float)
    hex_df['fall_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
    hex_df['fall_temp_anom'] = (hex_df['fall_temp'] - hex_df['fall_temp'].mean())/hex_df['fall_temp'].std().astype(float)
    hex_df['winter_temp'] = temp_seasonal[temp_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
    hex_df['winter_temp_anom'] = (hex_df['winter_temp'] - hex_df['winter_temp'].mean())/hex_df['winter_temp'].std().astype(float)
    pdsi_data = pdsi_data.dropna(axis=1, how='all')
    pdsi_data.index = pdsi_data['Unnamed: 0'].values
    pdsi_data.index = pd.to_datetime(pdsi_data.index)
    pdsi_data = pdsi_data.drop('Unnamed: 0', axis=1)
    pdsi_data = pdsi_data[pdsi_data.index >= '2003-01-01']
    pdsi_data['days_in_month'] = pdsi_data.index.days_in_month
    pdsi_data['year'] = pdsi_data.index.year
    # temp_data[temp_data.sum(axis=0) > 0]
    pdsi_data['date'] = pdsi_data.index
    pdsi_data.index = np.arange(np.shape(pdsi_data)[0])
    seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
    pdsi_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
    pdsi_seasonal['seasonal_date'] = seasonal_date
    for index, row in pdsi_seasonal.iterrows():
        if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
            i = pdsi_data[pdsi_data['date'] == row['seasonal_date']].index[0]
            seasonal_data = pdsi_data.iloc[i:i+3]
            test = seasonal_mean(seasonal_data, cols)
            for c in np.arange(len(cols)):  # this is awful but whatever it works
                pdsi_seasonal.loc[pdsi_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]
    hex_df['spring_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
    hex_df['spring_pdsi_anom'] = (hex_df['spring_pdsi'] - hex_df['spring_pdsi'].mean())/hex_df['spring_pdsi'].std().astype(float)
    hex_df['summer_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
    hex_df['summer_pdsi_anom'] = (hex_df['summer_pdsi'] - hex_df['summer_pdsi'].mean())/hex_df['summer_pdsi'].std().astype(float)
    hex_df['fall_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
    hex_df['fall_pdsi_anom'] = (hex_df['fall_pdsi'] - hex_df['fall_pdsi'].mean())/hex_df['fall_pdsi'].std().astype(float)
    hex_df['winter_pdsi'] = pdsi_seasonal[pdsi_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
    hex_df['winter_pdsi_anom'] = (hex_df['winter_pdsi'] - hex_df['winter_pdsi'].mean())/hex_df['winter_pdsi'].std().astype(float)
    precip_data = precip_data.dropna(axis=1, how='all')
    precip_data.index = precip_data['Unnamed: 0'].values
    precip_data.index = pd.to_datetime(precip_data.index)
    precip_data = precip_data.drop('Unnamed: 0', axis=1)
    precip_data = precip_data[precip_data.index >= '2003-01-01']
    precip_data['days_in_month'] = precip_data.index.days_in_month
    precip_data['year'] = precip_data.index.year
    # temp_data[temp_data.sum(axis=0) > 0]
    precip_data['date'] = precip_data.index
    precip_data.index = np.arange(np.shape(precip_data)[0])
    seasonal_date = pd.date_range(start='2003-03', end='2022-11', freq='3MS')
    precip_seasonal = pd.DataFrame(columns=cols, index=np.arange(len(seasonal_date)))
    precip_seasonal['seasonal_date'] = seasonal_date
    for index, row in precip_seasonal.iterrows():
        if np.isin(row['seasonal_date'].month, [3, 6, 9, 12]):
            i = precip_data[precip_data['date'] == row['seasonal_date']].index[0]
            seasonal_data = precip_data.iloc[i:i+3]
            test = seasonal_mean(seasonal_data, cols)
            for c in np.arange(len(cols)):  # this is awful but whatever it works
                precip_seasonal.loc[precip_seasonal['seasonal_date'] == row['seasonal_date'], cols[c]] = test[c]
    hex_df['spring_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 3][col_label].values[1:].astype(float)
    hex_df['spring_precip_anom'] = (hex_df['spring_precip'] - hex_df['spring_precip'].mean())/hex_df['spring_precip'].std().astype(float)
    hex_df['summer_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 6][col_label].values[1:].astype(float)
    hex_df['summer_precip_anom'] = (hex_df['summer_precip'] - hex_df['summer_precip'].mean())/hex_df['summer_precip'].std().astype(float)
    hex_df['fall_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 9][col_label].values[:-1].astype(float)
    hex_df['fall_precip_anom'] = (hex_df['fall_precip'] - hex_df['fall_precip'].mean())/hex_df['fall_precip'].std().astype(float)
    hex_df['winter_precip'] = precip_seasonal[precip_seasonal['seasonal_date'].dt.month == 12][col_label].values.astype(float)
    hex_df['winter_precip_anom'] = (hex_df['winter_precip'] - hex_df['winter_precip'].mean())/hex_df['winter_precip'].std().astype(float)
    wnv_threshold = hex_df['wnv'].nlargest(3).values[-1]
    if counter == 0:
        print('Check: ' + str(counter))
        outbreak_anoms = np.round(hex_df[hex_df['wnv'] >= wnv_threshold].mean(axis=0)[1:].values, 2)
        pearson_corrs = np.round(hex_df.corr(method='pearson')['wnv'].values[1::2], 2)
        spearman_corrs = np.round(hex_df.corr(method='spearman')['wnv'].values[1::2], 2)
        counter += 1
        continue
    outbreak_anoms = np.vstack((outbreak_anoms, np.round(hex_df[hex_df['wnv'] >= wnv_threshold].mean(axis=0)[1:].values, 2)))
    pearson_corrs = np.vstack((pearson_corrs, np.round(hex_df.corr(method='pearson')['wnv'].values[1::2], 2)))
    spearman_corrs = np.vstack((spearman_corrs, np.round(hex_df.corr(method='spearman')['wnv'].values[1::2], 2)))
    counter += 1


df_outbreak_anoms = pd.DataFrame(data = outbreak_anoms, index = cols, columns = hex_df.columns[1:])
df_pearson_corrs = pd.DataFrame(data = pearson_corrs, index = cols, columns = hex_df.columns[1::2])
df_spearman_corrs = pd.DataFrame(data = spearman_corrs, index = cols, columns = hex_df.columns[1::2])

df_outbreak_anoms.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/outbreak_anoms.csv')
df_pearson_corrs.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/pearson_corrs.csv')
df_spearman_corrs.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/spearman_corrs.csv')

test_anoms = hex_grid.merge(df_outbreak_anoms, left_on='grid_id', right_on=df_outbreak_anoms.index.astype(int), how='left')
test20_anoms = test[test.index.isin(wnv_sums[:20].index.astype(int))]

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

col_name = 'spring_temp_anom'

for col_name in test20_anoms.columns[7:]:
    v_lim = np.max([np.abs(np.quantile(test20_anoms[col_name], q=.95)), np.abs(np.quantile(test20_anoms[col_name], q=.05))])
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
    ax.set_aspect('equal')
    test20_anoms.plot(column=col_name, ax=ax, legend=True, cmap='bwr_r', vmin=-1*v_lim, vmax=v_lim)
    plt.xlim([-129, -62])
    plt.ylim([23, 52])
    plt.xticks(np.arange(-125, -55, 10), fontsize=18)
    plt.yticks(np.arange(30, 60, 10), fontsize=18)
    plt.title(col_name)
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/anom_map_' + col_name + '.png')



test_pearson = hex_grid.merge(df_pearson_corrs, left_on='grid_id', right_on=df_pearson_corrs.index.astype(int), how='left')
test20_pearson = test_pearson[test_pearson.index.isin(wnv_sums[:20].index.astype(int))]

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

for col_name in test20_pearson.columns[7:]:
    v_lim = np.max([np.abs(np.quantile(test20_pearson[col_name], q=.95)), np.abs(np.quantile(test20_pearson[col_name], q=.05))])
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
    ax.set_aspect('equal')
    test20_pearson.plot(column=col_name, ax=ax, legend=True, cmap='bwr_r', vmin=-1*v_lim, vmax=v_lim)
    plt.xlim([-129, -62])
    plt.ylim([23, 52])
    plt.xticks(np.arange(-125, -55, 10), fontsize=18)
    plt.yticks(np.arange(30, 60, 10), fontsize=18)
    plt.title(col_name)
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/pearson_map_' + col_name + '.png')



test_spearman = hex_grid.merge(df_spearman_corrs, left_on='grid_id', right_on=df_spearman_corrs.index.astype(int), how='left')
test20_spearman = test_spearman[test_spearman.index.isin(wnv_sums[:20].index.astype(int))]

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

for col_name in test20_spearman.columns[7:]:
    v_lim = np.max([np.abs(np.quantile(test20_spearman[col_name], q=.95)), np.abs(np.quantile(test20_spearman[col_name], q=.05))])
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
    ax.set_aspect('equal')
    test20_spearman.plot(column=col_name, ax=ax, legend=True, cmap='bwr_r', vmin=-1*v_lim, vmax=v_lim)
    plt.xlim([-129, -62])
    plt.ylim([23, 52])
    plt.xticks(np.arange(-125, -55, 10), fontsize=18)
    plt.yticks(np.arange(30, 60, 10), fontsize=18)
    plt.title(col_name)
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/spearman_map_' + col_name + '.png')





num_hexes = 20
cols = wnv_sums[:num_hexes].index 
# pearson_corrs = np.zeros([num_hexes, 12]); pearson_corrs[:] = np.nan
# spearman_corrs = np.zeros([num_hexes, 12]); spearman_corrs[:] = np.nan
# outbreak_anoms = np.zeros([num_hexes, 24]); outbreak_anoms[:] = np.nan

for col in cols:
    col_data = wnv_data[col]
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.set_style(style='darkgrid')
    fig = sns.lineplot(x=col_data.index.month, y=col_data.values, hue=col_data.index.year)
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/seasonal_cycle_' + col + '.png')



# gdf.plot(column='corr_map', ax=ax, legend=True, vmin=-1, vmax=1, cmap='bwr_r')







# creating some composites
hex_df['wnv'].nlargest(3)  # pulling three largest outbreaks
wnv_threshold = hex_df['wnv'].nlargest(3).values[-1]

hex_df[hex_df['wnv'] >= wnv_threshold].mean(axis=0)
hex_df.corr()['wnv']





#%% Maps of regional correlations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
corr_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/monthly_log_wnv_monthly_clim_corr.csv')
corr_df = corr_df.drop(labels='Unnamed: 0', axis=1)
clim_vars = corr_df['climate_var'].unique()


#--------------------------------------------------------------------------------------
#
#%% Making maps of correlations between climate variables and monthly log WNV caseload
#
#--------------------------------------------------------------------------------------

matplotlib.use('AGG')
for wnv_month in np.arange(1, 13, 1):
    for clim_month in np.arange(1, 13, 1):
        for clim_var in clim_vars:
        # for clim_var in ['Temperature', 'Drought', 'Precipitation']:
            hex_grid = gpd.read_file(
                '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
            test = corr_df[(corr_df['WNV_month']==wnv_month) & (corr_df['climate_var']==clim_var) & (corr_df['climate_month']==clim_month)]
            test.loc[:, 'hex'] = test['hex'].astype(int).values
            test = test[['hex', 'correlation']]
            gdf = hex_grid.merge(test, how='right', left_on='grid_id', right_on='hex')

            fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS)
            ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
            ax.set_aspect('equal')
            gdf.plot(column='correlation', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'correlation'},vmin=-0.5, vmax=0.5, cmap='bwr')
            plt.xlim([-129, -62])
            plt.ylim([23, 52])
            plt.xticks(np.arange(-125, -55, 10), fontsize=14)
            plt.yticks(np.arange(30, 60, 10), fontsize=14)
            plt.title(month_names[wnv_month] + ' WNND and ' + month_names[clim_month] + ' ' + clim_var, fontsize=18)
            plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/compiled_correlations/monthly_log_wnv_'  + str(wnv_month) + '_monthly_' + clim_var.lower() + '_' + str(clim_month))
            plt.close(fig)




#--------------------------------------------------------------------------------------
#
#%% Making maps of correlations between climate variables and annual log WNV caseload
#
#--------------------------------------------------------------------------------------

hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
corr_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/annual_log_wnv_monthly_clim_corr.csv')
clim_vars = corr_df.columns
clim_vars = clim_vars.drop(labels=['hex', 'month'])

matplotlib.use('AGG')
for clim_month in np.arange(1, 13, 1):
    for clim_var in clim_vars:
        hex_grid = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
        test = corr_df.loc[corr_df['month']==clim_month][['hex', clim_var]]
        gdf = hex_grid.merge(test, how='left', left_on='grid_id', right_on='hex')

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        ax.set_aspect('equal')
        gdf.plot(column=clim_var, ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'correlation'}, vmin=-0.5, vmax=0.5, cmap='bwr')
        plt.xlim([-129, -62])
        plt.ylim([23, 52])
        plt.xticks(np.arange(-125, -55, 10), fontsize=14)
        plt.yticks(np.arange(30, 60, 10), fontsize=14)
        # plt.title('Annual WNND and ' + month_names[clim_month] + ' ' + clim_var, fontsize=18)
        plt.title('Annual WNND (log) and March Drought', fontsize=24)
        cb_ax = ax.figure.axes[1]
        cb_ax.tick_params(labelsize=14)
        # plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/compiled_correlations/annual_log_wnv_monthly_' + clim_var.lower() + '_' + str(clim_month))
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/compiled_correlations/annual_log_wnv_March_drought.svg')
        plt.close(fig)



#--------------------------------------------------------------------------------------
#
#%% Making maps of correlations between monthly log WNV caseload and monthly log WNV caseload
#
#--------------------------------------------------------------------------------------

hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
corr_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/monthly_log_wnv_monthly_log_wnv_corr.csv')
corr_df = corr_df.drop(labels='Unnamed: 0', axis=1)

matplotlib.use('AGG')
clim_var = 'wnv'  # not bothering with a refactoring here
for wnv_month in np.arange(1, 13, 1):
    for clim_month in np.arange(1, 13, 1):
        hex_grid = gpd.read_file(
            '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
        test = corr_df[(corr_df['WNV_month'] == wnv_month) & (corr_df['climate_var'] == clim_var) & (
                    corr_df['climate_month'] == clim_month)]
        test.loc[:, 'hex'] = test['hex'].astype(int).values
        test = test[['hex', 'correlation']]
        gdf = hex_grid.merge(test, how='right', left_on='grid_id', right_on='hex')

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        ax.set_aspect('equal')
        gdf.plot(column='correlation', ax=ax, legend=True, legend_kwds={'shrink': .7, 'label': 'correlation'},
                 vmin=-0.5, vmax=0.5, cmap='bwr')
        plt.xlim([-129, -62])
        plt.ylim([23, 52])
        plt.xticks(np.arange(-125, -55, 10), fontsize=14)
        plt.yticks(np.arange(30, 60, 10), fontsize=14)
        plt.title(month_names[wnv_month] + ' WNND and ' + month_names[clim_month] + ' ' + clim_var, fontsize=18)
        plt.savefig(
            '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/compiled_correlations/wnv_to_wnv/monthly_log_wnv_' + str(
                wnv_month) + '_monthly_log_' + clim_var.lower() + '_' + str(clim_month))
        plt.close(fig)




#--------------------------------------------------------------------------------------
#
#%% Making maps of correlations between climate variables and seasonal timing of WNV caseload
#
#--------------------------------------------------------------------------------------

hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
timing_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_season_center.csv')
timing_df = timing_df[timing_df['num_cases']>=5]
timing_std = timing_df.groupby('hex')['timing_center'].std()
timing_std = timing_std * 30.5  # roughly translating to days
timing_mean = timing_df.groupby('hex')['timing_center'].mean()
matplotlib.use('macosx')

# mean of seasonal center
matplotlib.use('AGG')
hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
gdf = hex_grid.merge(timing_mean, how='left', left_on='grid_id', right_on=timing_mean.index)
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column='timing_center', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'timing of seasonal center'},cmap='plasma_r')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
plt.title('Mean Timing of Seasonal Center of WNND', fontsize=18)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/mean_seasonal_center')
plt.close()

# standard deviation of seasonal center
hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
gdf = hex_grid.merge(timing_std, how='right', left_on='grid_id', right_on=timing_std.index)
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column='timing_center', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'standard deviation of seasonal center (days)'},cmap='RdPu')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
cb_ax = ax.figure.axes[1]
cb_ax.tick_params(labelsize=14)
plt.title('Standard Deviation of Seasonal Center of WNND', fontsize=22)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/std_seasonal_center.svg')

# calculating mean seasonality index
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv')  # for hex-grid analysis
wnv_data = wnv_data.rename(columns={'Unnamed: 0': 'date'})
wnv_data = wnv_data[wnv_data['date'] >= '2004-01-01']

seasonality_index = pd.Series(index=timing_df.hex.unique())
for hex in seasonality_index.index:
    hex = str(hex)
    hex_wnv = wnv_data[['date', hex]]
    mon_counts = hex_wnv.groupby(pd.to_datetime(hex_wnv.date).dt.month)[hex].sum()
    seasonality_index.loc[int(hex)] = np.sum(np.abs(mon_counts - mon_counts.mean())) / mon_counts.sum()

# plotting mean seasonality index
hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
seasonality_index = pd.DataFrame(seasonality_index)
seasonality_index = seasonality_index.rename(columns={0: 'seasonality'})
gdf = hex_grid.merge(seasonality_index, how='right', left_on='grid_id', right_on=seasonality_index.index)
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column='seasonality', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'seasonality index (unitless)'},cmap='RdPu')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
plt.title('Seasonality of WNND Burden', fontsize=18)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/test')






## trying to identify which hexes are in which climate divisions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

hex_grid = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
# clim_regions = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/CONUS_CLIMATE_DIVISIONS/GIS.OFFICIAL_CLIM_DIVISIONS.shp')
hex_grid['coords'] = hex_grid['geometry'].apply(lambda x: x.representative_point().coords[:])
hex_grid['coords'] = [coords[0] for coords in hex_grid['coords']]
# gdf = hex_grid.merge(seasonality_index, how='right', left_on='grid_id', right_on=seasonality_index.index)
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
# hex_grid.plot(column='seasonality', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'seasonality index (unitless)'},cmap='RdPu')
hex_grid.plot(column='grid_id', ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'seasonality index (unitless)'},cmap='RdPu')
for idx, row in hex_grid.iterrows():
    plt.annotate(text=str(row['grid_id']), xy=row['coords'],
                 horizontalalignment='center')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
plt.title('Seasonality of WNND Burden', fontsize=18)


# ok, trying to aggregate by climate regions for AGU example
# south_hex_list = ['84', '86', '85', '87', '89', '102', '104', '106', '108', '110', '105', '107', '124', '126', '128', '125', '127']
south_hex_list = ['108', '106', '105', '125', '126', '127']
west_central_hex_list = ['112', '114']
ohio_valley_hex_list = ['111', '131', '']

fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
# clim_vars = ['temp_mean', 'temp_max', 'temp_min', 'temp_dew', 'drought', 'precipitation', 'soil_moisture_50', 'snow_depth', 'VPD max']  # add in something with moisture?
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv', index_col=0)  # for hex-grid analysis
temp_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_tmean_198101-202211_gridded-2.csv', index_col=0)
precip_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_precip_198101-202211_gridded-2.csv')
dew_point_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_td_198101-202211_gridded-2.csv')
pdsi_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/gridMET_PDSI_199901-202212_gridded-2.csv')
vpd_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_vpdmax_198101-202303_gridded-2.csv')
soil_moisture_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')

south_wnv = wnv_data[south_hex_list].sum(axis=1)
south_wnv.index = pd.to_datetime(south_wnv.index)
south_wnv_ann = south_wnv.groupby(south_wnv.index.year).sum()
south_temp = temp_data[south_hex_list].mean(axis=1)
south_temp.index = pd.to_datetime(south_temp.index)
south_precip = precip_data[south_hex_list].mean(axis=1)
south_precip.index = pd.date_range(start='1/1/1981', end='11/2/2022', freq='MS')
south_dew_point = dew_point_data[south_hex_list].mean(axis=1)
south_dew_point.index = pd.date_range(start='1/1/1981', end='11/2/2022', freq='MS')
south_pdsi = pdsi_data[south_hex_list].mean(axis=1)
south_pdsi.index = pd.date_range(start='1/1/1999', end='12/2/2022', freq='MS')
south_vpd = vpd_data[south_hex_list].mean(axis=1)
south_vpd.index = pd.date_range(start='1/1/1981', end='3/1/2023', freq='MS')
south_soil_moisture = soil_moisture_data[south_hex_list].mean(axis=1)
south_soil_moisture.index = pd.date_range(start='1/1/1999', end='12/2/2022', freq='MS')

for clim_data in [temp_data, precip_data, dew_point_data, pdsi_data, vpd_data, soil_moisture_data]:
    clim_data = clim_data[south_hex_list].mean(axis=1)
    clim_data.index = pd.to_datetime(clim_data.index)
    for m in np.arange(1, 13):
        if m <= 9:
            corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), south_soil_moisture[(south_soil_moisture.index.year > 2004) & (south_soil_moisture.index.year < 2023) & (south_soil_moisture.index.month == m)])
        elif m > 9:
            corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), south_soil_moisture[(south_soil_moisture.index.year > 2003) & (south_soil_moisture.index.year < 2022) & (south_soil_moisture.index.month == m)])
        print(np.round(corr[0,1], 2))

clim_data = south_precip
spring_clim = clim_data[(clim_data.index.year > 2004) & (clim_data.index.year < 2023) & (clim_data.index.month.isin([3, 4, 5]))]
spring_clim = spring_clim.groupby(spring_clim.index.year).mean()
spring_corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), spring_clim)
#
summer_clim = clim_data[(clim_data.index.year > 2004) & (clim_data.index.year < 2023) & (clim_data.index.month.isin([6, 7, 8]))]
summer_clim = summer_clim.groupby(summer_clim.index.year).mean()
summer_corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), summer_clim)
#
winter_clim = clim_data[(clim_data.index.year > 2003) & (clim_data.index.year < 2023) & (clim_data.index.month.isin([12, 1, 2]))]
winter_clim_year = winter_clim.index.year.values
winter_clim_year[winter_clim.index.month == 12] = np.add(winter_clim_year[winter_clim.index.month == 12], 1)
winter_clim = winter_clim[winter_clim_year > 2004].groupby(winter_clim_year[winter_clim_year > 2004]).mean()
winter_clim = winter_clim[winter_clim.index < 2023]
winter_corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), winter_clim)
#
fall_clim = clim_data[(clim_data.index.year > 2003) & (clim_data.index.year < 2022) & (clim_data.index.month.isin([9, 10, 11]))]
fall_clim = fall_clim.groupby(fall_clim.index.year).mean()
fall_corr = np.corrcoef(np.log10(south_wnv_ann[south_wnv_ann.index > 2004]), fall_clim)
print(np.round(fall_corr[0, 1], 2))
print(np.round(winter_corr[0, 1], 2))
print(np.round(spring_corr[0, 1], 2))
print(np.round(summer_corr[0, 1], 2))

for i in winter_clim.values:
    print(i)






##########
#
# Assigning county-level WNV data to NEON domains
#
##########
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import Point
import time

# old list of counties (only included those with reported WNND)
# wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_county_time_series.csv', index_col=0)

# new complete list of counties
wnv_counties = pd.read_csv('/Users/ryan.harp/Documents/GitHub/WNV-forecast-data-2022/data-locations/locations.csv')
wnv_counties = wnv_counties.drop(columns='population')
wnv_counties['fips'] = wnv_counties['fips'].astype(str).str.rjust(5, '0')
counties_fips = wnv_counties['fips']
county_map = pd.DataFrame(index=counties_fips, columns=['NEON_domain'])

# loading NEON domains pulled from: https://www.neonscience.org/data-samples/data/spatial-data-maps
neon_domains = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NEON_Domains/NEON_Domains.shp')
neon_domains = neon_domains[neon_domains['DomainID'] <= 17]  # subsetting to contiguous U.S.
neon_domains['min_lon'] = neon_domains.bounds['minx']
neon_domains['max_lon'] = neon_domains.bounds['maxx']
neon_domains['min_lat'] = neon_domains.bounds['miny']
neon_domains['max_lat'] = neon_domains.bounds['maxy']


# loading county population centers pulled from: https://www.census.gov/geographies/reference-files/time-series/geo/centers-population.html
county_pop_centers = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_centers_of_pop.txt', sep=',', header=0, converters={'STATEFP': str, 'COUNTYFP': str})
county_pop_centers['FIPS'] = county_pop_centers['STATEFP'] + county_pop_centers['COUNTYFP']

def find_which_shapefile(lon, lat, shapefile_list):
    pt = Point([lon, lat])
    for i in np.arange(len(shapefile_list)):
        if pt.within(shapefile_list[i]):
            break
    return i


c = 0
start_time = time.time()
for fip in counties_fips:
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(fip) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = county_pop_centers.loc[county_pop_centers['FIPS'] == fip, 'LONGITUDE'].values[0]
    lat = county_pop_centers.loc[county_pop_centers['FIPS'] == fip, 'LATITUDE'].values[0]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        county_map.loc[fip] = neon_domains[possible_domains]['DomainID'].values
    elif np.sum(possible_domains) == 0:
        if fip == '12087':
            county_map.loc[fip] = 4
        print('Failed: ' + str(fip))
    else:
        domain_list = neon_domains[possible_domains]['DomainID'].values
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        county_map.loc[fip] = domain_list[shapefile_ind]

wnv_counties = wnv_counties.merge(county_map, how='left', left_on='fips', right_on = county_map.index)
wnv_counties.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv', index_label='FIPS')



##########
#
# Creating PRISM, gridMET, NLDAS NEON domains grid mask
#
##########
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import Point
import time
import xarray as xr


# loading NEON domains pulled from: https://www.neonscience.org/data-samples/data/spatial-data-maps
neon_domains = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NEON_Domains/NEON_Domains.shp')
neon_domains = neon_domains[neon_domains['DomainID'] <= 17]  # subsetting to contiguous U.S.
neon_domains['min_lon'] = neon_domains.bounds['minx']
neon_domains['max_lon'] = neon_domains.bounds['maxx']
neon_domains['min_lat'] = neon_domains.bounds['miny']
neon_domains['max_lat'] = neon_domains.bounds['maxy']

neon_domains = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/cb_2022_us_county_500k/cb_2022_us_county_500k.shp')
# neon_domains = neon_domains[neon_domains['DomainID'] <= 17]  # subsetting to contiguous U.S.
neon_domains['min_lon'] = neon_domains.bounds['minx']
neon_domains['max_lon'] = neon_domains.bounds['maxx']
neon_domains['min_lat'] = neon_domains.bounds['miny']
neon_domains['max_lat'] = neon_domains.bounds['maxy']
neon_domains['FIPS'] = neon_domains['STATEFP'] + neon_domains['COUNTYFP']


# opening PRISM file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_tmean_198101-202211.nc')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

def find_which_shapefile(lon, lat, shapefile_list):
    pt = Point([lon, lat])
    for i in np.arange(len(shapefile_list)):
        if pt.within(shapefile_list[i]):
            break
    return i

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties.npy', grid_mask)


# opening gridMET file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/gridMET_PDSI_199901-202212.nc')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/gridMET_to_counties.npy', grid_mask)



# opening NLDAS file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/snow_depth/NLDAS_NOAH0125_M.A199901.002.grb.SUB.nc4')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties.npy', grid_mask)





## aggregating WNV counts onto NEON domains
import pandas as pd
import numpy as np

# loading data and conversion file
wnv_county_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_county_time_series.csv', index_col=0)
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')

# creating NEON-based dataframe
neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
wnv_neon_data = pd.DataFrame(0, index=wnv_county_data.index, columns=neon_domains)
for i, row in county_to_neon.iterrows():
    wnv_neon_data[row['NEON_domain']] = wnv_neon_data[row['NEON_domain']] + wnv_county_data[str(row['FIPS'])]

wnv_neon_data.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_neon_monthly_time_series.csv')






# creating county-averages of climatic conditions
import pandas as pd
import numpy as np
import xarray as xr
import time

# starting with PRISM
prism_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties.npy')
clim_var = 'vpdmin'
prism_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_' + clim_var + '_198101-202303.nc')
prism_data = prism_data.sel(time=slice('1999-01-01', '2022-11-01'))

# getting mean of all points in a given grid cell and transforming to a dataframe for saving
county_array = np.zeros([len(prism_data.time), len(np.unique(prism_to_county_grid)[:-1])])  # [:-1] to remove nan
county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
start_time = time.time()
for i, fip in np.ndenumerate(np.unique(prism_to_county_grid)[:-1]):
    if i[0] % 10 == 0:
        print('FIP Index: ' + str(i[0]) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    county_array[:, i[0]] = np.squeeze(prism_data[list(prism_data.keys())[0]].where(np.transpose(prism_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

fips = np.unique(prism_to_county_grid)[:-1]
ds_processed = pd.DataFrame(data = county_array, index = prism_data.time.values, columns = fips)
ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_' + clim_var + '_199901-202211.csv')




# starting with PRISM
prism_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties.npy')
clim_vars = ['tmin', 'tmax', 'td']
for clim_var in clim_vars:
    prism_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_' + clim_var + '_198101-202211.nc')
    prism_data = prism_data.sel(time=slice('1999-01-01', '2022-11-01'))

    # getting mean of all points in a given grid cell and transforming to a dataframe for saving
    county_array = np.zeros([len(prism_data.time), len(np.unique(prism_to_county_grid)[:-1])])  # [:-1] to remove nan
    county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
    # grid_array = np.zeros([len(ds.day), np.shape(grid)[0]]); grid_array[:,:] = np.nan
    start_time = time.time()
    for i, fip in np.ndenumerate(np.unique(prism_to_county_grid)[:-1]):
        if i[0] % 10 == 0:
            print('FIP Index: ' + str(i[0]) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
        county_array[:, i[0]] = np.squeeze(prism_data[list(prism_data.keys())[0]].where(np.transpose(prism_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

    fips = np.unique(prism_to_county_grid)[:-1]
    ds_processed = pd.DataFrame(data = county_array, index = prism_data.time.values, columns = fips)
    ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_' + clim_var + '_199901-202211.csv')





# gridMET for PDSI
gridMET_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/gridMET_to_counties.npy')
gridMET_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/gridMET_PDSI_199901-202212.nc')
gridMET_data = gridMET_data.sel(day=slice('1999-01-01', '2022-11-01'))

# getting mean of all points in a given grid cell and transforming to a dataframe for saving
county_array = np.zeros([len(gridMET_data.day), len(np.unique(gridMET_to_county_grid)[:-1])])  # [:-1] to remove nan
county_array[:,:] = np.nan
start_time = time.time()
for i, fip in np.ndenumerate(np.unique(gridMET_to_county_grid)[:-1]):
    if (i[0] % 10 == 0) & (i[0] > 0):
        current_time = time.time()
        print('Index: ' + str(i[0]) + ' at ' + str(np.round((current_time - start_time)/60, 1)) + ' minutes.')
        print('Projected remaining time: ' + str(np.round(((current_time - start_time)*3105/i[0]-(current_time - start_time))/60, 1)) + ' minutes.')
    county_array[:, i[0]] = np.squeeze(gridMET_data[list(gridMET_data.keys())[1]].where(np.transpose(gridMET_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

fips = np.unique(gridMET_to_county_grid)[:-1]
ds_processed = pd.DataFrame(data = county_array, index = gridMET_data.day.values, columns = fips)
ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/gridMET_PDSI_199901-202211.csv')





# finishing the set with NLDAS
NLDAS_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties.npy')
clim_var = 'snow_melt'
NLDAS_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/snow_melt/*.nc4')
NLDAS_data = NLDAS_data.sel(time=slice('1999-01-01', '2022-11-01'))

county_array = np.zeros([len(NLDAS_data.time), len(np.unique(NLDAS_to_county_grid)[:-1])])  # [:-1] to remove nan
county_array[:,:] = np.nan
start_time = time.time()
for i, fip in np.ndenumerate(np.unique(NLDAS_to_county_grid)[:-1]):
    if (i[0] % 10 == 0) & (i[0] > 0):
        current_time = time.time()
        print('Index: ' + str(i[0]) + ' at ' + str(np.round((current_time - start_time)/60, 1)) + ' minutes.')
        print('Projected remaining time: ' + str(np.round(((current_time - start_time)*3105/i[0]-(current_time - start_time))/60, 1)) + ' minutes.')
    county_array[:, i[0]] = np.squeeze(NLDAS_data[list(NLDAS_data.keys())[0]].where(np.transpose(NLDAS_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

fips = np.unique(NLDAS_to_county_grid)[:-1]
ds_processed = pd.DataFrame(data=county_array, index=NLDAS_data.time.values, columns=fips)
ds_processed.to_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/NLDAS_' + clim_var + '_199901-202211.csv')





# finishing the set with NLDAS multi-depth variables
NLDAS_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties.npy')
clim_var = 'soil_temp'
NLDAS_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/soil_temp/*.nc4')
NLDAS_data = NLDAS_data.sel(time=slice('1999-01-01', '2022-11-01'))

fips = np.unique(NLDAS_to_county_grid)[:-1]
for d in NLDAS_data.depth.values:
    county_array = np.zeros([len(NLDAS_data.time), len(np.unique(NLDAS_to_county_grid)[:-1])])  # [:-1] to remove nan
    county_array[:, :] = np.nan
    start_time = time.time()
    if ~np.isin(d, [70]):  # limiting to depths of interest (5, 50, 100)
        continue
    ds_level = NLDAS_data.sel(depth=d)
    for i, fip in np.ndenumerate(np.unique(NLDAS_to_county_grid)[:-1]):
        if (i[0] % 10 == 0) & (i[0] > 0):
            current_time = time.time()
            print('Index: ' + str(i[0]) + ' at ' + str(np.round((current_time - start_time)/60, 1)) + ' minutes.')
            print('Projected remaining time: ' + str(np.round(((current_time - start_time)*3105/i[0]-(current_time - start_time))/60, 1)) + ' minutes.')
        county_array[:, i[0]] = np.squeeze(ds_level[list(ds_level.keys())[1]].where(np.transpose(NLDAS_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())
    ds_processed = pd.DataFrame(data=county_array, index=ds_level.time.values, columns = fips)
    ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/NLDAS_' + clim_var + '_' + str(int(d)) + '_199901-202212.csv')




## rerunning the grid masks and aggregation over Connecticut because they changed their admin districts...
import geopandas as gpd
import xarray as xr
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.geometry import Point
import time


neon_domains = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/cb_2018_us_county_500k/cb_2018_us_county_500k.shp')
# neon_domains = neon_domains[neon_domains['DomainID'] <= 17]  # subsetting to contiguous U.S.
neon_domains['min_lon'] = neon_domains.bounds['minx']
neon_domains['max_lon'] = neon_domains.bounds['maxx']
neon_domains['min_lat'] = neon_domains.bounds['miny']
neon_domains['max_lat'] = neon_domains.bounds['maxy']
neon_domains['FIPS'] = neon_domains['STATEFP'] + neon_domains['COUNTYFP']


# opening PRISM file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_tmean_198101-202211.nc')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

def find_which_shapefile(lon, lat, shapefile_list):
    pt = Point([lon, lat])
    for i in np.arange(len(shapefile_list)):
        if pt.within(shapefile_list[i]):
            break
    return i

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties_CT.npy', grid_mask)


# opening gridMET file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/gridMET_PDSI_199901-202212.nc')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/gridMET_to_counties_CT.npy', grid_mask)



# opening NLDAS file
ds = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/snow_depth/NLDAS_NOAH0125_M.A199901.002.grb.SUB.nc4')
ds_lat = ds.lat.values
ds_lon = ds.lon.values

grid_mask = np.zeros([len(ds_lon), len(ds_lat)]); grid_mask[:] = np.nan

c = 0
start_time = time.time()
for ind, dummy in np.ndenumerate(grid_mask):
    c += 1
    if c % 100 == 0:
        print('Index: ' + str(ind) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lon = ds_lon[ind[0]]
    lat = ds_lat[ind[1]]
    possible_domains = (lat > neon_domains['min_lat']) & (lat < neon_domains['max_lat']) & (lon > neon_domains['min_lon']) & (lon < neon_domains['max_lon'])
    if np.sum(possible_domains) == 1:
        grid_mask[ind] = neon_domains['FIPS'][possible_domains]
        # grid_mask[ind] = neon_domains['DomainID'][possible_domains]
    elif np.sum(possible_domains) == 0:
        continue
        # print('Failed: ' + str(ind))
    else:
        domain_list = neon_domains['FIPS'][possible_domains]
        # domain_list = neon_domains['DomainID'][possible_domains]
        domain_shapefiles = list()
        for domain_ind in domain_list:
            domain_shapefiles.append(neon_domains[neon_domains['FIPS'] == domain_ind].unary_union)
            # domain_shapefiles.append(neon_domains[neon_domains['DomainID'] == domain_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, domain_shapefiles)
        grid_mask[ind] = domain_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time, 1)) + ' seconds.')

np.save('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties_CT.npy', grid_mask)


# starting with PRISM
prism_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties_CT.npy')
clim_var = 'tmean'
prism_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_' + clim_var + '_198101-202211.nc')
prism_data = prism_data.sel(time=slice('1999-01-01', '2022-11-01'))

# getting mean of all points in a given grid cell and transforming to a dataframe for saving
CT_fips = [9001, 9003, 9005, 9007, 9009, 9011, 9013, 9015]
county_array = np.zeros([len(prism_data.time), len(CT_fips)])  # [:-1] to remove nan
county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
start_time = time.time()
for i, fip in np.ndenumerate(CT_fips):
    if i[0] % 10 == 0:
        print('FIP Index: ' + str(i[0]) + ' at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    county_array[:, i[0]] = np.squeeze(prism_data[list(prism_data.keys())[0]].where(np.transpose(prism_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

fips = np.unique(prism_to_county_grid)[:-1]
ds_processed = pd.DataFrame(data = county_array, index = prism_data.time.values, columns = CT_fips)
ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_' + clim_var + '_199901-202211_CT.csv')


# starting with PRISM
prism_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/PRISM_to_counties_CT.npy')
clim_vars = ['vpdmin', 'vpdmax']
for clim_var in clim_vars:
    prism_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_' + clim_var + '_198101-202303.nc')
    prism_data = prism_data.sel(time=slice('1999-01-01', '2022-11-01'))

    # getting mean of all points in a given grid cell and transforming to a dataframe for saving
    CT_fips = [9001, 9003, 9005, 9007, 9009, 9011, 9013, 9015]
    county_array = np.zeros([len(prism_data.time), len(CT_fips)])  # [:-1] to remove nan
    county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
    start_time = time.time()
    for i, fip in np.ndenumerate(CT_fips):
        county_array[:, i[0]] = np.squeeze(prism_data[list(prism_data.keys())[0]].where(np.transpose(prism_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

    ds_processed = pd.DataFrame(data = county_array, index = prism_data.time.values, columns = CT_fips)
    ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_' + clim_var + '_199901-202211_CT.csv')


# gridMET for PDSI
gridMET_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/gridMET_to_counties_CT.npy')
gridMET_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/gridMET_PDSI_199901-202212.nc')
gridMET_data = gridMET_data.sel(day=slice('1999-01-01', '2022-11-01'))

# getting mean of all points in a given grid cell and transforming to a dataframe for saving
CT_fips = [9001, 9003, 9005, 9007, 9009, 9011, 9013, 9015]
county_array = np.zeros([len(gridMET_data.day), len(CT_fips)])  # [:-1] to remove nan
county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
start_time = time.time()
for i, fip in np.ndenumerate(CT_fips):
    county_array[:, i[0]] = np.squeeze(gridMET_data[list(gridMET_data.keys())[1]].where(np.transpose(gridMET_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

ds_processed = pd.DataFrame(data = county_array, index = gridMET_data.day.values, columns = CT_fips)
ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/gridMET_PDSI_199901-202211_CT.csv')


# finishing the set with NLDAS
NLDAS_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties_CT.npy')
clim_var = 'water_equivalent_snow_depth'
NLDAS_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/snow_depth/*.nc4')
NLDAS_data = NLDAS_data.sel(time=slice('1999-01-01', '2022-11-01'))

CT_fips = [9001, 9003, 9005, 9007, 9009, 9011, 9013, 9015]
county_array = np.zeros([len(NLDAS_data.time), len(CT_fips)])  # [:-1] to remove nan
county_array[:,:] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
start_time = time.time()
for i, fip in np.ndenumerate(CT_fips):
    county_array[:, i[0]] = np.squeeze(NLDAS_data[list(NLDAS_data.keys())[0]].where(np.transpose(NLDAS_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())

ds_processed = pd.DataFrame(data=county_array, index=NLDAS_data.time.values, columns=CT_fips)
ds_processed.to_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/NLDAS_' + clim_var + '_199901-202211_CT.csv')


# finishing the set with NLDAS multi-depth variables
NLDAS_to_county_grid = np.load('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NLDAS_to_counties_CT.npy')
clim_var = 'soil_moisture'
NLDAS_data = xr.open_mfdataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/NLDAS/soil_moisture/*.nc4')
NLDAS_data = NLDAS_data.sel(time=slice('1999-01-01', '2022-11-01'))

for d in NLDAS_data.depth.values:
    CT_fips = [9001, 9003, 9005, 9007, 9009, 9011, 9013, 9015]
    county_array = np.zeros([len(NLDAS_data.time), len(CT_fips)])  # [:-1] to remove nan
    county_array[:, :] = np.nan  # np.shape(grid)[1] for PRISM, [0] for gridMET, NLDAS
    start_time = time.time()
    if ~np.isin(d, [5, 50, 100]):  # limiting to depths of interest (5, 50, 100)
        continue
    ds_level = NLDAS_data.sel(depth=d)
    for i, fip in np.ndenumerate(CT_fips):
        county_array[:, i[0]] = np.squeeze(ds_level[list(ds_level.keys())[1]].where(np.transpose(NLDAS_to_county_grid) == fip).mean(dim=['lat', 'lon'], skipna=True).compute())
    ds_processed = pd.DataFrame(data=county_array, index=ds_level.time.values, columns = CT_fips)
    ds_processed.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/NLDAS_' + clim_var + '_' + str(int(d)) + '_199901-202212_CT.csv')





## creating county caseload-weighted averages of climatic conditions for NEON domains
import pandas as pd
import numpy as np
from glob import glob

clim_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_tmean_199901-202211.csv', index_col = 0)
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv', index_col = 0)
missing_clim_counties = np.setdiff1d(wnv_data.columns, clim_data.columns)
wnv_data = wnv_data.drop(columns=missing_clim_counties)
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')

neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
weighted_clim = np.zeros([len(neon_domains), np.shape(clim_data)[0]])
for domain in neon_domains:
    neon_fips = county_to_neon.loc[county_to_neon['NEON_domain'] == domain]['FIPS'].values.astype(str)
    neon_fips = np.delete(neon_fips, np.where(neon_fips == missing_clim_counties))
    wnv_sub = wnv_data[neon_fips]
    clim_sub = clim_data[wnv_sub.columns]  # only taking climate data from counties with historical WNV caseload
    weighted_clim[domain-1, :] = np.sum(clim_sub * np.log10(wnv_sub.sum()), axis=1) / np.log10(wnv_sub.sum()).sum()
weighted_clim_data = pd.DataFrame(data=np.transpose(weighted_clim), index=clim_data.index, columns=neon_domains)


# building NEON domain WNV case counts
wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv', index_col = 0)
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')

neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
wnv_data_neon = np.zeros([len(neon_domains), np.shape(wnv_data)[0]])
for domain in neon_domains:
    neon_fips = county_to_neon.loc[county_to_neon['NEON_domain'] == domain]['FIPS'].values.astype(str)
    wnv_data_neon[domain-1, :] = wnv_data[neon_fips].sum(axis=1)
wnv_data_neon_df = pd.DataFrame(data=np.transpose(wnv_data_neon), index=wnv_data.index, columns=neon_domains)
wnv_data_neon_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv')


# looking at a quick annual comparison
wnv_data_ann = wnv_data_neon_df.groupby(pd.to_datetime(wnv_data_neon_df.index).year).sum()
domain_wnv = wnv_data_ann[9][wnv_data_ann.index >=2005]
weighted_clim_data.index = pd.to_datetime(weighted_clim_data.index)
clim_data = weighted_clim_data[9][(weighted_clim_data.index.year >= 2005) & (weighted_clim_data.index.month == 7)]
print(np.corrcoef(np.log10(domain_wnv), clim_data))





path = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/*1.csv'
all_files = sorted(glob(path, recursive=True))
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')
neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
for file in all_files:
    name = file.split('/')[-1].split('1999')[0]
    clim_data = pd.read_csv(file, index_col=0)
    wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv', index_col=0)
    missing_clim_counties = np.setdiff1d(wnv_data.columns, clim_data.columns)
    wnv_data = wnv_data.drop(columns=missing_clim_counties)

    weighted_clim = np.zeros([len(neon_domains), np.shape(clim_data)[0]])
    weighted_clim[:] = np.nan
    for domain in neon_domains:
        neon_fips = county_to_neon.loc[county_to_neon['NEON_domain'] == domain]['FIPS'].values.astype(str)
        neon_fips = np.delete(neon_fips, np.where(np.isin(neon_fips, missing_clim_counties)))
        wnv_sub = wnv_data[neon_fips]
        clim_sub = clim_data[wnv_sub.columns]  # only taking climate data from counties with historical WNV caseload
        weighted_clim[domain - 1, :] = np.sum(clim_sub * np.log10(wnv_sub.sum()), axis=1) / np.log10(wnv_sub.sum()).sum()
    weighted_clim_data = pd.DataFrame(data=np.transpose(weighted_clim), index=clim_data.index, columns=neon_domains)
    weighted_clim_data.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/' + name + '199901-202211.csv')





#%% Creating hex- and county-level datasets of climate factors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')
import geopandas as gpd
import xarray as xr
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_clim_data(fname):
    clim_df = pd.read_csv(fname, index_col=0)
    clim_df.index = pd.to_datetime(clim_df.index)
    clim_df = clim_df[clim_df.index >= '2003-01-01']
    return clim_df

# prepping results df
wnv_columns = ['climate_month', 'climate_variable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']  # not sure why I couldn't pull out the numbers programatically but whatever
wnv_clim_mon_corrs = pd.DataFrame(columns = wnv_columns)

# prepping WNND data
wnv_data_mon = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv',
    index_col=0)
wnv_data_mon.index = pd.to_datetime(wnv_data_mon.index)
wnv_data_mon = wnv_data_mon[wnv_data_mon.index.year >= 2005]
wnv_data_ann = wnv_data_mon.groupby(wnv_data_mon.index.year).sum()
wnv_data_log = np.log10(wnv_data_ann + 1)
wnv_data_ln = np.log(wnv_data_ann + 1)

# calculating seasonal center
wnv_seasonal_center = pd.DataFrame().reindex_like(wnv_data_log)
wnv_seasonal_center = wnv_data_mon.mul(wnv_data_mon.index.month, axis=0).groupby(wnv_data_mon.index.year).sum()/wnv_data_ann
wnv_seasonal_center = wnv_seasonal_center.round(3)
wnv_seasonal_center.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_seasonal_center.csv')

# plotting standard deviation of seasonal center
wnv_seasonal_center_std = wnv_seasonal_center.std()*30.5  # loosely translated to days

gdf_neon = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NEON_domains/NEON_domains.shp')
gdf_neon = gdf_neon[gdf_neon['DomainID'] <= 17]
gdf_neon = gdf_neon[gdf_neon['OBJECTID'] != 0]

gdf = gdf_neon.merge(pd.DataFrame(wnv_seasonal_center_std), how='left', left_on='DomainID', right_on=wnv_seasonal_center_std.index.astype(int))

fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column=0, ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'standard deviation of seasonal center (~days)'},cmap='RdPu')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
plt.title('Standard Deviation of Seasonal Center of WNND Burden', fontsize=18)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/standard_deviation_of_seasonal_center')
plt.close(fig)

# calculating seasonality index
seasonality_index = pd.Series(index=wnv_seasonal_center.columns)
for domain in seasonality_index.index:
    mon_counts = wnv_data_mon.groupby(wnv_data_mon.index.month).sum()
    seasonality_index = np.sum(np.abs(mon_counts - mon_counts.mean())) / mon_counts.sum()

gdf_neon = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NEON_domains/NEON_domains.shp')
gdf_neon = gdf_neon[gdf_neon['DomainID'] <= 17]
gdf_neon = gdf_neon[gdf_neon['OBJECTID'] != 0]

gdf = gdf_neon.merge(pd.DataFrame(seasonality_index), how='left', left_on='DomainID', right_on=seasonality_index.index.astype(int))

import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 9
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
gdf.plot(column=0, ax=ax, legend=True, legend_kwds={'shrink':.7, 'label': 'seasonality index (unitless)'},cmap='RdPu', edgecolor='k', linewidth=1)
# gdf['mask'] = np.array(
#     [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0])  # masking out domains with fewer than 1000 historical cases
# gdf.loc[gdf['mask'] == 1].plot(ax=ax, color='none', edgecolor='lightgrey', hatch='/')
# gdf.loc[gdf['mask'] == 1].plot(ax=ax, color='none', edgecolor='k')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)
plt.title('Seasonality of WNND Burden', fontsize=18)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/seasonality_index.svg')
plt.close(fig)

# running over each climate variable and month
clim_vars = ['tmean', 'tmax', 'tmin', 'td', 'vpdmax', 'vpdmin', 'precip', 'PDSI', 'snow_depth', 'snow_melt', 'soil_moisture_5',
             'soil_moisture_50', 'soil_moisture_100', 'soil_temp_5', 'soil_temp_25', 'soil_temp_70', 'soil_temp_150', 'water_equivalent_snow_depth']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
for clim_var in clim_vars:
    if np.isin(clim_var, clim_vars[:7]):  # prism variables
        clim_data = load_clim_data(fpath + 'prism_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[7]):  # gridMET variable
        clim_data = load_clim_data(fpath + 'gridMET_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[8:]):  # NLDAS variables
        clim_data = load_clim_data(fpath + 'NLDAS_' + clim_var + '_199901-202211.csv')
    clim_data.index = pd.to_datetime(clim_data.index)
    clim_data = clim_data[clim_data.index.year >= 2004]

    for mon in np.arange(1, 13, 1):
        if mon <= 9:
            mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year > 2004)]
        elif mon >= 10:
            mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year < 2022)]
        mon_corrs = wnv_data_log.reset_index().drop(columns = 'index').corrwith(mon_clim_data.reset_index().drop(columns = 'index'))
        wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index), 'climate_month'] = mon
        wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index)-1, 'climate_variable'] = clim_var
        wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index)-1, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']] = np.round(mon_corrs, 2)

wnv_clim_mon_corrs.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_log_annual_wnv_corrs.csv')


# making monthly maps
wnv_clim_mon_corrs = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_log_annual_wnv_corrs.csv')
wnv_clim_mon_corrs = wnv_clim_mon_corrs.drop(columns='Unnamed: 0')
month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
clim_var_name = {'tmean': 'Mean Temp', 'tmax': 'Max Temp', 'tmin': 'Min Temp', 'td': 'Dew Point',
                 'vpdmax': 'Max Vapor Pressure Deficit', 'vpdmin': 'Min Vapor Pressure Deficit',
                 'precip': 'Precipitation', 'PDSI': 'PDSI', 'snow_depth': 'Snow Depth', 'snow_melt': 'Snow Melt',
                 'soil_moisture_5': '5 mm Soil Moisture', 'soil_moisture_50': '50 mm Soil Moisture',
                 'soil_moisture_100': '100 mm Soil Moisture', 'soil_temp_5': '5 mm Soil Temp',
                 'soil_temp_25': '25 mm Soil Temp', 'soil_temp_70': '70 mm Soil Temp', 'soil_temp_150': '150 mm Soil Temp',
                 'water_equivalent_snow_depth': 'Water Equivalent Snow Depth'}

matplotlib.use('AGG')
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

# prepping NEON domain shapefiles
gdf_neon = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/NEON_domains/NEON_domains.shp')
gdf_neon = gdf_neon[gdf_neon['DomainID'] <= 17]
gdf_neon = gdf_neon[gdf_neon['OBJECTID'] != 0]

# looping over each month/climate variable
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 9
for mon in np.arange(1, 13):
    for clim_var in wnv_clim_mon_corrs['climate_variable'].unique():
        mon_clim_corrs = wnv_clim_mon_corrs.loc[(wnv_clim_mon_corrs['climate_month']==mon) & (wnv_clim_mon_corrs['climate_variable']==clim_var)]
        mon_clim_corrs = mon_clim_corrs.drop(columns=['climate_month', 'climate_variable']).transpose().astype(float)
        gdf = gdf_neon.merge(mon_clim_corrs, how='left', left_on='DomainID', right_on=mon_clim_corrs.index.astype(int))

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='grey', linewidth=0.5)
        ax.set_aspect('equal')
        gdf.plot(column=gdf.columns[-1], ax=ax, vmin=-0.47, vmax=0.47, legend=True, cmap='bwr', edgecolor='k', linewidth=1)
        gdf['mask'] = np.array([0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0])  # masking out domains with fewer than 1000 historical cases
        gdf.loc[gdf['mask']==1].plot(ax=ax, color='none', edgecolor='lightgrey', hatch='/')
        gdf.loc[gdf['mask']==1].plot(ax=ax, color='none', edgecolor='k')
        plt.xlim([-129, -62])
        plt.ylim([23, 52])
        plt.xticks(np.arange(-125, -55, 10), fontsize=14)
        plt.yticks(np.arange(30, 60, 10), fontsize=14)
        # plt.title('WNND (log) and ' + month_names[mon] + ' ' + clim_var_name[clim_var], fontsize=18)
        plt.title('Annual WNND (log) and June Temperature', fontsize = 18)
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/correlations/annual_log_wnv_monthly_' + clim_var.lower() + '_' + str(mon) + '.svg')
        plt.close(fig)

# running over each climate variable and month
wnv_columns = ['climate_month', 'climate_variable', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']  # not sure why I couldn't pull out the numbers programatically but whatever
seasonal_center_mon_corrs = pd.DataFrame(columns = wnv_columns)
clim_vars = ['tmean', 'tmax', 'tmin', 'td', 'vpdmax', 'vpdmin', 'precip', 'PDSI', 'snow_depth', 'snow_melt', 'soil_moisture_5',
             'soil_moisture_50', 'soil_moisture_100', 'soil_temp_5', 'soil_temp_25', 'soil_temp_70', 'soil_temp_150', 'water_equivalent_snow_depth']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
for clim_var in clim_vars:
    if np.isin(clim_var, clim_vars[:7]):  # prism variables
        clim_data = load_clim_data(fpath + 'prism_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[7]):  # gridMET variable
        clim_data = load_clim_data(fpath + 'gridMET_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[8:]):  # NLDAS variables
        clim_data = load_clim_data(fpath + 'NLDAS_' + clim_var + '_199901-202211.csv')
    clim_data.index = pd.to_datetime(clim_data.index)
    clim_data = clim_data[clim_data.index.year >= 2004]

    for mon in np.arange(1, 13, 1):
        if mon <= 9:
            mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year > 2004)]
        elif mon >= 10:
            mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year < 2022)]
        mon_corrs = wnv_seasonal_center.reset_index().drop(columns = 'index').corrwith(mon_clim_data.reset_index().drop(columns = 'index'))
        seasonal_center_mon_corrs.loc[len(seasonal_center_mon_corrs.index), 'climate_month'] = mon
        seasonal_center_mon_corrs.loc[len(seasonal_center_mon_corrs.index)-1, 'climate_variable'] = clim_var
        seasonal_center_mon_corrs.loc[len(seasonal_center_mon_corrs.index)-1, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']] = np.round(mon_corrs, 2)

seasonal_center_mon_corrs.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_wnv_seasonal_center_corrs.csv')


# making heat maps of correlations
wnv_clim_mon_corrs = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_log_annual_wnv_corrs.csv', index_col=0)
month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
clim_var_name = {'tmean': 'Mean Temp', 'tmax': 'Max Temp', 'tmin': 'Min Temp', 'td': 'Dew Point',
                 'vpdmax': 'Max Vapor Pressure Deficit', 'vpdmin': 'Min Vapor Pressure Deficit',
                 'precip': 'Precipitation', 'PDSI': 'PDSI', 'snow_depth': 'Snow Depth', 'snow_melt': 'Snow Melt',
                 'soil_moisture_5': '5 mm Soil Moisture', 'soil_moisture_50': '50 mm Soil Moisture',
                 'soil_moisture_100': '100 mm Soil Moisture', 'soil_temp_5': '5 mm Soil Temp',
                 'soil_temp_25': '25 mm Soil Temp', 'soil_temp_70': '70 mm Soil Temp', 'soil_temp_150': '150 mm Soil Temp',
                 'water_equivalent_snow_depth': 'Water Equivalent Snow Depth'}
domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                '16': 'Pacific Northwest', '17': 'Pacific Southwest'}

matplotlib.use('AGG')
for domain in np.arange(1, 18):
    domain_corrs = wnv_clim_mon_corrs[['climate_month', 'climate_variable', str(domain)]]  # pulling out domain-specific values
    domain_corrs = domain_corrs.loc[domain_corrs['climate_variable'].isin(['tmean', 'td', 'vpdmax', 'precip', 'PDSI', 'soil_moisture_50'])]
    domain_corrs = domain_corrs.reset_index().pivot(columns='climate_variable',index='climate_month',values=str(domain))
    fig, ax = plt.subplots(1, 1, figsize=(6, 8), layout='tight')
    sns.heatmap(domain_corrs, ax=ax, annot=True, cmap='bwr_r', square=True, mask = abs(domain_corrs) < 0.32, cbar_kws={'shrink': 0.6},
                annot_kws={'fontsize':11}, vmax=0.601, vmin=-0.600)
    plt.xlabel('')
    plt.ylabel('')
    plt.title(domain_names[str(domain)] + ' Correlation Heat Map', fontsize=20)
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/correlations/test_annual_log_wnv_heat_map_domain_' + str(domain) + '.svg')
    plt.close(fig)






##########
#
# Calculating Correlations with Monthly WNV/Climate
#
##########

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
clim_var_name = {'tmean': 'Mean Temp', 'tmax': 'Max Temp', 'tmin': 'Min Temp', 'td': 'Dew Point',
                 'vpdmax': 'Max Vapor Pressure Deficit', 'vpdmin': 'Min Vapor Pressure Deficit',
                 'precip': 'Precipitation', 'PDSI': 'PDSI', 'snow_depth': 'Snow Depth', 'snow_melt': 'Snow Melt',
                 'soil_moisture_5': '5 mm Soil Moisture', 'soil_moisture_50': '50 mm Soil Moisture',
                 'soil_moisture_100': '100 mm Soil Moisture', 'soil_temp_5': '5 mm Soil Temp',
                 'soil_temp_25': '25 mm Soil Temp', 'soil_temp_70': '70 mm Soil Temp', 'soil_temp_150': '150 mm Soil Temp',
                 'water_equivalent_snow_depth': 'Water Equivalent Snow Depth'}
domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                '16': 'Pacific Northwest', '17': 'Pacific Southwest'}

import warnings
warnings.filterwarnings('ignore', category=UserWarning)  # tread lightly but enacting this to ignore palette length warnings

# prepping WNND data
wnv_data_mon = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv',
    index_col=0)
wnv_data_mon.index = pd.to_datetime(wnv_data_mon.index)
wnv_data_mon = wnv_data_mon[wnv_data_mon.index.year >= 2005]
wnv_data_mon_log = np.log10(wnv_data_mon + 1)

mon_corr_df = pd.DataFrame(columns = ['wnv_month', 'climate_month', 'climate_variable',
                                      '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'])
clim_vars = ['tmean', 'tmax', 'tmin', 'td', 'vpdmax', 'vpdmin', 'precip', 'PDSI', 'snow_depth', 'snow_melt', 'soil_moisture_5',
             'soil_moisture_50', 'soil_moisture_100', 'soil_temp_5', 'soil_temp_25', 'soil_temp_70', 'soil_temp_150', 'water_equivalent_snow_depth']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
for clim_var in clim_vars:
    if np.isin(clim_var, clim_vars[:7]):  # prism variables
        clim_data = load_clim_data(fpath + 'prism_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[7]):  # gridMET variable
        clim_data = load_clim_data(fpath + 'gridMET_' + clim_var + '_199901-202211.csv')
    elif np.isin(clim_var, clim_vars[8:]):  # NLDAS variables
        clim_data = load_clim_data(fpath + 'NLDAS_' + clim_var + '_199901-202211.csv')
    clim_data.index = pd.to_datetime(clim_data.index)
    clim_data = clim_data[clim_data.index.year >= 2004]

    for wnv_month in np.arange(1, 13, 1):
        mon_wnv_data = wnv_data_mon_log[wnv_data_mon_log.index.month == wnv_month]
        for clim_month in np.arange(wnv_month, wnv_month - 12, -1):
            if clim_month < 1:
                clim_month += 12
            if clim_month < wnv_month:
                mon_clim_data = clim_data[(clim_data.index.month == clim_month) & (clim_data.index.year >= mon_wnv_data.index.year.min())]
            elif clim_month >= wnv_month:
                mon_clim_data = clim_data[(clim_data.index.month == clim_month) & (clim_data.index.year >= mon_wnv_data.index.year.min() - 1) & (clim_data.index.year < mon_wnv_data.index.year.max())]
            mon_corrs = mon_wnv_data.reset_index().drop(columns='index').corrwith(mon_clim_data.reset_index().drop(columns='index'))
            row_num = len(mon_corr_df.index)
            mon_corr_df.loc[row_num, 'wnv_month'] = wnv_month
            mon_corr_df.loc[row_num, 'climate_month'] = clim_month
            mon_corr_df.loc[row_num, 'climate_variable'] = clim_var
            mon_corr_df.loc[row_num, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']] = np.round(mon_corrs, 2)


mon_corr_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_monthly_log_wnv_monthly_clim_corr.csv', na_rep='nan')
warnings.filterwarnings('default', category=UserWarning)  # reactivating warnings




# Looking at annual autocorrelation of WNV cases
# prepping WNND data
wnv_data_mon = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv',
    index_col=0)
wnv_data_mon.index = pd.to_datetime(wnv_data_mon.index)
wnv_data_mon = wnv_data_mon[wnv_data_mon.index.year >= 2005]
wnv_data_ann = wnv_data_mon.groupby(wnv_data_mon.index.year).sum()
wnv_data_log = np.log10(wnv_data_ann + 1)

auto_corr_df = pd.DataFrame(columns=['lag', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'])
for lag in np.arange(1, 6):
    auto_corrs = wnv_data_log.iloc[lag:].reset_index().drop(columns='index').corrwith(wnv_data_log[:(-lag)].reset_index().drop(columns='index'))
    auto_corr_df.loc[lag - 1, 'lag'] = lag
    auto_corr_df.loc[lag - 1, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']] = np.round(auto_corrs, 2)
auto_corr_df.index = auto_corr_df.lag
auto_corr_df = auto_corr_df.drop(columns='lag')

fig, ax = plt.subplots(1, 1, figsize=(8, 4), layout='constrained')
sns.heatmap(data=auto_corr_df.values.astype(float), ax=ax, annot=True, cmap='bwr_r', square=True, vmin=-0.5, vmax=0.5, cbar_kws={'shrink': 0.6},
            annot_kws={'fontsize':8},
            xticklabels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'],
            yticklabels = ['1', '2', '3', '4', '5'])
plt.title(domain_names[str(domain)] + ' Correlation Heat Map', fontsize=20)
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/correlations/annual_log_wnv_heat_map_domain_' + str(domain))
plt.close(fig)






mon_corrs = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_monthly_log_wnv_monthly_clim_corr.csv', index_col=0)
wnv_data_mon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv', index_col=0)
wnv_data_mon.index = pd.to_datetime(wnv_data_mon.index)
wnv_data_mon = wnv_data_mon[wnv_data_mon.index.year >= 2005]

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}
clim_var_name = {'tmean': 'Mean Temp', 'tmax': 'Max Temp', 'tmin': 'Min Temp', 'td': 'Dew Point',
                 'vpdmax': 'Max Vapor Pressure Deficit', 'vpdmin': 'Min Vapor Pressure Deficit',
                 'precip': 'Precipitation', 'PDSI': 'PDSI', 'snow_depth': 'Snow Depth', 'snow_melt': 'Snow Melt',
                 'soil_moisture_5': '5 mm Soil Moisture', 'soil_moisture_50': '50 mm Soil Moisture',
                 'soil_moisture_100': '100 mm Soil Moisture', 'soil_temp_5': '5 mm Soil Temp',
                 'soil_temp_25': '25 mm Soil Temp', 'soil_temp_70': '70 mm Soil Temp', 'soil_temp_150': '150 mm Soil Temp',
                 'water_equivalent_snow_depth': 'Water Equivalent Snow Depth'}
domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                '16': 'Pacific Northwest', '17': 'Pacific Southwest'}

matplotlib.use('AGG')
for domain in np.arange(1, 18):
    wnv_domain_mon = wnv_data_mon[str(domain)].groupby(wnv_data_mon.index.month).sum()
    for wnv_month in np.arange(1, 13):
        if wnv_domain_mon[wnv_month] == 0:
            continue
        mon_domain_corrs = mon_corrs.loc[mon_corrs['wnv_month']==wnv_month, ['climate_month', 'climate_variable', str(domain)]]
        mon_domain_corrs = mon_domain_corrs.loc[mon_domain_corrs['climate_variable'].isin(['tmean', 'td', 'vpdmax', 'precip', 'PDSI', 'soil_moisture_50'])]
        mon_domain_corrs = mon_domain_corrs.reset_index().pivot(columns='climate_variable',index='climate_month',values=str(domain))
        fig, ax = plt.subplots(1, 1, figsize=(6, 8), layout='constrained')
        sns.heatmap(mon_domain_corrs, ax=ax, annot=True, cmap='bwr_r', square=True, mask = abs(mon_domain_corrs) < 0.32, cbar_kws={'shrink': 0.6}, annot_kws={'fontsize':8})  # .32 is p=.2
        plt.title(domain_names[str(domain)] + ' Correlation Heat Map with '+ month_names[wnv_month] + ' WNV', fontsize=14)
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/neon_domain_exploratory/correlations/domain_' + str(domain) + '_mon_' + str(wnv_month) + '_corr_heat_map')
        plt.close(fig)








#% --------------------------------------------------------------------------------------------------------------------
# Testing out a LOOCV (Leave One Out Cross-Validation) to test effectiveness of regressions
# ----------------------------------------------------------------------------------------------------------------------

# importing libraries
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error as mae

# loading WNV data
wnv_data_mon = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv',
    index_col=0)
wnv_data_mon.index = pd.to_datetime(wnv_data_mon.index)
wnv_data_mon = wnv_data_mon[wnv_data_mon.index.year >= 2005]
wnv_data_ann = wnv_data_mon.groupby(wnv_data_mon.index.year).sum()
wnv_data_log = np.log10(wnv_data_ann + 1)

# loading climate data function
def load_clim_data(fname):
    clim_df = pd.read_csv(fname, index_col=0)
    clim_df.index = pd.to_datetime(clim_df.index)
    clim_df = clim_df[clim_df.index >= '2003-01-01']
    return clim_df

# loading climate data
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
# test_clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202211.csv')
# test_clim_data = load_clim_data(fpath + 'prism_tmean_199901-202211.csv')
test_clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202211.csv')
test_clim_data.index = pd.to_datetime(test_clim_data.index)


# clim_vars = ['tmean', 'tmax', 'tmin', 'td', 'vpdmax', 'vpdmin', 'precip', 'PDSI', 'snow_depth', 'snow_melt', 'soil_moisture_5',
#              'soil_moisture_50', 'soil_moisture_100', 'soil_temp_5', 'soil_temp_25', 'soil_temp_70', 'soil_temp_150', 'water_equivalent_snow_depth']
# for clim_var in clim_vars:
#     if np.isin(clim_var, clim_vars[:7]):  # prism variables
#         clim_data = load_clim_data(fpath + 'prism_' + clim_var + '_199901-202211.csv')
#     elif np.isin(clim_var, clim_vars[7]):  # gridMET variable
#         clim_data = load_clim_data(fpath + 'gridMET_' + clim_var + '_199901-202211.csv')
#     elif np.isin(clim_var, clim_vars[8:]):  # NLDAS variables
#         clim_data = load_clim_data(fpath + 'NLDAS_' + clim_var + '_199901-202211.csv')
#     clim_data.index = pd.to_datetime(clim_data.index)
#     clim_data = clim_data[clim_data.index.year >= 2004]
#
#     for mon in np.arange(1, 13, 1):
#         if mon <= 9:
#             mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year > 2004)]
#         elif mon >= 10:
#             mon_clim_data = clim_data[(clim_data.index.month == mon) & (clim_data.index.year < 2022)]
#         mon_corrs = wnv_data_log.reset_index().drop(columns = 'index').corrwith(mon_clim_data.reset_index().drop(columns = 'index'))
#         wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index), 'climate_month'] = mon
#         wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index)-1, 'climate_variable'] = clim_var
#         wnv_clim_mon_corrs.loc[len(wnv_clim_mon_corrs.index)-1, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17']] = np.round(mon_corrs, 2)



# testing out a LOOCV

# bootstrapping relationships
bootstrap_num = 1000
neon_domain = 1
mon = 6
if np.isin(mon, [10, 11, 12]):
    test_clim_data = test_clim_data[(test_clim_data.index.year >= 2004) & (test_clim_data.index.year < 2022)]
else:
    test_clim_data = test_clim_data[test_clim_data.index.year >= 2005]

test_wnv = wnv_data_log[str(neon_domain)]
test_clim = test_clim_data[str(neon_domain)]
test_clim = test_clim_data[test_clim_data.index.month == mon][str(neon_domain)]
array_len = len(test_wnv)
array_index = np.arange(array_len)
fit_results_slope = np.zeros([array_len, bootstrap_num])
fit_results_slope[:] = np.nan
fit_results_intercept = np.zeros([array_len, bootstrap_num])
fit_results_intercept[:] = np.nan
test_median = np.zeros(array_len)
test_median[:] = np.nan
test_expected = np.zeros(array_len)
test_expected[:] = np.nan
test_05 = np.zeros(array_len)
test_05[:] = np.nan
test_95 = np.zeros(array_len)
test_95[:] = np.nan
all_predictions = np.zeros([array_len, bootstrap_num])
all_predictions[:] = np.nan

for i in np.arange(array_len):
    wnv_sub = np.take(test_wnv.values, array_index[array_index != i])
    clim_sub = np.take(test_clim.values, array_index[array_index != i])
    test_int = np.random.randint(len(wnv_sub), size=[len(wnv_sub), bootstrap_num])
    for b in np.arange(bootstrap_num):
        fit_iteration = np.polyfit(clim_sub[test_int[:,b]], wnv_sub[test_int[:,b]], 1)
        fit_results_slope[i, b] = fit_iteration[0]
        fit_results_intercept[i, b] = fit_iteration[1]
        all_predictions[i, b] = fit_iteration[0] * test_clim[i] + fit_iteration[1]

for i in np.arange(array_len):
    test_median[i] = np.median(all_predictions[i, :])
    test_05[i] = np.quantile(all_predictions[i, :], q=0.05)
    test_95[i] = np.quantile(all_predictions[i, :], q=0.95)
    # test_expected[i] = np.median(fit_results_slope[i, :])*test_clim[i] + np.median(fit_results_intercept[i, :])
    # test_05_ind = np.where(fit_results_slope[i,:] == np.quantile(fit_results_slope[i, :], q=0.05, interpolation='nearest'))[0][0]
    # test_05[i] = fit_results_slope[i, test_05_ind] * test_clim[i] + fit_results_intercept[i, test_05_ind]
    # test_95_ind = np.where(fit_results_slope[i,:] == np.quantile(fit_results_slope[i, :], q=0.95, interpolation='nearest'))[0][0]
    # test_95[i] = fit_results_slope[i, test_95_ind] * test_clim[i] + fit_results_intercept[i, test_95_ind]

mean_abs_error = mae(test_wnv, test_median)
correlation = np.corrcoef(test_wnv, test_median)[1][0]
r_squared = np.square(correlation)
slope_median = np.median(fit_results_slope)
slope_05 = np.quantile(fit_results_slope, q=0.05)
slope_95 = np.quantile(fit_results_slope, q=0.95)
clim_std = np.std(test_clim)
wnv_median = np.median(test_wnv)
delta_median = 10 ** ((slope_median * clim_std)/2 + wnv_median) - 10 ** (wnv_median - (slope_median * clim_std)/2)
delta_median_per = delta_median / 10 ** wnv_median * 100  # percentage change of median attributable to regression over one standard deviation (hopefully that makes sense)
delta_05 = 10 ** ((slope_05 * clim_std)/2 + wnv_median) - 10 ** (wnv_median - (slope_05 * clim_std)/2)
delta_05_per = delta_05 / 10 ** wnv_median * 100
delta_95 = 10 ** ((slope_95 * clim_std)/2 + wnv_median) - 10 ** (wnv_median - (slope_95 * clim_std)/2)
delta_95_per = delta_95 / 10 ** wnv_median * 100

print('Correlation: ' + str(np.round(correlation, 3)))
print('R-Squared: ' + str(np.round(r_squared, 3)))
print('MAE: ' + str(np.round(mean_abs_error, 3)))
print('Median Slope: ' + str(np.round(slope_median, 3)))
print('Median Delta: ' + str(np.round(delta_median, 3)))
print('Median Delta (%): ' + str(np.round(delta_median_per, 3)))
print('5th Slope: ' + str(np.round(slope_05, 3)))
print('5th Delta: ' + str(np.round(delta_05, 3)))
print('5th Delta (%): ' + str(np.round(delta_05_per, 3)))
print('95th Slope: ' + str(np.round(slope_95, 3)))
print('95th Delta: ' + str(np.round(delta_95, 3)))
print('95th Delta (%): ' + str(np.round(delta_95_per, 3)))

plt.scatter(test_expected, test_wnv)
plt.xlabel('Modeled Cases (log 10)')
plt.ylabel('Reported Cases (log 10)')
plt.title('Forecast vs Actual WNND Caseload (Pacific Southwest)')
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/NEON_poc_' + str(neon_domain) + '.png')






year_ind_array = np.random.randint(array_len, size=bootstrap_num)
for i in np.arange(bootstrap_num):
    year_ind = year_ind_array[i]
    wnv_sub = np.take(test_wnv.values, array_index[array_index != i])
    clim_sub = np.take(test_clim.values, array_index[array_index != i])
    fit_results[i] = np.polyfit(clim_sub, wnv_sub, 1)[0]








## Converting R output to format needed for WIS
for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    print(domain_num)

    clim_forecast_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_loo_clim_' + domain_num + '.csv', index_col=0)

    quantile_per = {'1%': '0.01', '2.5%': '0.025', '5%': '0.05', '10%': '0.10', '15%': '0.15', '20%': '0.20', '25%': '0.25',
                    '30%': '0.30', '35%': '0.35', '40%': '0.40', '45%': '0.45', '50%': '0.5', '55%': '0.55', '60%': '0.6',
                    '65%': '0.65', '70%': '0.70', '75%': '0.75', '80%': '0.80', '85%': '0.85', '90%': '0.90',
                    '95%': '0.95', '97.5%': '0.975', '99%': '0.99'}
    domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                    '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                    '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                    '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                    '16': 'Pacific Northwest', '17': 'Pacific Southwest'}

    clim_forecast_df_exploded = pd.DataFrame(columns=['domain_name', 'domain_num', 'year', 'quantile', 'value'])
    # domain_num = 9
    for i in clim_forecast_df.columns:
        for j, row in clim_forecast_df[i].items():
            row_num = len(clim_forecast_df_exploded)
            clim_forecast_df_exploded.loc[row_num, 'domain_name'] = domain_names[str(domain_num)]
            clim_forecast_df_exploded.loc[row_num, 'domain_num'] = domain_num
            clim_forecast_df_exploded.loc[row_num, 'year'] = int(i)
            clim_forecast_df_exploded.loc[row_num, 'quantile'] = quantile_per[j]
            clim_forecast_df_exploded.loc[row_num, 'value'] = np.round(row, 2)

    clim_forecast_df_exploded.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_multi_' + domain_num + '_exploded.csv')


## Attempting to compare negative binomial vs climate-informed forecast
import scoring_2  # see GitHub for code

# first getting the score for the negative binomial baseline
file_name = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/neon_domains_nb_exploded_loo.csv'

wnv_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv', index_col=0)
wnv_df.index = pd.to_datetime(wnv_df.index)
wnv_df = wnv_df[wnv_df.index.year >= 2005]
wnv_df_ann = wnv_df.groupby(wnv_df.index.year).sum()

domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                '16': 'Pacific Northwest', '17': 'Pacific Southwest'}

print(file_name)

for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    print(domain_num)
    df_forecast = pd.read_csv(file_name, index_col=0)
    df_forecast = df_forecast[df_forecast['domain_num'] == int(domain_num)]

    # looping through all years
    array_p010 = np.array([])
    array_p025 = np.array([])
    array_p050 = np.array([])
    array_p100 = np.array([])
    array_p150 = np.array([])
    array_p200 = np.array([])
    array_p250 = np.array([])
    array_p300 = np.array([])
    array_p350 = np.array([])
    array_p400 = np.array([])
    array_p450 = np.array([])
    array_p500 = np.array([])
    array_p550 = np.array([])
    array_p600 = np.array([])
    array_p650 = np.array([])
    array_p700 = np.array([])
    array_p750 = np.array([])
    array_p800 = np.array([])
    array_p850 = np.array([])
    array_p900 = np.array([])
    array_p950 = np.array([])
    array_p975 = np.array([])
    array_p990 = np.array([])
    obs_cases = np.array([])
    year_list = df_forecast['year'].unique()
    for i, year in enumerate(year_list):
        df_forecast_sub = df_forecast[df_forecast['year'] == year]
        array_p010 = np.append(array_p010, df_forecast_sub['value'].iloc[0])
        array_p025 = np.append(array_p025, df_forecast_sub['value'].iloc[1])
        array_p050 = np.append(array_p050, df_forecast_sub['value'].iloc[2])
        array_p100 = np.append(array_p100, df_forecast_sub['value'].iloc[3])
        array_p150 = np.append(array_p150, df_forecast_sub['value'].iloc[4])
        array_p200 = np.append(array_p200, df_forecast_sub['value'].iloc[5])
        array_p250 = np.append(array_p250, df_forecast_sub['value'].iloc[6])
        array_p300 = np.append(array_p300, df_forecast_sub['value'].iloc[7])
        array_p350 = np.append(array_p350, df_forecast_sub['value'].iloc[8])
        array_p400 = np.append(array_p400, df_forecast_sub['value'].iloc[9])
        array_p450 = np.append(array_p450, df_forecast_sub['value'].iloc[10])
        array_p500 = np.append(array_p500, df_forecast_sub['value'].iloc[11])
        array_p550 = np.append(array_p550, df_forecast_sub['value'].iloc[12])
        array_p600 = np.append(array_p600, df_forecast_sub['value'].iloc[13])
        array_p650 = np.append(array_p650, df_forecast_sub['value'].iloc[14])
        array_p700 = np.append(array_p700, df_forecast_sub['value'].iloc[15])
        array_p750 = np.append(array_p750, df_forecast_sub['value'].iloc[16])
        array_p800 = np.append(array_p800, df_forecast_sub['value'].iloc[17])
        array_p850 = np.append(array_p850, df_forecast_sub['value'].iloc[18])
        array_p900 = np.append(array_p900, df_forecast_sub['value'].iloc[19])
        array_p950 = np.append(array_p950, df_forecast_sub['value'].iloc[20])
        array_p975 = np.append(array_p975, df_forecast_sub['value'].iloc[21])
        array_p990 = np.append(array_p990, df_forecast_sub['value'].iloc[22])
    # obs_cases = np.array([72, 85, 168, 16, 12, 9, 3, 106, 157, 30, 31, 70, 53, 142, 9, 10, 31, 40])
    obs_cases = np.array(wnv_df_ann[domain_num])
    quantile_dict_test = {
        0.01: np.log(array_p010 + 1),
        0.025: np.log(array_p025 + 1),
        0.05: np.log(array_p050 + 1),
        0.1: np.log(array_p100 + 1),
        0.15: np.log(array_p150 + 1),
        0.2: np.log(array_p200 + 1),
        0.25: np.log(array_p250 + 1),
        0.3: np.log(array_p300 + 1),
        0.35: np.log(array_p350 + 1),
        0.4: np.log(array_p400 + 1),
        0.45: np.log(array_p450 + 1),
        # 0.5: np.log(array_p500 + 1),
        0.55: np.log(array_p550 + 1),
        0.6: np.log(array_p600 + 1),
        0.65: np.log(array_p650 + 1),
        0.7: np.log(array_p700 + 1),
        0.75: np.log(array_p750 + 1),
        0.8: np.log(array_p800 + 1),
        0.85: np.log(array_p850 + 1),
        0.9: np.log(array_p900 + 1),
        0.95: np.log(array_p950 + 1),
        0.975: np.log(array_p975 + 1),
        0.99: np.log(array_p990 + 1),
    }
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    scores, sharpness, calibrations, overprediction, underprediction = scoring_2.weighted_interval_score_fast(np.log(obs_cases + 1), alphas=alphas, medians=np.log(array_p500 + 1), weights=None, q_dict=quantile_dict_test)
    print(file_name)
    print('Overall score:')
    print(np.round(scores.mean(), 4))
    print('Calibration score:')
    print(np.round(calibrations.mean(), 4))
    print('Sharpness score:')
    print(np.round(sharpness.mean(), 4))
    d = {'scores': scores, 'sharpness': sharpness, 'calibration': calibrations, 'overprediction': overprediction,
         'underprediction': underprediction}  # creating data dictionary
    df_scoring = pd.DataFrame(data=d, index=np.arange(2005, 2023))  # creating dataframe for saving
    df_scoring['domain_num'] = domain_num
    df_scoring['domain_name'] = domain_names[str(domain_num)]
    df_scoring = df_scoring[
        ['domain_name', 'domain_num', 'scores', 'sharpness', 'calibration', 'overprediction',
         'underprediction']]
    df_scoring.to_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_hist_nb_wis_' + domain_num + '.csv')


# now getting the score for the climate-informed forecast
import glob
all_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_quantiles_*.csv')

# executing scoring code
file_list = all_csv_files.copy()
for file in file_list:
    file_name = file
    domain_num = int(file.split('/')[-1].split('_')[4][0])
    print(file)
    df_forecast = pd.read_csv(file_name)
    df_forecast = df_forecast[df_forecast['domain_num']==domain_num]
    # looping through all years
    array_p010 = np.array([])
    array_p025 = np.array([])
    array_p050 = np.array([])
    array_p100 = np.array([])
    array_p150 = np.array([])
    array_p200 = np.array([])
    array_p250 = np.array([])
    array_p300 = np.array([])
    array_p350 = np.array([])
    array_p400 = np.array([])
    array_p450 = np.array([])
    array_p500 = np.array([])
    array_p550 = np.array([])
    array_p600 = np.array([])
    array_p650 = np.array([])
    array_p700 = np.array([])
    array_p750 = np.array([])
    array_p800 = np.array([])
    array_p850 = np.array([])
    array_p900 = np.array([])
    array_p950 = np.array([])
    array_p975 = np.array([])
    array_p990 = np.array([])
    obs_cases = np.array([])
    year_list = df_forecast['year'].unique()
    for i, year in enumerate(year_list):
        df_forecast_sub = df_forecast[df_forecast['year'] == year]
        array_p010 = np.append(array_p010, df_forecast_sub['value'].iloc[0])
        array_p025 = np.append(array_p025, df_forecast_sub['value'].iloc[1])
        array_p050 = np.append(array_p050, df_forecast_sub['value'].iloc[2])
        array_p100 = np.append(array_p100, df_forecast_sub['value'].iloc[3])
        array_p150 = np.append(array_p150, df_forecast_sub['value'].iloc[4])
        array_p200 = np.append(array_p200, df_forecast_sub['value'].iloc[5])
        array_p250 = np.append(array_p250, df_forecast_sub['value'].iloc[6])
        array_p300 = np.append(array_p300, df_forecast_sub['value'].iloc[7])
        array_p350 = np.append(array_p350, df_forecast_sub['value'].iloc[8])
        array_p400 = np.append(array_p400, df_forecast_sub['value'].iloc[9])
        array_p450 = np.append(array_p450, df_forecast_sub['value'].iloc[10])
        array_p500 = np.append(array_p500, df_forecast_sub['value'].iloc[11])
        array_p550 = np.append(array_p550, df_forecast_sub['value'].iloc[12])
        array_p600 = np.append(array_p600, df_forecast_sub['value'].iloc[13])
        array_p650 = np.append(array_p650, df_forecast_sub['value'].iloc[14])
        array_p700 = np.append(array_p700, df_forecast_sub['value'].iloc[15])
        array_p750 = np.append(array_p750, df_forecast_sub['value'].iloc[16])
        array_p800 = np.append(array_p800, df_forecast_sub['value'].iloc[17])
        array_p850 = np.append(array_p850, df_forecast_sub['value'].iloc[18])
        array_p900 = np.append(array_p900, df_forecast_sub['value'].iloc[19])
        array_p950 = np.append(array_p950, df_forecast_sub['value'].iloc[20])
        array_p975 = np.append(array_p975, df_forecast_sub['value'].iloc[21])
        array_p990 = np.append(array_p990, df_forecast_sub['value'].iloc[22])
        # obs_cases = np.append(obs_cases, df_county_wnv['COUNT'][df_county_wnv['location'] == location].values[0])
        obs_cases = np.array(wnv_df_ann[str(domain_num)])
    quantile_dict_test = {
        0.01: np.log(array_p010 + 1),
        0.025: np.log(array_p025 + 1),
        0.05: np.log(array_p050 + 1),
        0.1: np.log(array_p100 + 1),
        0.15: np.log(array_p150 + 1),
        0.2: np.log(array_p200 + 1),
        0.25: np.log(array_p250 + 1),
        0.3: np.log(array_p300 + 1),
        0.35: np.log(array_p350 + 1),
        0.4: np.log(array_p400 + 1),
        0.45: np.log(array_p450 + 1),
        # 0.5: np.log(array_p500 + 1),
        0.55: np.log(array_p550 + 1),
        0.6: np.log(array_p600 + 1),
        0.65: np.log(array_p650 + 1),
        0.7: np.log(array_p700 + 1),
        0.75: np.log(array_p750 + 1),
        0.8: np.log(array_p800 + 1),
        0.85: np.log(array_p850 + 1),
        0.9: np.log(array_p900 + 1),
        0.95: np.log(array_p950 + 1),
        0.975: np.log(array_p975 + 1),
        0.99: np.log(array_p990 + 1),
    }
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    scores, sharpness, calibrations, overprediction, underprediction = scoring_2.weighted_interval_score_fast(np.log(obs_cases + 1), alphas=alphas, medians=np.log(array_p500 + 1), weights=None, q_dict=quantile_dict_test)
    print(file_name)
    print('Overall score:')
    print(np.round(scores.mean(), 4))
    print('Calibration score:')
    print(np.round(calibrations.mean(), 4))
    print('Sharpness score:')
    print(np.round(sharpness.mean(), 4))
    d = {'scores': scores, 'sharpness': sharpness, 'calibration': calibrations, 'overprediction': overprediction, 'underprediction': underprediction}  # creating data dictionary
    df_scoring = pd.DataFrame(data=d, index=np.arange(2005, 2023))  # creating dataframe for saving
    df_scoring['domain_num'] = domain_num
    df_scoring['domain_name'] = domain_names[str(domain_num)]
    df_scoring = df_scoring[
        ['domain_name', 'domain_num', 'scores', 'sharpness', 'calibration', 'overprediction',
         'underprediction']]
    df_scoring.to_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_' + str(domain_num) + '.csv')


plot_df = pd.DataFrame()
# testing for statistical significance
for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    hist_nb = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_hist_nb_wis_' + domain_num + '.csv')
    hist_nb = hist_nb.rename(columns={'Unnamed: 0': 'year'})
    # clim_fore = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_loo_clim_' + domain_num + '.csv')
    clim_fore = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_' + domain_num + '.csv')
    clim_fore = clim_fore.rename(columns={'Unnamed: 0': 'year'})

    bootstrap_num = 1000
    year_ind = hist_nb.index.values
    arr_len = len(year_ind)
    hist_nb_vals = np.zeros([bootstrap_num, arr_len])
    hist_nb_vals[:] = np.nan
    clim_fore_vals = np.zeros([bootstrap_num, arr_len])
    clim_fore_vals[:] = np.nan
    for i in np.arange(bootstrap_num):
        bootstrap_arr = np.random.choice(hist_nb.index, arr_len, replace=True)
        hist_nb_vals[i, :] = hist_nb['scores'].iloc[bootstrap_arr]
        clim_fore_vals[i, :] = clim_fore['scores'].iloc[bootstrap_arr]
    # for i in np.arange(bootstrap_num):
    #     hist_nb_vals[i, :] = hist_nb['calibration'].sample(n=arr_len, replace=True)
    #     clim_fore_vals[i, :] = clim_fore['calibration'].sample(n=arr_len, replace=True)
    plot_df[domain_num]= np.mean(clim_fore_vals, axis=1) - np.mean(hist_nb_vals, axis=1)
    num_diff = np.sum(np.mean(clim_fore_vals, axis=1) < np.mean(hist_nb_vals, axis=1))
    print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_vals), 3)) + ',' + str(np.round(np.mean(hist_nb_vals), 3)) + ',' + str((bootstrap_num - num_diff)/bootstrap_num))



# testing an all-domain forecast for statistical significance
import glob

# gathering all climate-informed regression scores
all_clim_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_*.csv')
df_clim_list = []
for filename in all_clim_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_clim_list.append(df)

clim_forecast_wis_df = pd.concat(df_clim_list, axis=0, ignore_index=True)
clim_forecast_wis_df = clim_forecast_wis_df.rename(columns={'Unnamed: 0': 'year'})

# gathering all historical negative binomial scores
all_hist_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_hist_nb_wis_*.csv')
df_hist_list = []
for filename in all_hist_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_hist_list.append(df)

hist_nb_wis_df = pd.concat(df_hist_list, axis=0, ignore_index=True)
hist_nb_wis_df = hist_nb_wis_df.rename(columns={'Unnamed: 0': 'year'})

# sorting columns to ensure they match
clim_forecast_sorted = clim_forecast_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)
hist_nb_sorted = hist_nb_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)

# running bootstrapping over the two
bootstrap_num = 1000
arr_len = len(clim_forecast_sorted)
hist_nb_wis_vals = np.zeros([bootstrap_num, arr_len])
hist_nb_wis_vals[:] = np.nan
clim_fore_wis_vals = np.zeros([bootstrap_num, arr_len])
clim_fore_wis_vals[:] = np.nan
for i in np.arange(bootstrap_num):
    bootstrap_arr = np.random.choice(clim_forecast_sorted.index, arr_len, replace=True)
    hist_nb_wis_vals[i, :] = hist_nb_sorted['scores'].iloc[bootstrap_arr]
    clim_fore_wis_vals[i, :] = clim_forecast_sorted['scores'].iloc[bootstrap_arr]
plot_df['national'] = np.mean(clim_fore_wis_vals, axis=1) - np.mean(hist_nb_wis_vals, axis=1)
num_diff = np.sum(np.mean(clim_fore_wis_vals, axis=1) < np.mean(hist_nb_wis_vals, axis=1))
print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_wis_vals), 3)) + ',' + str(
    np.round(np.mean(hist_nb_wis_vals), 3)) + ',' + str((bootstrap_num - num_diff) / bootstrap_num))

import seaborn as sns
sns.set_style('whitegrid')
fig, ax = plt.subplots(layout='constrained')
plt.axhline(0, c='k')
ax = sns.boxplot(data=plot_df, ax=ax, whis=(5, 95), width=0.6, showfliers=False, flierprops={"marker": ".", 'markersize': 1}, boxprops=dict(alpha=0.7))
plt.xlabel('NEON Region', fontsize=10, weight='bold')
plt.ylabel('Difference in Forecast Skill (WIS)', fontsize=10, weight='bold')
x_tick_labels = ['Pacific Southwest', 'Desert Southwest', 'Great Lakes', 'Prairie Peninsula', 'Southeast', 'Southern Plains',
                 'Central Plains', 'Ozarks Complex', 'Northern Plains', 'Mid Atlantic', 'Northeast', 'National']
ax.set_xticklabels(x_tick_labels, rotation=60, ha='right')
plt.title('Climate-Informed Forecast vs Historical Benchmark', fontsize=14, weight='bold')
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/diff_in_forecast_skill_multi_loo_clim.svg')


# now getting the score for the upscaled 2022 ensemble forecasts
import glob

# executing scoring code
for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    df_forecast = pd.read_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_forecast_spatially_dependent.csv',
        index_col=0)
    df_forecast = df_forecast[df_forecast['neon_domain']==int(domain_num)]
    # looping through all domains
    array_p010 = np.array([])
    array_p025 = np.array([])
    array_p050 = np.array([])
    array_p100 = np.array([])
    array_p150 = np.array([])
    array_p200 = np.array([])
    array_p250 = np.array([])
    array_p300 = np.array([])
    array_p350 = np.array([])
    array_p400 = np.array([])
    array_p450 = np.array([])
    array_p500 = np.array([])
    array_p550 = np.array([])
    array_p600 = np.array([])
    array_p650 = np.array([])
    array_p700 = np.array([])
    array_p750 = np.array([])
    array_p800 = np.array([])
    array_p850 = np.array([])
    array_p900 = np.array([])
    array_p950 = np.array([])
    array_p975 = np.array([])
    array_p990 = np.array([])
    obs_cases = np.array([])
    year_list = 2022

    df_forecast_sub = df_forecast
    array_p010 = np.append(array_p010, df_forecast_sub['value'].iloc[0])
    array_p025 = np.append(array_p025, df_forecast_sub['value'].iloc[1])
    array_p050 = np.append(array_p050, df_forecast_sub['value'].iloc[2])
    array_p100 = np.append(array_p100, df_forecast_sub['value'].iloc[3])
    array_p150 = np.append(array_p150, df_forecast_sub['value'].iloc[4])
    array_p200 = np.append(array_p200, df_forecast_sub['value'].iloc[5])
    array_p250 = np.append(array_p250, df_forecast_sub['value'].iloc[6])
    array_p300 = np.append(array_p300, df_forecast_sub['value'].iloc[7])
    array_p350 = np.append(array_p350, df_forecast_sub['value'].iloc[8])
    array_p400 = np.append(array_p400, df_forecast_sub['value'].iloc[9])
    array_p450 = np.append(array_p450, df_forecast_sub['value'].iloc[10])
    array_p500 = np.append(array_p500, df_forecast_sub['value'].iloc[11])
    array_p550 = np.append(array_p550, df_forecast_sub['value'].iloc[12])
    array_p600 = np.append(array_p600, df_forecast_sub['value'].iloc[13])
    array_p650 = np.append(array_p650, df_forecast_sub['value'].iloc[14])
    array_p700 = np.append(array_p700, df_forecast_sub['value'].iloc[15])
    array_p750 = np.append(array_p750, df_forecast_sub['value'].iloc[16])
    array_p800 = np.append(array_p800, df_forecast_sub['value'].iloc[17])
    array_p850 = np.append(array_p850, df_forecast_sub['value'].iloc[18])
    array_p900 = np.append(array_p900, df_forecast_sub['value'].iloc[19])
    array_p950 = np.append(array_p950, df_forecast_sub['value'].iloc[20])
    array_p975 = np.append(array_p975, df_forecast_sub['value'].iloc[21])
    array_p990 = np.append(array_p990, df_forecast_sub['value'].iloc[22])

    obs_cases = wnv_df_ann[str(domain_num)].loc[2022]
    quantile_dict_test = {
        0.01: np.log(array_p010 + 1),
        0.025: np.log(array_p025 + 1),
        0.05: np.log(array_p050 + 1),
        0.1: np.log(array_p100 + 1),
        0.15: np.log(array_p150 + 1),
        0.2: np.log(array_p200 + 1),
        0.25: np.log(array_p250 + 1),
        0.3: np.log(array_p300 + 1),
        0.35: np.log(array_p350 + 1),
        0.4: np.log(array_p400 + 1),
        0.45: np.log(array_p450 + 1),
        # 0.5: np.log(array_p500 + 1),
        0.55: np.log(array_p550 + 1),
        0.6: np.log(array_p600 + 1),
        0.65: np.log(array_p650 + 1),
        0.7: np.log(array_p700 + 1),
        0.75: np.log(array_p750 + 1),
        0.8: np.log(array_p800 + 1),
        0.85: np.log(array_p850 + 1),
        0.9: np.log(array_p900 + 1),
        0.95: np.log(array_p950 + 1),
        0.975: np.log(array_p975 + 1),
        0.99: np.log(array_p990 + 1),
    }
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    scores, sharpness, calibrations, overprediction, underprediction = scoring_2.weighted_interval_score_fast(np.log(obs_cases + 1), alphas=alphas, medians=np.log(array_p500 + 1), weights=None, q_dict=quantile_dict_test)
    print(domain_num)
    print('Overall score:')
    print(np.round(scores.mean(), 4))
    # print('Calibration score:')
    # print(np.round(calibrations.mean(), 4))
    # print('Sharpness score:')
    # print(np.round(sharpness.mean(), 4))
    d = {'scores': scores, 'sharpness': sharpness, 'calibration': calibrations, 'overprediction': overprediction, 'underprediction': underprediction}  # creating data dictionary
    df_scoring = pd.DataFrame(data=d, index=np.arange(2005, 2023))  # creating dataframe for saving
    df_scoring['domain_num'] = domain_num
    df_scoring['domain_name'] = domain_names[str(domain_num)]
    df_scoring = df_scoring[
        ['domain_name', 'domain_num', 'scores', 'sharpness', 'calibration', 'overprediction',
         'underprediction']]
    df_scoring.to_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_' + str(domain_num) + '.csv')














## creating county caseload-weighted averages of climatic conditions for NEON domains
import pandas as pd
import numpy as np
from glob import glob

path = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/*1.csv'
all_files = sorted(glob(path, recursive=True))
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')
neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
for file in all_files:
    name = file.split('/')[-1].split('1999')[0]
    clim_data = pd.read_csv(file, index_col=0)
    wnv_data = pd.read_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv',
        index_col=0)
    missing_clim_counties = np.setdiff1d(wnv_data.columns, clim_data.columns)
    wnv_data = wnv_data.drop(columns=missing_clim_counties)
    # clim_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/prism_tmean_199901-202211.csv', index_col = 0)    neon_domains = np.sort(county_to_neon['NEON_domain'].unique())
    weighted_clim = np.zeros([len(neon_domains), np.shape(clim_data)[0]])
    for domain in neon_domains:
        neon_fips = county_to_neon.loc[county_to_neon['NEON_domain'] == domain]['FIPS'].values.astype(str)
        neon_fips = np.delete(neon_fips, np.isin(neon_fips, missing_clim_counties))
        # neon_fips = np.delete(neon_fips, np.where(neon_fips == missing_clim_counties))
        wnv_sub = wnv_data[neon_fips]
        wnv_year = pd.to_datetime(wnv_sub.index).year
        clim_sub = clim_data[wnv_sub.columns]  # only taking climate data from counties with historical WNV caseload
        clim_year = pd.to_datetime(clim_sub.index).year
        for yr in clim_year.unique():
            wnv_fips_yr_sub = wnv_sub.loc[:, wnv_sub[wnv_year != yr].sum() > 0].columns
            weighted_clim[domain-1, clim_year == yr] = np.sum(clim_sub.loc[clim_year == yr, wnv_fips_yr_sub] * np.log(wnv_sub.loc[wnv_year != yr, wnv_fips_yr_sub].sum()), axis=1) / np.log(wnv_sub.loc[wnv_year != yr, wnv_fips_yr_sub].sum()).sum()
    weighted_clim_data = pd.DataFrame(data=np.transpose(weighted_clim), index=clim_data.index, columns=neon_domains)
    weighted_clim_data.to_csv(
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/' + name + '199901-202211_loo.csv')

# 2005 onward?
# what log to use? 10 or natural



# starting with upscaling
df_wnv = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_county_annual_time_series.csv')
county_to_neon = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')
county_to_neon = county_to_neon.drop(columns='FIPS')

df_ensemble = pd.read_csv('/Users/ryan.harp/Documents/GitHub/WNV-forecast-data-2022/data-forecasts/ensemble-ensemble/2022-04-30-ensemble-ensemble.csv')

from itertools import product

# first attempt assuming there's complete spatial dependence within NEON regions
neon_domain_forecasts = pd.DataFrame(list(product(county_to_neon['NEON_domain'].unique(), df_ensemble['quantile'].unique())), columns = ['neon_domain', 'quantile'])
neon_domain_forecasts['value'] = np.nan
for domain in county_to_neon['NEON_domain'].unique():
    for q in df_ensemble['quantile'].unique():
        fips_list = county_to_neon[county_to_neon['NEON_domain']==domain]['location'].values
        neon_domain_forecasts.loc[(neon_domain_forecasts['neon_domain'] == domain) & (neon_domain_forecasts['quantile'] == q), 'value']= df_ensemble.loc[(df_ensemble['location'].isin(fips_list)) & (df_ensemble['quantile'] == q)]['value'].sum()
neon_domain_forecasts.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_forecast_spatially_dependent.csv')

# second attempt assuming there's complete spatial dependence within NEON regions
neon_domain_forecasts = pd.DataFrame(list(product(county_to_neon['NEON_domain'].unique(), df_ensemble['quantile'].unique())), columns = ['neon_domain', 'quantile'])
neon_domain_forecasts['value'] = np.nan

import random
import time
bootstrap_num = 1000
quantiles = df_ensemble['quantile'].unique()
start = time.time()
for domain in county_to_neon['NEON_domain'].unique():
    print(domain)
    fips_list = county_to_neon[county_to_neon['NEON_domain']==domain]['location'].values
    bootstrap_arr = np.zeros([len(fips_list), bootstrap_num])
    bootstrap_arr[:] = np.nan
    for ind, location in np.ndenumerate(fips_list):
        print(np.round(ind[0]/len(fips_list)*100, 2))
        for b in np.arange(bootstrap_num):
            # if b%100 == 0:
            #     print('bootstrap_num: ' + str(b) + ' at ' + str(time.time() - start))
            bootstrap_arr[ind[0], b] = df_ensemble.loc[(df_ensemble['location'] == location) & (df_ensemble['quantile'] == random.choice(quantiles))]['value'].values[0]
            # neon_domain_forecasts.loc[(neon_domain_forecasts['neon_domain'] == domain) & (neon_domain_forecasts['quantile'] == q), 'value']= df_ensemble.loc[(df_ensemble['location'] == location) & (df_ensemble['quantile'] == random.choice(quantiles))]['value'].values[0]
    for q in quantiles:
        neon_domain_forecasts.loc[(neon_domain_forecasts['neon_domain'] == domain) & (neon_domain_forecasts['quantile'] == q), 'value'] = np.quantile(bootstrap_arr.sum(axis=0), q)
neon_domain_forecasts.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_forecast_spatially_independent.csv')




# loading WIS score dataframe to get full complement of counties
df_obs = df_obs.merge(df_median, how='left', left_on='location', right_on='location')






# pulling all columns
from scipy.stats import spearmanr

wnv_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv', index_col=0)
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data[wnv_data.index.year >= 2005]
wnv_data_ann = wnv_data.groupby(wnv_data.index.year).sum()
wnv_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv', index_col=0)
wnv_df.index = pd.to_datetime(wnv_df.index)
wnv_df = wnv_df[wnv_df.index.year >= 2005]
wnv_df_ann = wnv_df.groupby(wnv_df.index.year).sum()

df_proportion = pd.Series(index = wnv_data.columns)
df_proportion_ann = pd.DataFrame(index = wnv_df.index.year.unique(), columns = wnv_data.columns)

for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    neon_domain = int(domain_num)
    fips_list = county_to_neon.loc[county_to_neon['NEON_domain'] == neon_domain, 'fips'].astype(str)
    test = wnv_data[np.intersect1d(wnv_data.columns, fips_list)]
    test.index = test.index.astype('datetime64[ns]')
    test = test[test.index.year >= 2005]
    test_sum = test.resample('Y').sum()
    test_sum.index = test_sum.index.year
    for fip in test_sum.columns:
        df_proportion_ann[fip] = np.round(test_sum[fip]/wnv_df_ann[domain_num],3) * 100
        df_proportion.loc[fip] = np.round(test_sum[fip].sum()/wnv_df_ann[domain_num].sum(), 3)

df_ks_test = pd.DataFrame(columns=['fips_1', 'fips_2', 'prop_p_value', 'raw_p_value'])
fips1 = []
fips2 = []
prop_p = []
raw_p = []
non_zero_fips = wnv_data_ann.sum(axis=0) > 0
fips_all = df_proportion.sort_values(ascending=False).index
fips_sub = wnv_data_ann[non_zero_fips[non_zero_fips == True].index]
fips = fips_sub.columns
fips_len = len(fips)
for i in np.arange(fips_len):
    print(i)
    if i == fips_len:
        continue
    for j in np.arange(i+1, fips_len):
        # df_len = len(df_ks_test)
        fips1.append(fips[i])
        fips2.append(fips[j])
        prop_p.append(np.round(ss.ks_2samp(df_proportion_ann[fips[i]], df_proportion_ann[fips[j]])[1], 2))
        raw_p.append(np.round(ss.ks_2samp(wnv_data_ann[fips[i]], wnv_data_ann[fips[j]])[1], 2))
df_ks_test['fips_1'] = fips1
df_ks_test['fips_2'] = fips2
df_ks_test['prop_p_value'] = prop_p
df_ks_test['raw_p_value'] = raw_p
df_ks_test.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/ks_tests.csv')

domain_num = '17'
test_fip = '6037'
p_threshold = .1
test_fip_sub = df_ks_test[(df_ks_test['fips_1'] == test_fip) | (df_ks_test['fips_2'] == test_fip)].sort_values('prop_p_value', ascending=False)
test_fip_sub = test_fip_sub[(test_fip_sub['prop_p_value'] >= p_threshold) & (test_fip_sub['raw_p_value'] >= p_threshold)]
test_fip_arr = np.concatenate((test_fip_sub['fips_1'].values, test_fip_sub['fips_2'].values))
test_fip_arr = test_fip_arr[test_fip_arr != test_fip]
test_fip_arr = np.append(test_fip_arr, test_fip)
prop_arr = np.sort(df_proportion_ann[test_fip_arr].values.flatten())
spearman_corr = spearmanr(a=wnv_df_ann[domain_num], b=df_proportion_ann[test_fip])
test_df = pd.DataFrame([wnv_df_ann[domain_num], df_proportion_ann[test_fip]]).transpose()

domain_forecast_pdf = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_loo_clim_pdf_2022_17.csv', index_col=0)
domain_forecast = domain_forecast_pdf.iloc[100].values[0]

test_arr = test_df[domain_num]
forecast_per = sum(np.abs(test_arr) < domain_forecast) / float(len(test_arr))

from scipy.stats import norm
p = norm.rvs(size=len(prop_arr), loc=forecast_per)
np.random.choice(prop_arr, 1, p=[])


for i in test_sum.sum(axis=0).sort_values(ascending=False).index:
    print(np.round(spearmanr(a=test_sum[i], b=test_sum.sum(axis=1))[0],3))

for i in test_sum.sum(axis=0).sort_values(ascending=False).index:
    print(np.round(np.corrcoef(test_sum[i], test_sum.sum(axis=1))[0][1],3))


test = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/random_effects_test_pdf_9.csv', index_col=0)
test_sum = test.sum(axis=1)
quant_array = np.array([0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
                       0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99])
for quant in np.arange(23):
    print(np.quantile(test_sum, q = quant_array[quant]))

test = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/random_effects_test_9.csv', index_col=0)









# creating county-level forecast quantiles
for neon_region in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    print(neon_region)

    reg_forecast = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/county_level_re_forecast_' + neon_region + '.csv', index_col = 0)
    quantile_array = np.array([0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
                          0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99])

    neon_region = int(neon_region)
    forecast_quantiles_county_df = pd.DataFrame(columns=['year', 'region', 'fip', 'quantile', 'value'])
    for yr in reg_forecast['year'].unique():
        print(yr)
        for county in reg_forecast['county'].unique():
            county_forecast_sub = reg_forecast.loc[(reg_forecast['year'] == yr) & (reg_forecast['county'] == county)]
            # forecast_quants = np.quantile(forecast_sub['.prediction'], quantile_array)
            county_arr = np.zeros([23, 5])
            county_arr[:] = np.nan
            county_arr[:, 0] = yr
            county_arr[:, 1] = neon_region
            county_arr[:, 2] = county
            county_arr[:, 3] = quantile_array
            county_arr[:, 4] = np.quantile(county_forecast_sub['.prediction'], quantile_array)
            county_df = pd.DataFrame(data=county_arr, columns=['year', 'region', 'fip', 'quantile', 'value'])
            forecast_quantiles_county_df = pd.concat([forecast_quantiles_county_df, county_df], ignore_index=True)
    forecast_quantiles_county_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_county_quantiles_' + str(neon_region) + '.csv')

    # creating regional-level forecast quantiles
    forecast_quantiles_region_df = pd.DataFrame(columns=['year', 'region', 'quantile', 'value'])
    for yr in reg_forecast['year'].unique():
        print(yr)
        region_forecast_sub = reg_forecast.loc[reg_forecast['year'] == yr]
        forecast_sum = region_forecast_sub.groupby('.draw')['.prediction'].sum()
        region_arr = np.zeros([23, 4])
        region_arr[:] = np.nan
        region_arr[:, 0] = yr
        region_arr[:, 1] = neon_region
        region_arr[:, 2] = quantile_array
        region_arr[:, 3] = np.quantile(forecast_sum, quantile_array)
        region_df = pd.DataFrame(data=region_arr, columns=['year', 'region', 'quantile', 'value'])
        forecast_quantiles_region_df = pd.concat([forecast_quantiles_region_df, region_df], ignore_index=True)
    forecast_quantiles_region_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_quantiles_' + str(neon_region) + '.csv')









# executing scoring code
all_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_county_quantiles_*.csv')
file_list = all_csv_files.copy()
for file in file_list:
    file_name = file
    num_region = file_name.split('/')[-1].split('_')[-1].split('.')[0]
    # model = file_name.split('/')[-1].split('.')[0][11:]
    # submission_date = file_name.split('/')[-1].split('.')[0][:10]
    # print(file)
#    df_forecast = pd.read_csv('/Users/ryan.harp/Documents/GitHub/WNV-forecast-data-2022/data-forecasts/' + file_name)
    df_forecast = pd.read_csv(file_name, index_col=0)
    df_forecast = df_forecast[df_forecast['year']==2022]
    df_county_fips = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/2021_Gaz_counties_national.txt', sep='\t')  # using "old" county file because 2022 version has weird Connecticut planning disticts
    df_state_abbr = pd.read_csv('/Users/ryan.harp/Documents/2023_WNV_Forecasting_Challenge/supplemental_files/state_fips_codes.txt', sep='|')
    df_wnv = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/raw_data/WNV_NeurobyCounty_Final2022.csv')
    df_wnv = df_wnv[['County', 'COUNT']]
    df_wnv = df_wnv.groupby('County').sum()
    df_locations = pd.read_csv('/Users/ryan.harp/Documents/2022_WNV_Forecasting_Challenge/supplemental_files/suppled-locations.csv')
    df_county_wnv = df_locations.merge(df_wnv, how='left', left_on='fips', right_on=df_wnv.index.values)
    df_county_wnv.loc[df_county_wnv['COUNT'].isna(),'COUNT'] = 0
    # looping through all locations
    array_p010 = np.array([])
    array_p025 = np.array([])
    array_p050 = np.array([])
    array_p100 = np.array([])
    array_p150 = np.array([])
    array_p200 = np.array([])
    array_p250 = np.array([])
    array_p300 = np.array([])
    array_p350 = np.array([])
    array_p400 = np.array([])
    array_p450 = np.array([])
    array_p500 = np.array([])
    array_p550 = np.array([])
    array_p600 = np.array([])
    array_p650 = np.array([])
    array_p700 = np.array([])
    array_p750 = np.array([])
    array_p800 = np.array([])
    array_p850 = np.array([])
    array_p900 = np.array([])
    array_p950 = np.array([])
    array_p975 = np.array([])
    array_p990 = np.array([])
    obs_cases = np.array([])
    location_list = df_forecast['fip'].unique()
    for i, location in enumerate(location_list):
        # print(location)
        # 74 = Maricopa County
        # i = 74
        # location = location_list[74]
        county_forecast = df_forecast[df_forecast['fip'] == location]
        array_p010 = np.append(array_p010, county_forecast['value'].iloc[0])
        array_p025 = np.append(array_p025, county_forecast['value'].iloc[1])
        array_p050 = np.append(array_p050, county_forecast['value'].iloc[2])
        array_p100 = np.append(array_p100, county_forecast['value'].iloc[3])
        array_p150 = np.append(array_p150, county_forecast['value'].iloc[4])
        array_p200 = np.append(array_p200, county_forecast['value'].iloc[5])
        array_p250 = np.append(array_p250, county_forecast['value'].iloc[6])
        array_p300 = np.append(array_p300, county_forecast['value'].iloc[7])
        array_p350 = np.append(array_p350, county_forecast['value'].iloc[8])
        array_p400 = np.append(array_p400, county_forecast['value'].iloc[9])
        array_p450 = np.append(array_p450, county_forecast['value'].iloc[10])
        array_p500 = np.append(array_p500, county_forecast['value'].iloc[11])
        array_p550 = np.append(array_p550, county_forecast['value'].iloc[12])
        array_p600 = np.append(array_p600, county_forecast['value'].iloc[13])
        array_p650 = np.append(array_p650, county_forecast['value'].iloc[14])
        array_p700 = np.append(array_p700, county_forecast['value'].iloc[15])
        array_p750 = np.append(array_p750, county_forecast['value'].iloc[16])
        array_p800 = np.append(array_p800, county_forecast['value'].iloc[17])
        array_p850 = np.append(array_p850, county_forecast['value'].iloc[18])
        array_p900 = np.append(array_p900, county_forecast['value'].iloc[19])
        array_p950 = np.append(array_p950, county_forecast['value'].iloc[20])
        array_p975 = np.append(array_p975, county_forecast['value'].iloc[21])
        array_p990 = np.append(array_p990, county_forecast['value'].iloc[22])
        obs_cases = np.append(obs_cases, df_county_wnv['COUNT'][df_county_wnv['fips'].astype(int) == int(location)].values[0])
    quantile_dict_test = {
        0.01: np.log(array_p010 + 1),
        0.025: np.log(array_p025 + 1),
        0.05: np.log(array_p050 + 1),
        0.1: np.log(array_p100 + 1),
        0.15: np.log(array_p150 + 1),
        0.2: np.log(array_p200 + 1),
        0.25: np.log(array_p250 + 1),
        0.3: np.log(array_p300 + 1),
        0.35: np.log(array_p350 + 1),
        0.4: np.log(array_p400 + 1),
        0.45: np.log(array_p450 + 1),
        # 0.5: np.log(array_p500 + 1),
        0.55: np.log(array_p550 + 1),
        0.6: np.log(array_p600 + 1),
        0.65: np.log(array_p650 + 1),
        0.7: np.log(array_p700 + 1),
        0.75: np.log(array_p750 + 1),
        0.8: np.log(array_p800 + 1),
        0.85: np.log(array_p850 + 1),
        0.9: np.log(array_p900 + 1),
        0.95: np.log(array_p950 + 1),
        0.975: np.log(array_p975 + 1),
        0.99: np.log(array_p990 + 1),
    }
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    scores, sharpness, calibrations, overprediction, underprediction = scoring_2.weighted_interval_score_fast(np.log(obs_cases + 1), alphas=alphas, medians=np.log(array_p500 + 1), weights=None, q_dict=quantile_dict_test)
    print(file_name)
    print('Overall score:')
    print(np.round(scores.mean(), 4))
    print('Calibration score:')
    print(np.round(calibrations.mean(), 4))
    print('Sharpness score:')
    print(np.round(sharpness.mean(), 4))
    d = {'scores': scores, 'sharpness': sharpness, 'calibration': calibrations, 'overprediction': overprediction, 'underprediction': underprediction}  # creating data dictionary
    df_scoring = pd.DataFrame(data = d, index = location_list.astype(int).astype(str))  # creating dataframe for saving
    df_county_sub = df_county_wnv[['location', 'fips']]
    df_county_sub['fips'] = df_county_sub['fips'].astype(str)
    df_scoring = df_scoring.merge(df_county_sub, how='left', left_on=df_scoring.index, right_on='fips')
    df_scoring = df_scoring[['location', 'fips', 'scores', 'sharpness', 'calibration', 'overprediction', 'underprediction']]
    df_scoring.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/re_county_wis_scores_' + num_region + '.csv')





import glob
all_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_quantiles_*.csv')

# executing scoring code
file_list = all_csv_files.copy()
for file in file_list:
    file_name = file
    domain_num = int(file.split('/')[-1].split('_')[-1].split('.')[0])
    print(file)
    df_forecast = pd.read_csv(file_name, index_col = 0)
    df_forecast = df_forecast[df_forecast['region']==domain_num]
    df_forecast = df_forecast[df_forecast['year']==2022]
    # looping through all years
    array_p010 = np.array([])
    array_p025 = np.array([])
    array_p050 = np.array([])
    array_p100 = np.array([])
    array_p150 = np.array([])
    array_p200 = np.array([])
    array_p250 = np.array([])
    array_p300 = np.array([])
    array_p350 = np.array([])
    array_p400 = np.array([])
    array_p450 = np.array([])
    array_p500 = np.array([])
    array_p550 = np.array([])
    array_p600 = np.array([])
    array_p650 = np.array([])
    array_p700 = np.array([])
    array_p750 = np.array([])
    array_p800 = np.array([])
    array_p850 = np.array([])
    array_p900 = np.array([])
    array_p950 = np.array([])
    array_p975 = np.array([])
    array_p990 = np.array([])
    obs_cases = np.array([])
    year_list = df_forecast['year'].unique()
    for i, year in enumerate(year_list):
        df_forecast_sub = df_forecast[df_forecast['year'] == year]
        array_p010 = np.append(array_p010, df_forecast_sub['value'].iloc[0])
        array_p025 = np.append(array_p025, df_forecast_sub['value'].iloc[1])
        array_p050 = np.append(array_p050, df_forecast_sub['value'].iloc[2])
        array_p100 = np.append(array_p100, df_forecast_sub['value'].iloc[3])
        array_p150 = np.append(array_p150, df_forecast_sub['value'].iloc[4])
        array_p200 = np.append(array_p200, df_forecast_sub['value'].iloc[5])
        array_p250 = np.append(array_p250, df_forecast_sub['value'].iloc[6])
        array_p300 = np.append(array_p300, df_forecast_sub['value'].iloc[7])
        array_p350 = np.append(array_p350, df_forecast_sub['value'].iloc[8])
        array_p400 = np.append(array_p400, df_forecast_sub['value'].iloc[9])
        array_p450 = np.append(array_p450, df_forecast_sub['value'].iloc[10])
        array_p500 = np.append(array_p500, df_forecast_sub['value'].iloc[11])
        array_p550 = np.append(array_p550, df_forecast_sub['value'].iloc[12])
        array_p600 = np.append(array_p600, df_forecast_sub['value'].iloc[13])
        array_p650 = np.append(array_p650, df_forecast_sub['value'].iloc[14])
        array_p700 = np.append(array_p700, df_forecast_sub['value'].iloc[15])
        array_p750 = np.append(array_p750, df_forecast_sub['value'].iloc[16])
        array_p800 = np.append(array_p800, df_forecast_sub['value'].iloc[17])
        array_p850 = np.append(array_p850, df_forecast_sub['value'].iloc[18])
        array_p900 = np.append(array_p900, df_forecast_sub['value'].iloc[19])
        array_p950 = np.append(array_p950, df_forecast_sub['value'].iloc[20])
        array_p975 = np.append(array_p975, df_forecast_sub['value'].iloc[21])
        array_p990 = np.append(array_p990, df_forecast_sub['value'].iloc[22])
        # obs_cases = np.append(obs_cases, df_county_wnv['COUNT'][df_county_wnv['location'] == location].values[0])
        obs_cases = np.array(wnv_df_ann[str(domain_num)])
    quantile_dict_test = {
        0.01: np.log(array_p010 + 1),
        0.025: np.log(array_p025 + 1),
        0.05: np.log(array_p050 + 1),
        0.1: np.log(array_p100 + 1),
        0.15: np.log(array_p150 + 1),
        0.2: np.log(array_p200 + 1),
        0.25: np.log(array_p250 + 1),
        0.3: np.log(array_p300 + 1),
        0.35: np.log(array_p350 + 1),
        0.4: np.log(array_p400 + 1),
        0.45: np.log(array_p450 + 1),
        # 0.5: np.log(array_p500 + 1),
        0.55: np.log(array_p550 + 1),
        0.6: np.log(array_p600 + 1),
        0.65: np.log(array_p650 + 1),
        0.7: np.log(array_p700 + 1),
        0.75: np.log(array_p750 + 1),
        0.8: np.log(array_p800 + 1),
        0.85: np.log(array_p850 + 1),
        0.9: np.log(array_p900 + 1),
        0.95: np.log(array_p950 + 1),
        0.975: np.log(array_p975 + 1),
        0.99: np.log(array_p990 + 1),
    }
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    scores, sharpness, calibrations, overprediction, underprediction = scoring_2.weighted_interval_score_fast(np.log(obs_cases + 1), alphas=alphas, medians=np.log(array_p500 + 1), weights=None, q_dict=quantile_dict_test)
    print(file_name)
    print('Overall score:')
    print(np.round(scores.mean(), 4))
    print('Calibration score:')
    print(np.round(calibrations.mean(), 4))
    print('Sharpness score:')
    print(np.round(sharpness.mean(), 4))
    d = {'scores': scores, 'sharpness': sharpness, 'calibration': calibrations, 'overprediction': overprediction, 'underprediction': underprediction}  # creating data dictionary
    df_scoring = pd.DataFrame(data=d, index=np.arange(2005, 2023))  # creating dataframe for saving
    df_scoring['domain_num'] = domain_num
    df_scoring['domain_name'] = domain_names[str(domain_num)]
    df_scoring = df_scoring[
        ['domain_name', 'domain_num', 'scores', 'sharpness', 'calibration', 'overprediction',
         'underprediction']]
    df_scoring.to_csv(
        # '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_wis_' + str(domain_num) + '.csv')
        '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_county_wis_' + str(
            domain_num) + '.csv')

plot_df = pd.DataFrame()
# testing for statistical significance
for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    re = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_wis_' + domain_num + '.csv')
    re = re.rename(columns={'Unnamed: 0': 'year'})
    clim_fore = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_' + domain_num + '.csv')
    clim_fore = clim_fore.rename(columns={'Unnamed: 0': 'year'})

    bootstrap_num = 1000
    year_ind = re.index.values
    arr_len = len(year_ind)
    re_vals = np.zeros([bootstrap_num, arr_len])
    re_vals[:] = np.nan
    clim_fore_vals = np.zeros([bootstrap_num, arr_len])
    clim_fore_vals[:] = np.nan
    for i in np.arange(bootstrap_num):
        bootstrap_arr = np.random.choice(re.index, arr_len, replace=True)
        re_vals[i, :] = re['scores'].iloc[bootstrap_arr]
        clim_fore_vals[i, :] = clim_fore['scores'].iloc[bootstrap_arr]
    # for i in np.arange(bootstrap_num):
    #     hist_nb_vals[i, :] = hist_nb['calibration'].sample(n=arr_len, replace=True)
    #     clim_fore_vals[i, :] = clim_fore['calibration'].sample(n=arr_len, replace=True)
    plot_df[domain_num]= np.mean(clim_fore_vals, axis=1) - np.mean(re_vals, axis=1)
    num_diff = np.sum(np.mean(clim_fore_vals, axis=1) < np.mean(re_vals, axis=1))
    print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_vals), 3)) + ',' + str(np.round(np.mean(re_vals), 3)) + ',' + str((bootstrap_num - num_diff)/bootstrap_num))

# testing an all-domain forecast for statistical significance
import glob

# gathering all climate-informed regression scores
all_clim_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_multi_*.csv')
df_clim_list = []
for filename in all_clim_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_clim_list.append(df)

clim_forecast_wis_df = pd.concat(df_clim_list, axis=0, ignore_index=True)
clim_forecast_wis_df = clim_forecast_wis_df.rename(columns={'Unnamed: 0': 'year'})

# gathering all random effects scores
all_re_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_wis_*.csv')
df_re_list = []
for filename in all_re_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_re_list.append(df)

re_wis_df = pd.concat(df_re_list, axis=0, ignore_index=True)
re_wis_df = re_wis_df.rename(columns={'Unnamed: 0': 'year'})

# sorting columns to ensure they match
clim_forecast_sorted = clim_forecast_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)
re_sorted = re_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)

# running bootstrapping over the two
bootstrap_num = 1000
arr_len = len(clim_forecast_sorted)
re_wis_vals = np.zeros([bootstrap_num, arr_len])
re_wis_vals[:] = np.nan
clim_fore_wis_vals = np.zeros([bootstrap_num, arr_len])
clim_fore_wis_vals[:] = np.nan
for i in np.arange(bootstrap_num):
    bootstrap_arr = np.random.choice(clim_forecast_sorted.index, arr_len, replace=True)
    re_wis_vals[i, :] = re_sorted['scores'].iloc[bootstrap_arr]
    clim_fore_wis_vals[i, :] = clim_forecast_sorted['scores'].iloc[bootstrap_arr]
plot_df['national'] = np.mean(clim_fore_wis_vals, axis=1) - np.mean(re_wis_vals, axis=1)
num_diff = np.sum(np.mean(clim_fore_wis_vals, axis=1) < np.mean(re_wis_vals, axis=1))
print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_wis_vals), 3)) + ',' + str(
    np.round(np.mean(re_wis_vals), 3)) + ',' + str((bootstrap_num - num_diff) / bootstrap_num))

import seaborn as sns
sns.set_style('whitegrid')
fig, ax = plt.subplots(layout='constrained')
plt.axhline(0, c='k')
ax = sns.boxplot(data=plot_df, ax=ax, whis=(5, 95), width=0.6, showfliers=False, flierprops={"marker": ".", 'markersize': 1}, boxprops=dict(alpha=0.7))
plt.xlabel('NEON Region', fontsize=10, weight='bold')
plt.ylabel('Difference in Forecast Skill (WIS)', fontsize=10, weight='bold')
x_tick_labels = ['Pacific Southwest', 'Desert Southwest', 'Great Lakes', 'Prairie Peninsula', 'Southeast', 'Southern Plains',
                 'Central Plains', 'Ozarks Complex', 'Northern Plains', 'Mid Atlantic', 'Northeast', 'National']
ax.set_xticklabels(x_tick_labels, rotation=60, ha='right')
plt.title('County-Level vs Regional Forecast', fontsize=14, weight='bold')
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/diff_in_forecast_skill_clim_re.svg')









plot_df = pd.DataFrame()
# testing for statistical significance
for domain_num in ['17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1']:
    hist_nb = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_hist_nb_wis_' + domain_num + '.csv')
    hist_nb = hist_nb.rename(columns={'Unnamed: 0': 'year'})
    # clim_fore = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_clim_regression_wis_loo_clim_' + domain_num + '.csv')
    clim_fore = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_wis_' + domain_num + '.csv')
    clim_fore = clim_fore.rename(columns={'Unnamed: 0': 'year'})

    bootstrap_num = 1000
    year_ind = hist_nb.index.values
    arr_len = len(year_ind)
    hist_nb_vals = np.zeros([bootstrap_num, arr_len])
    hist_nb_vals[:] = np.nan
    clim_fore_vals = np.zeros([bootstrap_num, arr_len])
    clim_fore_vals[:] = np.nan
    for i in np.arange(bootstrap_num):
        bootstrap_arr = np.random.choice(hist_nb.index, arr_len, replace=True)
        hist_nb_vals[i, :] = hist_nb['scores'].iloc[bootstrap_arr]
        clim_fore_vals[i, :] = clim_fore['scores'].iloc[bootstrap_arr]
    # for i in np.arange(bootstrap_num):
    #     hist_nb_vals[i, :] = hist_nb['calibration'].sample(n=arr_len, replace=True)
    #     clim_fore_vals[i, :] = clim_fore['calibration'].sample(n=arr_len, replace=True)
    plot_df[domain_num]= np.mean(clim_fore_vals, axis=1) - np.mean(hist_nb_vals, axis=1)
    num_diff = np.sum(np.mean(clim_fore_vals, axis=1) < np.mean(hist_nb_vals, axis=1))
    print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_vals), 3)) + ',' + str(np.round(np.mean(hist_nb_vals), 3)) + ',' + str((bootstrap_num - num_diff)/bootstrap_num))



# testing an all-domain forecast for statistical significance
import glob

# gathering all climate-informed regression scores
all_clim_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_region_wis_*.csv')
df_clim_list = []
for filename in all_clim_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_clim_list.append(df)

clim_forecast_wis_df = pd.concat(df_clim_list, axis=0, ignore_index=True)
clim_forecast_wis_df = clim_forecast_wis_df.rename(columns={'Unnamed: 0': 'year'})

# gathering all historical negative binomial scores
all_hist_csv_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/neon_domains_hist_nb_wis_*.csv')
df_hist_list = []
for filename in all_hist_csv_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_hist_list.append(df)

hist_nb_wis_df = pd.concat(df_hist_list, axis=0, ignore_index=True)
hist_nb_wis_df = hist_nb_wis_df.rename(columns={'Unnamed: 0': 'year'})

# sorting columns to ensure they match
clim_forecast_sorted = clim_forecast_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)
hist_nb_sorted = hist_nb_wis_df.sort_values(['domain_name', 'year'], ascending=[True, True]).reset_index(drop=True)

# running bootstrapping over the two
bootstrap_num = 1000
arr_len = len(clim_forecast_sorted)
hist_nb_wis_vals = np.zeros([bootstrap_num, arr_len])
hist_nb_wis_vals[:] = np.nan
clim_fore_wis_vals = np.zeros([bootstrap_num, arr_len])
clim_fore_wis_vals[:] = np.nan
for i in np.arange(bootstrap_num):
    bootstrap_arr = np.random.choice(clim_forecast_sorted.index, arr_len, replace=True)
    hist_nb_wis_vals[i, :] = hist_nb_sorted['scores'].iloc[bootstrap_arr]
    clim_fore_wis_vals[i, :] = clim_forecast_sorted['scores'].iloc[bootstrap_arr]
plot_df['national'] = np.mean(clim_fore_wis_vals, axis=1) - np.mean(hist_nb_wis_vals, axis=1)
num_diff = np.sum(np.mean(clim_fore_wis_vals, axis=1) < np.mean(hist_nb_wis_vals, axis=1))
print(domain_num + ',' + domain_names[domain_num] + ',' + str(np.round(np.mean(clim_fore_wis_vals), 3)) + ',' + str(
    np.round(np.mean(hist_nb_wis_vals), 3)) + ',' + str((bootstrap_num - num_diff) / bootstrap_num))

import seaborn as sns
sns.set_style('whitegrid')
fig, ax = plt.subplots(layout='constrained')
plt.axhline(0, c='k')
ax = sns.boxplot(data=plot_df, ax=ax, whis=(5, 95), width=0.6, showfliers=False, flierprops={"marker": ".", 'markersize': 1}, boxprops=dict(alpha=0.7))
plt.xlabel('NEON Region', fontsize=10, weight='bold')
plt.ylabel('Difference in Forecast Skill (WIS)', fontsize=10, weight='bold')
x_tick_labels = ['Pacific Southwest', 'Desert Southwest', 'Great Lakes', 'Prairie Peninsula', 'Southeast', 'Southern Plains',
                 'Central Plains', 'Ozarks Complex', 'Northern Plains', 'Mid Atlantic', 'Northeast', 'National']
ax.set_xticklabels(x_tick_labels, rotation=60, ha='right')
plt.title('Random Effects Climate-Informed Forecast vs Historical Benchmark', fontsize=14, weight='bold')
plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/diff_in_forecast_skill_re_clim.svg')




# ok, let's score county-level forecast vs historical negative binomial
def bootstrap_means_wis(a1, a2, bootstrap_num = 1000):
    a_diff = np.zeros([bootstrap_num]); a_diff[:] = np.nan
    all_locs = a1.location.unique()
    num_locs = np.shape(all_locs)[0]
    for k in np.arange(bootstrap_num):
        samp_locs = np.random.choice(all_locs, num_locs)
        a1_samp = a1[a1['location'].isin(samp_locs)]['wis_score'].values
        a2_samp = a2[a2['location'].isin(samp_locs)]['wis_score'].values
        a_diff[k] = a1_samp.mean() - a2_samp.mean()
    return np.sum(a_diff > 0)

## starting with WIS

# loading and prepping WIS results
all_re_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_re_county_wis_*.csv')
df_re_list = []
for filename in all_re_files:
    df = pd.read_csv(filename, index_col=0, header=0)
    df_re_list.append(df)

re_wis_df = pd.concat(df_re_list, axis=0)
re_wis_df = hist_nb_wis_df.rename(columns={'Unnamed: 0': 'year'})

random_effects_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal/results/merged_WIS_log_scores.csv',index_col=0)
merged_score_df = pd.read_csv('/Users/ryan.harp/Documents/2022_WNV_Forecasting_Challenge/results/merged_WIS_log_scores.csv',index_col=0)
scores_sub = merged_score_df[merged_score_df['forecast_date']=='2022-04-30'].copy()
df_wis_results = pd.DataFrame(columns=['model1', 'model2', 'p-value-sub'])
df_wis = scores_sub.rename(columns={'model': 'team'})  # when using janky wis conversion
wis_models = df_wis['team'].unique()
num_wis_models = len(wis_models)
c = 0

# calculating differences between each pair of models
for i in np.arange(num_wis_models):
    for j in np.arange(i+1, num_wis_models):
        c += 1
        print(str(c) + ': ' + wis_models[i] + ' + ' + wis_models[j])
        df_wis_results.loc[len(df_wis_results.index)] = [wis_models[i], wis_models[j], bootstrap_means_wis(df_wis[df_wis['team'] == wis_models[i]], df_wis[df_wis['team'] == wis_models[j]])]

# testing for statistical significance using Holm's method of adjusting for family-wise errors
df_wis_results.loc[df_wis_results['p-value-sub'] > 500, 'p-value-sub'] = 1000 - df_wis_results[df_wis_results['p-value-sub'] > 500]['p-value-sub']
df_wis_results['p-value-sub'] = df_wis_results['p-value-sub']/1000
df_wis_results = df_wis_results.sort_values('p-value-sub', ascending=False)
num_comparisons = np.shape(df_wis_results)[0]
alpha = np.array([0.1, 0.05, 0.01])



all_re_files = glob.glob('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/re_county_wis_scores_*.csv')
df_re_list = []
for filename in all_re_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df_re_list.append(df)

re_wis_df = pd.concat(df_re_list, axis=0, ignore_index=True)
re_wis_df = re_wis_df.rename(columns={'Unnamed: 0': 'year'})


hist_nb_df = pd.read_csv('/Users/ryan.harp/Documents/2022_WNV_Forecasting_Challenge/results/merged_wis_log_scores.csv', index_col=0)
hist_nb_df = hist_nb_df[hist_nb_df['model']=='CDC-HistNB']
# hist_nb_df = hist_nb_df[hist_nb_df['fips'].isin(re_wis_df['fips'])]
hist_nb_df = hist_nb_df[hist_nb_df['fips'].isin(fips.astype(int))]


wnv_ann = wnv_data.resample('Y').sum()

# identifying highest case load counties
total_threshold = 50  # set percentage of caseload to be considered
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_sums = wnv_data.sum(axis=0)
wnv_sums = wnv_sums.sort_values(ascending=False)
wnv_percent = wnv_sums/wnv_sums.sum()*100
wnv_cum_percent = np.divide(wnv_sums.cumsum(), wnv_sums.sum())*100

# different ways to subset results
fips = wnv_sums[wnv_cum_percent <= total_threshold].index  # returns highest count counties which collectively contribute the total_threshold amount of WNV
# fips = wnv_sums[wnv_sums >= 1].index
# fips = wnv_sums[(wnv_sums >= 19) & ~(wnv_cum_percent <= total_threshold)].index
# fips = wnv_sums[(wnv_sums <= 19) & (wnv_sums > 0)].index
# fips = wnv_sums[fips_2022].index
