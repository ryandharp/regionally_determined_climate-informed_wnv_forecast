# This script reads in a previously created hexagonal grid and assigns counts of WNV to the respective hexes.

import numpy as np
import pandas as pd
import geopandas as gpd


# prepping for WNV data assignment
# loading in hex grid
gdf = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_' + str(side_len) + '_deg.shp')
bounds = gdf.bounds
gdf = gdf.join(bounds)


# loading in WNV data and supporting county info
df_wnv = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/NeuroWNV_by_county_2000-2022_monthly.csv')
df_pop_cens = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_centers_of_pop.txt', dtype={'STATEFP': str, 'COUNTYFP': str})
df_pop_cens['FIPS'] = df_pop_cens['STATEFP'] + df_pop_cens['COUNTYFP']
excluded_states = ['02', '15', '72']
df_pop_cens = df_pop_cens[~df_pop_cens['STATEFP'].isin(excluded_states)]
df_pop_cens = df_pop_cens.reset_index()
df_pop_cens = df_pop_cens[['LATITUDE', 'LONGITUDE', 'FIPS']]
df_pop_cens['hex'] = np.nan

# identifying which hex a given county falls within
def find_which_shapefile(lon, lat, shapefile_list):
    pt = Point([lon, lat])
    for i in np.arange(len(shapefile_list)):
        if pt.within(shapefile_list[i]):
            break
    return i

start_time = time.time()
for i, row in df_pop_cens.iterrows():
    if i % 100 == 0:
        print('Index: ' + str(i) +  'at ' + str(np.round(time.time() - start_time, 1)) + ' seconds.')
    lat = df_pop_cens['LATITUDE'].iloc[i]
    lon = df_pop_cens['LONGITUDE'].iloc[i]
    possible_hexes = (lat > gdf['miny']) & (lat < gdf['maxy']) & (lon > gdf['minx']) & (lon < gdf['maxx'])
    if np.sum(possible_hexes) == 1:
        df_pop_cens.loc[df_pop_cens.index==i, 'hex'] = possible_hexes[possible_hexes == True].index.values[0]
    elif np.sum(possible_hexes) == 0:
        print('Failed: ' + str(i))
        print(row)
    else:
        hex_list = gdf[possible_hexes]['grid_id']
        hex_shapefiles = list()
        for hex_ind in hex_list:
            hex_shapefiles.append(gdf[gdf['grid_id'] == hex_ind].unary_union)
        shapefile_ind = find_which_shapefile(lon, lat, hex_shapefiles)
        df_pop_cens.loc[df_pop_cens.index==i, 'hex'] = hex_list.iloc[shapefile_ind]

print('Finished after ' + str(np.round(time.time()-start_time,1)) + ' seconds.')
df_pop_cens.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/fips_to_hex_grid_2.csv')

# merging data sets
df_pop_cens['FIPS'] = df_pop_cens['FIPS'].astype(float)
df_wnv = df_wnv.merge(df_pop_cens, how='left', left_on='County', right_on='FIPS')
df_wnv = df_wnv[['County', 'month', 'Year', 'COUNT', 'hex']]
wnv_hex_counts = df_wnv['COUNT'].groupby(df_wnv['hex']).sum()
gdf = gdf.merge(wnv_hex_counts, how='left', left_on='grid_id', right_on='hex')

## a test plot
import matplotlib
fig, ax = plt.subplots(1, 1)
gdf.plot(column='COUNT', ax=ax, legend=True, norm=matplotlib.colors.LogNorm(vmin=1, vmax=gdf.COUNT.max()))

# combining
df_wnv['date'] = pd.to_datetime(dict(year=df_wnv['Year'], month=df_wnv['month'], day=1))
time_index = pd.date_range(df_wnv['date'].min(), df_wnv['date'].max(), freq='MS')

df_hexes_wnv = pd.DataFrame(index = time_index, columns=np.sort(df_wnv['hex'].unique().astype(int)))
df_hexes_wnv = df_hexes_wnv.fillna(value=0)
for i, row in df_wnv.iterrows():
    if (np.isnan(row['month']) | np.isnan(row['hex'])):
        print(row)
        continue
    df_hexes_wnv.loc[row['date'], int(row['hex'])] = df_hexes_wnv.loc[row['date']][int(row['hex'])] + row['COUNT']

df_hexes_wnv.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_time_series_' + str(side_len) + '_deg.csv')


