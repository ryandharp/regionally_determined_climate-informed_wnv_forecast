import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')
import geopandas as gpd
import xarray as xr
import seaborn as sns

# reading in data
# temp_grid = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_tmean_198101-202211.nc')
wnv_data = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv')
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_ann = wnv_data.resample('YS').sum()

hex_locations = {27: 'Los Angeles', 47: 'Maricopa', 131: 'Chicago', 106: 'Dallas', 30: 'Sacramento-Central Valley',
                 105: 'Houston', 90: 'Front Range', 190: 'New York Metro', 126: 'Louisiana', 151: 'NW Ohio',
                 125: 'New Orleans-Gulf Shore', 29: 'Fresno-Central Valley', 111: 'Kansas City', 112: 'W Nebraska-SD',
                 127: 'N Mississippi-Alabama', 192: 'Albany', 171: 'Central Pennsylvania', 108: 'OKC',
                 114: 'W ND-SD', 152: 'Grand Rapids', 87: 'Texas Panhandle', 109: 'Tulsa-Wichita-Springfield',
                 65: 'El Paso', 66: 'Tucson', 107: 'Texarkana', 69: 'S Colorado', 130: 'Columbia-MO',
                 150: 'Indianapolis', 104: 'San Antonio', 169: 'Richmond-VA', 52: 'SW Idaho'}

df_wnv_timing_center = pd.DataFrame(columns=['hex', 'year', 'timing_center', 'num_cases'])
for hex in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
    for year in wnv_data.index.year.unique():
        hex_wnv_data = wnv_data.loc[wnv_data.index.year == year, hex]
        df_wnv_timing_center.loc[len(df_wnv_timing_center)] = [hex, year, np.sum(hex_wnv_data.index.month*hex_wnv_data)/hex_wnv_data.sum(), hex_wnv_data.sum()]

df_wnv_timing_center.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_season_center.csv')

# creating histograms of timing centers for all hexes
for hex in df_wnv_timing_center['hex'].unique():
    plt.hist(df_wnv_timing_center[df_wnv_timing_center['hex']==hex].timing_center, np.arange(5, 11, 0.25))
    plt.title('Hex ' + hex + ' (' + hex_locations[int(hex)] + ') Timing Center')
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/' + hex + '_timing_center')
    plt.close()
    plt.scatter(df_wnv_timing_center[df_wnv_timing_center['hex']==hex].timing_center.values, wnv_ann[hex])
    plt.title('Hex ' + hex + ' (' + hex_locations[int(hex)] + ') Timing vs Cases')
    plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/' + hex + '_timing_vs_cases')
    plt.close()

df_wnv[df_wnv['hex']==int(hex)].COUNT.groupby(df_wnv['Year']).sum()
