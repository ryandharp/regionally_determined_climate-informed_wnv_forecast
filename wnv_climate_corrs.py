#%% Creating hex- and county-level datasets of climate factors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')
import geopandas as gpd
import xarray as xr
import seaborn as sns


# reading in shapefiles
hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
# county_shapefiles = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_shapefiles/tl_2022_us_county.shp')

# reading in data
# temp_grid = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_tmean_198101-202211.nc')
wnv_data = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv')
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_ann = wnv_data.resample('YS').sum()

temp_data = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_tmean_198101-202211_gridded-2.csv')
temp_data.index = temp_data['Unnamed: 0'].values
temp_data.index = pd.to_datetime(temp_data.index)
temp_data = temp_data.drop('Unnamed: 0', axis=1)
temp_data = temp_data[temp_data.index >= '2003-01-01']

pdsi_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/gridMET_PDSI_199901-202212_gridded-2.csv')
pdsi_data.index = pdsi_data['Unnamed: 0'].values
pdsi_data.index = pd.to_datetime(pdsi_data.index)
pdsi_data = pdsi_data.drop('Unnamed: 0', axis=1)
pdsi_data = pdsi_data[pdsi_data.index >= '2003-01-01']

precip_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_precip_198101-202211_gridded-2.csv')
precip_data.index = precip_data['Unnamed: 0'].values
precip_data.index = pd.to_datetime(precip_data.index)
precip_data = precip_data.drop('Unnamed: 0', axis=1)
precip_data = precip_data[precip_data.index >= '2003-01-01']

def load_clim_data(fname):
    clim_df = pd.read_csv(fname)
    clim_df.index = clim_df['Unnamed: 0'].values
    clim_df.index = pd.to_datetime(clim_df.index)
    clim_df = clim_df.drop('Unnamed: 0', axis=1)
    clim_df = clim_df[clim_df.index >= '2003-01-01']
    return clim_df

# sorting hexes by highest caseload
wnv_ann_sum  = wnv_ann.sum().sort_values(ascending=False)
wnv_per = wnv_ann_sum/wnv_ann_sum.sum()
wnv_cum_per = wnv_per.cumsum()
# 9 hexes to get to 50%
# 19 hexes to get to 67%
# 26 hexes to get to 75%
# 29 hexes contribute at least 1% of cases
# 31 hexes to get to 80%
# 45 hexes to get to 90%


hex = 27  # Colorado hex for example. 7th highest caseload of all hexes at 758
month = 3
lag = 0
num_hexes = 19
hex_locations = {27: 'Los Angeles', 47: 'Maricopa', 131: 'Chicago', 106: 'Dallas', 30: 'Sacramento-Central Valley',
                 105: 'Houston', 90: 'Front Range', 190: 'New York Metro', 126: 'Louisiana', 151: 'NW Ohio',
                 125: 'New Orleans-Gulf Shore', 29: 'Fresno-Central Valley', 111: 'Kansas City', 112: 'W Nebraska-SD',
                 127: 'N Mississippi-Alabama', 192: 'Albany', 171: 'Central Pennsylvania', 108: 'OKC',
                 114: 'W ND-SD', 152: 'Grand Rapids', 87: 'Texas Panhandle', 109: 'Tulsa-Wichita-Springfield',
                 65: 'El Paso', 66: 'Tucson', 107: 'Texarkana', 69: 'S Colorado', 130: 'Columbia-MO',
                 150: 'Indianapolis', 104: 'San Antonio', 169: 'Richmond-VA', 52: 'SW Idaho'}

def climate_comparison(wnv_df, clim_df, hex, month):
    wnv_ts = wnv_df[str(hex)]
    clim_ts = clim_df[str(hex)]
    mon_clim_ts = clim_ts[clim_ts.index.month == month]
    if month < 10:
        lagged_mon_clim_ts = mon_clim_ts[1:]
    elif month >= 10:
        lagged_mon_clim_ts = mon_clim_ts[:-1]
    elif month == 12:
        lagged_mon_clim_ts = mon_clim_ts[:]
    wnv_year = wnv_ts.index
    wnv_ts = wnv_ts.reset_index(drop=True)
    lagged_mon_clim_ts = lagged_mon_clim_ts.reset_index(drop=True)
    comparison_df = pd.DataFrame({'wnv': wnv_ts, 'lagged_clim': lagged_mon_clim_ts, 'wnv_year': wnv_year})
    return comparison_df

fig, axes = plt.subplots(ncols=4, nrows=3, layout='constrained')
for m, ax in zip(np.arange(1, 13, 1), axes.flat):
    # comparison_df = climate_comparison(wnv_ann, clim_data, hex, m, 0)
    comparison_df = climate_comparison(wnv_ann, wnv_data, hex, m)
    sns.scatterplot(comparison_df, x='lagged_clim', y='wnv', hue = 'wnv_year', ax=ax, legend=False, palette=sns.color_palette('mako_r', n_colors=19))
    ax.title.set_text('mon: ' + str(m) + ', R=' + str(np.round(np.corrcoef(comparison_df['lagged_clim'], comparison_df['wnv'])[0,1], 2)))
    plt.show()


log = False
index_iterables = [
    wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index,
    np.arange(1, 13, 1)
]
multiindex = pd.MultiIndex.from_product(index_iterables, names=["hex", "month"])

clim_vars = ['Temperature', 'Drought', 'Precipitation', 'Snow Depth', 'Snow Melt', 'SMWE', 'Soil Moisture 5', 'Soil Moisture 50', 'Soil Moisture 100', 'Soil Temp 5', 'VPD Max', 'VPD Min']
results_df = pd.DataFrame(columns = clim_vars, index=multiindex)
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
matplotlib.use('AGG')
for clim_var in clim_vars:
    if clim_var == 'Temperature':
        clim_data = load_clim_data(fpath + 'PRISM_tmean_198101-202211_gridded-2.csv')
    elif clim_var == 'Drought':
        clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202212_gridded-2.csv')
    elif clim_var == 'Precipitation':
        clim_data = load_clim_data(fpath + 'PRISM_precip_198101-202211_gridded-2.csv')
    elif clim_var == 'Snow Depth':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Snow Melt':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_melt_199901-202212_gridded-2.csv')
    elif clim_var == 'SMWE':
        clim_data = load_clim_data(fpath + 'NLDAS_water_equivalent_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 50':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 100':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_100_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Temp 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_temp_5_199901-202212_gridded-2.csv')
    elif clim_var == 'VPD Max':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmax_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    elif clim_var == 'VPD Min':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmin_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    for hex in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
        fig, axes = plt.subplots(ncols=4, nrows=3, layout='constrained', figsize=(10, 8))
        hex_corrs = np.zeros([12]); hex_corrs[:] = np.nan
        i = 0
        for m, ax in zip(np.arange(1, 13, 1), axes.flat):
            comparison_df = climate_comparison(wnv_ann, clim_data, int(hex), m)
            mask = ~comparison_df.isna().any(axis=1)
            if log:
                sns.scatterplot(comparison_df, x='lagged_clim', y=np.log10(comparison_df['wnv']+1), hue='wnv_year', ax=ax,
                                legend=False, palette=sns.color_palette('mako_r', n_colors=19))
                # sns.scatterplot(comparison_df, x=np.log10(comparison_df['lagged_clim']+1), y=np.log10(comparison_df['wnv']+1), hue='wnv_year', ax=ax,
                #                 legend=False, palette=sns.color_palette('mako_r', n_colors=19))
                ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                    np.round(np.corrcoef(comparison_df['lagged_clim'][mask], np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)))
                # ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                #     np.round(np.corrcoef(np.log10(comparison_df['lagged_clim'][mask]+1), np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)))
                fig.suptitle('WNV (log) and Annual ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                             fontsize=20)
                hex_corrs[i] = np.round(np.corrcoef(comparison_df['lagged_clim'][mask], np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)
            else:
                sns.scatterplot(comparison_df, x='lagged_clim', y='wnv', hue='wnv_year', ax=ax, legend=False,
                            palette=sns.color_palette('mako_r', n_colors=19))
                ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                    np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask])[0, 1], 2)))
                fig.suptitle('WNV and Annual ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                             fontsize=20)
                hex_corrs[i] = np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask])[0, 1], 2)
            i += 1
            results_df.loc[hex, clim_var] = hex_corrs
        if log:
            plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/annual_log_wnv_monthly_' + clim_var.lower().replace(' ', '-') + '_' + hex + '_' + hex_locations[int(hex)].replace(' ', '-').lower())
        else:
            plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/annual_wnv_monthly_' + clim_var.lower().replace(' ', '-') + '_' + hex + '_' + hex_locations[int(hex)].replace(' ', '-').lower())
        plt.close(fig)
results_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/annual_wnv_monthly_clim_corr.csv')



#
#
#%% Looking at relationships between climatic variables and timing
#
#
index_iterables = [
    wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index,
    np.arange(1, 13, 1)
]
multiindex = pd.MultiIndex.from_product(index_iterables, names=["hex", "month"])

clim_vars = ['Temperature', 'Drought', 'Precipitation', 'Snow Depth', 'Snow Melt', 'SMWE', 'Soil Moisture 5', 'Soil Moisture 50', 'Soil Moisture 100', 'Soil Temp 5', 'VPD Max', 'VPD Min']
results_df = pd.DataFrame(columns = clim_vars, index=multiindex)
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
matplotlib.use('AGG')
wnv_timing = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/wnv_season_center.csv')
for clim_var in clim_vars:
    if clim_var == 'Temperature':
        clim_data = load_clim_data(fpath + 'PRISM_tmean_198101-202211_gridded-2.csv')
    elif clim_var == 'Drought':
        clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202212_gridded-2.csv')
    elif clim_var == 'Precipitation':
        clim_data = load_clim_data(fpath + 'PRISM_precip_198101-202211_gridded-2.csv')
    elif clim_var == 'Snow Depth':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Snow Melt':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_melt_199901-202212_gridded-2.csv')
    elif clim_var == 'SMWE':
        clim_data = load_clim_data(fpath + 'NLDAS_water_equivalent_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 50':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 100':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_100_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Temp 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_temp_5_199901-202212_gridded-2.csv')
    elif clim_var == 'VPD Max':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmax_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    elif clim_var == 'VPD Min':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmin_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    for hex in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
        fig, axes = plt.subplots(ncols=4, nrows=3, layout='constrained', figsize=(10, 8))
        hex_corrs = np.zeros([12]); hex_corrs[:] = np.nan
        i = 0
        for m, ax in zip(np.arange(1, 13, 1), axes.flat):
            comparison_df = climate_comparison(wnv_ann, clim_data, int(hex), m)
            comparison_df['wnv_timing'] = wnv_timing.loc[wnv_timing['hex']==int(hex), 'timing_center'].values
            comparison_df = comparison_df.drop(columns='wnv')
            mask = ~comparison_df.isna().any(axis=1)
            hex_corrs[i] = np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv_timing'][mask])[0, 1], 2)
            sns.scatterplot(comparison_df, x='lagged_clim', y='wnv_timing', hue='wnv_year', ax=ax, legend=False,
                            palette=sns.color_palette('mako_r', n_colors=19))
            ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv_timing'][mask])[0, 1], 2)))
            fig.suptitle('WNV Timing and Monthly ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                         fontsize=20)
            i += 1
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/timing/annual_wnv_timing_monthly_' + clim_var.lower().replace(
                ' ', '-') + '_' + hex + '_' + hex_locations[int(hex)].replace(' ', '-').lower())
        results_df.loc[hex, clim_var] = hex_corrs
results_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/annual_wnv_timing_monthly_clim_corr.csv')




#
#
#%% Testing out 3D scatter plots to look for clustering?
#
#

test = climate_comparison(wnv_ann, pdsi_data, hex, 4)
test2 = climate_comparison(wnv_ann, precip_data, hex, 4)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(test['lagged_clim'], test2['lagged_clim'], test['wnv'])


#
#
#%% Plotting monthly climatic influences on monthly WNV cases
#
#

def climate_comparison_monthly(wnv_df, clim_df, hex, month, lag):
    wnv_ts = wnv_df[str(hex)]
    clim_ts = clim_df[str(hex)]
    mon_clim_ts = clim_ts[clim_ts.index.month == month]
    wnv_year = wnv_ts.index
    wnv_ts = wnv_ts.reset_index(drop=True)
    if lag:
        lagged_mon_clim_ts = mon_clim_ts[:-lag]
        # lagged_mon_clim_ts = mon_clim_ts[:-1]
        # wnv_ts = wnv_ts[1:].reset_index(drop=True)
        # wnv_year = wnv_year[1:]
    else:
        lagged_mon_clim_ts = mon_clim_ts[1:]
        # lagged_mon_clim_ts = mon_clim_ts
    lagged_mon_clim_ts = lagged_mon_clim_ts.reset_index(drop=True)
    comparison_df = pd.DataFrame({'wnv': wnv_ts, 'lagged_clim': lagged_mon_clim_ts, 'wnv_year': wnv_year})
    return comparison_df


month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}


log = True
results_df = pd.DataFrame(columns = ['hex', 'WNV_month', 'climate_month', 'climate_var', 'correlation'])
clim_vars = ['Temperature', 'Drought', 'Precipitation', 'Snow Depth', 'Snow Melt', 'SMWE', 'Soil Moisture 5', 'Soil Moisture 50', 'Soil Moisture 100', 'Soil Temp 5', 'VPD Max', 'VPD Min']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
matplotlib.use('AGG')
start_time = time.time()
for clim_var in clim_vars:
    if clim_var == 'Temperature':
        clim_data = load_clim_data(fpath + 'PRISM_tmean_198101-202211_gridded-2.csv')
    elif clim_var == 'Drought':
        clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202212_gridded-2.csv')
    elif clim_var == 'Precipitation':
        clim_data = load_clim_data(fpath + 'PRISM_precip_198101-202211_gridded-2.csv')
    elif clim_var == 'Snow Depth':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Snow Melt':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_melt_199901-202212_gridded-2.csv')
    elif clim_var == 'SMWE':
        clim_data = load_clim_data(fpath + 'NLDAS_water_equivalent_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 50':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 100':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_100_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Temp 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_temp_5_199901-202212_gridded-2.csv')
    elif clim_var == 'VPD Max':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmax_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    elif clim_var == 'VPD Min':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmin_198101-202303_gridded-2.csv')
        clim_data = clim_data[clim_data.index <= '2022-12-01']
    for hex in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
        print(hex + ' ' + clim_var + ': ' + str(time.time() - start_time))
        for wnv_month in np.arange(1, 13, 1):
            print(month_names[wnv_month])
            # fig, axes = plt.subplots(ncols=4, nrows=3, layout='constrained', figsize=(10, 8))
            mon_wnv_data = wnv_data[wnv_data.index.month == wnv_month]
            lag = 0
            for m, ax in zip(np.arange(wnv_month, wnv_month - 12, -1), axes.flat):
                if m < 1:
                    m += 12
                    lag = 1
                if mon_wnv_data[str(hex)].sum() == 0:
                    results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.nan]
                    # if log:
                        # plt.savefig(
                        #     '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_' +
                        #         str(wnv_month) + '_log_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' +
                        #     hex_locations[int(hex)].replace(' ', '-').lower())
                    # else:
                    #     plt.savefig(
                    #         '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_' +
                    #             str(wnv_month) + '_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' +
                    #         hex_locations[int(hex)].replace(' ', '-').lower())
                    # plt.close(fig)
                    # continue
                comparison_df = climate_comparison_monthly(mon_wnv_data, clim_data, int(hex), m, lag)
                mask = ~comparison_df.isna().any(axis=1)
                if log:
                    # sns.scatterplot(comparison_df, x='lagged_clim', y=np.log10(comparison_df['wnv']+1), hue='wnv_year', ax=ax,
                    #                 legend=False, palette=sns.color_palette('mako_r', n_colors=19))
                    # ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                    #     np.round(np.corrcoef(comparison_df['lagged_clim'][mask], np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)))
                    # fig.suptitle(month_names[wnv_month] + ' WNV (log) and Monthly ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                    #              fontsize=20)
                    results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.round(np.corrcoef(comparison_df['lagged_clim'][mask], np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)]
                # else:
                #     sns.scatterplot(comparison_df, x='lagged_clim', y='wnv', hue='wnv_year', ax=ax, legend=False,
                #                 palette=sns.color_palette('mako_r', n_colors=19))
                #     ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                #         np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask])[0, 1], 2)))
                #     fig.suptitle(month_names[wnv_month] + ' WNV and Monthly ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                #                  fontsize=20)
                #     results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask]+1)[0, 1], 2)]
                # if log:
                    # plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_' + str(wnv_month) + '_log_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' + hex_locations[int(hex)].replace(' ', '-').lower())
                # else:
                #     plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_'  + str(wnv_month) + '_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' + hex_locations[int(hex)].replace(' ', '-').lower())
            # plt.close(fig)

results_df.to_csv(
            '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/monthly_log_wnv_monthly_clim_corr.csv', na_rep='nan')



#
#
#%% Plotting monthly climatic influences on monthly WNV cases
#
#

month_names = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August',
               9: 'September', 10: 'October', 11: 'November', 12: 'December'}

import warnings
warnings.filterwarnings('ignore', category=UserWarning)  # tread lightly but enacting this to ignore palette length warnings

hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')
# county_shapefiles = gpd.read_file('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_shapefiles/tl_2022_us_county.shp')

# reading in data
# temp_grid = xr.open_dataset('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/gridded_climate/PRISM_tmean_198101-202211.nc')
wnv_data = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv')
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_ann = wnv_data.resample('YS').sum()

log = 'True'
results_df = pd.DataFrame(columns = ['hex', 'WNV_month', 'climate_month', 'climate_var', 'correlation'])
clim_vars = ['Temperature', 'Drought', 'Precipitation', 'Snow Depth', 'Snow Melt', 'SMWE', 'Soil Moisture 5', 'Soil Moisture 50', 'Soil Moisture 100', 'Soil Temp 5']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
matplotlib.use('AGG')
clim_data = wnv_data
start_time = time.time()
for hex in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
    print(hex + ' ' + clim_var + ': ' + str(time.time() - start_time))
    for wnv_month in np.arange(1, 13, 1):
        print(month_names[wnv_month])
        # fig, axes = plt.subplots(ncols=4, nrows=3, layout='constrained', figsize=(10, 8))
        mon_wnv_data = wnv_data[wnv_data.index.month == wnv_month]
        lag = 0
        for m, ax in zip(np.arange(wnv_month, wnv_month - 12, -1), axes.flat):
            if m < 1:
                m += 12
                lag = 1
            if mon_wnv_data[str(hex)].sum() == 0:
                results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.nan]
                continue
            comparison_df = climate_comparison_monthly(mon_wnv_data, clim_data, int(hex), m, lag)
            mask = ~comparison_df.isna().any(axis=1)
            if log:
                # sns.scatterplot(comparison_df, x=np.log10(comparison_df['lagged_clim']+1), y=np.log10(comparison_df['wnv']+1), hue='wnv_year', ax=ax,
                #                 legend=False, palette=sns.color_palette('mako_r', n_colors=19))
                # ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                #     np.round(np.corrcoef(np.log10(comparison_df['lagged_clim'][mask]+1), np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)))
                # fig.suptitle(month_names[wnv_month] + ' WNV (log) and Monthly ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                #              fontsize=20)
                # results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.round(np.corrcoef(np.log10(comparison_df['lagged_clim'][mask]+1), np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)]
                results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.round(np.corrcoef(comparison_df['lagged_clim'][mask], np.log10(comparison_df['wnv'][mask]+1))[0, 1], 2)]
            else:
                # sns.scatterplot(comparison_df, x='lagged_clim', y='wnv', hue='wnv_year', ax=ax, legend=False,
                #             palette=sns.color_palette('mako_r', n_colors=19))
                # ax.title.set_text('mon: ' + str(m) + ', R=' + str(
                #     np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask])[0, 1], 2)))
                # fig.suptitle(month_names[wnv_month] + ' WNV and Monthly ' + clim_var + ': ' + hex_locations[int(hex)] + ' (hex ' + str(hex) + ')',
                #              fontsize=20)
                results_df.loc[len(results_df)] = [hex, wnv_month, m, clim_var, np.round(np.corrcoef(comparison_df['lagged_clim'][mask], comparison_df['wnv'][mask]+1)[0, 1], 2)]
        #     if log:
        #         plt.savefig(
        #             '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_' + str(
        #                 wnv_month) + '_log_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' +
        #             hex_locations[int(hex)].replace(' ', '-').lower())
        #     else:
        #         plt.savefig(
        #             '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/exploratory/mon_' +
        #                 str(wnv_month) + '_wnv_monthly_' + clim_var.lower() + '_' + hex + '_' + hex_locations[
        #                 int(hex)].replace(' ', '-').lower())
        # plt.close(fig)

results_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/monthly_log_wnv_monthly_clim_corr.csv', na_rep='nan')
warnings.filterwarnings('default', category=UserWarning)  # reactivating warnings

