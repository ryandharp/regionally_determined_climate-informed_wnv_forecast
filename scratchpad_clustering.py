#%% Clustering with SOMs
import os
import math
# Essential Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Preprocessing
from sklearn.preprocessing import MinMaxScaler
# Algorithms
from minisom import MiniSom
from tslearn.barycenters import dtw_barycenter_averaging
from tslearn.clustering import TimeSeriesKMeans
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import geopandas as gpd

# gathering input data
hex_grid = gpd.read_file(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/updated_hex_grid_2_deg.shp')

wnv_data = pd.read_csv(
    '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/updated_hex_time_series_2_deg.csv')
wnv_data.index = wnv_data['Unnamed: 0'].values
wnv_data.index = pd.to_datetime(wnv_data.index)
wnv_data = wnv_data.drop('Unnamed: 0', axis=1)
wnv_data = wnv_data[wnv_data.index >= '2004-01-01']
wnv_ann = wnv_data.resample('YS').sum()

test = wnv_ann.sum().sort_values(ascending=False)
test_cum = test.cumsum()/test.sum()*100

num_hexes = 9
i = 0
mySeries = []
for col in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
    print(col)
    scaler = MinMaxScaler()
    mySeries.append(MinMaxScaler().fit_transform(wnv_ann[col].values.reshape(-1, 1)))
    mySeries[i]= mySeries[i].reshape(len(mySeries[i]))
    i += 1


som_x = som_y = math.ceil(math.sqrt(math.sqrt(len(mySeries))))
# I didn't see its significance but to make the map square,
# I calculated square root of map size which is
# the square root of the number of series
# for the row and column counts of som

som = MiniSom(som_x, som_y,len(mySeries[0]), sigma=0.3, learning_rate = 0.1)

som.random_weights_init(mySeries)
som.train(mySeries, 50000)

def plot_som_series_averaged_center(som_x, som_y, win_map):
    fig, axs = plt.subplots(som_x,som_y,figsize=(25,25), layout='constrained')
    fig.suptitle('Clusters')
    for x in range(som_x):
        for y in range(som_y):
            cluster = (x,y)
            if cluster in win_map.keys():
                for series in win_map[cluster]:
                    axs[cluster].plot(series,c="gray",alpha=0.5)
                axs[cluster].plot(np.average(np.vstack(win_map[cluster]),axis=0),c="red")
            cluster_number = x*som_y+y+1
            axs[cluster].set_title(f"Cluster {cluster_number}")

    plt.show()

win_map = som.win_map(mySeries)
# Returns the mapping of the winner nodes and inputs

plot_som_series_averaged_center(som_x, som_y, win_map)

for series in mySeries:
    print(som.winner(series))

cluster_map = []
for idx in range(len(mySeries)):
    winner_node = som.winner(mySeries[idx])
    cluster_map.append((mySeries[idx],f"Cluster {winner_node[0]*som_y+winner_node[1]+1}"))

# pd.DataFrame(cluster_map,columns=["Series","Cluster"]).sort_values(by="Cluster").set_index("Series")

test_clusters = pd.DataFrame(cluster_map,columns=["Series","Cluster"])
test_clusters.index = wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index

hex_locations = {27: 'Los Angeles', 47: 'Maricopa', 131: 'Chicago', 106: 'Dallas', 30: 'Sacramento-Central Valley',
                 105: 'Houston', 90: 'Front Range', 190: 'New York Metro', 126: 'Louisiana', 151: 'NW Ohio',
                 125: 'New Orleans-Gulf Shore', 29: 'Fresno-Central Valley', 111: 'Kansas City', 112: 'W Nebraska-SD',
                 127: 'N Mississippi-Alabama', 192: 'Albany', 171: 'Central Pennsylvania', 108: 'OKC',
                 114: 'W ND-SD', 152: 'Grand Rapids', 87: 'Texas Panhandle', 109: 'Tulsa-Wichita-Springfield',
                 65: 'El Paso', 66: 'Tucson', 107: 'Texarkana', 69: 'S Colorado', 130: 'Columbia-MO',
                 150: 'Indianapolis', 104: 'San Antonio', 169: 'Richmond-VA', 52: 'SW Idaho'}


# now trying k-means clustering
cluster_count = math.ceil(math.sqrt(len(mySeries)))
cluster_count = 9
# A good rule of thumb is choosing k as the square root of the number of points in the training data set in kNN

km = TimeSeriesKMeans(n_clusters=cluster_count, metric="dtw")

labels = km.fit_predict(mySeries)

plot_count = math.ceil(math.sqrt(cluster_count))

fig, axs = plt.subplots(plot_count, plot_count, figsize=(10, 10))
fig.suptitle('Clusters')
row_i = 0
column_j = 0
# For each label there is,
# plots every series with that label
for label in set(labels):
    cluster = []
    for i in range(len(labels)):
        if (labels[i] == label):
            axs[row_i, column_j].plot(mySeries[i], c="gray", alpha=0.4)
            cluster.append(mySeries[i])
    if len(cluster) > 0:
        axs[row_i, column_j].plot(np.average(np.vstack(cluster), axis=0), c="red")
    axs[row_i, column_j].set_title("Cluster " + str(row_i * som_y + column_j))
    column_j += 1
    if column_j % plot_count == 0:
        row_i += 1
        column_j = 0

plt.show()

cluster_labels = pd.DataFrame(zip(labels, wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index), columns=['Cluster', 'hex']).sort_values(by='Cluster')




# trying some fuzzy clustering methods
from sklearn import datasets
import skfuzzy as fuzz

# collecting/preprocessing data
num_hexes = 9
i = 0
mySeries = []
for col in wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index:
    print(col)
    scaler = MinMaxScaler()
    mySeries.append(MinMaxScaler().fit_transform(wnv_ann[col].values.reshape(-1, 1)))
    mySeries[i]= mySeries[i].reshape(len(mySeries[i]))
    i += 1

# number of clusters
n_clusters = np.sqrt(num_hexes).round(0).astype(int)

# Fuzzy c-means algorithm
cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(mySeries), n_clusters, 2, error=0.005, maxiter=1000)

# Assign each data point to the cluster with the highest membership value
cluster_membership = np.argmax(u, axis=0)

fpc_arr = np.zeros([9]); fpc_arr[:] = np.nan
num_hexes = np.shape(wnv_ann)[1]
# ok, let's do this with some climate data to create regions for aggregation
temp_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_tmean_198101-202211_gridded-2.csv')
temp_data.index = temp_data['Unnamed: 0'].values
temp_data = temp_data.drop(columns=['Unnamed: 0'])
temp_data.index = pd.to_datetime(temp_data.index)
temp_data_sub = temp_data[temp_data.index.year <= 2021]
# temp_mon_anom = temp_data_sub.groupby(temp_data_sub.index.month).transform(lambda x: x-x.mean())
# temp_data_sub = temp_data_sub.groupby(temp_data_sub.index.month).values() - temp_data_sub.groupby(temp_data_sub.index.month).mean()
test = temp_data_sub[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index].transpose()
# test_ann_means = temp_data_sub.groupby(temp_data_sub.index.year).mean()
# test_ann_stds = temp_data_sub.groupby(temp_data_sub.index.year).std()
cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(test), n_clusters, 2, error=0.005, maxiter=1000)
cluster_membership = np.argmax(u, axis=0)
test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')

import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
crs_new = ccrs.PlateCarree()

cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(test), n_clusters, 2, error=0.005, maxiter=1000)
cluster_membership = np.argmax(u, axis=0)
test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')

fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
ax.set_aspect('equal')
test_grid.plot(column='cluster', ax=ax, cmap='tab20c')
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=14)
plt.yticks(np.arange(30, 60, 10), fontsize=14)



def load_clim_data(fname):
    clim_df = pd.read_csv(fname)
    clim_df.index = clim_df['Unnamed: 0'].values
    clim_df.index = pd.to_datetime(clim_df.index)
    clim_df = clim_df.drop('Unnamed: 0', axis=1)
    # clim_df = clim_df[clim_df.index <= '2022-01-01']
    return clim_df

# clim_vars = ['Temperature', 'Drought', 'Precipitation', 'Snow Depth', 'Snow Melt', 'SMWE', 'Soil Moisture 5', 'Soil Moisture 50', 'Soil Moisture 100', 'Soil Temp 5', 'VPD Max', 'VPD Min']
fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
clim_vars = ['temp_mean', 'temp_max', 'temp_min', 'temp_dew', 'drought', 'precipitation', 'soil_moisture_50', 'snow_depth', 'VPD max']  # add in something with moisture?
summary_arr = np.zeros([2*len(clim_vars), num_hexes])
for index, clim_var in np.ndenumerate(clim_vars):
    if clim_var == 'temp_mean':
        clim_data = load_clim_data(fpath + 'PRISM_tmean_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_max':
        clim_data = load_clim_data(fpath + 'PRISM_tmax_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_min':
        clim_data = load_clim_data(fpath + 'PRISM_tmin_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_dew':
        clim_data = load_clim_data(fpath + 'PRISM_td_198101-202211_gridded-2.csv')
    elif clim_var == 'drought':
        clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202212_gridded-2.csv')
    elif clim_var == 'precipitation':
        clim_data = load_clim_data(fpath + 'PRISM_precip_198101-202211_gridded-2.csv')
    elif clim_var == 'snow_depth':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Snow Melt':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_melt_199901-202212_gridded-2.csv')
    elif clim_var == 'SMWE':
        clim_data = load_clim_data(fpath + 'NLDAS_water_equivalent_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')
    elif clim_var == 'soil_moisture_50':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 100':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_100_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Temp 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_temp_5_199901-202212_gridded-2.csv')
    elif clim_var == 'VPD_max':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmax_198101-202303_gridded-2.csv')
    elif clim_var == 'VPD_min':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmin_198101-202303_gridded-2.csv')

    clim_df = clim_data[clim_data.index < '2022-01-01']
    clim_df = clim_df.dropna(axis=1, how='all')
    test_ann_means = clim_df.groupby(clim_df.index.year).mean().mean()
    test_ann_stds = clim_df.groupby(clim_df.index.year).std().mean()
    summary_arr[index[0]*2,:] = test_ann_means[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index]
    summary_arr[index[0]*2 + 1, :] = test_ann_stds[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index]

cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(summary_arr, n_clusters, 2, error=0.005, maxiter=1000)
cluster_membership = np.argmax(u, axis=0)
test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')

    test = clim_df[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index].transpose()




max_n_clusters = 10  # square root of num_hexes
fpc_arr = np.zeros([9]); fpc_arr[:] = np.nan
fpc_arr_weighted = np.zeros([9]); fpc_arr_weighted[:] = np.nan
num_hexes = np.shape(wnv_ann)[1]
# ok, let's do this with some climate data to create regions for aggregation
temp_data = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/PRISM_tmean_198101-202211_gridded-2.csv')
temp_data.index = temp_data['Unnamed: 0'].values
temp_data = temp_data.drop(columns=['Unnamed: 0'])
temp_data.index = pd.to_datetime(temp_data.index)
temp_data_sub = temp_data[temp_data.index.year <= 2021]
temp_mon_anom = temp_data_sub.groupby(temp_data_sub.index.month).transform(lambda x: x-x.mean())
# temp_data_sub = temp_data_sub.groupby(temp_data_sub.index.month).values() - temp_data_sub.groupby(temp_data_sub.index.month).mean()
test = temp_mon_anom[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index].transpose()
# test_ann_means = temp_data_sub.groupby(temp_data_sub.index.year).mean()
# test_ann_stds = temp_data_sub.groupby(temp_data_sub.index.year).std()
for n_clusters in np.arange(2, max_n_clusters+1):
    cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(test), n_clusters, 2, error=0.005, maxiter=1000)
    fpc_arr[n_clusters - 2] = fpc
    fpc_arr_weighted[n_clusters - 2] = fpc * n_clusters
fig, ax1 = plt.subplots()
ax1.plot(np.arange(2, max_n_clusters+1), fpc_arr)
ax2 = ax1.twinx()
ax2.plot(np.arange(2, max_n_clusters+1), fpc_arr_weighted, color='r')
plt.show()
cluster_membership = np.argmax(u, axis=0)
test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')



fpath = '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/hex_grid/'
clim_vars = ['temp_mean', 'temp_max', 'temp_min', 'temp_dew', 'drought', 'precipitation', 'soil_moisture_50', 'snow_depth', 'VPD max']  # add in something with moisture?
summary_arr = np.zeros([2*len(clim_vars), num_hexes])
matplotlib.use('AGG')
for index, clim_var in np.ndenumerate(clim_vars):
    if clim_var == 'temp_mean':
        clim_data = load_clim_data(fpath + 'PRISM_tmean_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_max':
        clim_data = load_clim_data(fpath + 'PRISM_tmax_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_min':
        clim_data = load_clim_data(fpath + 'PRISM_tmin_198101-202211_gridded-2.csv')
    elif clim_var == 'temp_dew':
        clim_data = load_clim_data(fpath + 'PRISM_td_198101-202211_gridded-2.csv')
    elif clim_var == 'drought':
        clim_data = load_clim_data(fpath + 'gridMET_PDSI_199901-202212_gridded-2.csv')
    elif clim_var == 'precipitation':
        clim_data = load_clim_data(fpath + 'PRISM_precip_198101-202211_gridded-2.csv')
    elif clim_var == 'snow_depth':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Snow Melt':
        clim_data = load_clim_data(fpath + 'NLDAS_snow_melt_199901-202212_gridded-2.csv')
    elif clim_var == 'SMWE':
        clim_data = load_clim_data(fpath + 'NLDAS_water_equivalent_snow_depth_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_5_199901-202212_gridded-2.csv')
    elif clim_var == 'soil_moisture_50':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_50_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Moisture 100':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_moisture_100_199901-202212_gridded-2.csv')
    elif clim_var == 'Soil Temp 5':
        clim_data = load_clim_data(fpath + 'NLDAS_soil_temp_5_199901-202212_gridded-2.csv')
    elif clim_var == 'VPD_max':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmax_198101-202303_gridded-2.csv')
    elif clim_var == 'VPD_min':
        clim_data = load_clim_data(fpath + 'PRISM_vpdmin_198101-202303_gridded-2.csv')

    temp_mon_anom = clim_data.groupby(clim_data.index.month).transform(lambda x: x-x.mean())
    test_mon = temp_mon_anom[wnv_ann.sum().sort_values(ascending=False).index.astype(int)][:num_hexes].transpose()
    test_raw = clim_data[wnv_ann.sum().sort_values(ascending=False)[:num_hexes].index].transpose()

    for n_clusters in np.arange(2, 11):
        cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(test_mon), n_clusters, 2, error=0.005, maxiter=1000)
        cluster_membership = np.argmax(u, axis=0)
        test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
        test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        ax.set_aspect('equal')
        test_grid.plot(column='cluster', ax=ax, cmap='tab20c')
        plt.xlim([-129, -62])
        plt.ylim([23, 52])
        plt.xticks(np.arange(-125, -55, 10), fontsize=14)
        plt.yticks(np.arange(30, 60, 10), fontsize=14)
        plt.title(clim_var + ' monthly anom: ' + str(n_clusters) + ' clusters')
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/cluster_analysis_' + clim_var + '_mon_anoms_' + str(n_clusters) + '_clusters')
        plt.close()

        cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(np.transpose(test_raw), n_clusters, 2, error=0.005, maxiter=1000)
        cluster_membership = np.argmax(u, axis=0)
        test_df = pd.DataFrame(cluster_membership, columns=['cluster'], index=test.index)
        test_grid = hex_grid.merge(test_df, left_on='grid_id', right_on=test_df.index.astype(int), how='right')

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, layout='constrained', figsize=(10, 5))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        ax.set_aspect('equal')
        test_grid.plot(column='cluster', ax=ax, cmap='tab20c')
        plt.xlim([-129, -62])
        plt.ylim([23, 52])
        plt.xticks(np.arange(-125, -55, 10), fontsize=14)
        plt.yticks(np.arange(30, 60, 10), fontsize=14)
        plt.title(clim_var + ' raw: ' + str(n_clusters) + ' clusters')
        plt.savefig('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/figures/cluster_analysis_' + clim_var + '_raw_' + str(n_clusters) + '_clusters')
        plt.close()

