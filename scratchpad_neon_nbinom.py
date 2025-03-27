# This script creates a negative binomial distribution baseline from historical WNV data for NEON region regression.
#
# input: annual-NEON region level neuroinvasive WNV case counts (csv)
# output: annual-NEON region level negative binomial forecasts of 2023 neuroinvasive WNV reported cases
#
# note: use binom_param_exact
#
# Initialized by RDH on 5/11/2023

## Importing libraries
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import nbinom
from itertools import product

## Opening WNV data
wnv_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv', index_col=0)
wnv_df.index = pd.to_datetime(wnv_df.index)
wnv_df = wnv_df[wnv_df.index.year >= 2005]
wnv_df_ann = wnv_df.groupby(wnv_df.index.year).sum()
neon_regions = wnv_df_ann.columns.values

def counts_to_neg_binom_quantiles(count_array):
    X = np.ones_like(count_array)
    res = sm.NegativeBinomial(count_array, X).fit(start_params=[1,1], disp=False)
    mu = np.exp(res.params[0])
    p = 1/(1+np.exp(res.params[0])*res.params[1])
    n = np.exp(res.params[0])*p/(1-p)
    count_dist = np.random.negative_binomial(n, p, size=np.size(count_array))
    quantiles = [1, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 97.5, 99]
    quantile_values = np.percentile(count_dist, quantiles)
    return quantile_values

def neg_binom_param_calculator(count_array):
    X = np.ones_like(count_array)
    res = sm.NegativeBinomial(count_array, X).fit(start_params=[1,1], disp=False)
    mu = np.exp(res.params[0])
    p = 1/(1+np.exp(res.params[0])*res.params[1])
    n = np.exp(res.params[0])*p/(1-p)
    return n, p

def binom_param_exact(count_array):
    n, p = neg_binom_param_calculator(count_array)
    quantiles = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99]
    quantile_1 = nbinom.ppf(quantiles[0], n, p)
    quantile_2p5 = nbinom.ppf(quantiles[1], n, p)
    quantile_5 = nbinom.ppf(quantiles[2], n, p)
    quantile_10 = nbinom.ppf(quantiles[3], n, p)
    quantile_15 = nbinom.ppf(quantiles[4], n, p)
    quantile_20 = nbinom.ppf(quantiles[5], n, p)
    quantile_25 = nbinom.ppf(quantiles[6], n, p)
    quantile_30 = nbinom.ppf(quantiles[7], n, p)
    quantile_35 = nbinom.ppf(quantiles[8], n, p)
    quantile_40 = nbinom.ppf(quantiles[9], n, p)
    quantile_45 = nbinom.ppf(quantiles[10], n, p)
    quantile_50 = nbinom.ppf(quantiles[11], n, p)
    quantile_55 = nbinom.ppf(quantiles[12], n, p)
    quantile_60 = nbinom.ppf(quantiles[13], n, p)
    quantile_65 = nbinom.ppf(quantiles[14], n, p)
    quantile_70 = nbinom.ppf(quantiles[15], n, p)
    quantile_75 = nbinom.ppf(quantiles[16], n, p)
    quantile_80 = nbinom.ppf(quantiles[17], n, p)
    quantile_85 = nbinom.ppf(quantiles[18], n, p)
    quantile_90 = nbinom.ppf(quantiles[19], n, p)
    quantile_95 = nbinom.ppf(quantiles[20], n, p)
    quantile_97p5 = nbinom.ppf(quantiles[21], n, p)
    quantile_99 = nbinom.ppf(quantiles[22], n, p)
    median_quantiles = np.array([quantile_1, quantile_2p5, quantile_5, quantile_10, quantile_15, quantile_20, quantile_25, quantile_30, quantile_35, quantile_40, quantile_45, quantile_50, quantile_55, quantile_60, quantile_65, quantile_70, quantile_75, quantile_80, quantile_85, quantile_90, quantile_95, quantile_97p5, quantile_99])
    return median_quantiles

## Running state-month combination loop
# multiindex = [states, months]
# ind = pd.MultiIndex.from_product(multiindex, names=['state', 'month'])
df = pd.DataFrame(index = neon_regions, columns= ['quantiles'])
# df_exact = pd.DataFrame(index = neon_regions, columns = ['quantiles'])

# creating dataframe with all available
years = wnv_df_ann.index
new_columns = ['q001', 'q025', 'q050', 'q100', 'q150', 'q200', 'q250', 'q300', 'q350', 'q400', 'q450',
           'q500', 'q550', 'q600', 'q650', 'q700', 'q750', 'q800', 'q850', 'q900', 'q950', 'q975', 'q990']
df_exact = pd.DataFrame(list(product(neon_regions, years)), columns=['neon_region', 'year'])
df_exact[new_columns] = np.nan

for index, row in df_exact.iterrows():
    count_array = np.take(wnv_df_ann[row['neon_region']].values, array_index[array_index != np.where(years == row['year'])[0]])
    if np.sum(count_array) == 0:
       df_exact.loc[(df_exact['neon_region'] == row['neon_region']) & (df_exact['year'] == row['year']), new_columns] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1])  # using 0 count Bayesian model from Mike
        continue
    df_exact.loc[(df_exact['neon_region'] == row['neon_region']) & (df_exact['year'] == row['year']), new_columns] = binom_param_exact(count_array)

df_exact.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/neon_domains_nb_loo.csv')

# exploding/reformatting to something I can input into the wis code
domain_names = {'1': 'Northeast', '2': 'Mid Atlantic', '3': 'Southeast', '4': 'Atlantic Neotropical', '5': 'Great Lakes',
                '6': 'Prairie Peninsula', '7': 'Appalachians & Cumberland Plateau', '8': 'Ozarks Complex',
                '9': 'Northern Plains', '10': 'Central Plains', '11': 'Southern Plains', '12': 'Northern Rockies',
                '13': 'Southern Rockies and Colorado Plateau', '14': 'Desert Southwest', '15': 'Great Basin',
                '16': 'Pacific Northwest', '17': 'Pacific Southwest'}
# quantile = {0: '0.01', 1: '0.025', 2: '0.05', 3: '0.10', 4: '0.15', 5: '0.20', 6: '0.25', 7: '0.30', 8: '0.35',
#             9: '0.40', 10: '0.45', 11: '0.5', 12: '0.55', 13: '0.6', 14: '0.65', 15: '0.70', 16: '0.75',
#             17: '0.80', 18: '0.85', 19: '0.90', 20: '0.95', 21: '0.975', 22: '0.99'}
quantile = {'q001': '0.01', 'q025': '0.025', 'q050': '0.05', 'q100': '0.10', 'q150': '0.15', 'q200': '0.20',
            'q250': '0.25', 'q300': '0.30', 'q350': '0.35', 'q400': '0.40', 'q450': '0.45', 'q500': '0.5',
            'q550': '0.55', 'q600': '0.6', 'q650': '0.65', 'q700': '0.70', 'q750': '0.75', 'q800': '0.80',
            'q850': '0.85', 'q900': '0.90', 'q950': '0.95', 'q975': '0.975', 'q990': '0.99'}

nb_df = pd.DataFrame(columns=['domain_name', 'domain_num', 'year', 'quantile', 'value'])
for i, row in df_exact.iterrows():
    for j in new_columns:
        row_num = len(nb_df)
        nb_df.loc[row_num, 'domain_name'] = domain_names[row['neon_region']]
        nb_df.loc[row_num, 'domain_num'] = row['neon_region']
        nb_df.loc[row_num, 'year'] = row['year']
        nb_df.loc[row_num, 'quantile'] = quantile[j]
        nb_df.loc[row_num, 'value'] = np.round(row[j], 2)

nb_df.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/neon_domains_nb_exploded_loo.csv')


## Converting R output to same wis format as above
clim_forecast_df = pd.read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/test_9_forecast2.csv', index_col=0)

quantile_per = {'1%': '0.01', '2.5%': '0.025', '5%': '0.05', '10%': '0.10', '15%': '0.15', '20%': '0.20', '25%': '0.25',
                '30%': '0.30', '35%': '0.35', '40%': '0.40', '45%': '0.45', '50%': '0.5', '55%': '0.55', '60%': '0.6',
                '65%': '0.65', '70%': '0.70', '75%': '0.75', '80%': '0.80', '85%': '0.85', '90%': '0.90',
                '95%': '0.95', '97.5%': '0.975', '99%': '0.99'}

clim_forecast_df_exploded = pd.DataFrame(columns=['domain_name', 'domain_num', 'year', 'quantile', 'value'])
domain_num = 9
for i in clim_forecast_df.columns:
    for j, row in clim_forecast_df[i].items():
        row_num = len(clim_forecast_df_exploded)
        clim_forecast_df_exploded.loc[row_num, 'domain_name'] = domain_names[str(domain_num)]
        clim_forecast_df_exploded.loc[row_num, 'domain_num'] = domain_num
        clim_forecast_df_exploded.loc[row_num, 'year'] = int(i)
        clim_forecast_df_exploded.loc[row_num, 'quantile'] = quantile_per[j]
        clim_forecast_df_exploded.loc[row_num, 'value'] = np.round(row, 2)

clim_forecast_df_exploded.to_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/test_9_forecast_exploded.csv')

    for j in row['quantiles']:
        row_num = len(nb_df)
        clim_forecast_df_exploded.loc[row_num, 'domain_name'] = domain_names[i]
        clim_forecast_df_exploded.loc[row_num, 'domain_num'] = i
        clim_forecast_df_exploded.loc[row_num, 'year'] = i
        clim_forecast_df_exploded.loc[row_num, 'quantile'] = quantile[np.where(row['quantiles'] == j)[0][0]]
        clim_forecast_df_exploded.loc[row_num, 'value'] = np.round(j, 2)



## Setting up csv file (filling in the template creator code with calculated values)
from itertools import product

def expand_grid(dictionary):
   return pd.DataFrame([row for row in product(*dictionary.values())],
                       columns=dictionary.keys())

dictionary = {'location': ['Alabama', 'Arizona', 'Arkansas', 'California',
                           'Colorado', 'Connecticut', 'Delaware', 'District of Columbia',
                           'Florida', 'Georgia', 'Idaho', 'Illinois','Indiana',
                           'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine',
                           'Maryland', 'Massachusetts', 'Michigan', 'Minnesota',
                           'Mississippi', 'Missouri', 'Montana', 'Nebraska',
                           'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico',
                           'New York', 'North Carolina', 'North Dakota', 'Ohio',
                           'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island',
                           'South Carolina', 'South Dakota', 'Tennessee', 'Texas',
                           'Utah', 'Vermont', 'Virginia', 'Washington',
                           'West Virginia', 'Wisconsin', 'Wyoming'],
              'forecast_date': ['2023-04-30'],
              'target_end_date': ['2023-05-31', '2023-06-30', '2023-07-31', '2023-08-31',
                                  '2023-09-30', '2023-10-31', '2023-11-30', '2023-12-31'],
              'target': ['Monthly WNV neuroinvasive disease cases'],
              'type': ['quantile'],
              'quantile': ['0.01', '0.025', '0.05', '0.1', '0.15', '0.20', '0.25',
                           '0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6',
                           '0.65', '0.7', '0.75', '0.8', '0.85', '0.9', '0.95',
                           '0.975', '0.99']}

template_df = expand_grid(dictionary)
template_df['value'] = np.nan
template_df[['target_end_year', 'target_end_month', 'target_end_day']] = template_df['target_end_date'].str.split('-', expand=True).astype(int)

end_date_str = ['May WNV neuroinvasive disease cases', 'June WNV neuroinvasive disease cases', 'July WNV neuroinvasive disease cases',
                'August WNV neuroinvasive disease cases', 'September WNV neuroinvasive disease cases', 'October WNV neuroinvasive disease cases',
                'November WNV neuroinvasive disease cases', 'December WNV neuroinvasive disease cases']
for i in np.arange(5, 13):
    template_df['target'][template_df['target_end_month'] == i] = end_date_str[i-5]

df_fips = pd.read_csv('../supplemental_files/state_fips_codes.txt', delimiter='|')
df_fips = df_fips.drop(['STATE', 'STATENS'], axis=1)
df_fips = df_fips.rename(columns={'STATE_NAME': 'location', 'STUSAB': 'state_abbrev'})

template_df = template_df.merge(df_fips)

df_exact = df_exact[np.in1d(df.index.get_level_values(1), [np.arange(5, 13)])]
for index, row in df_exact.iterrows():
    template_df['value'][((template_df['state_abbrev'] == index[0]) & (template_df['target_end_month'] == index[1]))] = row['quantiles']
    # template_df['value'][((template_df['state_abbrev'] == index[0]) & (template_df['target_end_month'] == index[1]))] = np.round(row['quantiles'],0).astype(int)
template_df = template_df.drop(['target_end_year', 'target_end_month', 'target_end_day', 'state_abbrev'], axis=1)
template_df['value'] = template_df['value'].round(0).astype(int)
template_df.to_csv('/Users/ryan.harp/Documents/2023_WNV_Forecasting_Challenge/supplemental_files/WNV_neg_binom_forecast_template_test_exact_no_bootstrap_2010_with_bayesian.csv', index=False)

