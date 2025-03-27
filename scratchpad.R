library(tidyverse)
library(dplyr)
library(rstanarm)
library(lubridate)
library(MASS)

options(mc.cores = parallel::detectCores())  # to parallelize
neon_domain_month <- c('17'=3, '14'=7, '5'=7, '6'=7, '3'=1, '11'=8, '10'=10,
                       '8'=11, '9'=3, '2'=5, '1'=6)
neon_domain_var <- c('17'='temp', '14'='temp', '5'='temp', '6'='PDSI',
                     '3'='soil_moisture', '11'='temp', '10'='PDSI', 
                     '8'='PDSI', '9'='PDSI', '2'='temp', '1'='soil_moisture')
var_file <- c('temp'='prism_tmean_199901-202211_loo.csv', 
              'PDSI'='gridMET_PDSI_199901-202211_loo.csv',
              'soil_moisture'='NLDAS_soil_moisture_50_199901-202211_loo.csv')

# prepping WNV data for analysis
wnv_df <- read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv')
colnames(wnv_df)[1] <- "date"
wnv_df$date <-  as.Date(wnv_df$date,'%m/%d/%Y')
wnv_df$year <- as.numeric(format(wnv_df$date,'%Y'))

# getting annual WNV counts
wnv_ann <- wnv_df %>% 
  group_by(wnv_df$year) %>%
  summarise(across(where(is.numeric), sum))
colnames(wnv_ann)[1] <- "year"

# getting climate data
neon_domains <- c('17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1')
for (neon_domain in neon_domains) {
  print(neon_domain)

  fpath <- '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
  clim_df <- read_csv(paste(fpath, var_file[neon_domain_var[neon_domain]], sep=''))
  colnames(clim_df)[1] <- "date"
  clim_df$date <-  as.Date(clim_df$date,'%m/%d/%Y')
  clim_df$year <- as.numeric(format(clim_df$date,'%Y'))
  clim_df$year <- ifelse(clim_df$year>90, clim_df$year+1900, clim_df$year+2000)
  clim_df$month <- as.numeric(format(clim_df$date,'%m'))
  
  # subsetting climate data
  if (neon_domain_month[neon_domain] < 10) {
    clim_sub <- clim_df[(clim_df$year >= 2005) & (clim_df$month == neon_domain_month[neon_domain]),][neon_domain]
  }
  if (neon_domain_month[neon_domain] >= 10) {
    clim_sub <- clim_df[(clim_df$year >=2004) & (clim_df$year < 2022) & (clim_df$month == neon_domain_month[neon_domain]),][neon_domain]
  }
  wnv_sub <- wnv_ann[(wnv_ann$year >= 2005),][neon_domain]
  sub_df <- do.call(rbind, Map(data.frame, 'clim'=clim_sub, 'wnv'=wnv_sub))
  
  # mod <- glm.nb(wnv ~ clim, data=sub_df)
  # wnv_forecast <- predict(mod, sub_df, type='response', se.fit = T)
  # cor(log(wnv_forecast), log(wnv_sub), method='pearson')
  
  quantile_func <- function(x) {
    value <- quantile(x, c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                           0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                           0.95, 0.975, 0.99))
    return(value)
  }
  
  mod2_quantiles_all <- data.frame(matrix(NA, nrow = 23, ncol = nrow(sub_df)))
  colnames(mod2_quantiles_all) <- c(2005:2022)
  # rownames(mod2_quantiles_all)
  for (i in 1:nrow(sub_df)) {
    mod2 <- stan_glm.nb(wnv ~ clim, data=sub_df[-i,], chains=4, iter=1000)
    clim_test_input <- data.frame(clim = sub_df[i,]$clim)
    mod2_pdf <- posterior_predict(mod2, clim_test_input)
    mod2_quantiles <- apply(X = mod2_pdf, MARGIN = 2, FUN = quantile_func)
    mod2_quantiles_all[,i] <- mod2_quantiles
  }
  rownames(mod2_quantiles_all) <- rownames(mod2_quantiles)
  
  # mod2_pdf <- posterior_predict(mod2)
  # mod2_quantiles <- apply(X = mod2_pdf, MARGIN = 2, FUN = quantile_func)
  # colnames(mod2_quantiles) <- c(2005:2022)
  write.csv(mod2_quantiles_all, file=paste('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_loo_clim_', neon_domain, '.csv', sep=''))
}



stan_summary <- round(mod$stan_summary, 3) #param estimates
stan_summary
rmse <- round(sqrt(mean(mod$residuals^2)),3) #rmse
mae <- round(mean(abs(mod$residuals)),3) #MAE
error_df <- data.frame(rmse, mae)
error_df
loo_summary <- loo(mod) #ELPD
loo_summary

# saveRDS(mod, file='/Users/ryan.harp/Documents/2022_WNV_Forecasting_Challenge/results/all_county_regression.rds')








## Trying some multivariate regressions

options(mc.cores = parallel::detectCores())  # to parallelize
neon_domain_month1 <- c('17'=3, '14'=7, '5'=7, '6'=3, '3'=1, '11'=10, '10'=10,
                       '8'=11, '9'=3, '2'=5, '1'=7)
neon_domain_var1 <- c('17'='temp', '14'='temp', '5'='temp', '6'='PDSI',
                     '3'='soil_moisture', '11'='soil_moisture', '10'='PDSI',
                    '8'='PDSI', '9'='PDSI', '2'='temp', '1'='soil_moisture')
neon_domain_month2 <- c('17'=6, '14'=7, '5'=3, '6'=5, '3'=5, '11'=3, '10'=6,
                        '8'=6, '9'=5, '2'=7, '1'=5)
neon_domain_var2 <- c('17'='soil_moisture', '14'='precip', '5'='soil_moisture', '6'='temp',
                      '3'='vpd_max', '11'='precip', '10'='dew_point',
                      '8'='dew_point', '9'='temp', '2'='precip', '1'='temp')
var_file <- c('temp'='prism_tmean_199901-202211_loo.csv',
             'PDSI'='gridMET_PDSI_199901-202211_loo.csv',
             'soil_moisture'='NLDAS_soil_moisture_50_199901-202211_loo.csv',
             'precip' = 'prism_precip_199901-202211_loo.csv',
             'vpd_max' = 'prism_vpdmax_199901-202211_loo.csv',
             'dew_point' = 'prism_td_199901-202211_loo.csv')

# prepping WNV data for analysis
wnv_df <- read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv')
colnames(wnv_df)[1] <- "date"
wnv_df$date <-  as.Date(wnv_df$date,'%m/%d/%Y')
wnv_df$year <- as.numeric(format(wnv_df$date,'%Y'))

# getting annual WNV counts
wnv_ann <- wnv_df %>% 
  group_by(wnv_df$year) %>%
  summarise(across(where(is.numeric), sum))
colnames(wnv_ann)[1] <- "year"

# getting climate data
neon_domains <- c('17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1')
for (neon_domain in neon_domains) {
  print(neon_domain)
  
  fpath <- '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
  clim_df1 <- read_csv(paste(fpath, var_file[neon_domain_var1[neon_domain]], sep=''))
  colnames(clim_df1)[1] <- "date"
  clim_df1$date <-  as.Date(clim_df1$date,'%m/%d/%Y')
  clim_df1$year <- as.numeric(format(clim_df1$date,'%Y'))
  clim_df1$year <- ifelse(clim_df1$year>90, clim_df1$year+1900, clim_df1$year+2000)
  clim_df1$month <- as.numeric(format(clim_df1$date,'%m'))
  
  
  # subsetting climate data
  if (neon_domain_month1[neon_domain] < 10) {
    clim_sub1 <- clim_df1[(clim_df1$year >= 2005) & (clim_df1$month == neon_domain_month1[neon_domain]),][neon_domain]
  }
  if (neon_domain_month1[neon_domain] >= 10) {
    clim_sub1 <- clim_df1[(clim_df1$year >=2004) & (clim_df1$year < 2022) & (clim_df1$month == neon_domain_month1[neon_domain]),][neon_domain]
  }
  
  clim_df2 <- read_csv(paste(fpath, var_file[neon_domain_var2[neon_domain]], sep=''))
  colnames(clim_df2)[1] <- "date"
  clim_df2$date <-  as.Date(clim_df2$date,'%m/%d/%Y')
  clim_df2$year <- as.numeric(format(clim_df2$date,'%Y'))
  clim_df2$year <- ifelse(clim_df2$year>90, clim_df2$year+1900, clim_df2$year+2000)
  clim_df2$month <- as.numeric(format(clim_df2$date,'%m'))
  
  # subsetting climate data
  if (neon_domain_month2[neon_domain] < 10) {
    clim_sub2 <- clim_df2[(clim_df2$year >= 2005) & (clim_df2$month == neon_domain_month2[neon_domain]),][neon_domain]
  }
  if (neon_domain_month2[neon_domain] >= 10) {
    clim_sub2 <- clim_df2[(clim_df2$year >=2004) & (clim_df2$year < 2022) & (clim_df2$month == neon_domain_month2[neon_domain]),][neon_domain]
  }
  
  wnv_sub <- wnv_ann[(wnv_ann$year >= 2005),][neon_domain]
  
  sub_df <- do.call(rbind, Map(data.frame, 'clim1'=clim_sub1, 'clim2'=clim_sub2, 'wnv'=wnv_sub))
  
  quantile_func <- function(x) {
    value <- quantile(x, c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                           0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                           0.95, 0.975, 0.99))
    return(value)
  }
  
  mod2_quantiles_all <- data.frame(matrix(NA, nrow = 23, ncol = nrow(sub_df)))
  colnames(mod2_quantiles_all) <- c(2005:2022)
  # rownames(mod2_quantiles_all)
  for (i in 1:nrow(sub_df)) {
    mod2 <- stan_glm.nb(wnv ~ clim1 + clim2, data=sub_df[-i,], chains=4, iter=1000)
    clim_test_input <- data.frame(clim = sub_df[i,c('clim1', 'clim2')])
    mod2_pdf <- posterior_predict(mod2, sub_df[i,c('clim1', 'clim2')])
    # mod2_pdf <- posterior_predict(mod2, )
    mod2_quantiles <- apply(X = mod2_pdf, MARGIN = 2, FUN = quantile_func)
    mod2_quantiles_all[,i] <- mod2_quantiles
  }
  rownames(mod2_quantiles_all) <- rownames(mod2_quantiles)
  
  write.csv(mod2_pdf, file=paste('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_loo_clim_pdf_2022_', neon_domain, '.csv', sep=''))
  write.csv(mod2_quantiles_all, file=paste('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_loo_clim_', neon_domain, '.csv', sep=''))
}



region_wnv <- c(46, 25, 23, 41, 6, 101, 47, 118, 31, 29, 65, 41, 62, 196, 18, 34, 65, 93)
county_wnv <- c(11, 4, 2, 17, 0, 38, 7, 5, 5, 2, 5, 6, 7, 13, 2, 6, 4, 22)
county_wnv_per <- county_wnv/region_wnv
prop_mod <- stan_glm(county_wnv_per ~ region_wnv, family=poisson(), chains=4, iter=1000)
prop_test <- posterior_predict(prop_mod)


# getting climate data for Northern Plains

options(mc.cores = parallel::detectCores())  # to parallelize
neon_domain_month1 <- c('17'=3, '14'=7, '5'=7, '6'=3, '3'=1, '11'=10, '10'=10,
                        '8'=11, '9'=3, '2'=5, '1'=7)
neon_domain_var1 <- c('17'='temp', '14'='temp', '5'='temp', '6'='PDSI',
                      '3'='soil_moisture', '11'='soil_moisture', '10'='PDSI',
                      '8'='PDSI', '9'='PDSI', '2'='temp', '1'='soil_moisture')
neon_domain_month2 <- c('17'=6, '14'=7, '5'=3, '6'=5, '3'=5, '11'=3, '10'=6,
                        '8'=6, '9'=5, '2'=7, '1'=5)
neon_domain_var2 <- c('17'='soil_moisture', '14'='precip', '5'='soil_moisture', '6'='temp',
                      '3'='vpd_max', '11'='precip', '10'='dew_point',
                      '8'='dew_point', '9'='temp', '2'='precip', '1'='temp')
var_file <- c('temp'='prism_tmean_199901-202211_loo.csv',
              'PDSI'='gridMET_PDSI_199901-202211_loo.csv',
              'soil_moisture'='NLDAS_soil_moisture_50_199901-202211_loo.csv',
              'precip' = 'prism_precip_199901-202211_loo.csv',
              'vpd_max' = 'prism_vpdmax_199901-202211_loo.csv',
              'dew_point' = 'prism_td_199901-202211_loo.csv')

# prepping WNV data for analysis
wnv_df <- read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/wnv_neon_monthly_time_series.csv')
colnames(wnv_df)[1] <- "date"
wnv_df$date <-  as.Date(wnv_df$date,'%m/%d/%Y')
wnv_df$year <- as.numeric(format(wnv_df$date,'%Y'))

# getting annual WNV counts
wnv_ann <- wnv_df %>% 
  group_by(wnv_df$year) %>%
  summarise(across(where(is.numeric), sum))
colnames(wnv_ann)[1] <- "year"

# getting climate data
neon_domains <- c('17', '14', '5', '6', '3', '11', '10', '8', '9', '2', '1')
for (neon_domain in neon_domains) {
  print(neon_domain)
  
  fpath <- '/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/neon_domains/'
  clim_df1 <- read_csv(paste(fpath, var_file[neon_domain_var1[neon_domain]], sep=''))
  colnames(clim_df1)[1] <- "date"
  clim_df1$date <-  as.Date(clim_df1$date,'%m/%d/%Y')
  clim_df1$year <- as.numeric(format(clim_df1$date,'%Y'))
  clim_df1$year <- ifelse(clim_df1$year>90, clim_df1$year+1900, clim_df1$year+2000)
  clim_df1$month <- as.numeric(format(clim_df1$date,'%m'))
  
  # subsetting climate data
  if (neon_domain_month1[neon_domain] < 10) {
    clim_sub1 <- clim_df1[(clim_df1$year >= 2005) & (clim_df1$month == neon_domain_month1[neon_domain]),][neon_domain]
  }
  if (neon_domain_month1[neon_domain] >= 10) {
    clim_sub1 <- clim_df1[(clim_df1$year >=2004) & (clim_df1$year < 2022) & (clim_df1$month == neon_domain_month1[neon_domain]),][neon_domain]
  }
  
  clim_df2 <- read_csv(paste(fpath, var_file[neon_domain_var2[neon_domain]], sep=''))
  colnames(clim_df2)[1] <- "date"
  clim_df2$date <-  as.Date(clim_df2$date,'%m/%d/%Y')
  clim_df2$year <- as.numeric(format(clim_df2$date,'%Y'))
  clim_df2$year <- ifelse(clim_df2$year>90, clim_df2$year+1900, clim_df2$year+2000)
  clim_df2$month <- as.numeric(format(clim_df2$date,'%m'))
  
  # subsetting climate data
  if (neon_domain_month2[neon_domain] < 10) {
    clim_sub2 <- clim_df2[(clim_df2$year >= 2005) & (clim_df2$month == neon_domain_month2[neon_domain]),][neon_domain]
  }
  if (neon_domain_month2[neon_domain] >= 10) {
    clim_sub2 <- clim_df2[(clim_df2$year >=2004) & (clim_df2$year < 2022) & (clim_df2$month == neon_domain_month2[neon_domain]),][neon_domain]
  }
  
  wnv_sub <- wnv_ann[(wnv_ann$year >= 2005),][neon_domain]
  
  # pulling WNV data for all counties
  county_wnv <- read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/cleaned_data/counties/wnv_county_time_series.csv')
  colnames(county_wnv)[1] <- "date"
  county_wnv_ann <- county_wnv %>%
    mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
    group_by(year) %>%
    filter(year >= 2005) %>%
    summarize_if(is.numeric, sum)

  county_to_neon <- read_csv('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/supplemental_files/county_to_neon_domain.csv')
  county_to_neon$fips <- str_remove(county_to_neon$fips, '^0+')
  fip_list = county_to_neon$fips[county_to_neon['NEON_domain'] == as.numeric(neon_domain)]
  names.use <- names(county_wnv_ann)[(names(county_wnv_ann) %in% fip_list)]
  county_wnv_ann_sub <- county_wnv_ann[, c('year', names.use)]
  zero_counties <- setdiff(fip_list, names.use)
  for (col in zero_counties) {
    county_wnv_ann_sub[col] = 0
  }

  reg_df <- do.call(rbind, Map(data.frame, 'year' = county_wnv_ann['year'], 
                               'clim1'=clim_sub1, 'clim2'=clim_sub2, 
                               'wnv'=wnv_sub))
  
  reg_county_df <- merge(reg_df, county_wnv_ann_sub, by='year')
  reg_county_df <- pivot_longer(reg_county_df, cols=-c('year', 'clim1', 'clim2', 'wnv'), names_to='county', values_to = 'wnv_county') 

  rm(stacked_county_pdf)
  for (yr in c(2005:2022)) {
    print(yr)
    reg_county_df_sub <- reg_county_df %>% filter(year != yr)
    reg_county_df_sub_new <- reg_county_df %>% filter(year == yr)
    county_mod <- stan_glmer.nb(wnv_county ~ clim1 + clim2 + (1 | county), data=reg_county_df_sub)
    # test_pdf <- posterior_predict(test_mod, test_df_test)
    county_pdf <- predicted_draws(county_mod, reg_county_df_sub_new)
    if (yr == 2005) {
      stacked_county_pdf <- county_pdf
    }
    if (yr != 2005) {
      stacked_county_pdf <- rbind(stacked_county_pdf, county_pdf)
    }
    # test_quantiles <- apply(X = test_pdf, MARGIN = 2, FUN = quantile_func)
    # test_quantiles <- as.data.frame(test_quantiles)
    # test_pivot <- pivot_longer(test_quantiles, cols=where(is.numeric), names_to='county', values_to='quantile_values' )
  }
  stacked_county_pdf_sub <- stacked_county_pdf[c('year', 'county', '.draw', '.prediction')]
  write.csv(stacked_county_pdf_sub, file=paste('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/county_level_re_forecast_', neon_domain, '.csv', sep=''))
}

# testing some stuff
library(tidybayes)
new_test <- predicted_draws(test_mod, test_df_test)
test_wide <- pivot_wider(new_test[c('year', 'county', '.prediction', '.draw')] , names_from = county, values_from = .prediction)
test_test <- subset(test_wide, select=-c(year, .draw))
apply(X = test_test, MARGIN=2, FUN=quantile_func)



m <- test_quantiles
idk <- m %>%
  tibble::rownames_to_column(var = "Row") %>%
  pivot_longer(cols=where(is.numeric), names_to='county', values_to='quantile_values')
  # pivot_longer(-Row, names_to = "Name", values_to = "Prob") %>%
  # group_by(Name) %>%
  # top_n(n = 3) %>%
  # select(c(2, 3, 1))

# test_pdf <- posterior_predict(test_mod, test_df[i,c('clim1', 'clim2')])
mod2_quantiles <- apply(X = mod2_pdf, MARGIN = 2, FUN = quantile_func)
mod2_quantiles_all[,i] <- mod2_quantiles

mod2_quantiles_all <- data.frame(matrix(NA, nrow = 23, ncol = nrow(sub_df)))
colnames(mod2_quantiles_all) <- c(2005:2022)
# rownames(mod2_quantiles_all)
for (i in 1:nrow(sub_df)) {
  mod2 <- stan_glmer.nb(wnv ~ clim1 + clim2, data=sub_df[-i,], chains=4, iter=1000)
  # clim_test_input <- data.frame(clim = sub_df[i,c('clim1', 'clim2')])
  mod2_pdf <- posterior_predict(mod2, sub_df[i,c('clim1', 'clim2')])
  mod2_quantiles <- apply(X = mod2_pdf, MARGIN = 2, FUN = quantile_func)
  mod2_quantiles_all[,i] <- mod2_quantiles
}
rownames(mod2_quantiles_all) <- rownames(mod2_quantiles)

write.csv(mod2_quantiles_all, file=paste('/Users/ryan.harp/Documents/WNV_Spatiotemporal_Patterns/results/forecast_9_multi_', neon_domain, '.csv', sep=''))


