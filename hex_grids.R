
#making hexagonal grid across US
# library(sp)
# counties_df <- usmap::us_map(regions = "counties", exclude = c("HI","AK")) #all counties in CONUS
# spsample
library(tidyr)
library(dplyr)
library(sf)
library(ggplot2)
library(cowplot)

options(tigris_use_cache = TRUE)

#counties shapefile
library(tigris)
cts <- counties(state = unique(fips_codes$state[!fips_codes$state %in% c("AK","HI","PR","AS","GU","MP","UM","VI")]),
                        resolution = "20m", class="sf") %>% st_transform(., "EPSG:5070")

#states shapefile
sts <- states() %>% filter(., !(STUSPS %in% c("AK","HI","PR","AS","GU","MP","UM","VI"))) %>% st_transform(., "EPSG:5070")

#label states as climate region #https://www.ncei.noaa.gov/access/monitoring/reference-maps/us-climate-regions
sts$clim.regn <- ifelse(sts$STUSPS %in% c("IA","MI","MN","WI"), "UpperMidwest", 
                            ifelse(sts$STUSPS %in% c("IL","IN","KY","MO","OH","TN","WV"), "OhioValley",
                                   ifelse(sts$STUSPS %in% c("AL","FL","GA","NC","SC","VA"), "Southeast",
                                          ifelse(sts$STUSPS %in% c("MT","NE","ND","SD","WY"), "NRockiesandPlains",
                                                 ifelse(sts$STUSPS %in% c("AR","KS","LA","MS","OK","TX"),"South",
                                                        ifelse(sts$STUSPS %in% c("AZ","CO","NM","UT"), "Southwest",
                                                               ifelse(sts$STUSPS %in% c("ID","OR","WA"), "Northwest",
                                                                      ifelse(sts$STUSPS %in% c("CA","NV"), "West", 
                                                                             "Northeast"))))))))

#plot of climate regions by state
ggsave2("clim-region_states.tiff",
        ggplot(sts) + 
          geom_polygon(aes(x=x, y=y, group=group, fill=clim.regn), color="darkgrey", show.legend = F) + 
          theme_map(), width=6, height=4, units = "in", dpi = 300, compression="lzw")

#2000-2021 WNV counts
df_all <- read.csv("../NeuroWNV_by_county_2000-2021FULL.csv", header=T) %>%
  mutate(fips = sprintf("%05d", fips))

totpop <- read.csv("../data/census_data.csv", header=T) %>% #census data for 2010
  select(geoid, tot_pop, pop65) %>% mutate(geoid = sprintf("%05d", geoid)) %>% 
  left_join(., cts %>% as.data.frame() %>% select(GEOID,ALAND), by = c("geoid" = "GEOID")) #still has AL, HI


# cts <- merge(cts, df_all, by.x = "GEOID", by.y = "fips")
cts_cents <- st_centroid(cts) #centroids for counties

#make hex grida cross US
h <- st_make_grid(cts, cellsize = 2e5, square=FALSE) #cell size is diameter of hex (m)
h_cts <- h[cts] #crop to US outline
# h_cts <- st_transform(h_cts)

#which county centroids in which hex
b <- st_intersects(h_cts, cts_cents) #list of which county centroid(s) in which hexagon

#which climate region each hex associated with
library(parallel)
h2 <- mclapply(1:length(h_cts), function(x) { #for each hex, for computation
  print(x) #monitor progress
  w_sts <- st_intersection(sts, h_cts[x]) %>%
    mutate(
      a_ovr = st_area(.), #area of overalp with state
      hnum = x) #which hex this intersection for

  w_sts %>% group_by(clim.regn) %>%
    summarise(tot = sum(a_ovr)) %>% slice_max(tot) %>%
    as.data.frame() %>% select(clim.regn) %>% #which clim region with largest area overlap (across all overlapping states in that region)
    st_set_geometry(., h_cts[x]) #replace geometry with hex geometry

}, mc.cores= 30)

h_clim <- do.call(rbind, h2)

h_clim <- readRDS("../data/h_clim_proj.rds") #R object of hexagonal grid


# calculate total counts, pop per hex
sums_df <-lapply(seq_along(b), function(x) {
  ovr_locs <- cts_cents$GEOID[b[[x]]] #GEOID/fips for any overlapping county centroids for that hexagon
  if(length(ovr_locs) > 0) {
    wnv = df_all[df_all$fips %in% ovr_locs,] #WNV data for overlap
    pop = totpop[totpop$geoid %in% ovr_locs,] #demographics data for overlap
    wnv %>% group_by(year) %>% 
      summarise(tot_count = sum(count)) %>% #sum of cases from overlapping counties per year
      mutate(totpop = sum(pop$tot_pop), #sum of total population from overlapping counties
             pop65 = sum(pop$pop65), #sum of pop > 65 from overlapping counties
             pdens = sum(pop$tot_pop)/(sum(pop$ALAND)/1e6), #total pop density (/km^2) (tot pop/area of overlapping counties)
             geometry = list(h_cts[[x]])) %>%
      st_as_sf(sf_column_name = "geometry", crs = st_crs(h_cts)) #add sfc polygon, sf object
    } else {
    data.frame(year = 2000:2021, tot_count = NA, totpop = NA, pop65 = NA, pdens = NA) %>% #NA counts if no centroid 
      mutate(geometry = list(h_cts[[x]])) %>%
      st_as_sf(sf_column_name = "geometry", crs = st_crs(h_cts)) #add sfc polygon, sf object
  }
})


sd_all <- do.call(rbind, sums_df) #all hexagons with data together

saveRDS(sd_all, "../data/WNND_hex_2000-2021.rds") #r object with cases, pop per hex

# ggplot(sd_all[sd_all$year == 2010, ]) +
#   geom_sf(aes(fill = tot_count)) + #facet_wrap(~year) +
#   theme_map()
# 
# #z-score of cases (case - mean)/SD; using mean and SD calc from previous years (not including current year for calc)
# zdf <- lapply(sums_df, function(x) {
#   # x$zscore = sapply(1:nrow(x), function(i) {
#   #   if(i == 1) { #first year so no mean or SD of previous year's cases
#   #     NA
#   #   } else {
#   #     rmean = mean(x$tot_count[1:(i-1)]) #running mean # cases of all previous years
#   #     rSD = sd(x$tot_count[1:(i-1)]) #SD of running # cases of all previous years
#   #     (x$tot_count[i] - rmean)/rSD
#   #   }
#   # })
#   # x
#   x$zscore = (x$tot_count - mean(x$tot_count)) / sd(x$tot_count) #zscore based on mean/SD of whole range (2000-2021)
#   if(sum(x$tot_count) == 0 & !is.na(sum(x$tot_count))) {
#     x$zscore = replace(x$zscore, is.na(x$zscore), 0)} #if never cases -> zscore is 0 vs NA (NAs stay for no overlap hexes)
#   x
# })
# 
# sd_all2 <- do.call(rbind, zdf)
# saveRDS(sd_all2, "../data/WNND_hex_2000-2021_zscores.rds")
