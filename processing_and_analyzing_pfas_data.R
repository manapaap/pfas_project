# Data analysis portion of project
# Created: 7/7/21
# Updated: 7/15/21
# Name: Aakash Manapat
# Contact: aakash.p.manapat@vanderbilt.edu


# -----------------------
# Packages (more than I will actually use)
# -----------------------

library(rgdal)
library(leaflet)
library(dplyr)
library(tidyverse)
library(tidymodels)
library(stringr)

library(sf)
library(raster)
library(RSAGA)
library(rgrass7)

library(rgeos)
library(mapview)
library(ggplot2)

library(concaveman)

library(Hmisc)
library(progress)
library(future)

library(furrr)

library(ISLR)
library(corrplot)
library(caret) # https://topepo.github.io/caret/index.html

library(imputeTS)

library(doParallel)
library(randomForest)
library(inTrees)
library(caTools)
library(gbm)
library(naivebayes)
library(rFerns)

library(RANN)
# -----------------------

setwd('C:/Users/Aakas/Desktop/Stuff for DWJ/PFAS Project/Data/')


# -----------------------
# Load information as sf-type objects
# -----------------------


# Shapefiles

state_grid <- rgdal::readOGR('Shapefiles/plant_elev_temp_precp_gwreach_STATE/hex_raster_stuff_STATE_7_15_21.shp') %>%
  st_as_sf %>%
  st_transform(3857)

state_soils <- rgdal::readOGR('Shapefiles/state_soils_JOINED/state_soils_JOINED_7_15_21.shp') %>%
  st_as_sf %>%
  st_transform(3857)

tn_esab <- rgdal::readOGR('Shapefiles/tn_esab/TN_Standardized_1_6_21.shp') %>%
  st_as_sf %>%
  st_transform(3857)

ky_lines <- rgdal::readOGR('Shapefiles/ky_lines/ky_lines_FUSED_PWSID.shp') %>%
  st_as_sf %>%
  st_transform(3857) %>%     # Using the fused file to begin with as it is too large- fusing was done in QGIS
  st_simplify        # as the file is very large and not in polygon form natively

ky_pfas <- rgdal::readOGR('Shapefiles/ky_pfas_contam_map_6_23_21/ky_pfas_contam_map_6_23_21.shp') %>%
  st_as_sf %>%
  st_transform(3857)

naics_pfas <- rgdal::readOGR('Shapefiles/ECHO_NAICS/ECHO_NAICS_data.shp') %>%
  st_as_sf %>%
  st_transform(3857)

us_military_bases <- rgdal::readOGR('Shapefiles/US_military_bases/Military_Bases.shp') %>%
  st_as_sf %>%
  st_transform(3857) %>%
  st_zm

# Rasters

na_elevation <- raster('Shapefiles/USGS_elevation/na_elevation.tif') %>%
  projectRaster(crs = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs')

na_rain <- raster('Shapefiles/weather_RAW/Ppt_annual_historical.tif') %>%
  projectRaster(crs = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs')

na_temps <- raster('Shapefiles/weather_RAW/Tavg_annual_historical.tif') %>%
  projectRaster(crs = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs')

na_ground_rech <- raster('Shapefiles/effective_groundwater_recharge/RC_eff_2013.tif') %>%
  projectRaster(crs = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs')

# CSVs

sdwis_query_ky <- read.csv('CSVs/water_system_detail_ky_sdwis_Q42020.csv')


# -----------------------
# Generic setup stuff
# -----------------------

my.cluster <- parallel::makeCluster(
  parallel::detectCores(logical = FALSE) - 1, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

prog_bar <- function(count, total) {
  # Workable basic progress bar for longer loops
  # NOT WORKING RIGHT NOW
  if ((count / total) %% 10) {
    cat('XX10XX')
  }
}

workspace.size <- function() {
  # Size o all objects in environment
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object.size(get(x))))
  class(ws) <- "object_size"
  ws
}

# -----------------------
# Raster processing
# -----------------------

# Subset the full NA elevation map to be more manageable

na_elevation <- na_elevation[raster(xmn = as.numeric(st_bbox(naics_pfas_rel)[1]),
                                    xmx = as.numeric(st_bbox(naics_pfas_rel)[3]),
                                    ymn = as.numeric(st_bbox(naics_pfas_rel)[2]),
                                    ymx = as.numeric(st_bbox(naics_pfas_rel)[4])),
                             drop = FALSE]

na_rain <- na_rain %>% crop(state_soils)

na_temps <- na_temps %>% crop(state_soils)

na_ground_rech <- na_ground_rech %>% crop(state_soils)

# Rasterize relevant parameters from shapefiles so average values can later be assigned to the 
# ESABs

template <- raster(extent(state_soils), resolution = 1000,
                   crs = st_crs(state_soils)$proj4string)

template_mode <- raster(extent(state_soils), resolution = 1000)


state_soils$H_ion_conc <- 10 ^ -(state_soils$pH)


clay_prc_raster <- rasterize(state_soils, template,
                             field = 'cly_prc', fun = 'mean')

H_ion_conc_raster <- rasterize(state_soils, template,
                               field = 'H_ion_conc', fun = 'mean')

na_soc <- rasterize(state_grid, template,
                    field = 'soc_mean', fun = 'mean')

na_bel_carb <- rasterize(state_grid, template,
                         field = 'bel_carbme', fun = 'mean')

na_abv_carb <- rasterize(state_grid, template,
                         field = 'abv_carbme', fun = 'mean')

# https://gis.stackexchange.com/questions/325586/r-rasterize-spatialpolygonsdataframe-and-keep-factor-field
# For the rasterizing of categorical values

state_soils$cations <- as.factor(state_soils$cations)

nam <- unique(state_soils$cations)
nam_df <- data.frame(ID = 1:length(nam), nam = nam)
state_soils$ID <- nam_df$ID[match(state_soils$cations, nam_df$nam)]

na_cations_raster <- rasterize(state_soils, template_mode,
                               field = 'cations')

na_cations_raster <- ratify(na_cations)

rat <- levels(na_cations)[[1]]
rat$names <- nam_df$nam
rat$IDs <- nam_df$ID
levels(na_cations_raster) <- rat

# -----------------------
# NAICS data processing
# -----------------------

# Preprocessing before distance measures can be accommodated- filter to TN, KY, states around them
# Have to separate NAICS code as some have multiple codes that that is an inconvenience 

naics_pfas_rel <- naics_pfas %>% dplyr::filter(FacState %in% c('TN', 'KY', 'AL', 'GA', 'NC', 'VA',
                                                               'WV', 'OH', 'IN', 'IL', 'MO', 'AR', 'MS')) %>%
  separate(col = FacNAICSCo, sep = ' ', into = c('NAICS1', 'NAICS2', 'NAICS3', 'NAICS4', 'NAICS5', 'NAICS6', 
                                                 'NAICS7', 'NAICS8', 'NAICS9', 'NAICS10', 'NAICS11', 'NAICS12',
                                                 'NAICS13', 'NAICS14', 'NAICS15', 'NAICS16'),
           remove = FALSE)

# Give it height information

for (n in 1: nrow(naics_pfas_rel)) {
  height <- st_coordinates(naics_pfas_rel[n,])
  height <- raster::extract(na_elevation, height) * (.3048 ^ 2)
  naics_pfas_rel$height[n] <- height
}


# This is the ugliest possible way to do this filtering properly, but I don't care anymore. It works
# Goes one-by-one and filters matching values in each NAICS parameter. Yay

dumb_filter <- function(x) {
  filtered <- naics_pfas_rel %>% dplyr::filter(NAICS1 %in% x | NAICS2 %in% x | NAICS3 %in% x | NAICS4 %in% x |
                                                 NAICS5 %in% x | NAICS6 %in% x | NAICS7 %in% x | NAICS8 %in% x |
                                                 NAICS9 %in% x | NAICS10 %in% x | NAICS11 %in% x | NAICS12 %in% x |
                                                 NAICS13 %in% x | NAICS14 %in% x | NAICS15 %in% x | NAICS16 %in% x)
  return(filtered)
}


naics_pfas_transport <- dumb_filter(c(4811, 4812, 481111, 481112, 481211, 481212, 481219))

naics_pfas_firefight <- dumb_filter(c(922160, 611519, 561990))

naics_pfas_waste <- dumb_filter(c(562111, 562112, 562119, 562211, 562212, 562213, 562219, 562991))

naics_pfas_industry <- dumb_filter(c(238320, 238330, 
                                     313110, 313210, 313220, 313320, 314910,
                                     315210, 315280, 315990, 316110, 316210, 
                                     316998, 322110, 322121, 322130, 322212, 
                                     322220, 322230, 324110, 325199, 325510, 
                                     325520, 325611, 325612, 325613, 325620, 
                                     325998, 326111, 326112, 326113, 326119,
                                     326150, 32619, 32629, 332215, 33281,
                                     332812, 332813, 333241, 333242, 333318, 
                                     33351, 333517, 334413, 334419, 334515, 
                                     335210, 335220, 335999, 336412, 339114,
                                     339920))

# Similar processing for US military bases

us_military_bases <- us_military_bases %>% filter(STATE_TERR %in% c('Tennessee', 'Kentucky', 'Illinois', 'Ohio',
                                                                    'West Virginia', 'Virginia', 'North Carolina', 
                                                                    'Georgia', 'Alabama', 'Mississippi', 'Missouri', 
                                                                    'Indiana', 'Arkansas')) %>%
  
  st_geometry(us_military_bases) <- st_point_on_surface(us_military_bases)

us_military_bases <- us_military_bases %>% st_cast('POINT')

us_military_bases <- us_military_bases[!duplicated(us_military_bases$SHAPE_Leng), ]

for (n in 1: nrow(us_military_bases)) {
  height <- st_coordinates(us_military_bases$geometry[n])
  height <- raster::extract(na_elevation, height) * (.3048 ^ 2)
  us_military_bases$height[n] <- height
}

# -----------------------
# Tennessee map processing
# -----------------------


# Fixing population served for TN which was being annoying
tn_esab$Pop_Served <- tn_esab$Pop_Served %>%
  str_replace_all(',', '') %>%
  as.numeric

# Calculate centroids of each service area boundary for sake of calculations

tn_esab$centroids <- tn_esab %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  # use st_point_on_surface to ensure centroid is actually inside service area boundary?
  st_geometry

# plots as a sanity check

mapview(list(tn_esab, tn_esab$centroids))

# time to approximate areas for the same of this

tn_esab$area <- tn_esab %>% st_area()
tn_esab$radius <- sqrt(tn_esab$area / pi)

# create circle polygons representing this simplified area now

tn_esab$circle_area <- st_buffer(tn_esab$centroids, tn_esab$radius)

mapview(list(tn_esab, tn_esab$circle_area))



# -----------------------
# Kentucky map processing
# -----------------------

# Will follow to a large extent the same steps as TN but its easier to segregate it

# TODO: spatial join of Louiville data to this data
# TODO: left join of population parameters to this data


# Change to TN style conventions for ease of use. 
# Remove irrelevant columns (can always add back in later)

ky_lines <- ky_lines %>%
  dplyr::rename(PWS_ID = PWSID) %>%
  relocate(PWS_ID, .before = OBJECTID) %>%
  dplyr::select(-c(MATERIAL, XY_SOURCE, XY_ISSUES, COMMENTS, Shape_len))


# Using concaveman package to generate concave hulls to approximate ESABs
# Need to run a loop to iterate one-by-one as concaveman will otherwise try to generate
# a hull across all of Kentucky

ky_esab <- ky_lines 

for (n in 1:length(ky_esab$geometry)) {
  
  points <- ky_esab$geometry[n] %>% 
    st_sf %>% 
    concaveman(2) # This parameter can be reduced to get more concave shapes- ~3 appears to be the 
  # best compromise between getting an accurate polygon and overfitting
  
  ky_esab$conc[n] <- points %>%
    st_as_sf %>%
    st_cast('POLYGON') %>%
    st_zm %>%
    st_simplify %>%
    st_geometry
}

ky_esab$conc <- ky_esab$conc %>% st_sfc %>% st_set_crs(3857)

st_geometry(ky_esab) <- ky_esab$conc
ky_esab <- ky_esab %>% dplyr::select(-c(conc))

mapview(ky_esab)


# Time to perform the same centroid operations as TN. Gonna do it in fewer steps as we know it works now
# Comments are preserved in the TN data

ky_esab$centroids <- ky_esab %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  st_geometry

ky_esab$area <- ky_esab %>% st_area()
ky_esab$radius <- sqrt(ky_esab$area / pi)
ky_esab$circle_area <- st_buffer(ky_esab$centroids, ky_esab$radius)

mapview(list(ky_esab, ky_esab$circle_area))

# Write to file to test loading

st_write(ky_esab, dsn = 'Shapefiles/ky_esab_test', 
         layer ='ky_esab', driver = "ESRI Shapefile")


# Time to add in PFAS data by merging with the PFAS data
# Joining by name and not spatially as many of the points overlap with multiple SAs,
# and some lie outside the SA boundary

any(ky_pfas$Location %in% ky_esab$SYS_NAME)

# Clearly, the location names do not match up perfectly
# Need to do some sort of 'fuzzy merge' then

ky_esab$matched_names <- vector(mode = 'character', length = nrow(ky_esab))

# Preparing all the names by removing redundant phrases that will complicate the algorithm
# and converting all to lowercase

prepare_names <- function(names_list) {
  simple_names <- names_list %>%
    tolower %>%
    str_replace_all('city', '') %>%
    str_replace_all('water', '') %>%
    str_replace_all('department', '') %>%
    str_replace_all('district', '') %>%
    str_replace_all('works', '') %>%
    str_replace_all('utilities', '') %>%
    str_replace_all('association', '') %>%
    str_replace_all('school', '') %>%
    str_replace_all('park', '') %>%
    str_replace_all('sewer', '') %>%
    str_replace_all(' of', '') %>%
    str_replace_all('comission', '') %>%
    str_replace_all('system', '') %>%
    str_replace_all('trailer', '') %>%
    str_replace_all('municipal', '') %>%
    str_replace_all('inc', '') %>%
    str_replace_all('wtp', '') %>%
    str_replace_all('"', '')
  return(simple_names)
}

ky_esab$prep_names <- prepare_names(ky_esab$SYS_NAME)

ky_pfas$prep_names <- prepare_names(ky_pfas$Location)

# Fuzzy match time

for(n in 1:length(ky_esab$prep_names)) {
  x <- agrep(ky_esab$prep_names[n], ky_pfas$prep_names,
             max.distance = 0.1,
             value = TRUE,
             ignore.case = TRUE,
             useBytes = TRUE)
  
  if (length(x) > 0) {
    ky_esab$matched_names[n] <- x
  } else {
    ky_esab$matched_names[n] <- NA
  }
} 

# Self-evident that fuzzy match is not functioning well. Do it manually?

for (n in 1:length(ky_esab$matched_names)) {
  if (is.na(ky_esab$matched_names[n]) == FALSE) {
    cat(ky_esab$SYS_NAME[n],'      ||      ')
    print(ky_esab$matched_names[n])
  }
}

# Manually matching names from PFAS to ESAB as the fuzzy match is terrible

write.csv(as.data.frame(ky_esab$SYS_NAME, ky_esab$matched_names), 'CSVs/name_merge.csv')
write.csv(as.data.frame(ky_pfas$Location), 'CSVs/ref_names.csv')


# Sum- ky-american (KAWC), owensboro, nkwd, Gallatin Co WD
# ofc, louiville still an issue

matched_names <- read_csv('CSVs/name_merge_MERGED.csv') %>% 
  as_tibble

ky_esab <- ky_esab %>% left_join(matched_names %>% dplyr::select(SYS_NAME, Location), by = 'SYS_NAME')

# Having to extract the relevant PFAS information as sf is being annoying

ky_pfas_info <- ky_pfas %>% 
  dplyr::select(Location, Type, 
                PFBS, HFPO..DA, PFHpA, PFHxS, ADONA, PFOA, PFOS, PFNA,
                net_PFAS, Units) %>%
  as_tibble

ky_esab <- ky_esab %>% left_join(ky_pfas_info %>%
                                   dplyr::select(Location, Type, 
                                                 PFBS, HFPO..DA, PFHpA, PFHxS, ADONA, PFOA, PFOS, PFNA,
                                                 net_PFAS, Units),
                                 by = 'Location')

# Time to fuse values with multiple entries. Only bothering for locations where PFAS was found, and only for net PFAS
# Without this transformation, set to highest value

n <- match('KY-American KY River Sta II', ky_esab$Location)
ky_esab$net_PFAS[n] <- ky_esab$net_PFAS[n] + ky_pfas$net_PFAS[23] + ky_pfas$net_PFAS[61]

n <- match('Owensboro Mun Utilities', ky_esab$Location)
ky_esab$net_PFAS[n] <- ky_esab$net_PFAS[n] + ky_pfas$net_PFAS[35] 

n <- match('NKWD FT Thomas WTP', ky_esab$Location)
ky_esab$net_PFAS[n] <- ky_esab$net_PFAS[n] + ky_pfas$net_PFAS[33]

# Fantastic, PFAS information is now coded in! Time to visualize this

mapview(ky_esab, legend = TRUE, at = c(0, 1, 10, 70), zcol = 'net_PFAS')

# Time to add elevation information to these points

# I might need to find a different elevation basemap- something is seriously wrong with this one
# in terms of absolute heights. The scale factor should help for now

# Height assignment process for ky_esab as that was not available before

for (n in 1:nrow(ky_esab)) {
  height <-  st_coordinates(ky_esab$centroids[n])
  height <- raster::extract(na_elevation, height) * (.3048 ^ 2)
  ky_esab$height[n] <- height
}


# Time to add distance information for this guy

# Exponential decay function for addressing relative impact. Scaled to 30 km = 1 unit (arbitrary)

exp_decay <- function(x) {
  # Scaled exponential decay function to scale polluter impact by distance
  
  x <- as.numeric(x)
  
  if (x < 0) {
    
    return(1)
  } else {
    
    scaled_value <- exp(-x / 30000)
    
    return(scaled_value)
  }
}


calc_impact_score <- function(esab, relevant_naics, robust) {
  # Calculates impact score by service area boundary and polluter type
  
  impact_column <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)){
    
    # Filter for those points with higher elevation first
    flow_heights <- dplyr::filter(relevant_naics, 
                                  relevant_naics$height > esab$height[n])
    
    
    # Additional filter for those polluters where there may be a valley between
    # the polluter and the CWS
    if (robust == TRUE) {
      
      # WIll index all columns to be removed
      index_col = vector(mode = 'numeric', length = nrow(flow_heights)) 
      
      for (x in 1:nrow(flow_heights)) {
        
        # Create a long connecting the two points
        point_1 <- st_coordinates(esab$centroids[n]) %>% as_vector
        point_2 <- st_coordinates(relevant_naics[x, ]) %>% as_vector
        
        line <- cbind(c(point_1[1], point_2[2]), c(point_1[2], point_2[2])) %>%
          st_linestring %>%
          st_sfc(crs = st_crs(ky_esab)) %>%
          st_sf %>%
          st_zm
        
        
        # Calculate minimum height along the valley
        valley <- raster::extract(na_elevation, line, 
                                  along = TRUE, cellnumbers = TRUE) %>%
          as_vector
        
        valley <- valley[!is.na(valley)]
        
        
        if (esab$height[n] > (min(valley) * (.3048 ^ 2)) | (max(valley)* (.3048 ^ 2)) > relevant_naics$height[x]) {
          
          index_col[x] <- x
          # If the minimum height is lesser than the height of the esab, by definition
          # there is a ridge between the esab and polluter, so no groundwater transport can occur
          
        }
      }
      
      index_col <- index_col[!is.na(index_col)]
      
      flow_heights <- flow_heights[-index_col, ]
    }
    
    
    # Filter for within the service area boundary, 
    # as these may be lower elevation but still affect transport
    flow_within <- relevant_naics[which(st_intersects(esab[n, ], relevant_naics, 
                                                      sparse = FALSE)), ]
    
    flow_heights <- rbind(flow_heights, flow_within) %>%
      unique
    
    
    if (nrow(flow_heights) > 0) {
      
      # Obtain distance matrix for esab to all polluters
      distance_mat <- st_distance(ky_esab[n,], flow_heights)
      
      for (x in 1:length(distance_mat)) {
        
        # Create exponential decay fit for distance matrix
        distance_mat[x] <- distance_mat[x] - esab$radius[n]
        
        distance_mat[x] <- exp_decay(distance_mat[x])
        
      }
      
      # Net impact from all polluters to a given esab
      eff_impact <- sum(distance_mat)
      
      impact_column[n] <- eff_impact
      
    } else {
      impact_column[n] <- 0
    }
  }
  
  return(impact_column)
}


ky_esab$transport_impact <- calc_impact_score(ky_esab, naics_pfas_transport, FALSE)
ky_esab$firefight_impact <- calc_impact_score(ky_esab, naics_pfas_firefight, FALSE)
ky_esab$waste_impact <- calc_impact_score(ky_esab, naics_pfas_waste, FALSE)
ky_esab$industry_impact <- calc_impact_score(ky_esab, naics_pfas_industry, FALSE)
ky_esab$military_impact <- calc_impact_score(ky_esab, us_military_bases, FALSE)


# Maps load fine, can change the parameter in question to confirm
mapview(ky_esab, zcol = 'military_impact')


# TODO- Find out why there are more pfas positive detections than there should be
# TODO- column for atmospheric deposition (within 20 km)
# TODO- averages for other criteria


# Atmospheric deposition column

ky_esab_atmosph <- ky_esab %>%
  st_buffer(dist = 20000)


atmospheric_dep_count <- function(esab, relevant_naics) {
  
  fac_count_list <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)) {
    # Loop over every esab in the region and count how many polluters are 
    # in this buffered region around the esab
    
    fac_count <- st_intersection(esab[n, ], relevant_naics) %>%
      nrow %>%
      as.numeric
    
    fac_count_list[n] <- fac_count
    
  }
  
  return(fac_count_list)
}

ky_esab$industry_atmos <- atmospheric_dep_count(ky_esab_atmosph, naics_pfas_industry)

rm(ky_esab_atmosph)

# Time to take averages for other relevant parameters

raster_average <- function(esab, relevant_raster) {
  
  raster_avg <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)) {
    # Loop over each esab and get the average raster value for it
    
    raster_values <- raster::extract(relevant_raster, esab[n, ], FUN = mean) %>%
      as_vector
    
    polygon_mean <- mean(raster_values)
    
    raster_avg[n] <- polygon_mean
  }
  
  return(raster_avg)  
  
}

getmode <- function(v) {
  # Function to calculate mode (most frequent value)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

raster_mode <- function(esab, relevant_raster) {
  
  raster_avg <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)) {
    # Loop over each esab and get the value corresponding to the type for it
    
    cation_str <- raster::extract(relevant_raster, esab[n, ], FUN = mean) %>%
      as_vector() %>%
      getmode()
    
    # Use this index value and go back to the matrix originally created to assign 
    # the categorical value for the variable
    
    cation_str <- which(nam_df == cation_str, arr.ind = T)
    cation_str <- nam_df$nam[cation_str[1]]
    
    raster_avg[n] <- as.character(cation_str)
  }
  
  raster_avg <- as.factor(raster_avg)
  
  return(raster_avg)  
  
}

ky_esab$rain <- raster_average(ky_esab, na_rain)
ky_esab$temp <- raster_average(ky_esab, na_temps)
ky_esab$soc <- raster_average(ky_esab, na_soc)
ky_esab$gw_reach <- raster_average(ky_esab, na_ground_rech)
ky_esab$abv_carb <- raster_average(ky_esab, na_abv_carb)
ky_esab$bel_carb <- raster_average(ky_esab, na_bel_carb)

ky_esab$clay_prc <- raster_average(ky_esab, clay_prc_raster)
ky_esab$H_ion_conc <- raster_average(ky_esab, H_ion_conc_raster)
ky_esab$cations <- raster_mode(ky_esab, na_cations_raster) %>% as.character()

# Grouping as numeric values are easier to handle and as the number of active/subactive  cations is very small

for (n in 1:length(ky_esab$cations)) {
  if (ky_esab$cations[n] == 'active') {
    ky_esab$cations[n] <- 1
  } else if (ky_esab$cations[n] == 'superactive'){
    ky_esab$cations[n] <- 1
  } else if (ky_esab$cations[n] == 'subactive'){
    ky_esab$cations[n] <- 0
  } else if (ky_esab$cations[n] == 'semiactive'){
    ky_esab$cations[n] <- 0
  } else {
    ky_esab$cations[n] <- 0
  }
}

ky_esab$cations <- ky_esab$cations %>% as.factor() %>% as.numeric() %>% -1

# TODO: Groundwater/ Surface water source information
# TODO: Alternate source information as this clearly is not working

sdwis_query_ky <- sdwis_query_ky %>% 
  dplyr::rename(PWS_ID = PWS.ID) %>% 
  dplyr::rename(source = Primary.Source) %>%
  dplyr::rename(population = Population.Served.Count) %>%
  dplyr::select(PWS_ID, source, population)

# Fiing esoteric categories and changing it into a dummy variable where
# 0 = surface water and 1 = ground water

sdwis_query_ky$source <- gsub("Surface water purchased", 0, sdwis_query_ky$source)
sdwis_query_ky$source <- gsub("Surface water", 0, sdwis_query_ky$source)

sdwis_query_ky$source <- gsub("Ground water purchased", 1, sdwis_query_ky$source)
sdwis_query_ky$source <- gsub("Groundwater under influence of surface water", 
                              1, sdwis_query_ky$source)
sdwis_query_ky$source <- gsub("Purchased ground water under influence of surface water source", 
                              1, sdwis_query_ky$source)
sdwis_query_ky$source <- gsub("Ground water", 1, sdwis_query_ky$source)

sdwis_query_ky$source <- as.numeric(sdwis_query_ky$source)

ky_esab <- ky_esab %>% left_join(sdwis_query_ky, by = 'PWS_ID')


# Filters for the rest

ky_esab$pH <- -log10(ky_esab$H_ion_conc)

ky_esab_PFAS_info <- ky_esab %>% filter(is.na(net_PFAS) == FALSE)

ky_esab_PFAS_info$PFAS_detect <- cut(ky_esab_PFAS_info$net_PFAS, breaks = c(-100, 1, 10000), labels = c(0, 1))


ky_esab_PFAS_info$PFAS_category <- cut(ky_esab_PFAS_info$net_PFAS, breaks = c(-100, 1, 10, 100000), 
                                       labels = c('none', 'low', 'high'))



# Subset only to relevant variables for ease of analysis

ky_esab_PFAS_info <- ky_esab_PFAS_info %>% dplyr::select('PWS_ID', 'transport_impact',
                                                         'firefight_impact', 'waste_impact', 'industry_impact',
                                                         'military_impact', 
                                                         'industry_atmos', 'rain', 'temp', 'soc', 'gw_reach', 
                                                         'bel_carb', 'clay_prc', 'pH', 'cations',
                                                         'PFAS_detect', 'PFAS_category', 'geometry', 'source') %>%
  relocate(source, .before = PFAS_detect) %>%
  as_tibble() # In data frame form so it doesn't confuse any algorithms. WIll put back later

# Impute any annoying missing values by median (because of categorical data)

ky_esab_PFAS_info$cations <- as.factor(ky_esab_PFAS_info$cations)

impute_model <- preProcess(ky_esab_PFAS_info[2:17] %>% as.data.frame(),
                           method = 'knnImpute')

ky_esab_PFAS_normal <- predict(impute_model, newdata = ky_esab_PFAS_info[2:17])

ky_esab_PFAS_normal$cations <- as.numeric(ky_esab_PFAS_normal$cations) %>% -1

ky_esab_PFAS_info <- na_mean(ky_esab_PFAS_info, option = 'median')

#-----------------
# Data analysis time
#-----------------


# Correlation matrices to determine some basic effects

ky_esab_PFAS_info$PFAS_detect <- as.numeric(as.character(ky_esab_PFAS_info$PFAS_detect))

ky_esab_PFAS_info$cations <- as.numeric(ky_esab_PFAS_info$cations) %>% -1

correlations <- cor(ky_esab_PFAS_info[2:17])
corrplot(correlations, method="circle")

ky_esab_PFAS_info$PFAS_detect <- as.factor(ky_esab_PFAS_info$PFAS_detect)


# Observe distributions with/without PFAS

x <- ky_esab_PFAS_normal %>% dplyr::select('transport_impact',
                                         'firefight_impact', 'waste_impact', 'industry_impact',
                                         'military_impact',
                                         'industry_atmos', 'rain', 'temp', 'soc', 'gw_reach', 
                                         'bel_carb', 'clay_prc', 'pH', 'cations', 'source') 


# Density plot

scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=ky_esab_PFAS_normal$PFAS_detect, plot="density", 
            scales=scales, strip=strip.custom(par.strip.text=list(cex=.7)),
            auto.key = list(columns = 2))

# Box plot

featurePlot(x = x, y = ky_esab_PFAS_normal$PFAS_detect,
            plot = 'box', strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = scales, auto.key = list(columns = 2))


# Logistic regression trials

logreg <- glm(PFAS_detect ~ transport_impact + firefight_impact + waste_impact + industry_impact + 
                military_impact +
                industry_atmos + rain + temp + soc + gw_reach +  
                bel_carb + clay_prc + pH + cations + source, data = ky_esab_PFAS_info, family = binomial)
summary(logreg)


# Predictions to assess accuracy (probably not good)

PFAS_predict <- predict(logreg, newdata = x, type = 'response')

pred_50 <- ifelse(PFAS_predict > 0.5, "1", "0")

nrow(dplyr::filter(ky_esab_PFAS_info, PFAS_detect == pred_50)) / nrow(ky_esab_PFAS_info)

# Confusion matrix- in left-right then bottom-down
# True positive, false negative, false positive, True negative

table(ky_esab_PFAS_info$PFAS_detect, pred_50)


# Around 85% accuracy with a threshold of .5, and around 75% accuracy with a threshold of .8. Not bad!


#-----------------
# Machine learning time
#-----------------


# Following guide at https://towardsdatascience.com/random-forest-in-r-f66adf80ec9
# Simple state to the problem

# First, put factors back in factor form

ky_esab_PFAS_info$PFAS_detect <- as.factor(ky_esab_PFAS_info$PFAS_detect)


# Caret guide- Gradient Boosting

set.seed(420)
inTraining <- createDataPartition(ky_esab_PFAS_info$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_info[ inTraining,]
pfas_testing  <- ky_esab_PFAS_info[-inTraining,]

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 10)

set.seed(2000)

gb_model <- caret::train(PFAS_detect ~ . ,
                         data = pfas_training[2:17],
                         method = 'gbm',
                         trControl = trainfolds,
                         verbose = FALSE)

PFAS_predict_gb <- predict(gb_model, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_gb)


# Random Forest


set.seed(9000)

rf_model <- caret::train(PFAS_detect ~ . ,
                         data = pfas_training[2:17],
                         method = 'rfRules',
                         trControl = trainfolds,
                         verbose = FALSE)

PFAS_predict_rf <- predict(rf_model, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_rf)


# Bayesian model

set.seed(6000)

by_model <- caret::train(PFAS_detect ~ . ,
                         data = pfas_training[2:17],
                         method = 'nb',
                         trControl = trainfolds,
                         verbose = FALSE)

PFAS_predict_by <- predict(by_model, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_by)



#-----------------
# Same analysis but with normalized data (and more models)
#-----------------

# Using the information generated by preProces

# Logistic regression

logreg_norm <- glm(PFAS_detect ~ transport_impact + firefight_impact + waste_impact + industry_impact + 
                     military_impact +  
                     industry_atmos + rain + temp + soc + gw_reach +  
                     bel_carb + clay_prc + pH + cations + source, data = ky_esab_PFAS_normal, family = binomial)
summary(logreg_norm)

PFAS_predict_norm <- predict(logreg_norm, newdata = ky_esab_PFAS_normal, type = 'response')

pred_50_norm <- ifelse(PFAS_predict > 0.5, "1", "0")

table(ky_esab_PFAS_normal$PFAS_detect, pred_50_norm) #Same predictive accuracy as previous


# Gradient boosting

set.seed(420)
inTraining <- createDataPartition(ky_esab_PFAS_normal$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_normal[ inTraining,]
pfas_testing  <- ky_esab_PFAS_normal[-inTraining,]

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 10)

set.seed(2000)

gb_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'gbm',
                              trControl = trainfolds,
                              verbose = FALSE)

PFAS_predict_gb_norm <- predict(gb_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_gb_norm)

# Random forest

set.seed(9000)

rf_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'rfRules',
                              trControl = trainfolds,
                              verbose = FALSE)

PFAS_predict_rf_norm <- predict(rf_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_rf_norm)

# Baynesian network

set.seed(6000)

by_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'nb',
                              trControl = trainfolds,
                              verbose = FALSE)

PFAS_predict_by_norm <- predict(by_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_by_norm)

# Boosted logistic regression

set.seed(790)

lg_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'LogitBoost',
                              trControl = trainfolds,
                              verbose = FALSE)

PFAS_predict_lg_norm <- predict(lg_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_lg_norm)

# Juan recommendation

set.seed(190)

Fe_model_norm <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'rFerns',
                              metric = 'Accuracy',
                              tuneLength = 20,
                              trControl = trainfolds)

PFAS_predict_Fe_norm <- predict(Fe_model_norm, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_Fe_norm)

# Comparing these models

resamps <- resamples(list(GBM = gb_model_norm,
                          RF = rf_model_norm,
                          BAY = by_model_norm,
                          LOG = lg_model_norm,
                          FER = Fe_model_norm))

difValues <- diff(resamps) # T-test for differences between models

trellis.par.set(caretTheme())
dotplot(difValues)       # Differences in accuracy 

trellis.par.set(caretTheme())
dotplot(resamps, metric = "Accuracy")

dotplot(resamps, metric = "Kappa")

# TODO: Parallel processing for impact factor loop

# TODO: Spear treasurer role stuff, birds project


# Plots of variable importance- can only do on certain types

varImp(gb_model_norm, scale = FALSE) %>%
  plot()

varImp(rf_model_norm, scale = FALSE) %>%
  plot()


# Time to use recursive feature selection on naive bayes and forests

# Bayes

set.seed(6970)


ctrl <- rfeControl(functions = nbFuncs,
                   method = "repeatedcv",
                   repeats = 10,
                   number = 5,
                   verbose = FALSE)

nbProfile <- rfe(x=pfas_training[1:15], y=pfas_training$PFAS_detect,
                 sizes = c(1:15),
                 rfeControl = ctrl)


# Forest

set.seed(6813)


ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 10,
                   number = 5,
                   verbose = FALSE)

rfProfile <- rfe(x=pfas_training[1:15], y=pfas_training$PFAS_detect,
                 sizes = c(1:15),
                 rfeControl = ctrl)











