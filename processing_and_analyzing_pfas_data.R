# Data analysis portion of project
# Created: 7/7/21
# Updated: 8/12/21
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
library(ggspatial)

library(concaveman)

library(Hmisc)
library(progress)
library(future)

library(furrr)
library(car)
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
library(MLeval)

library(RANN)
library(caretEnsemble)

library(imputeTS)
library(comprehenr)
library(sjPlot)
library(tmap)
# -----------------------

setwd('/data/water_lab/aakash_files/Data/')


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

co_contaminants <- rgdal::readOGR('Shapefiles/co_contaminants/WQD_co_contaminants.shp') %>%
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
  projectRaster(crs = '+init=epsg:3857')

# CSVs

sdwis_query_ky <- read.csv('CSVs/water_system_detail_ky_sdwis_Q42020.csv')

sdwis_query_tn <- read.csv('CSVs/TN_water_system_detail_Q4_2020_7_28_2021.csv')


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

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

RhpcBLASctl::omp_set_num_threads(1L)

# -----------------------
# Raster processing
# -----------------------

# Subset the full NA elevation map to be more manageable

na_elevation <- na_elevation[raster(xmn = as.numeric(st_bbox(naics_pfas)[1]),
                                    xmx = as.numeric(st_bbox(naics_pfas)[3]),
                                    ymn = as.numeric(st_bbox(naics_pfas)[2]),
                                    ymx = as.numeric(st_bbox(naics_pfas)[4])),
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

naics_pfas_rel_tn <- naics_pfas %>% dplyr::filter(FacState %in% c('TN', 'KY', 'AL', 'GA', 'NC', 'VA', 'AR', 'MS')) %>%
  separate(col = FacNAICSCo, sep = ' ', into = c('NAICS1', 'NAICS2', 'NAICS3', 'NAICS4', 'NAICS5', 'NAICS6', 
                                                 'NAICS7', 'NAICS8', 'NAICS9', 'NAICS10', 'NAICS11', 'NAICS12',
                                                 'NAICS13', 'NAICS14', 'NAICS15', 'NAICS16'),
           remove = FALSE)

naics_pfas_rel_ky <- naics_pfas %>% dplyr::filter(FacState %in% c( 'KY','VA','WV', 'OH', 'IN', 'IL', 'MO', 'TN')) %>%
  separate(col = FacNAICSCo, sep = ' ', into = c('NAICS1', 'NAICS2', 'NAICS3', 'NAICS4', 'NAICS5', 'NAICS6', 
                                                 'NAICS7', 'NAICS8', 'NAICS9', 'NAICS10', 'NAICS11', 'NAICS12',
                                                 'NAICS13', 'NAICS14', 'NAICS15', 'NAICS16'),
           remove = FALSE)

rm(naics_pfas)

# Give it height information

for (n in 1: nrow(naics_pfas_rel_ky)) {
  height <- st_coordinates(naics_pfas_rel_ky[n,])
  height <- raster::extract(na_elevation, height) 
  naics_pfas_rel_ky$height[n] <- height
}

for (n in 1: nrow(naics_pfas_rel_tn)) {
  height <- st_coordinates(naics_pfas_rel_tn[n,])
  height <- raster::extract(na_elevation, height) 
  naics_pfas_rel_tn$height[n] <- height
}

# This is the ugliest possible way to do this filtering properly, but I don't care anymore. It works
# Goes one-by-one and filters matching values in each NAICS parameter. Yay

dumb_filter <- function(state_naics, x) {
  filtered <- state_naics %>% dplyr::filter(NAICS1 %in% x | NAICS2 %in% x | NAICS3 %in% x | NAICS4 %in% x |
                                                 NAICS5 %in% x | NAICS6 %in% x | NAICS7 %in% x | NAICS8 %in% x |
                                                 NAICS9 %in% x | NAICS10 %in% x | NAICS11 %in% x | NAICS12 %in% x |
                                                 NAICS13 %in% x | NAICS14 %in% x | NAICS15 %in% x | NAICS16 %in% x)
  return(filtered)
}

# https://www.asdwa.org/wp-content/uploads/2020/05/ASDWA-PFAS-SWP-Technical-Appendix_FINAL3.pdf

naics_pfas_transport_ky <- dumb_filter(naics_pfas_rel_ky, c(4811, 4812, 481111, 481112, 481211, 481212, 481219))

naics_pfas_firefight_ky <- dumb_filter(naics_pfas_rel_ky, c(922160, 611519, 561990))

naics_pfas_waste_ky <- dumb_filter(naics_pfas_rel_ky, c(562111, 562112, 562119, 562211, 562212, 562213, 562219, 562991))

naics_pfas_industry_ky <- dumb_filter(naics_pfas_rel_ky, c(238320, 238330, 
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

naics_pfas_transport_tn <- dumb_filter(naics_pfas_rel_tn, c(4811, 4812, 481111, 481112, 481211, 481212, 481219))

naics_pfas_firefight_tn <- dumb_filter(naics_pfas_rel_tn, c(922160, 611519, 561990))

naics_pfas_waste_tn <- dumb_filter(naics_pfas_rel_tn, c(562111, 562112, 562119, 562211, 562212, 562213, 562219, 562991))

naics_pfas_industry_tn <- dumb_filter(naics_pfas_rel_tn, c(238320, 238330, 
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

us_military_bases <- us_military_bases[!duplicated(us_military_bases$SHAPE_Leng), ]

us_military_bases_ky <- us_military_bases %>% filter(STATE_TERR %in% c('Tennessee', 'Kentucky', 'Illinois', 'Ohio',
                                                                    'West Virginia', 'Virginia', 'Missouri', 
                                                                    'Indiana')) %>% 
  st_cast('POINT')

us_military_bases_ky <- us_military_bases_ky[!duplicated(us_military_bases_ky$SHAPE_Leng), ]

us_military_bases_tn <- us_military_bases %>% filter(STATE_TERR %in% c('Tennessee', 'Kentucky', 'Virginia', 'North Carolina', 
                                                                       'Georgia', 'Alabama', 'Mississippi', 'Missouri', 
                                                                       'Arkansas')) %>% 
  st_cast('POINT')

us_military_bases_tn <- us_military_bases_tn[!duplicated(us_military_bases_tn$SHAPE_Leng), ]

for (n in 1: nrow(us_military_bases_ky)) {
  height <- st_coordinates(us_military_bases_ky$geometry[n])
  height <- raster::extract(na_elevation, height) 
  us_military_bases_ky$height[n] <- height
}

for (n in 1: nrow(us_military_bases_tn)) {
  height <- st_coordinates(us_military_bases_tn$geometry[n])
  height <- raster::extract(na_elevation, height) 
  us_military_bases_tn$height[n] <- height
}


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
                                   dplyr::select(Location,  
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

n <- match("Louisville Water Co Payne Plant", ky_esab$Location)
ky_esab$net_PFAS[n] <- ky_esab$net_PFAS[n] + ky_pfas$net_PFAS[25]

# Fantastic, PFAS information is now coded in! Time to visualize this

mapview(ky_esab, legend = TRUE, at = c(0, 1, 10, 70), zcol = 'net_PFAS')

# Time to add elevation information to these points
# Height assignment process for ky_esab as that was not available before

raster_average <- function(esab, relevant_raster) {
  
  raster_avg <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)) {
    # Loop over each esab and get the average raster value for it
    
    raster_values <- raster::extract(relevant_raster, esab[n, ], FUN = mean, 
                                     small = TRUE, na.rm = TRUE) %>%
      as_vector
    
    polygon_mean <- mean(raster_values)
    
    raster_avg[n] <- polygon_mean
  }
  
  return(raster_avg)  
  
}


ky_esab$height <- raster_average(ky_esab, na_elevation)

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


impact_for_esab <- function(SA, relevant_naics, na_elevation, n_points){
  
  
  flow_heights <- dplyr::filter(relevant_naics, 
                                relevant_naics$height > SA$height)
  
  if (nrow(flow_heights) > 0) {
    
    index_col <-  vector(mode = 'numeric', length = nrow(flow_heights))
    
    for (x in 1:nrow(flow_heights)) {
      
      # Create a long connecting the two points
      point_1 <- st_coordinates(SA$centroids) %>% as_vector
      point_2 <- st_coordinates(relevant_naics[x, ]) %>% as_vector
      
      # Rather than calculate heights along a line, which is computationally intense
      # we will calculate heights along 10 points along the line and compare
      # This will significantly speed up runtime without sacrificing too much accuracy
      
      # Vector pointing from SA to polluter, with 1 / X magnitude (so X points can be generated)
      # Parameter to be increased if more points are to used for more accuracy is n_points
      
      esab_vector <- c(point_2[1] - point_1[1], point_2[2] - point_1[2]) / n_points
      
      
      valley <- vector(mode = 'numeric', length = (n_points - 1))
      
      for (m in 1:(n_points - 1)) {
        test_point <- (point_1 + (m * esab_vector)) %>%
          st_point %>%
          st_sfc(crs = st_crs(SA)) %>%
          st_sf %>%
          st_zm
        
        valley[m] <- raster::extract(na_elevation, test_point)
        
      }
      
      valley <- valley[!is.na(valley)]
      
      if (length(valley) == 0) {
        min <- SA$height
        max <- relevant_naics$height[x]
        
      } else {
        min <- min(valley) 
        max <- max(valley) 
      }
      
      
      if (SA$height > min | max > relevant_naics$height[x]) { 
        
        index_col[x] <- x
        # If the minimum height is lesser than the height of the esab, by definition
        # there is a ridge between the esab and polluter, so no groundwater transport can occur
        # Conversely, if the maximum height is greater than the height of the polluter, 
        # there must be a ridge between the polluter and esab, so no transport
        # this holds as we have already filtered only for those polluters who are higher
        # than the esabs
      }
    }
    # Removing those polluters that fail the transport criteria
    index_col <- index_col[!is.na(index_col)]
    
    flow_heights <- flow_heights[-index_col, ]
  }
  
  # Adds back in those polluters who are within the ESAB,
  # as they can still pollute even if 'below' the ESAB
  # TODO: Substitude to st_near and then get rid of atmospheric deposition column
  
  flow_within <- relevant_naics[which(st_intersects(SA, relevant_naics, 
                                                    sparse = FALSE)), ]
  
  flow_heights <- rbind(flow_heights, flow_within) %>%
    unique
  
  
  if (nrow(flow_heights) > 0) {
    
    # Obtain distance matrix for esab to all polluters
    distance_mat <- st_distance(SA, flow_heights)
    
    
    for (x in 1:length(distance_mat)) {
      
      
      # Create exponential decay fit for distance matrix
      distance_mat[x] <- distance_mat[x] - SA$radius
      
      distance_mat[x] <- exp_decay(distance_mat[x])
      
    }
    
    # Net impact from all polluters to a given esab
    eff_impact <- sum(distance_mat) %>% as.numeric
    
    return(eff_impact)
    
  } else {
    return(0)
  }
  
}

n_points = 20
ky_esab$transport_impact <- foreach(x = 1:nrow(ky_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
  impact_for_esab(ky_esab[x, ], naics_pfas_transport_ky, na_elevation, n_points)
}

ky_esab$firefight_impact <- foreach(x = 1:nrow(ky_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
  impact_for_esab(ky_esab[x, ], naics_pfas_firefight_ky, na_elevation, n_points)
}

ky_esab$waste_impact <- foreach(x = 1:nrow(ky_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
  impact_for_esab(ky_esab[x, ], naics_pfas_waste_ky, na_elevation, n_points)
}

ky_esab$industry_impact <- foreach(x = 1:nrow(ky_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
  impact_for_esab(ky_esab[x, ], naics_pfas_industry_ky, na_elevation, n_points)
}

ky_esab$military_impact <- foreach(x = 1:nrow(ky_esab), 
                                   .combine = 'c', 
                                   .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
 impact_for_esab(ky_esab[x, ], us_military_bases_ky, na_elevation, n_points)
}


fix_impact_outliers <- function(transport_column) {
  # As a result of the calculation methodology, some values are too high and I think skewing
  # my ML algorithms
  
  avg_impact <- mean(transport_column) # or just use quantile(.95) with no multiplication
  
  # Exclude by either greater than 1...or fixed multiple of mean value
  # .95 quantile / mean average for all but industry - 3.2; max- 3.9
  
  for (n in 1:length(transport_column)) {
    if (transport_column[n] > (3.2 * avg_impact)) {
      transport_column[n] = 3.2 * avg_impact
    }
  }
  
  return(transport_column)
}

ky_esab$transport_impact <- ky_esab$transport_impact %>% fix_impact_outliers()
ky_esab$firefight_impact <- ky_esab$firefight_impact %>% fix_impact_outliers()
ky_esab$waste_impact <- ky_esab$waste_impact %>% fix_impact_outliers()
ky_esab$industry_impact <- ky_esab$industry_impact %>% fix_impact_outliers()
ky_esab$military_impact <- ky_esab$military_impact %>% fix_impact_outliers()

# Maps load fine, can change the parameter in question to confirm
mapview(ky_esab, zcol = 'military_impact')


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

ky_esab$industry_atmos <- atmospheric_dep_count(ky_esab_atmosph, naics_pfas_industry_ky)
ky_esab$co_contam <- atmospheric_dep_count(ky_esab, co_contaminants) %>%
  as.logical() %>% # To make this about presence/absence of contaminants
  as.numeric()     # Rather than quantitative amounts as that wouldnt make sense

rm(ky_esab_atmosph)


# Time to take averages for other relevant parameters


getmode <- function(v) {
  # Function to calculate mode (most frequent value)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

raster_mode <- function(esab, relevant_raster) {
  
  raster_avg <- vector(mode = 'numeric', length = nrow(esab))
  
  for (n in 1:nrow(esab)) {
    # Loop over each esab and get the value corresponding to the type for it
    
    cation_str <- raster::extract(relevant_raster, esab[n, ], FUN = mean,
                                  small = TRUE, na.rm = TRUE) %>%
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
                                                         'PFAS_detect', 'geometry', 'source',
                                                         'co_contam') %>%
  relocate(source, .before = PFAS_detect) %>%
  relocate(co_contam, .before = PFAS_detect) %>%
  as_tibble() # In data frame form so it doesn't confuse any algorithms. WIll put back later

# Impute any annoying missing values by median (because of categorical data)

ky_esab_PFAS_info$cations <- as.factor(ky_esab_PFAS_info$cations)
ky_esab_PFAS_info$source <- as.factor(ky_esab_PFAS_info$source)
ky_esab_PFAS_info$co_contam <- as.factor(ky_esab_PFAS_info$co_contam)

impute_model <- preProcess(ky_esab_PFAS_info[2:18] %>% as.data.frame(),
                           method = 'knnImpute')

ky_esab_PFAS_normal <- predict(impute_model, newdata = ky_esab_PFAS_info[2:18])

ky_esab_PFAS_normal$cations <- as.numeric(ky_esab_PFAS_normal$cations) %>% -1
ky_esab_PFAS_normal$source <- as.numeric(ky_esab_PFAS_normal$source) %>% -1
ky_esab_PFAS_normal$co_contam <- as.numeric(ky_esab_PFAS_info$co_contam) %>% -1

#-----------------
# Data analysis time
#-----------------


# Correlation matrices to determine some basic effects

ky_esab_PFAS_normal$PFAS_detect <- ky_esab_PFAS_normal$PFAS_detect %>% as.numeric() %>% -1

correlations <- cor(ky_esab_PFAS_normal[1:17])
corrplot(correlations, method="circle", type = 'upper',
         tl.col = 'black', tl.srt = 45)

ky_esab_PFAS_normal$PFAS_detect <- ky_esab_PFAS_normal$PFAS_detect %>% as.factor()

# Observe distributions with/without PFAS

x <- ky_esab_PFAS_normal %>% dplyr::select('transport_impact',
                                         'firefight_impact', 'waste_impact', 'industry_impact',
                                         'military_impact',
                                         'industry_atmos', 'rain', 'temp', 'soc', 'gw_reach', 
                                         'bel_carb', 'clay_prc', 'pH', 'cations', 'source',
                                         'co_contam') 


# Density plot

scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=ky_esab_PFAS_normal$PFAS_detect, plot="density", 
            scales=scales, strip=strip.custom(par.strip.text=list(cex=.7)),
            auto.key = list(columns = 2))

# Box plot

featurePlot(x = x, y = ky_esab_PFAS_normal$PFAS_detect,
            plot = 'box', strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = scales, auto.key = list(columns = 2))

# Time to create a correlation matrix to observe variables a little more closely

ky_esab_PFAS_normal$PFAS_detect <- ky_esab_PFAS_normal$PFAS_detect %>% as.numeric() %>% -1

corr_mat <- cor(ky_esab_PFAS_normal, method = c('spearman')) %>% as_tibble() 
# returns correlation matrix with the R values. corr_mat[16, ] 
# will return column for just PFAS_detect

ky_esab_PFAS_normal$PFAS_detect <- ky_esab_PFAS_normal$PFAS_detect %>% as.factor()

# Time to get those values for which R-squared is greater than .05 (Arbitrary)

r_sq_mat <- corr_mat[16, ] ^ 2 %>% as_tibble()
r_sq_mat <- r_sq_mat[1:16]

relevant_variables <- c()

for (n in 1:length(r_sq_mat)) {
  if (r_sq_mat[n] > 0.05) {
    relevant_variables <- c(relevant_variables, c(colnames(r_sq_mat[n])))
  }
}  

# Confusion matrix- in left-right then bottom-down
# True positive, false negative, false positive, True negative

# also buffer to select sites near the state rather than simply TODO

# Following guide at https://towardsdatascience.com/random-forest-in-r-f66adf80ec9

#-----------------
# Machine learning with normalized data 
#-----------------

# Using the information generated by preProces

# Logistic regression

logreg_norm <- glm(PFAS_detect ~ transport_impact + firefight_impact + waste_impact + industry_impact + 
                     military_impact +  
                     industry_atmos + rain + temp + soc + gw_reach +  
                     bel_carb + clay_prc + pH + co_contam + 
                     cations + source, data = ky_esab_PFAS_normal, family = binomial)
summary(logreg_norm)

tab_model(logreg_norm, terms = c('rain', 'pH','industry_impact',
                                                'source', 'co_contam', 'waste_impact',
                                                'clay_prc'),
          pred.labels = c('Industry Impact', 'Waste Impact', 'Precipitation', 'Clay Percentage',
                          'pH', 'Co Contaminant Presence', 'Source'),
          string.ci = "|Conf. Int (95%)",
          string.p = "P-Value",
          p.style = 'stars',
          title = 'PFAS Presence')

PFAS_predict_norm <- predict(logreg_norm, newdata = ky_esab_PFAS_normal, type = 'response')

pred_50_norm <- ifelse(PFAS_predict_norm > 0.5, "1", "0")

table(ky_esab_PFAS_normal$PFAS_detect, pred_50_norm) 


# Based on significance values, add another variable to the relevant variables list

relevant_variables <- c(relevant_variables, 'pH', 'rain')


# Calculate variance inflation factor 

vif_values <- vif(logreg_norm)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
abline(v = 5, lwd = 3, lty = 2)
# bel_carb high negative correlation with temperature don't include in relevant variables?


# Test/Train splits

set.seed(420)
inTraining <- createDataPartition(ky_esab_PFAS_normal$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_normal[ inTraining,]
pfas_testing  <- ky_esab_PFAS_normal[-inTraining,]
testing_norm <- pfas_testing

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 5,
                           savePredictions = T, classProbs = T)

pfas_training$PFAS_detect <- pfas_training$PFAS_detect %>%
  as.factor() %>% as.numeric() %>% -1 %>%
  cut(breaks = c(-100, 0.5, 100), labels = c('no', 'yes'))


# Logistic regression (very simple, mostly for paper reasons)

set.seed(2021)

logreg_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'glm',
                              trControl = trainfolds,
                              family = 'binomial')

PFAS_predict_logreg_norm <- predict(logreg_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_logreg_norm)

# Gradient boosting

set.seed(2000)

gb_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'gbm',
                              trControl = trainfolds,
                              tuneLength = 15,
                              verbose = FALSE)

PFAS_predict_gb_norm <- predict(gb_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_gb_norm)

# Random forest

set.seed(9000)

rf_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'rf',
                              trControl = trainfolds,
                              tuneLength = 5,
                              verbose = FALSE)

PFAS_predict_rf_norm <- predict(rf_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_rf_norm)

# Bayesian network

set.seed(6000)

by_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'nb',
                              trControl = trainfolds,
                              tuneLength = 15,
                              verbose = FALSE)

PFAS_predict_by_norm <- predict(by_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_by_norm)

# Boosted logistic regression

set.seed(790)

lg_model_norm <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'LogitBoost',
                              trControl = trainfolds,
                              tuneLength = 50,
                              verbose = FALSE)

PFAS_predict_lg_norm <- predict(lg_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_lg_norm)

# Juan recommendation

set.seed(190)

Fe_model_norm <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'rFerns',
                              metric = 'Accuracy',
                              tuneLength = 25,
                              trControl = trainfolds)

PFAS_predict_Fe_norm <- predict(Fe_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_Fe_norm)

# Exreme gradient boosting (because that worked well)

set.seed(794)

unregister_dopar()

Xb_model_norm <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'xgbDART',
                              tuneLength = 5,
                              trControl = trainfolds)

doParallel::registerDoParallel(cl = my.cluster)

PFAS_predict_Xb_norm <- predict(Xb_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_Xb_norm)


# Neural net

set.seed(366727)

nn_model_norm <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'pcaNNet',
                              tuneLength = 20,
                              trControl = trainfolds)

PFAS_predict_nn_norm <- predict(nn_model_norm, newdata = testing_norm)

table(pfas_testing$PFAS_detect, PFAS_predict_nn_norm)


# Comparing these models

resamps <- resamples(list(GBM = gb_model_norm,
                          RF = rf_model_norm,
                          BAY = by_model_norm,
                          LOG = lg_model_norm,
                          FER = Fe_model_norm,
                          XBM = Xb_model_norm,
                          NET = nn_model_norm))

difValues <- diff(resamps) # T-test for differences between models

trellis.par.set(caretTheme())
dotplot(difValues)       # Differences in accuracy 

dotplot(resamps, metric = "Accuracy")

dotplot(resamps, metric = "Kappa")

bwplot(resamps, scales=scales)

# Subset of variables for poster usage as can't fit everything

resamps2 <- resamples(list(`Random Forest` = rf_model_norm,
                           `Baynesian Network` = by_model_norm,
                           `Extreme Gradient Boosting` = Xb_model_norm,
                           `Neural Net` = nn_model_norm))

dotplot(resamps2, metric = "Kappa")

# Plots of variable importance- can only do on certain types

varImp(gb_model_norm, scale = FALSE) %>%
  plot()

varImp(rf_model_norm, scale = FALSE) %>%
  plot()

varImp(Xb_model_norm, scale = FALSE) %>%
  plot()

# using top 5 from each varImp plot to add to relevant_variables, if not included already

relevant_variables <- c(relevant_variables, 'waste_impact')

# Time to use recursive feature selection on naive bayes and forests

# Bayes

set.seed(6970)


ctrl <- rfeControl(functions = nbFuncs,
                   method = "repeatedcv",
                   repeats = 10,
                   number = 5,
                   verbose = FALSE)

nbProfile <- rfe(x=ky_esab_PFAS_normal[1:15], y=ky_esab_PFAS_normal$PFAS_detect,
                 sizes = c(1:15),
                 rfeControl = ctrl)


# Forest

set.seed(6813)


ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 10,
                   number = 5,
                   verbose = FALSE)

rfProfile <- rfe(x=ky_esab_PFAS_normal[1:15], y=ky_esab_PFAS_normal$PFAS_detect,
                 sizes = c(1:15),
                 rfeControl = ctrl)

relevant_variables <- c(relevant_variables, 'gw_reach', 'temp')

#--------------
# Working with variable subsets
#--------------

# See if accuracy can be improved by working with the relevant_variables subset
# Rather than the whole dataset

relevant_variables <- relevant_variables[!duplicated(relevant_variables)]

ky_esab_PFAS_subset <- ky_esab_PFAS_normal %>%
  dplyr::select(c(relevant_variables, PFAS_detect))

# Logistic regression

logreg_sub <- glm(PFAS_detect ~ waste_impact + industry_impact + 
                     rain + temp + industry_atmos + gw_reach +
                     pH + source + co_contam,
                     data = ky_esab_PFAS_subset, family = binomial)
summary(logreg_sub)

PFAS_predict_sub <- predict(logreg_sub, newdata = ky_esab_PFAS_subset, type = 'response')

pred_50_sub <- ifelse(PFAS_predict_sub > 0.5, "1", "0")

table(ky_esab_PFAS_subset$PFAS_detect, pred_50_sub) 


# Machine learning time baybee
# Lets only use the most successful models from last time

set.seed(427)

inTraining <- createDataPartition(ky_esab_PFAS_subset$PFAS_detect, p = .8, list = FALSE)
pfas_training <- ky_esab_PFAS_subset[ inTraining,]
pfas_testing  <- ky_esab_PFAS_subset[-inTraining,]
testing_sub <- pfas_testing

trainfolds <- trainControl(method = 'repeatedcv', number = 5, repeats = 5,
                           savePredictions = T, classProbs = T)

pfas_training$PFAS_detect <- pfas_training$PFAS_detect %>%
  as.factor() %>% as.numeric() %>% -1 %>%
  cut(breaks = c(-100, 0.5, 100), labels = c('no', 'yes'))

# Random forest

set.seed(9027)

rf_model_sub <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'rf',
                              trControl = trainfolds,
                              tuneLength = 5,
                              verbose = FALSE)

PFAS_predict_rf_sub <- predict(rf_model_sub, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_rf_sub)

# Naive bayes

set.seed(6100)

by_model_sub <- caret::train(PFAS_detect ~ . ,
                              data = pfas_training,
                              method = 'nb',
                              trControl = trainfolds,
                              tuneLength = 15,
                              verbose = FALSE)

PFAS_predict_by_sub <- predict(by_model_sub, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_by_sub)

# Neural net

set.seed(366797)

nn_model_sub <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'pcaNNet',
                              tuneLength = 20,
                              trControl = trainfolds)

PFAS_predict_nn_sub <- predict(nn_model_sub, newdata = pfas_testing)

table(pfas_testing$PFAS_detect, PFAS_predict_nn_sub)

# Extreme gradient boosting


set.seed(794)

unregister_dopar()

Xb_model_sub <- caret::train(PFAS_detect ~ .,
                              data = pfas_training,
                              method = 'xgbDART',
                              tuneLength = 5,
                              trControl = trainfolds)

PFAS_predict_Xb_sub <- predict(Xb_model_sub, newdata = pfas_testing)
PFAS_predict_Xb_sub_prob <- predict(Xb_model_sub, newdata = pfas_testing,
                                    type='prob')

doParallel::registerDoParallel(cl = my.cluster)

table(pfas_testing$PFAS_detect, PFAS_predict_Xb_sub)


# Plots of accuracy

resamps_sub <- resamples(list(`Random Forest` = rf_model_sub,
                           `Baynesian Network` = by_model_sub,
                           `Extreme Gradient Boosting` = Xb_model_sub,
                           `Neural Net` = nn_model_sub))

dotplot(resamps_sub, metric = "Accuracy")
dotplot(resamps_sub, metric = "Kappa")

# Comparison to non-subsetted models

resamps_xbm <- resamples(list(`XBM norm` = Xb_model_norm,
                               `XBM sub` = Xb_model_sub)) # No significant difference- norm advantage
resamps_bay <- resamples(list(`BAY norm` = by_model_norm,
                              `BAY sub` = by_model_sub)) # Large difference, but not significant- sub advantage
resamps_for <- resamples(list(`RF norm` = rf_model_norm,
                              `RF sub` = rf_model_sub)) # No significant difference- norm advantage
resamps_net <- resamples(list(`NET norm` = nn_model_norm,
                              `NET sub` = nn_model_sub)) # No significant difference- norm advantage

dotplot(resamps_xbm, metric = "Accuracy")
dotplot(resamps_xbm, metric = "Kappa")

diff(resamps_net) %>% dotplot()

varImp(Xb_model_sub, scale = F) %>%
  plot()

resamples(list(`NET Full` = nn_model_norm,
               `NET Sub` = nn_model_sub,
               `XGB Full` = Xb_model_norm)) %>% dotplot(metric = 'Accuracy')

#-------------------
# Compute model accuracy with more rigor
#-------------------

norm_test_eval <- evalm(list(nn_model_norm, 
                            rf_model_norm, Xb_model_norm), 
                       gnames=c("Neural net", 
                                'Random Forest', 'Extreme Gradient Boosting'))


sub_test_eval <- evalm(list(by_model_sub, nn_model_sub, 
                            rf_model_sub, Xb_model_sub), 
                     gnames=c('Bayes', "Neural net", 
                              'Random Forest', 'Extreme Gradient Boosting'))

best_test_eval <- evalm(list(Xb_model_norm, nn_model_sub, nn_model_norm),
                      gnames=c('XGB- Full Data', 'NET- Subset', 'NET- Full Data'))


#-------------------
# Applying the model to KY data
#-------------------

ky_esab_predict <- ky_esab %>% dplyr::select('PWS_ID', 'transport_impact',
                                             'firefight_impact', 'waste_impact', 'industry_impact',
                                             'military_impact', 
                                             'industry_atmos', 'rain', 'temp', 'soc', 'gw_reach', 
                                             'bel_carb', 'clay_prc', 'pH', 'cations',
                                             'source', 'co_contam',
                                             'population', 'geometry') %>%
  as_tibble()

# Doing some similar preprocessing for the data to be predicted as was for training

ky_esab_predict$cations <- as.factor(ky_esab_predict$cations)
ky_esab_predict$source <- as.factor(ky_esab_predict$source)
ky_esab_predict$co_contam <- as.factor(ky_esab_predict$co_contam)


impute_model <- preProcess(ky_esab_predict[2:18] %>% as.data.frame(),
                           method = 'knnImpute')

ky_esab_predict <- predict(impute_model, newdata = ky_esab_predict[2:18])

ky_esab_predict$cations <- as.numeric(ky_esab_predict$cations) %>% -1
ky_esab_predict$source <- as.numeric(ky_esab_predict$source) %>% -1
ky_esab_predict$co_contam <- as.numeric(ky_esab_predict$co_contam) %>% -1

ky_esab_predict$cations <- na_mean(ky_esab_predict$cations, option = 'median')
ky_esab_predict$source <- na_mean(ky_esab_predict$source, option = 'median')
ky_esab_predict$co_contam <- na_mean(ky_esab_predict$co_contam, option = 'median')


apply_model_predictions <- function(esab, model, parent_esab) {
  esab <- esab %>%
    mutate(PWS_ID = parent_esab$PWS_ID) %>%
    mutate(population = parent_esab$population) %>%
    mutate(geometry = parent_esab$geometry) %>%
    mutate(PFAS_predict_prob = predict(model, newdata = esab, type = 'prob')) %>% 
    # Here is where we generate predictions using the model parameters
    mutate(PFAS_predict = predict(model, newdata = esab)) %>%
    st_as_sf()
  
  esab$PFAS_predict_prob <- esab$PFAS_predict_prob[2] %>% 
    unlist() %>%
    as.numeric()
  
  return(esab)
}

ky_esab_predict <- apply_model_predictions(ky_esab_predict, Xb_model_norm, ky_esab)

mapview(ky_esab_predict, zcol = 'PFAS_predict')
mapview(ky_esab_predict, zcol = 'PFAS_predict_prob',  at = c(0, .25, .5, .75, 1))


at_risk <- to_vec(for (`prob, value` in zip_lists(ky_esab_predict$PFAS_predict, ky_esab_predict$population)) 
  if(prob == 'yes') value) %>%
  na.omit() %>%
  str_replace_all(',','') %>%
  as.numeric() %>%
  sum()

ky_pop_on_record <- ky_esab$population %>% na.omit() %>% str_replace_all(',','') %>%
  as.numeric() %>% sum() 

# Represents 3.228 million people at risk, or ~72% of KY population. That is no good
# Let's try another measure of people at risk- multiply the probability of contamination
# with the population of a given drinking water system


cws_population <- ky_esab_predict$population %>% str_replace_all(',', '') %>%
  as.numeric()

at_risk_prob <- to_vec(for (`prob, value` in zip_lists(ky_esab_predict$PFAS_predict_prob,
                                                       cws_population)) prob * value) %>%
  na.omit() %>%
  str_replace_all(',','') %>%
  as.numeric() %>%
  sum()

# Oop, represents slightly more at ~73% of KY population. Also no good



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

tn_esab$circle_area <- st_buffer(tn_esab$centroids, tn_esab$radius)

mapview(list(tn_esab, tn_esab$circle_area))

# create circle polygons representing this simplified area now

# Generate the same KY parameters for TN- starting with impact scores


tn_esab$height <- raster_average(tn_esab, na_elevation)

tn_esab$transport_impact <- foreach(x = 1:nrow(tn_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
                                      impact_for_esab(tn_esab[x, ], naics_pfas_transport_tn, na_elevation, n_points)
                                    }

tn_esab$firefight_impact <- foreach(x = 1:nrow(tn_esab), 
                                    .combine = 'c', 
                                    .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
                                      impact_for_esab(tn_esab[x, ], naics_pfas_firefight_tn, na_elevation, n_points)
                                    }

tn_esab$waste_impact <- foreach(x = 1:nrow(tn_esab), 
                                .combine = 'c', 
                                .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
                                  impact_for_esab(tn_esab[x, ], naics_pfas_waste_tn, na_elevation, n_points)
                                }

tn_esab$industry_impact <- foreach(x = 1:nrow(tn_esab), 
                                   .combine = 'c', 
                                   .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
                                     impact_for_esab(tn_esab[x, ], naics_pfas_industry_tn, na_elevation, n_points)
                                   }

tn_esab$military_impact <- foreach(x = 1:nrow(tn_esab), 
                                   .combine = 'c', 
                                   .packages = c('magrittr', 'sf', 'tidyverse')) %dopar% {
                                     impact_for_esab(tn_esab[x, ], us_military_bases_tn, na_elevation, n_points)
                                   }

tn_esab$transport_impact <- tn_esab$transport_impact %>% fix_impact_outliers()
tn_esab$firefight_impact <- tn_esab$firefight_impact %>% fix_impact_outliers()
tn_esab$waste_impact <- tn_esab$waste_impact %>% fix_impact_outliers()
tn_esab$industry_impact <- tn_esab$industry_impact %>% fix_impact_outliers()
tn_esab$military_impact <- tn_esab$military_impact %>% fix_impact_outliers()

# Atmospheric impact and co-contamination

tn_esab_atmosph <- tn_esab %>%
  st_buffer(dist = 20000)

tn_esab_temp <- tn_esab %>% st_buffer(0)

tn_esab$industry_atmos <- atmospheric_dep_count(tn_esab_atmosph, naics_pfas_industry_ky)
tn_esab$co_contam <- atmospheric_dep_count(tn_esab_temp, co_contaminants) %>%
  as.logical() %>% # To make this about presence/absence of contaminants
  as.numeric()     # Rather than quantitative amounts as that wouldnt make sense

rm(tn_esab_atmosph)
rm(tn_esab_temp)

# Soil parameters

tn_esab$rain <- raster_average(tn_esab, na_rain)
tn_esab$temp <- raster_average(tn_esab, na_temps)
tn_esab$soc <- raster_average(tn_esab, na_soc)
tn_esab$gw_reach <- raster_average(tn_esab, na_ground_rech)
tn_esab$abv_carb <- raster_average(tn_esab, na_abv_carb)
tn_esab$bel_carb <- raster_average(tn_esab, na_bel_carb)

tn_esab$clay_prc <- raster_average(tn_esab, clay_prc_raster)
tn_esab$H_ion_conc <- raster_average(tn_esab, H_ion_conc_raster)

tn_esab$cations <- raster_mode(tn_esab, na_cations_raster) %>% as.character()
for (n in 1:length(tn_esab$cations)) {
  if (is.na(tn_esab$cations[n])) {
    tn_esab$cations[n] = NA
  } else if (tn_esab$cations[n] == 'active') {
    tn_esab$cations[n] <- 1
  } else if (tn_esab$cations[n] == 'superactive'){
    tn_esab$cations[n] <- 1
  } else if (tn_esab$cations[n] == 'subactive'){
    tn_esab$cations[n] <- 0
  } else if (tn_esab$cations[n] == 'semiactive'){
    tn_esab$cations[n] <- 0
  } else {
    tn_esab$cations[n] <- 0
  }
}

tn_esab$cations <- tn_esab$cations %>% as.factor() %>% as.numeric() %>% -1
tn_esab$cations <- na_mean(tn_esab$cations, option = 'median')

# Source information

sdwis_query_tn <- sdwis_query_tn %>% 
  dplyr::rename(PWS_ID = PWS.ID) %>% 
  dplyr::rename(source = Primary.Source) %>%
  dplyr::rename(population = Population.Served.Count) %>%
  dplyr::select(PWS_ID, source, population)


sdwis_query_tn$source <- gsub("Surface water purchased", 0, sdwis_query_tn$source)
sdwis_query_tn$source <- gsub("Surface water", 0, sdwis_query_tn$source)

sdwis_query_tn$source <- gsub("Ground water purchased", 1, sdwis_query_tn$source)
sdwis_query_tn$source <- gsub("Groundwater under influence of surface water", 
                              1, sdwis_query_tn$source)
sdwis_query_tn$source <- gsub("Purchased ground water under influence of surface water source", 
                              1, sdwis_query_tn$source)
sdwis_query_tn$source <- gsub("Ground water", 1, sdwis_query_tn$source)

sdwis_query_tn$source <- as.numeric(sdwis_query_tn$source)

tn_esab <- tn_esab %>% left_join(sdwis_query_tn, by = 'PWS_ID')

# Processing

tn_esab$pH <- -log10(tn_esab$H_ion_conc)

tn_esab_PFAS_info <- tn_esab %>% dplyr::select('PWS_ID', 'transport_impact',
                                                         'firefight_impact', 'waste_impact', 'industry_impact',
                                                         'military_impact', 
                                                         'industry_atmos', 'rain', 'temp', 'soc', 'gw_reach', 
                                                         'bel_carb', 'clay_prc', 'pH', 'cations',
                                                         'geometry', 'source', 'co_contam') %>%
  as_tibble() 


tn_esab_PFAS_info$cations <- as.factor(tn_esab_PFAS_info$cations)
tn_esab_PFAS_info$source <- as.factor(tn_esab_PFAS_info$source)
tn_esab_PFAS_info$co_contam <- as.factor(tn_esab_PFAS_info$co_contam)

impute_model_tn <- preProcess(tn_esab_PFAS_info[2:17] %>% as.data.frame(),
                           method = 'knnImpute')

tn_esab_PFAS_normal <- predict(impute_model_tn, newdata = tn_esab_PFAS_info[2:17])

tn_esab_PFAS_normal$cations <- as.numeric(tn_esab_PFAS_normal$cations) %>% -1
tn_esab_PFAS_normal$source <- as.numeric(tn_esab_PFAS_normal$source) %>% -1
tn_esab_PFAS_normal$co_contam <- as.numeric(tn_esab_PFAS_normal$co_contam) %>% -1

tn_esab_PFAS_normal$source <- na_mean(tn_esab_PFAS_normal$source, option = 'median')
tn_esab_PFAS_normal$co_contam <- na_mean(tn_esab_PFAS_normal$co_contam, option = 'median')


# Applying models to view results

tn_esab_predict <- apply_model_predictions(tn_esab_PFAS_normal, Xb_model_norm, tn_esab)
tn_esab_predict_agg <- apply_model_predictions(tn_esab_PFAS_normal, nn_model_norm, tn_esab)

mapview(tn_esab_predict, zcol = 'PFAS_predict')
mapview(tn_esab_predict, zcol = 'PFAS_predict_prob',  at = c(0, .25, .5, .75, 1))


at_risk_tn <- to_vec(for (`prob, value` in zip_lists(tn_esab_predict$PFAS_predict, tn_esab_predict$population)) 
  if(prob == 'yes') value) %>%
  na.omit() %>%
  str_replace_all(',','') %>%
  as.numeric() %>%
  sum()

sdwis_tn_pop <- sdwis_query_tn$population %>% na.omit() %>% str_replace_all(',','') %>%
  as.numeric() %>% sum() 

tn_pop_on_record <- tn_esab$Pop_Served %>% na.omit() %>% str_replace_all(',','') %>%
  as.numeric() %>% sum() 

cws_population_tn <- tn_esab_predict$population %>% str_replace_all(',', '') %>%
  as.numeric()

at_risk_prob_tn <- to_vec(for (`prob, value` in zip_lists(tn_esab_predict$PFAS_predict_prob,
                                                       cws_population_tn)) prob * value) %>%
  na.omit() %>%
  str_replace_all(',','') %>%
  as.numeric() %>%
  sum()


# -----------------------
# Maps for publication purposes
# -----------------------


# Making a map for publication reasons. Just PFAS detection across KY

ky_esab$`PFAS Presence` <- ky_esab$net_PFAS %>% as.logical() %>% as.numeric() %>%
  cut(breaks = c(-100, 0.5, 100), labels = c('Absent', 'Present')) %>% as.factor()


ggplot() +
  annotation_map_tile(type='hikebike', zoom=7) +
  geom_sf(data = ky_esab, aes(fill = `PFAS Presence`), col = NA) +
  theme_void() +
  ggtitle('Kentucky DEP PFAS Testing', subtitle = 'KY CWSs') +
  scale_fill_manual(values=c('#440154', '#fde725')) 


ggplot() +
  annotation_map_tile(type='hikebike',zoom=7) +
  geom_sf(data = ky_esab_predict, aes(fill = PFAS_predict_prob), col = NA) +
  geom_sf(data = tn_esab_predict, aes(fill = PFAS_predict_prob), col = NA) +
  theme_void() +
  ggtitle('Predictions of PFAS Contamination', 
          subtitle = 'KY and TN CWSs') +
  labs(fill='Detection Prob.') +
  scale_fill_viridis_c()
  
# For some reason this switches the logistic regression and xgboost lines when plotting
# had to fix it in microsoft paint. goddamn this society
xgb_eval <- evalm(list(Xb_model_norm, logreg_model_norm),
                  gnames=c('XGBoost', 'Logistic Regression'))


# Bless https://stackoverflow.com/questions/52200095/how-to-customize-the-importance-plot-generated-by-package-randomforest

xgb_imp <- varImp(Xb_model_norm, scale = F)
xgb_imp <- xgb_imp$importance %>% as.data.frame()
xgb_imp$varnames <- rownames(xgb_imp)
rownames(xgb_imp) <- NULL
xgb_imp$var_categ <- c('Soil information', 'Water source', 'Soil information',
                       'Soil information', 'Weather information',
                       'Distance to potential polluters', 'Co-contaminant presence', 'Soil information',
                       'Soil information', 'Weather information', 'Distance to potential polluters',
                       'Distance to potential polluters', 'Soil information', 
                       'Atmospheric deposition potential', 'Distance to potential polluters',
                       'Distance to potential polluters')
# As Points 
ggplot(xgb_imp[1:9, ], aes(x=reorder(varnames, Overall), y=Overall, color=as.factor(var_categ)))+ 
  geom_point() +
  geom_segment(aes(x=varnames,xend=varnames,y=0,yend=Overall)) +
  ylab("Relative Importance") +
  xlab("Predictor") +
  scale_colour_viridis_d(name="Variable Group", direction = -1) +
  ggtitle(' ', subtitle='Variable Importance') +
  coord_flip() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = rev(c('Soil pH', 'Source', 'Clay Percentage',
                                  'Belowground Carbon', 'Precipitation',
                                  'Industry Impact', 'Co Contaminant Presence',
                                  'Soil Organic Carbon', 'Groundwater Recharge')))


# As box plot
ggplot(xgb_imp[1:9, ], aes(x=reorder(varnames, Overall), weight=Overall, fill = as.factor(var_categ))) + 
  geom_bar() +
  scale_fill_discrete(name="Variable Group") +
  scale_colour_viridis_d() +
  ylab("Relative Importance") +
  xlab("Predictor") +
  coord_flip()

# Predictions for TDH People

predictions = data.frame(tn_esab_predict$PWS_ID, 
                            tn_esab_predict$PFAS_predict_prob,
                            tn_esab_predict_agg$PFAS_predict_prob)
colnames(predictions) <- c('PWS_ID', 'Conservative', 'Aggressive')
write_csv(predictions, 'TN_PFAS_predictions.csv')

# Model accuracy assessment

assess_model <- function(ml_model) {
  mean_acc <- ml_model$resample$Accuracy %>% mean()
  acc_stdev <- ml_model$resample$Accuracy %>% sd()
  min_acc <- (mean_acc - 1.96 * acc_stdev / sqrt(ml_model$resample$Accuracy %>% length())) %>%
    signif(digits = 2)
  max_acc <- (mean_acc + 1.96 * acc_stdev / sqrt(ml_model$resample$Accuracy %>% length())) %>%
    signif(digits = 2)
  acc_range <- paste0(min_acc, '-', max_acc)
  
  mean_kapp <- ml_model$resample$Kappa %>% mean()
  kapp_stdev <- ml_model$resample$Kappa %>% sd()
  min_kapp <- (mean_kapp - 1.96 * kapp_stdev / sqrt(ml_model$resample$Kappa %>% length())) %>%
    signif(digits = 2)
  max_kapp <- (mean_kapp + 1.96 * kapp_stdev / sqrt(ml_model$resample$Kappa %>% length())) %>%
    signif(digits = 2)
  kapp_range <- paste0(min_kapp, '-', max_kapp)
  
  auroc <- evalm(ml_model)$stdres$`Group 1`$CI[13]
  
  print(paste0('Accuracy range: ', acc_range))
  print(paste0('Kappa range: ', kapp_range))
  print(paste0('AUROC range: ', auroc))
}

 