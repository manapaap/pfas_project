# Testing usage of foreach to create a parallelized version of our current function mess


# Idea- use foreach to modify each value but with a function that is not parallelized

impact_for_esab <- function(SA, relevant_naics, na_elevation){
  
  
  flow_heights <- dplyr::filter(relevant_naics, 
                                relevant_naics$height > SA$height)
  
  index_col = vector(mode = 'numeric', length = nrow(flow_heights)) 
  
  print(nrow(flow_heights))
  
  for (x in 1:nrow(flow_heights)) {
    
    # Create a long connecting the two points
    point_1 <- st_coordinates(SA$centroids) %>% as_vector
    point_2 <- st_coordinates(relevant_naics[x, ]) %>% as_vector
    
    line <- cbind(c(point_1[1], point_2[1]), c(point_1[2], point_2[2])) %>%
      st_linestring %>%
      st_sfc(crs = st_crs(SA)) %>%
      st_sf %>%
      st_zm
    
    # Rather than calculate heights along a line, which is computationally intense
    # we will calculate heights along 10 points along the line and compare
    # This will significantly speed up runtime without sacrificing too much accuracy
     
    # Vector pointing from SA to polluter, with 1 / X magnitude (so X points can be generated)
    n_points <- 10 # Parameter to be increased if more points are to used for more accuracy
    
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
    
    print(SA$height)
    print(min)
    print(max)
    print(relevant_naics$height[x])
    
    if (SA$height > min) { #  | max > relevant_naics$height[x] but need better elevation raster TODO
      
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
    
    print(index_col)
    
    flow_heights <- flow_heights[-index_col, ]
    
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

  
impact_column <- vector(mode = 'numeric', length = nrow(esab))
)

system.time(
  ky_esab$transport_impact <- foreach(x = 1:nrow(ky_esab), .combine = 'c', .packages = c('magrittr', 
                                                         'sf', 'tidyverse')) %dopar% {
      
    impact_for_esab(ky_esab[x, ], naics_pfas_transport, na_elevation)
    
  }
)

system.time(
  ky_esab$industry_impact <- foreach(x = 1:nrow(ky_esab), .combine = 'c', .packages = c('magrittr', 
                                                         'sf', 'tidyverse')) %dopar% {
                                                           
   impact_for_esab(ky_esab[x, ], naics_pfas_industry, na_elevation)
   
 }
)



  
# Testing by using individual item
impact_for_esab(ky_esab[4, ], naics_pfas_transport, na_elevation) 


